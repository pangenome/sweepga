use anyhow::Result;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Mutex;

use crate::compact_index::{CompactPafIndex, CompactRecord};
use crate::paf_filter::{FilterConfig, FilterMode};

pub struct EfficientPafFilter {
    config: FilterConfig,
    keep_self: bool,
    scaffolds_only: bool,
}

impl EfficientPafFilter {
    pub fn new(config: FilterConfig) -> Self {
        EfficientPafFilter {
            config,
            keep_self: false,
            scaffolds_only: false,
        }
    }

    pub fn with_keep_self(mut self, keep: bool) -> Self {
        self.keep_self = keep;
        self
    }

    pub fn with_scaffolds_only(mut self, scaffolds: bool) -> Self {
        self.scaffolds_only = scaffolds;
        self
    }

    /// Extract prefix from sequence name for grouping
    fn extract_prefix(&self, name: &str) -> String {
        if !self.config.skip_prefix {
            return name.to_string();
        }

        if let Some(pos) = name.rfind(self.config.prefix_delimiter) {
            name[..pos].to_string()
        } else {
            name.to_string()
        }
    }

    /// First pass: Build compact index from PAF file
    fn build_index<P: AsRef<Path>>(&self, input_path: P) -> Result<CompactPafIndex> {
        let file = File::open(input_path)?;
        let mut reader = BufReader::new(file);
        let mut index = CompactPafIndex::new();

        let mut line = String::new();
        let mut byte_offset: u64 = 0;

        while reader.read_line(&mut line)? > 0 {
            let line_bytes = line.len() as u32;

            // Skip comments
            if !line.starts_with('#') && !line.trim().is_empty() {
                let fields: Vec<&str> = line.split('\t').collect();

                if fields.len() >= 12 {
                    let query_name = fields[0];
                    let query_len: u32 = fields[1].parse().unwrap_or(0);
                    let query_start: u32 = fields[2].parse().unwrap_or(0);
                    let query_end: u32 = fields[3].parse().unwrap_or(0);
                    let strand = fields[4].chars().next().unwrap_or('+');
                    let target_name = fields[5];
                    let target_len: u32 = fields[6].parse().unwrap_or(0);
                    let target_start: u32 = fields[7].parse().unwrap_or(0);
                    let target_end: u32 = fields[8].parse().unwrap_or(0);
                    let num_matches: u32 = fields[9].parse().unwrap_or(0);
                    let block_length: u32 = fields[10].parse().unwrap_or(0);

                    // Apply minimum block length filter
                    if block_length >= self.config.min_block_length {
                        // Apply self-mapping filter
                        if self.keep_self || query_name != target_name {
                            index.add_record(
                                byte_offset,
                                line_bytes,
                                query_name,
                                target_name,
                                query_start,
                                query_end,
                                target_start,
                                target_end,
                                block_length,
                                strand,
                            );
                        }
                    }
                }
            }

            byte_offset += line_bytes as u64;
            line.clear();
        }

        Ok(index)
    }

    /// Apply filtering algorithm on compact records
    fn apply_filters(&self, index: &CompactPafIndex) -> Result<HashSet<usize>> {
        let mut passing: HashSet<usize> = HashSet::new();
        let num_records = index.len();

        // For now, simple example: keep all records
        // TODO: Implement actual filtering logic working on CompactRecord

        // Collect records for filtering
        let mut records: Vec<CompactRecord> = Vec::new();
        for i in 0..num_records {
            if let Some(record) = index.get_record(i) {
                records.push(record);
            }
        }

        // Apply primary mapping filter
        let filtered_indices = self.apply_mapping_filter(index, &records)?;

        // If no scaffold filtering, return filtered indices
        if self.config.min_scaffold_length == 0 {
            return Ok(filtered_indices);
        }

        // Apply scaffold filtering
        let scaffold_indices = self.apply_scaffold_filter(index, &filtered_indices, &records)?;

        Ok(scaffold_indices)
    }

    /// Apply mapping filter to compact records (parallelized per chromosome pair)
    fn apply_mapping_filter(
        &self,
        index: &CompactPafIndex,
        records: &[CompactRecord],
    ) -> Result<HashSet<usize>> {
        // Group by prefix pairs if needed
        if self.config.mapping_filter_mode == FilterMode::ManyToMany {
            // No filtering - keep all
            return Ok(records.iter().map(|r| r.idx).collect());
        }

        // Group by prefix pairs (chromosome pairs)
        let mut prefix_groups: HashMap<(String, String), Vec<CompactRecord>> = HashMap::new();

        for &record in records {
            let query_name = index.get_name(record.query_id).unwrap_or("");
            let target_name = index.get_name(record.target_id).unwrap_or("");
            let query_prefix = self.extract_prefix(query_name);
            let target_prefix = self.extract_prefix(target_name);

            prefix_groups
                .entry((query_prefix, target_prefix))
                .or_insert_with(Vec::new)
                .push(record);
        }

        // Process each chromosome pair in parallel
        let passing_mutex = Mutex::new(HashSet::new());

        prefix_groups
            .into_par_iter()
            .try_for_each(|(_prefix_pair, group)| -> Result<()> {
                let filtered = self.filter_compact_group(group, self.config.mapping_filter_mode)?;
                let mut passing = passing_mutex.lock().unwrap();
                for record in filtered {
                    passing.insert(record.idx);
                }
                Ok(())
            })?;

        Ok(passing_mutex.into_inner().unwrap())
    }

    /// Filter a group of compact records based on mode
    fn filter_compact_group(
        &self,
        mut group: Vec<CompactRecord>,
        mode: FilterMode,
    ) -> Result<Vec<CompactRecord>> {
        // Sort by block length (descending)
        group.sort_by_key(|r| std::cmp::Reverse(r.block_length));

        match mode {
            FilterMode::OneToOne => {
                let mut kept_queries = HashSet::new();
                let mut kept_targets = HashSet::new();
                let mut result = Vec::new();

                for record in group {
                    if !kept_queries.contains(&record.query_id) &&
                       !kept_targets.contains(&record.target_id) {
                        kept_queries.insert(record.query_id);
                        kept_targets.insert(record.target_id);
                        result.push(record);
                    }
                }
                Ok(result)
            }
            FilterMode::OneToMany => {
                // Keep best per query, allowing multiple per target
                let mut kept_queries = HashSet::new();
                let mut result = Vec::new();

                for record in group {
                    if !kept_queries.contains(&record.query_id) {
                        kept_queries.insert(record.query_id);
                        result.push(record);
                    }
                }
                Ok(result)
            }
            FilterMode::ManyToMany => Ok(group),
        }
    }

    /// Apply scaffold filter (placeholder - needs full implementation)
    fn apply_scaffold_filter(
        &self,
        index: &CompactPafIndex,
        input_indices: &HashSet<usize>,
        all_records: &[CompactRecord],
    ) -> Result<HashSet<usize>> {
        // TODO: Implement scaffold filtering on compact records
        // For now, return input as-is
        Ok(input_indices.clone())
    }

    /// Second pass: Copy passing records from input to output
    fn write_filtered_output<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: P,
        index: &CompactPafIndex,
        passing_indices: &HashSet<usize>,
    ) -> Result<()> {
        let mut input = BufReader::new(File::open(input_path)?);
        let mut output = BufWriter::new(File::create(output_path)?);

        // Sort passing indices by byte offset for sequential reading
        let mut sorted_indices: Vec<usize> = passing_indices.iter().copied().collect();
        sorted_indices.sort_by_key(|&idx| index.get_record(idx).unwrap().byte_offset);

        for idx in sorted_indices {
            if let Some(record) = index.get_record(idx) {
                // Seek to record start
                input.seek(SeekFrom::Start(record.byte_offset))?;

                // Read the exact number of bytes for this record
                let mut buffer = vec![0u8; record.byte_length as usize];
                input.read_exact(&mut buffer)?;

                // Write to output (could add st:Z: tag here if needed)
                output.write_all(&buffer)?;
            }
        }

        output.flush()?;
        Ok(())
    }

    /// Main entry point: filter PAF file efficiently
    pub fn filter_paf<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: P,
    ) -> Result<()> {
        eprintln!("[EFFICIENT] First pass: building compact index...");
        let index = self.build_index(&input_path)?;
        eprintln!("[EFFICIENT] Indexed {} records", index.len());

        eprintln!("[EFFICIENT] Applying filters...");
        let passing_indices = self.apply_filters(&index)?;
        eprintln!("[EFFICIENT] {} records passed filtering", passing_indices.len());

        eprintln!("[EFFICIENT] Second pass: writing output...");
        self.write_filtered_output(&input_path, &output_path, &index, &passing_indices)?;

        Ok(())
    }
}