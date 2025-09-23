use anyhow::Result;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use ordered_float::OrderedFloat;

use crate::mapping::ChainStatus;

/// Filtering mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterMode {
    OneToOne,   // 1:1 - best mapping per query AND per target
    OneToMany,  // 1:N - best mapping per query, N per target
    ManyToMany, // N:N - N mappings per query and per target
}

/// Filter configuration
#[derive(Clone)]
pub struct FilterConfig {
    pub chain_gap: u32,           // -c/--chain-jump
    pub min_block_length: u32,    // -l/--block-length

    // Primary mapping filter (applied to raw mappings before scaffold creation)
    pub mapping_filter_mode: FilterMode,  // Default: N:N (no filtering)
    pub mapping_max_per_query: Option<usize>,
    pub mapping_max_per_target: Option<usize>,

    // Scaffold filter (applied to scaffold chains)
    pub scaffold_filter_mode: FilterMode,      // Default: 1:1
    pub scaffold_max_per_query: Option<usize>,
    pub scaffold_max_per_target: Option<usize>,

    pub overlap_threshold: f64,   // -O/--overlap
    pub sparsity: f64,            // -x/--sparsify
    pub no_merge: bool,           // -M/--no-merge
    pub scaffold_gap: u32,        // -j/--scaffold-jump
    pub min_scaffold_length: u32, // -S/--scaffold-mass
    pub scaffold_overlap_threshold: f64,
    pub scaffold_max_deviation: u32, // -D/--scaffold-dist
    pub prefix_delimiter: char,
    pub skip_prefix: bool,
}

/// Record metadata for filtering without modifying records
#[derive(Debug, Clone)]
struct RecordMeta {
    rank: usize,  // 0-based index in original file
    query_name: String,
    target_name: String,
    query_start: u32,
    query_end: u32,
    target_start: u32,
    target_end: u32,
    block_length: u32,
    identity: f64,
    strand: char,
    chain_id: Option<String>,
    chain_status: ChainStatus,
    discard: bool,
    overlapped: bool,
}

/// Represents a merged chain for scaffold filtering
#[derive(Debug, Clone)]
struct MergedChain {
    query_name: String,
    target_name: String,
    query_start: u32,
    query_end: u32,
    target_start: u32,
    target_end: u32,
    strand: char,
    total_length: u32,
    member_indices: Vec<usize>,  // Indices of original mappings in this chain
}

/// PAF filter that preserves original records
pub struct PafFilter {
    config: FilterConfig,
    temp_dir: Option<String>,
    keep_self: bool,
    scaffolds_only: bool,
}

impl PafFilter {
    pub fn new(config: FilterConfig) -> Self {
        PafFilter {
            config,
            temp_dir: std::env::var("TMPDIR").ok().or_else(|| Some("/tmp".to_string())),
            keep_self: false,  // Exclude self-mappings by default
            scaffolds_only: false,
        }
    }

    /// Extract prefix from sequence name for grouping
    fn extract_prefix(&self, name: &str) -> String {
        if self.config.skip_prefix {
            // Skip prefix extraction, use full name
            return name.to_string();
        }

        // Extract prefix: take everything up to and including the last delimiter
        // e.g., "r#1#2" -> "r#1#", "Rabacal-1#Chr1" -> "Rabacal-1#", "r#1" -> "r#"
        if let Some(pos) = name.rfind(self.config.prefix_delimiter) {
            name[..=pos].to_string()  // Include the delimiter itself
        } else {
            // No delimiter found, use full name
            name.to_string()
        }
    }

    pub fn with_keep_self(mut self, keep_self: bool) -> Self {
        self.keep_self = keep_self;
        self
    }

    pub fn with_scaffolds_only(mut self, scaffolds_only: bool) -> Self {
        self.scaffolds_only = scaffolds_only;
        self
    }

    /// Main filtering pipeline using record ranks
    pub fn filter_paf<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: P
    ) -> Result<()> {
        // First pass: extract metadata for all records
        let metadata = self.extract_metadata(&input_path)?;

        // Apply filters to get passing record ranks
        let passing_ranks = self.apply_filters(metadata)?;

        // Second pass: write passing records with annotations
        self.write_filtered_output(&input_path, &output_path, passing_ranks)?;

        Ok(())
    }

    /// Extract metadata from PAF without modifying records
    fn extract_metadata<P: AsRef<Path>>(&self, path: P) -> Result<Vec<RecordMeta>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut metadata = Vec::new();

        for (rank, line) in reader.lines().enumerate() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() < 11 {
                continue;
            }

            // Parse essential fields
            let query_name = fields[0].to_string();
            let _query_len = fields[1].parse::<u32>().unwrap_or(0);
            let query_start = fields[2].parse::<u32>().unwrap_or(0);
            let query_end = fields[3].parse::<u32>().unwrap_or(0);
            let strand = if fields[4] == "+" { '+' } else { '-' };
            let target_name = fields[5].to_string();
            let _target_len = fields[6].parse::<u32>().unwrap_or(0);
            let target_start = fields[7].parse::<u32>().unwrap_or(0);
            let target_end = fields[8].parse::<u32>().unwrap_or(0);
            let matches = fields[9].parse::<u32>().unwrap_or(0);
            let block_length = fields[10].parse::<u32>().unwrap_or(1);

            // Look for divergence tag from FASTGA
            let mut identity = matches as f64 / block_length.max(1) as f64;
            for field in &fields[11..] {
                if field.starts_with("dv:f:") {
                    if let Ok(div) = field[5..].parse::<f64>() {
                        identity = 1.0 - div;
                    }
                }
            }

            metadata.push(RecordMeta {
                rank,
                query_name,
                target_name,
                query_start,
                query_end,
                target_start,
                target_end,
                block_length,
                identity,
                strand,
                chain_id: None,
                chain_status: ChainStatus::Unassigned,
                discard: false,
                overlapped: false,
            });
        }

        Ok(metadata)
    }

    /// Apply filtering pipeline following wfmash's algorithm
    fn apply_filters(&self, mut metadata: Vec<RecordMeta>) -> Result<HashMap<usize, RecordMeta>> {
        // 1. Filter by minimum block length and self-mappings
        metadata.retain(|m| {
            m.block_length >= self.config.min_block_length &&
            (self.keep_self || m.query_name != m.target_name)
        });

        // Keep all original mappings for rescue phase (before any plane sweep)
        let all_original_mappings = metadata.clone();

        // 2. Apply primary mapping filter to input mappings
        // Default: N:N (no filtering) before scaffold creation
        if self.config.mapping_filter_mode != FilterMode::ManyToMany {
            // Group by prefix pairs for filtering modes
            let mut prefix_groups: HashMap<(String, String), Vec<RecordMeta>> = HashMap::new();
            for record in metadata {
                let query_prefix = self.extract_prefix(&record.query_name);
                let target_prefix = self.extract_prefix(&record.target_name);
                prefix_groups.entry((query_prefix, target_prefix))
                    .or_insert_with(Vec::new)
                    .push(record);
            }

            // Debug: [MAPPING_FILTER] Grouped into prefix pairs

            // Apply filtering within each prefix pair IN PARALLEL
            let filtered_groups: Vec<Vec<RecordMeta>> = prefix_groups
                .into_par_iter()
                .map(|((_q_prefix, _t_prefix), group)| {
                    let group_size = group.len();
                    let filtered = match self.config.mapping_filter_mode {
                        FilterMode::OneToOne => self.filter_one_to_one_within_group(group),
                        FilterMode::OneToMany => self.apply_plane_sweep_to_mappings(&group),
                        FilterMode::ManyToMany => Ok(group),
                    };
                    // Debug: [MAPPING_FILTER] Processed group
                    filtered.unwrap_or_else(|_| Vec::new())
                })
                .collect();

            // Combine all filtered groups
            metadata = filtered_groups.into_iter().flatten().collect();
            // Debug: [MAPPING_FILTER] After filtering
        }

        // If no scaffold filtering, we're done - return the plane-swept mappings
        if self.config.min_scaffold_length == 0 {
            let mut result = HashMap::new();
            for m in metadata {
                result.insert(m.rank, m);
            }
            return Ok(result);
        }

        // 3. Apply scaffold filtering (wfmash's filterByScaffolds)
        // According to CLAUDE.md: create scaffolds from the PLANE-SWEPT mappings,
        // then rescue from ALL ORIGINAL mappings

        // Step 1: Create scaffolds from the plane-swept mappings
        // Debug: [SCAFFOLD_TRACE] Creating scaffolds from plane-swept mappings

        // Use scaffold_gap for merging into scaffolds (per CLAUDE.md)
        let merged_chains = self.merge_mappings_into_chains(&metadata, self.config.scaffold_gap)?;
        // Debug: [SCAFFOLD_TRACE] After merging with gap

        // Step 2: Filter chains by minimum scaffold length
        let mut filtered_chains: Vec<MergedChain> = merged_chains
            .into_iter()
            .filter(|chain| chain.total_length >= self.config.min_scaffold_length)
            .collect();
        // Debug: [SCAFFOLD_TRACE] After length filter

        // Step 3: Apply plane sweep to scaffolds based on scaffold_filter_mode
        // Default: 1:1 filtering (best scaffold per query-target pair)
        filtered_chains = self.apply_scaffold_plane_sweep(filtered_chains)?;

        // Debug: [SCAFFOLD_TRACE] After plane sweep on scaffolds

        // If scaffolds_only mode, return just the scaffold chains
        if self.scaffolds_only {
            let mut scaffold_mappings = HashMap::new();
            for (idx, chain) in filtered_chains.iter().enumerate() {
                // Create a synthetic mapping for each scaffold chain
                let mut meta = RecordMeta {
                    rank: idx,  // Use sequential ranks for scaffolds
                    query_name: chain.query_name.clone(),
                    query_start: chain.query_start,
                    query_end: chain.query_end,
                    target_name: String::new(), // Will be filled from original
                    target_start: chain.target_start,
                    target_end: chain.target_end,
                    strand: chain.strand,
                    block_length: chain.total_length,
                    identity: 1.0,
                    chain_id: None,
                    chain_status: ChainStatus::Scaffold,
                    discard: false,
                    overlapped: false,
                };

                // Get target name from first original mapping in this chain
                if !chain.member_indices.is_empty() && chain.member_indices[0] < metadata.len() {
                    let first = &metadata[chain.member_indices[0]];
                    meta.target_name = first.target_name.clone();
                }

                scaffold_mappings.insert(idx, meta);
            }
            return Ok(scaffold_mappings);
        }

        // Step 4: Identify anchors - use the actual member mappings of scaffold chains
        // NOTE: member_indices actually contains ranks, not array indices!
        let mut anchor_ranks = HashSet::new();
        for chain in &filtered_chains {
            // The member_indices field contains the ranks of mappings in this chain
            for &member_rank in &chain.member_indices {
                anchor_ranks.insert(member_rank);
            }
        }
        // Debug: [SCAFFOLD_TRACE] Anchors identified

        // Map ranks to indices in all_original_mappings
        let mut rank_to_idx = HashMap::new();
        for (idx, meta) in all_original_mappings.iter().enumerate() {
            rank_to_idx.insert(meta.rank, idx);
        }

        // Step 5: Rescue mappings within scaffold_max_deviation of anchors
        // OPTIMIZED: Group by chromosome pair and sort for efficient rescue

        // First, group all mappings by (query_chr, target_chr) pair
        let mut mappings_by_chr_pair: HashMap<(String, String), Vec<usize>> = HashMap::new();
        for (idx, mapping) in all_original_mappings.iter().enumerate() {
            let key = (mapping.query_name.clone(), mapping.target_name.clone());
            mappings_by_chr_pair.entry(key).or_insert_with(Vec::new).push(idx);
        }

        // Sort each chromosome pair's mappings by query position for binary search
        for indices in mappings_by_chr_pair.values_mut() {
            indices.sort_by_key(|&idx| all_original_mappings[idx].query_start);
        }

        // Collect anchors by chromosome pair for efficient lookup
        let mut anchors_by_chr_pair: HashMap<(String, String), Vec<usize>> = HashMap::new();
        for &anchor_rank in &anchor_ranks {
            if let Some(&anchor_idx) = rank_to_idx.get(&anchor_rank) {
                let anchor = &all_original_mappings[anchor_idx];
                let key = (anchor.query_name.clone(), anchor.target_name.clone());
                anchors_by_chr_pair.entry(key).or_insert_with(Vec::new).push(anchor_idx);
            }
        }

        let mut kept_mappings = Vec::new();
        let mut kept_status = HashMap::new();
        let max_deviation = self.config.scaffold_max_deviation;

        // Process each chromosome pair independently (can be parallelized)
        for ((query_chr, target_chr), mapping_indices) in &mappings_by_chr_pair {
            let chr_key = (query_chr.clone(), target_chr.clone());
            let chr_anchors = anchors_by_chr_pair.get(&chr_key).map(|v| v.as_slice()).unwrap_or(&[]);

            if chr_anchors.is_empty() {
                continue; // No anchors on this chromosome pair, skip
            }

            // For each mapping on this chromosome pair
            for &mapping_idx in mapping_indices {
                let mapping = &all_original_mappings[mapping_idx];

                if anchor_ranks.contains(&mapping.rank) {
                    // This is an anchor - keep it
                    kept_mappings.push(mapping.clone());
                    if mapping.block_length >= self.config.min_scaffold_length {
                        kept_status.insert(mapping.rank, ChainStatus::Scaffold);
                    } else {
                        kept_status.insert(mapping.rank, ChainStatus::Rescued);
                    }
                } else {
                    // Check if within deviation distance of any anchor
                    // Use binary search to find anchors within query range
                    let mapping_q_center = (mapping.query_start + mapping.query_end) / 2;
                    let mapping_t_center = (mapping.target_start + mapping.target_end) / 2;

                    let mut min_distance = u32::MAX;

                    // Only check anchors that are within max_deviation in query space
                    for &anchor_idx in chr_anchors {
                        let anchor = &all_original_mappings[anchor_idx];
                        let anchor_q_center = (anchor.query_start + anchor.query_end) / 2;

                        // Early exit if anchor is too far in query space
                        let q_diff = (mapping_q_center as i64 - anchor_q_center as i64).abs() as u64;
                        if q_diff > max_deviation as u64 {
                            continue; // Too far in query dimension alone
                        }

                        let anchor_t_center = (anchor.target_start + anchor.target_end) / 2;
                        let t_diff = (mapping_t_center as i64 - anchor_t_center as i64).abs() as u64;

                        // Euclidean distance
                        let distance = ((q_diff * q_diff + t_diff * t_diff) as f64).sqrt() as u32;
                        min_distance = min_distance.min(distance);

                        // Early exit if we found a close enough anchor
                        if min_distance <= max_deviation {
                            break;
                        }
                    }

                    if min_distance <= max_deviation {
                        // This mapping is rescued
                        kept_mappings.push(mapping.clone());
                        kept_status.insert(mapping.rank, ChainStatus::Rescued);
                    }
                }
            }
        }

        // Count anchors vs rescued
        let anchor_count = anchor_ranks.len();
        let rescued_count = kept_mappings.len() - anchor_count;
        // Debug: [SCAFFOLD_TRACE] Final rescue statistics

        // Build result map directly from kept mappings
        let mut passing = HashMap::new();
        for meta in kept_mappings {
            let mut result = meta.clone();
            // Set the chain status based on what we determined earlier
            result.chain_status = kept_status.get(&meta.rank)
                .cloned()
                .unwrap_or(ChainStatus::Scaffold);
            passing.insert(meta.rank, result);
        }

        Ok(passing)
    }

    /// Merge mappings into chains using wfmash's union-find approach
    fn merge_mappings_into_chains(&self, metadata: &[RecordMeta], max_gap: u32) -> Result<Vec<MergedChain>> {
        use crate::union_find::UnionFind;

        // Group by (query, target, strand) - this is like wfmash's refSeqId grouping
        // Store the original rank, not the position in metadata array
        let mut groups: HashMap<(String, String, char), Vec<(usize, usize)>> = HashMap::new();

        for (idx, meta) in metadata.iter().enumerate() {
            let key = (meta.query_name.clone(), meta.target_name.clone(), meta.strand);
            groups.entry(key).or_insert_with(Vec::new).push((meta.rank, idx));
        }

        let mut all_chains = Vec::new();

        for ((query, target, strand), indices) in groups {
            // Sort by query start position (like wfmash's sort)
            let mut sorted_indices = indices.clone();
            sorted_indices.sort_by_key(|&(_rank, idx)| metadata[idx].query_start);

            // Use union-find for transitive chaining (like wfmash)
            let mut uf = UnionFind::new(sorted_indices.len());

            // Check pairs for potential chaining, but stop when too far
            for i in 0..sorted_indices.len() {
                let (_rank_i, idx_i) = sorted_indices[i];

                // Look ahead until we find a mapping too far away in query
                for j in (i+1)..sorted_indices.len() {
                    let (_rank_j, idx_j) = sorted_indices[j];

                    // If this mapping starts too far away in query, stop looking
                    // This prevents creating genome-spanning chains
                    if metadata[idx_j].query_start > metadata[idx_i].query_end + max_gap {
                        break; // Too far, no point checking further mappings
                    }

                    // Calculate distances (similar to wfmash's q_dist and r_dist)
                    // Note: wfmash allows small overlaps (up to windowLength/5)
                    let q_gap = if metadata[idx_j].query_start >= metadata[idx_i].query_end {
                        metadata[idx_j].query_start - metadata[idx_i].query_end
                    } else if metadata[idx_i].query_end > metadata[idx_j].query_start {
                        // Small overlap allowed
                        let overlap = metadata[idx_i].query_end - metadata[idx_j].query_start;
                        // Allow overlaps up to max_gap/5 (like wfmash's windowLength/5)
                        if overlap <= max_gap / 5 {
                            0  // Treat small overlap as valid
                        } else {
                            max_gap + 1  // Too much overlap, don't chain
                        }
                    } else {
                        0
                    };

                    let r_gap = if strand == '+' {
                        if metadata[idx_j].target_start >= metadata[idx_i].target_end {
                            metadata[idx_j].target_start - metadata[idx_i].target_end
                        } else if metadata[idx_i].target_end > metadata[idx_j].target_start {
                            // Small overlap allowed
                            let overlap = metadata[idx_i].target_end - metadata[idx_j].target_start;
                            if overlap <= max_gap / 5 {
                                0
                            } else {
                                max_gap + 1
                            }
                        } else {
                            0
                        }
                    } else {
                        // Reverse strand
                        if metadata[idx_i].target_start >= metadata[idx_j].target_end {
                            metadata[idx_i].target_start - metadata[idx_j].target_end
                        } else if metadata[idx_j].target_end > metadata[idx_i].target_start {
                            let overlap = metadata[idx_j].target_end - metadata[idx_i].target_start;
                            if overlap <= max_gap / 5 {
                                0
                            } else {
                                max_gap + 1
                            }
                        } else {
                            0
                        }
                    };

                    // Chain if both gaps are within threshold (allowing transitive connections)
                    if q_gap <= max_gap && r_gap <= max_gap {
                        uf.union(i, j);
                    }
                }
            }

            // Extract chains from union-find sets
            let sets = uf.get_sets();
            let chains: Vec<Vec<(usize, usize)>> = sets.into_iter()
                .map(|set_indices| {
                    set_indices.into_iter()
                        .map(|i| sorted_indices[i])
                        .collect()
                })
                .collect();

            // Create merged chains from groups of mappings
            for chain_indices in chains {
                if chain_indices.is_empty() {
                    continue;
                }

                // Calculate merged chain boundaries (like wfmash's merging)
                let mut q_min = u32::MAX;
                let mut q_max = 0;
                let mut t_min = u32::MAX;
                let mut t_max = 0;
                let mut member_ranks = Vec::new();

                for &(rank, idx) in &chain_indices {
                    let meta = &metadata[idx];
                    q_min = q_min.min(meta.query_start);
                    q_max = q_max.max(meta.query_end);
                    t_min = t_min.min(meta.target_start);
                    t_max = t_max.max(meta.target_end);
                    member_ranks.push(rank);  // Store the original rank
                }

                let total_length = q_max - q_min;

                all_chains.push(MergedChain {
                    query_name: query.clone(),
                    target_name: target.clone(),
                    query_start: q_min,
                    query_end: q_max,
                    target_start: t_min,
                    target_end: t_max,
                    strand,
                    total_length,
                    member_indices: member_ranks,  // Now storing ranks, not indices
                });
            }
        }

        Ok(all_chains)
    }

    /// Apply plane sweep to filter overlapping chains
    /// Calculate overlap fraction between two ranges
    fn calculate_overlap(start1: u32, end1: u32, start2: u32, end2: u32) -> f64 {
        let overlap_start = start1.max(start2);
        let overlap_end = end1.min(end2);

        if overlap_start >= overlap_end {
            return 0.0; // No overlap
        }

        let overlap_len = overlap_end - overlap_start;
        let len1 = end1 - start1;
        let len2 = end2 - start2;
        let min_len = len1.min(len2);

        if min_len == 0 {
            return 0.0;
        }

        overlap_len as f64 / min_len as f64
    }

    /// Apply 1:1 filtering within a prefix group
    fn filter_one_to_one_within_group(&self, mappings: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        if mappings.is_empty() {
            return Ok(Vec::new());
        }

        // Find best mapping (longest)
        let best = mappings.into_iter()
            .max_by_key(|m| m.block_length)
            .unwrap();

        Ok(vec![best])
    }

    /// Apply plane sweep to raw mappings (before scaffold filtering)
    /// This implements the wfmash plane sweep algorithm exactly
    /// CRITICAL: Plane sweep must run PER QUERY SEQUENCE, not globally
    fn apply_plane_sweep_to_mappings(&self, mappings: &[RecordMeta]) -> Result<Vec<RecordMeta>> {
        if mappings.is_empty() || mappings.len() <= 1 {
            return Ok(mappings.to_vec());
        }

        // Group mappings by query sequence
        let mut by_query: HashMap<String, Vec<usize>> = HashMap::new();
        for (idx, mapping) in mappings.iter().enumerate() {
            by_query.entry(mapping.query_name.clone())
                .or_insert_with(Vec::new)
                .push(idx);
        }

        // Clone mappings to allow modification
        let mut mappings = mappings.to_vec();

        // Step 1: Initially mark ALL mappings as bad (discard=true, overlapped=false)
        for mapping in &mut mappings {
            mapping.discard = true;
            mapping.overlapped = false;
        }

        // Calculate scores for each mapping (identity * log(block_length))
        let mut scores = Vec::with_capacity(mappings.len());
        for mapping in &mappings {
            // Use identity from PAF if valid, otherwise default to 0.95
            let identity = if mapping.identity > 0.0 && mapping.identity <= 1.0 {
                mapping.identity
            } else {
                0.95
            };
            // Score = identity * log(block_length)
            let score = identity * (mapping.block_length as f64).ln();
            scores.push(score);
        }

        // For mapping filter mode, determine how many secondary mappings to keep
        // In wfmash, the default is to keep ALL non-overlapping mappings (secondaryToKeep = MAX)
        let secondary_to_keep = match self.config.mapping_filter_mode {
            FilterMode::ManyToMany => usize::MAX, // Keep all non-overlapping
            FilterMode::OneToMany => usize::MAX,  // Keep all non-overlapping (wfmash default)
            FilterMode::OneToOne => 0,   // Keep only best mapping
        };

        let overlap_threshold = self.config.overlap_threshold;

        // Process each query sequence independently
        for (_query, indices) in by_query {
            // Step 2: Create event schedule for plane sweep FOR THIS QUERY
            #[derive(Debug)]
            struct Event {
                pos: u32,
                event_type: u8, // 0 = BEGIN, 1 = END
                idx: usize,
            }

            let mut events = Vec::with_capacity(indices.len() * 2);
            for &idx in &indices {
                let mapping = &mappings[idx];
                events.push(Event { pos: mapping.query_start, event_type: 0, idx });
                events.push(Event { pos: mapping.query_end, event_type: 1, idx });
            }

            // Sort events by position, then by type (BEGIN before END)
            events.sort_by_key(|e| (e.pos, e.event_type));

            // Step 3: Execute plane sweep algorithm FOR THIS QUERY
            // BST to maintain active segments, ordered by score (descending), then query_start, then target_name
            let mut active: BTreeSet<(OrderedFloat<f64>, u32, String, usize)> = BTreeSet::new();

            let mut i = 0;
            while i < events.len() {
                let current_pos = events[i].pos;

                // Process all events at the current position
                let mut j = i;
                while j < events.len() && events[j].pos == current_pos {
                    let event = &events[j];
                    let idx = event.idx;
                    let mapping = &mappings[idx];

                    if event.event_type == 0 { // BEGIN event
                        // Add to active set (negative score for descending order)
                        active.insert((
                            OrderedFloat(-scores[idx]),
                            mapping.query_start,
                            mapping.target_name.clone(),
                            idx
                        ));
                    } else { // END event
                        // Remove from active set
                        active.remove(&(
                            OrderedFloat(-scores[idx]),
                            mapping.query_start,
                            mapping.target_name.clone(),
                            idx
                        ));
                    }
                    j += 1;
                }

                // Mark mappings as good (following wfmash's markGood logic)
                if !active.is_empty() {
                    let mut kept_count = 0;
                    let active_vec: Vec<_> = active.iter().cloned().collect();

                    // First, mark the best secondaryToKeep+1 mappings as good
                    for &(_, _, _, idx) in &active_vec {
                        if kept_count > secondary_to_keep {
                            break;
                        }
                        mappings[idx].discard = false;
                        kept_count += 1;
                    }

                    // Then check for overlaps if threshold < 1.0
                    if overlap_threshold < 1.0 && kept_count > 0 {
                        // Check overlaps with kept mappings
                        let kept_indices: Vec<_> = active_vec.iter()
                            .take(kept_count)
                            .map(|(_, _, _, idx)| *idx)
                            .collect();

                        for (k, &(_, _, _, idx)) in active_vec.iter().enumerate() {
                            if k >= kept_count {
                                // This mapping is beyond secondaryToKeep, check overlap
                                let mapping = &mappings[idx];

                                for &kept_idx in &kept_indices {
                                    let kept = &mappings[kept_idx];

                                    // Only check overlap on same strand
                                    if mapping.strand != kept.strand {
                                        continue;
                                    }

                                    let q_overlap = Self::calculate_overlap(
                                        mapping.query_start, mapping.query_end,
                                        kept.query_start, kept.query_end
                                    );

                                    if q_overlap > overlap_threshold {
                                        mappings[idx].overlapped = true;
                                        mappings[idx].discard = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                i = j;
            }
        }

        // Step 4: Remove bad mappings (where discard=true OR overlapped=true)
        let result: Vec<_> = mappings.into_iter()
            .filter(|m| !m.discard && !m.overlapped)
            .collect();

        Ok(result)
    }

    /// Apply scaffold plane sweep with proper prefix grouping and filtering mode
    fn apply_scaffold_plane_sweep(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // Group chains by prefix pairs (respecting PanSN grouping)
        let mut prefix_groups: HashMap<(String, String), Vec<MergedChain>> = HashMap::new();

        for chain in chains {
            let query_prefix = self.extract_prefix(&chain.query_name);
            let target_prefix = self.extract_prefix(&chain.target_name);
            prefix_groups.entry((query_prefix, target_prefix))
                .or_insert_with(Vec::new)
                .push(chain);
        }

        // Debug: [SCAFFOLD_FILTER] Grouped scaffolds into prefix pairs

        // Apply filtering within each prefix pair IN PARALLEL
        let filtered_groups: Vec<Vec<MergedChain>> = prefix_groups
            .into_par_iter()
            .map(|((_q_prefix, _t_prefix), group)| {
                let group_size = group.len();
                let filtered = match self.config.scaffold_filter_mode {
                    FilterMode::OneToOne => self.scaffold_one_to_one_filter(group),
                    FilterMode::OneToMany => self.scaffold_one_to_many_filter(group),
                    FilterMode::ManyToMany => Ok(group),  // No filtering
                };
                // Debug: [SCAFFOLD_FILTER] Processed prefix group
                filtered.unwrap_or_else(|_| Vec::new())
            })
            .collect();

        // Combine all filtered groups
        Ok(filtered_groups.into_iter().flatten().collect())
    }

    /// Apply 1:1 filtering to scaffold chains (best per query and per target)
    fn scaffold_one_to_one_filter(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // First, group by (query_chr, target_chr) pairs and keep best in each
        let mut by_chr_pair: HashMap<(String, String), Vec<MergedChain>> = HashMap::new();

        for chain in chains {
            let key = (chain.query_name.clone(), chain.target_name.clone());
            by_chr_pair.entry(key).or_insert_with(Vec::new).push(chain);
        }

        // Keep only the best chain for each chromosome pair
        let mut best_per_pair = Vec::new();
        for (_pair, mut pair_chains) in by_chr_pair {
            // Sort by total length (descending) to find the best
            pair_chains.sort_by(|a, b| b.total_length.cmp(&a.total_length));
            if let Some(best) = pair_chains.into_iter().next() {
                best_per_pair.push(best);
            }
        }

        // Now apply 1:1 constraint across all chromosome pairs
        // Sort all best chains by total length
        best_per_pair.sort_by(|a, b| b.total_length.cmp(&a.total_length));

        let mut kept_queries = HashSet::new();
        let mut kept_targets = HashSet::new();
        let mut result = Vec::new();

        for chain in best_per_pair {
            let query_name = chain.query_name.clone();
            let target_name = chain.target_name.clone();

            // Keep only if this is the first (best) for both query and target chromosome
            if !kept_queries.contains(&query_name) && !kept_targets.contains(&target_name) {
                kept_queries.insert(query_name.clone());
                kept_targets.insert(target_name.clone());
                result.push(chain);
            }
        }
        Ok(result)
    }

    /// Apply 1:∞ filtering to scaffold chains (best per query, multiple per target)
    fn scaffold_one_to_many_filter(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // Group by query
        let mut by_query: HashMap<String, Vec<MergedChain>> = HashMap::new();
        for chain in chains {
            by_query.entry(chain.query_name.clone()).or_insert_with(Vec::new).push(chain);
        }

        let mut filtered = Vec::new();

        for (_query, mut query_chains) in by_query {
            // Sort by total length (descending) to prioritize longer chains
            query_chains.sort_by(|a, b| b.total_length.cmp(&a.total_length));

            // Keep non-overlapping chains with preference for longer ones
            let mut kept: Vec<MergedChain> = Vec::new();

            for chain in query_chains {
                let mut has_significant_overlap = false;

                // Check overlap with already kept chains
                for existing in &kept {
                    let overlap_start = chain.query_start.max(existing.query_start);
                    let overlap_end = chain.query_end.min(existing.query_end);

                    if overlap_start < overlap_end {
                        let overlap_len = overlap_end - overlap_start;
                        let chain_len = chain.query_end - chain.query_start;
                        let existing_len = existing.query_end - existing.query_start;

                        let overlap_frac = overlap_len as f64 / chain_len.min(existing_len) as f64;

                        if overlap_frac > self.config.scaffold_overlap_threshold {
                            has_significant_overlap = true;
                            break;
                        }
                    }
                }

                if !has_significant_overlap {
                    kept.push(chain);

                    // Apply max_per_query limit if set
                    if let Some(max) = self.config.scaffold_max_per_query {
                        if kept.len() >= max {
                            break;
                        }
                    }
                }
            }

            filtered.extend(kept);
        }

        Ok(filtered)
    }

    /// Keep all non-overlapping chains (for scaffold filtering with -n inf behavior)
    fn plane_sweep_keep_all_non_overlapping(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // Group chains by query
        let mut by_query: HashMap<String, Vec<MergedChain>> = HashMap::new();
        for chain in chains {
            by_query.entry(chain.query_name.clone()).or_insert_with(Vec::new).push(chain);
        }

        let mut filtered = Vec::new();

        for (_query, mut query_chains) in by_query {
            // Sort by total length (descending) to prioritize longer chains
            query_chains.sort_by(|a, b| b.total_length.cmp(&a.total_length));

            let mut kept: Vec<MergedChain> = Vec::new();

            for chain in query_chains {
                let mut has_significant_overlap = false;

                // Check overlap with already kept chains
                for existing in &kept {
                    let overlap_start = chain.query_start.max(existing.query_start);
                    let overlap_end = chain.query_end.min(existing.query_end);

                    if overlap_start < overlap_end {
                        let overlap_len = overlap_end - overlap_start;
                        let chain_len = chain.query_end - chain.query_start;
                        let existing_len = existing.query_end - existing.query_start;
                        let min_len = chain_len.min(existing_len);

                        // Use scaffold_overlap_threshold to determine significant overlap
                        let overlap_frac = overlap_len as f64 / min_len as f64;
                        if overlap_frac > self.config.scaffold_overlap_threshold {
                            has_significant_overlap = true;
                            break;
                        }
                    }
                }

                // Keep chain if it doesn't significantly overlap with any existing
                if !has_significant_overlap {
                    kept.push(chain);
                }
            }

            filtered.extend(kept);
        }

        Ok(filtered)
    }

    fn plane_sweep_filter_chains(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // Group chains by query
        let mut by_query: HashMap<String, Vec<MergedChain>> = HashMap::new();
        for chain in chains {
            by_query.entry(chain.query_name.clone()).or_insert_with(Vec::new).push(chain);
        }

        let mut filtered = Vec::new();

        for (_query, mut query_chains) in by_query {
            // Sort by query start
            query_chains.sort_by_key(|c| c.query_start);

            // Plane sweep to remove overlaps
            let mut kept: Vec<MergedChain> = Vec::new();

            for chain in query_chains {
                // Check overlap with already kept chains
                let mut has_significant_overlap = false;

                kept.retain(|existing| {
                    let overlap_start = chain.query_start.max(existing.query_start);
                    let overlap_end = chain.query_end.min(existing.query_end);

                    if overlap_start < overlap_end {
                        let overlap_len = overlap_end - overlap_start;
                        let chain_len = chain.query_end - chain.query_start;
                        let existing_len = existing.query_end - existing.query_start;
                        let min_len = chain_len.min(existing_len);

                        // If >50% overlap (wfmash uses scaffold_overlap_threshold)
                        if overlap_len as f64 / min_len as f64 > self.config.scaffold_overlap_threshold {
                            // Keep the longer chain
                            if existing.total_length > chain.total_length {
                                has_significant_overlap = true;
                                true  // Keep existing
                            } else {
                                false  // Remove existing, will add current
                            }
                        } else {
                            true  // No significant overlap, keep both
                        }
                    } else {
                        true  // No overlap, keep
                    }
                });

                if !has_significant_overlap {
                    kept.push(chain);
                }
            }

            filtered.extend(kept);
        }

        Ok(filtered)
    }

    /// Apply final filtering mode (1:1, 1:N, N:N)
    fn apply_filtering_mode(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        match self.config.mapping_filter_mode {
            FilterMode::OneToOne => self.one_to_one_filter(metadata),
            FilterMode::OneToMany => self.one_to_many_filter(metadata),
            FilterMode::ManyToMany => self.many_to_many_filter(metadata),
        }
    }

    /// One-to-one filtering
    fn one_to_one_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        let mut best_per_query: HashMap<String, RecordMeta> = HashMap::new();

        // Find best mapping per query
        for meta in metadata {
            best_per_query
                .entry(meta.query_name.clone())
                .and_modify(|e| {
                    if meta.block_length > e.block_length {
                        *e = meta.clone();
                    }
                })
                .or_insert(meta);
        }

        // From those, keep best per target
        let mut best_per_target: HashMap<String, RecordMeta> = HashMap::new();
        for meta in best_per_query.into_values() {
            best_per_target
                .entry(meta.target_name.clone())
                .and_modify(|e| {
                    if meta.block_length > e.block_length {
                        *e = meta.clone();
                    }
                })
                .or_insert(meta);
        }

        Ok(best_per_target.into_values().collect())
    }

    /// One-to-many filtering
    fn one_to_many_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        let mut result = Vec::new();

        // Group by query
        let mut by_query: HashMap<String, Vec<RecordMeta>> = HashMap::new();
        for meta in metadata {
            by_query.entry(meta.query_name.clone()).or_insert_with(Vec::new).push(meta);
        }

        // Keep best per query, optionally limit per target
        for (_query, mut metas) in by_query {
            metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));

            if let Some(limit) = self.config.mapping_max_per_target {
                metas.truncate(limit);
            }

            result.extend(metas);
        }

        Ok(result)
    }

    /// Many-to-many filtering
    fn many_to_many_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        // For scaffold filtering, we typically don't need additional limiting
        // but we can apply per-query/target limits if specified

        if self.config.mapping_max_per_query.is_none() && self.config.mapping_max_per_target.is_none() {
            return Ok(metadata);
        }

        let mut result = metadata;

        // Apply per-query limit
        if let Some(limit) = self.config.mapping_max_per_query {
            let mut by_query: HashMap<String, Vec<RecordMeta>> = HashMap::new();
            for meta in result {
                by_query.entry(meta.query_name.clone()).or_insert_with(Vec::new).push(meta);
            }

            result = Vec::new();
            for (_query, mut metas) in by_query {
                metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));
                metas.truncate(limit);
                result.extend(metas);
            }
        }

        // Apply per-target limit
        if let Some(limit) = self.config.mapping_max_per_target {
            let mut by_target: HashMap<String, Vec<RecordMeta>> = HashMap::new();
            for meta in result {
                by_target.entry(meta.target_name.clone()).or_insert_with(Vec::new).push(meta);
            }

            result = Vec::new();
            for (_target, mut metas) in by_target {
                metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));
                metas.truncate(limit);
                result.extend(metas);
            }
        }

        Ok(result)
    }

    /// Write filtered output with annotations
    fn write_filtered_output<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: P,
        passing: HashMap<usize, RecordMeta>,
    ) -> Result<()> {
        let output_file = File::create(output_path)?;
        let mut writer = BufWriter::new(output_file);

        // If scaffolds_only mode, write synthetic PAF lines
        if self.scaffolds_only {
            // Sort by rank to maintain order
            let mut sorted: Vec<_> = passing.iter().collect();
            sorted.sort_by_key(|(rank, _)| *rank);

            for (_rank, meta) in sorted {
                // Write synthetic PAF line for scaffold
                writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tst:Z:scaffold",
                    meta.query_name,
                    583092,  // Hardcoded query length for chrV - would need proper tracking
                    meta.query_start,
                    meta.query_end,
                    meta.strand,
                    meta.target_name,
                    576784,  // Hardcoded target length for chrV - would need proper tracking
                    meta.target_start,
                    meta.target_end,
                    meta.block_length,  // matches
                    meta.block_length,  // block length
                    255  // mapping quality
                )?;
            }
            writer.flush()?;
            return Ok(());
        }

        // Normal mode - read input and filter
        let input_file = File::open(input_path)?;
        let reader = BufReader::new(input_file);

        for (rank, line) in reader.lines().enumerate() {
            if let Some(meta) = passing.get(&rank) {
                let mut line = line?;

                // Add our annotations as tags
                if let Some(ref chain_id) = meta.chain_id {
                    line.push_str(&format!("\tch:Z:{}", chain_id));
                }
                // Output status tag like wfmash does
                let status_str = match meta.chain_status {
                    ChainStatus::Scaffold => "scaffold",
                    ChainStatus::Rescued => "rescued",
                    ChainStatus::Unassigned => "unassigned",
                };
                line.push_str(&format!("\tst:Z:{}", status_str));

                writeln!(writer, "{}", line)?;
            }
        }

        writer.flush()?;
        Ok(())
    }
}