//! Batch alignment for large genome sets
//!
//! When aligning many genomes, the k-mer index can exceed available disk space.
//! This module partitions genomes into batches based on estimated index size,
//! runs FastGA on each batch pair, and aggregates results.
//!
//! Index size formula (empirically derived):
//!   GIX_size ≈ 0.1 GB + 12 bytes/bp
//!
//! For a 100 Mbp batch → ~1.3 GB index
//! For a 1 Gbp batch → ~12.1 GB index

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Estimated index overhead in bytes (fixed cost per index)
const INDEX_OVERHEAD_BYTES: u64 = 100_000_000; // 100 MB

/// Estimated bytes of index per basepair of genome
const INDEX_BYTES_PER_BP: f64 = 12.0;

/// Information about a genome (from PanSN prefix)
#[derive(Debug, Clone)]
pub struct GenomeInfo {
    /// PanSN prefix (e.g., "SGDref#1#")
    pub prefix: String,
    /// Total basepairs in this genome
    pub total_bp: u64,
    /// Source FASTA file
    pub source_file: PathBuf,
}

/// Information about a single sequence (chromosome/contig)
#[derive(Debug, Clone)]
pub struct SequenceInfo {
    /// Full sequence name (e.g., "SGDref#1#chrI")
    pub name: String,
    /// Basepairs in this sequence
    pub length: u64,
    /// Source FASTA file
    pub source_file: PathBuf,
}

/// A batch of sequences to align together (for sequence-level splitting)
#[derive(Debug, Clone)]
pub struct SequenceBatch {
    /// Sequences in this batch
    pub sequences: Vec<SequenceInfo>,
    /// Total basepairs in batch
    pub total_bp: u64,
    /// Estimated index size in bytes
    pub estimated_index_bytes: u64,
}

impl Default for SequenceBatch {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceBatch {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            total_bp: 0,
            estimated_index_bytes: INDEX_OVERHEAD_BYTES,
        }
    }

    /// Add a sequence to this batch
    pub fn add(&mut self, seq: SequenceInfo) {
        self.total_bp += seq.length;
        self.estimated_index_bytes = estimate_index_size(self.total_bp);
        self.sequences.push(seq);
    }

    /// Check if adding a sequence would exceed the byte limit
    pub fn would_exceed(&self, seq: &SequenceInfo, max_bytes: u64) -> bool {
        let new_total = self.total_bp + seq.length;
        estimate_index_size(new_total) > max_bytes
    }
}

/// A batch of genomes to align together
#[derive(Debug, Clone)]
pub struct GenomeBatch {
    /// Genomes in this batch
    pub genomes: Vec<GenomeInfo>,
    /// Total basepairs in batch
    pub total_bp: u64,
    /// Estimated index size in bytes
    pub estimated_index_bytes: u64,
}

impl Default for GenomeBatch {
    fn default() -> Self {
        Self::new()
    }
}

impl GenomeBatch {
    pub fn new() -> Self {
        Self {
            genomes: Vec::new(),
            total_bp: 0,
            estimated_index_bytes: INDEX_OVERHEAD_BYTES,
        }
    }

    /// Add a genome to this batch
    pub fn add(&mut self, genome: GenomeInfo) {
        self.total_bp += genome.total_bp;
        self.estimated_index_bytes = estimate_index_size(self.total_bp);
        self.genomes.push(genome);
    }

    /// Check if adding a genome would exceed the byte limit
    pub fn would_exceed(&self, genome: &GenomeInfo, max_bytes: u64) -> bool {
        let new_total = self.total_bp + genome.total_bp;
        estimate_index_size(new_total) > max_bytes
    }
}

/// Estimate index size in bytes for a given genome size
pub fn estimate_index_size(genome_bp: u64) -> u64 {
    INDEX_OVERHEAD_BYTES + (genome_bp as f64 * INDEX_BYTES_PER_BP) as u64
}

/// Format bytes as human-readable string
pub fn format_bytes(bytes: u64) -> String {
    if bytes >= 1_000_000_000 {
        format!("{:.1} GB", bytes as f64 / 1_000_000_000.0)
    } else if bytes >= 1_000_000 {
        format!("{:.1} MB", bytes as f64 / 1_000_000.0)
    } else if bytes >= 1_000 {
        format!("{:.1} KB", bytes as f64 / 1_000.0)
    } else {
        format!("{} bytes", bytes)
    }
}

/// Parse a single FASTA file and return genome info
fn parse_single_fasta_genomes(fasta_file: &str) -> Result<Vec<(String, u64, PathBuf)>> {
    let path = Path::new(fasta_file);
    let file =
        File::open(path).with_context(|| format!("Failed to open FASTA: {}", fasta_file))?;

    let reader: Box<dyn BufRead> =
        if fasta_file.ends_with(".gz") || fasta_file.ends_with(".bgz") {
            Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

    let mut results: Vec<(String, u64, PathBuf)> = Vec::new();
    let mut current_prefix: Option<String> = None;
    let mut current_bp: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with('>') {
            // Save previous genome's data
            if let Some(prefix) = current_prefix.take() {
                results.push((prefix, current_bp, path.to_path_buf()));
            }

            // Extract PanSN prefix (genome#haplotype#)
            let name = trimmed
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("");
            let prefix = extract_pansn_prefix(name);
            current_prefix = Some(prefix);
            current_bp = 0;
        } else if !trimmed.is_empty() {
            current_bp += trimmed.len() as u64;
        }
    }

    // Don't forget the last sequence
    if let Some(prefix) = current_prefix {
        results.push((prefix, current_bp, path.to_path_buf()));
    }

    Ok(results)
}

/// Parse genome sizes from FASTA file(s) in parallel
/// Returns map of genome prefix -> GenomeInfo
pub fn parse_genome_sizes(fasta_files: &[String]) -> Result<Vec<GenomeInfo>> {
    // Parse files in parallel
    let file_results: Vec<Result<Vec<(String, u64, PathBuf)>>> = fasta_files
        .par_iter()
        .map(|f| parse_single_fasta_genomes(f))
        .collect();

    // Merge results
    let mut genomes: HashMap<String, GenomeInfo> = HashMap::new();
    for result in file_results {
        for (prefix, bp, source_file) in result? {
            let entry = genomes.entry(prefix.clone()).or_insert_with(|| GenomeInfo {
                prefix: prefix.clone(),
                total_bp: 0,
                source_file,
            });
            entry.total_bp += bp;
        }
    }

    // Convert to sorted vector for deterministic ordering
    let mut genome_list: Vec<GenomeInfo> = genomes.into_values().collect();
    genome_list.sort_by(|a, b| a.prefix.cmp(&b.prefix));

    Ok(genome_list)
}

/// Extract PanSN genome prefix from sequence name
/// E.g., "SGDref#1#chrI" -> "SGDref#1#"
fn extract_pansn_prefix(name: &str) -> String {
    let parts: Vec<&str> = name.split('#').collect();
    if parts.len() >= 2 {
        format!("{}#{}#", parts[0], parts[1])
    } else {
        format!("{}#", name)
    }
}

/// Parse a single FASTA file and return sequence info
fn parse_single_fasta_sequences(fasta_file: &str) -> Result<Vec<SequenceInfo>> {
    let path = Path::new(fasta_file);
    let file =
        File::open(path).with_context(|| format!("Failed to open FASTA: {}", fasta_file))?;

    let reader: Box<dyn BufRead> =
        if fasta_file.ends_with(".gz") || fasta_file.ends_with(".bgz") {
            Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

    let mut sequences: Vec<SequenceInfo> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_bp: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with('>') {
            // Save previous sequence
            if let Some(name) = current_name.take() {
                sequences.push(SequenceInfo {
                    name,
                    length: current_bp,
                    source_file: path.to_path_buf(),
                });
            }

            // Start new sequence
            let name = trimmed
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_name = Some(name);
            current_bp = 0;
        } else if !trimmed.is_empty() {
            current_bp += trimmed.len() as u64;
        }
    }

    // Don't forget the last sequence
    if let Some(name) = current_name {
        sequences.push(SequenceInfo {
            name,
            length: current_bp,
            source_file: path.to_path_buf(),
        });
    }

    Ok(sequences)
}

/// Parse individual sequences from FASTA file(s) in parallel
/// Returns a list of all sequences with their sizes
pub fn parse_sequences(fasta_files: &[String]) -> Result<Vec<SequenceInfo>> {
    // Parse files in parallel
    let file_results: Vec<Result<Vec<SequenceInfo>>> = fasta_files
        .par_iter()
        .map(|f| parse_single_fasta_sequences(f))
        .collect();

    // Merge results
    let mut sequences: Vec<SequenceInfo> = Vec::new();
    for result in file_results {
        sequences.extend(result?);
    }

    // Sort by size (largest first) for better bin packing
    sequences.sort_by(|a, b| b.length.cmp(&a.length));

    Ok(sequences)
}

/// Partition sequences into batches based on max index size
/// Uses best-fit decreasing bin packing algorithm for optimal space utilization
pub fn partition_sequences_into_batches(
    sequences: Vec<SequenceInfo>,
    max_index_bytes: u64,
) -> Vec<SequenceBatch> {
    // Divide by 2 because during cross-batch alignment, both query and target
    // indexes exist simultaneously on disk
    let per_batch_limit = max_index_bytes / 2;

    let mut batches: Vec<SequenceBatch> = Vec::new();

    // sequences are already sorted by size (largest first) from parse_sequences()
    for seq in sequences {
        // Check if this sequence alone exceeds the limit
        let single_seq_size = estimate_index_size(seq.length);
        if single_seq_size > per_batch_limit {
            eprintln!(
                "[batch] WARNING: Sequence {} ({}) has estimated index {} which exceeds per-batch limit {}",
                seq.name,
                format_bytes(seq.length),
                format_bytes(single_seq_size),
                format_bytes(per_batch_limit)
            );
            eprintln!("[batch] Including it as a single-sequence batch anyway");

            // Add as its own batch
            let mut oversized = SequenceBatch::new();
            oversized.add(seq);
            batches.push(oversized);
            continue;
        }

        // Best-fit: find the batch with least remaining space that can still fit this sequence
        let mut best_fit_idx: Option<usize> = None;
        let mut best_fit_remaining = u64::MAX;

        for (idx, batch) in batches.iter().enumerate() {
            if !batch.would_exceed(&seq, per_batch_limit) {
                // Calculate remaining space after adding this sequence
                let new_size = estimate_index_size(batch.total_bp + seq.length);
                let remaining = per_batch_limit.saturating_sub(new_size);

                // Best-fit: prefer the batch with least remaining space (tightest fit)
                if remaining < best_fit_remaining {
                    best_fit_remaining = remaining;
                    best_fit_idx = Some(idx);
                }
            }
        }

        if let Some(idx) = best_fit_idx {
            batches[idx].add(seq);
        } else {
            // No batch has room, create a new one
            let mut new_batch = SequenceBatch::new();
            new_batch.add(seq);
            batches.push(new_batch);
        }
    }

    batches
}

/// Stream-write sequence batches from a single source file using Mutex-protected writers
/// This allows parallel processing of multiple source files
fn stream_write_seq_batches_from_source_parallel(
    source_file: &Path,
    name_to_batches: &HashMap<String, Vec<usize>>,
    batch_writers: &[std::sync::Mutex<std::io::BufWriter<File>>],
) -> Result<()> {
    let file = File::open(source_file)
        .with_context(|| format!("Failed to open source FASTA: {}", source_file.display()))?;

    let source_str = source_file.to_string_lossy();
    let reader: Box<dyn BufRead> = if source_str.ends_with(".gz") || source_str.ends_with(".bgz") {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            noodles::bgzf::io::reader::Reader::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(256 * 1024, file))
    };

    let mut current_batch_indices: Vec<usize> = Vec::new();
    // Buffer to accumulate a complete sequence before writing
    let mut seq_buffer: Vec<u8> = Vec::with_capacity(64 * 1024);

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // Flush previous sequence to its batches
            if !seq_buffer.is_empty() && !current_batch_indices.is_empty() {
                for &batch_idx in &current_batch_indices {
                    let mut writer = batch_writers[batch_idx].lock().unwrap();
                    writer.write_all(&seq_buffer)?;
                }
                seq_buffer.clear();
            }

            // Extract sequence name and find which batches need it
            let name = line
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("");
            current_batch_indices = name_to_batches
                .get(name)
                .cloned()
                .unwrap_or_default();
        }

        // Buffer the line
        if !current_batch_indices.is_empty() {
            seq_buffer.extend_from_slice(line.as_bytes());
            seq_buffer.push(b'\n');
        }
    }

    // Flush final sequence
    if !seq_buffer.is_empty() && !current_batch_indices.is_empty() {
        for &batch_idx in &current_batch_indices {
            let mut writer = batch_writers[batch_idx].lock().unwrap();
            writer.write_all(&seq_buffer)?;
        }
    }

    Ok(())
}

/// Write all sequence batches by processing source files in parallel
/// Each source file is processed by a separate thread, writing to Mutex-protected batch files
pub fn write_all_sequence_batches(
    batches: &[SequenceBatch],
    batch_dir: &Path,
) -> Result<Vec<PathBuf>> {
    use std::sync::Mutex;

    // Build sequence_name -> batch indices mapping
    let mut name_to_batches: HashMap<String, Vec<usize>> = HashMap::new();
    for (batch_idx, batch) in batches.iter().enumerate() {
        for seq in &batch.sequences {
            name_to_batches
                .entry(seq.name.clone())
                .or_default()
                .push(batch_idx);
        }
    }

    // Collect all unique source files
    let mut source_files: Vec<PathBuf> = batches
        .iter()
        .flat_map(|b| b.sequences.iter().map(|s| s.source_file.clone()))
        .collect();
    source_files.sort();
    source_files.dedup();

    // Create batch subdirectories and output files with Mutex protection
    let mut batch_paths: Vec<PathBuf> = Vec::with_capacity(batches.len());
    let mut batch_writers: Vec<Mutex<std::io::BufWriter<File>>> = Vec::with_capacity(batches.len());

    for i in 0..batches.len() {
        let batch_subdir = batch_dir.join(format!("batch_{}", i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("sequences.fa");
        let file = File::create(&batch_path)
            .with_context(|| format!("Failed to create batch FASTA: {}", batch_path.display()))?;
        batch_writers.push(Mutex::new(std::io::BufWriter::with_capacity(256 * 1024, file)));
        batch_paths.push(batch_path);
    }

    // Process source files in parallel
    let results: Vec<Result<()>> = source_files
        .par_iter()
        .map(|source_file| {
            stream_write_seq_batches_from_source_parallel(source_file, &name_to_batches, &batch_writers)
        })
        .collect();

    // Check for errors
    for result in results {
        result?;
    }

    // Flush all writers
    for writer in &batch_writers {
        writer.lock().unwrap().flush()?;
    }

    Ok(batch_paths)
}

/// Report sequence batch statistics
pub fn report_sequence_batch_plan(
    batches: &[SequenceBatch],
    total_bp: u64,
    max_index_bytes: u64,
    quiet: bool,
) {
    if quiet {
        return;
    }

    let total_sequences: usize = batches.iter().map(|b| b.sequences.len()).sum();
    let total_pairs = batches.len() * batches.len();
    let self_pairs = batches.len();
    let cross_pairs = total_pairs - self_pairs;

    eprintln!(
        "[batch] Partitioned {} sequences ({}) into {} batches (sequence-level splitting):",
        total_sequences,
        format_bytes(total_bp),
        batches.len()
    );

    // Find the two largest batches
    let mut batch_sizes: Vec<u64> = batches.iter().map(|b| b.estimated_index_bytes).collect();
    batch_sizes.sort_by(|a, b| b.cmp(a));

    let largest_batch = batch_sizes.first().copied().unwrap_or(0);
    let second_largest = batch_sizes.get(1).copied().unwrap_or(0);
    let peak_disk_estimate = largest_batch + second_largest;

    for (i, batch) in batches.iter().enumerate() {
        eprintln!(
            "[batch]   Batch {}: {} sequences, {}, est. index {}",
            i + 1,
            batch.sequences.len(),
            format_bytes(batch.total_bp),
            format_bytes(batch.estimated_index_bytes)
        );
    }

    eprintln!(
        "[batch] Will run {} batch alignments ({} self + {} cross-batch)",
        total_pairs, self_pairs, cross_pairs
    );

    if peak_disk_estimate > max_index_bytes {
        eprintln!(
            "[batch] WARNING: Peak disk usage (~{}) may exceed --batch-bytes limit ({})",
            format_bytes(peak_disk_estimate),
            format_bytes(max_index_bytes)
        );
    } else {
        eprintln!(
            "[batch] Peak disk usage estimate: {} (within {} limit)",
            format_bytes(peak_disk_estimate),
            format_bytes(max_index_bytes)
        );
    }
}

/// Partition genomes into batches based on max index size
/// Uses best-fit decreasing bin packing algorithm for optimal space utilization
/// Note: max_index_bytes is the TOTAL peak disk limit. Since FastGA creates indexes
/// for both query and target during cross-batch alignment, we use half the limit
/// per batch to ensure peak usage (2 batches) stays within the user's limit.
pub fn partition_into_batches(genomes: Vec<GenomeInfo>, max_index_bytes: u64) -> Vec<GenomeBatch> {
    // Divide by 2 because during cross-batch alignment, both query and target
    // indexes exist simultaneously on disk
    let per_batch_limit = max_index_bytes / 2;

    // Sort genomes by size (largest first) for best-fit decreasing
    let mut sorted_genomes = genomes;
    sorted_genomes.sort_by(|a, b| b.total_bp.cmp(&a.total_bp));

    let mut batches: Vec<GenomeBatch> = Vec::new();

    for genome in sorted_genomes {
        // Check if this genome alone exceeds the limit
        let single_genome_size = estimate_index_size(genome.total_bp);
        if single_genome_size > per_batch_limit {
            eprintln!(
                "[batch] WARNING: Genome {} ({}) has estimated index {} which exceeds per-batch limit {} (total limit {})",
                genome.prefix,
                format_bytes(genome.total_bp),
                format_bytes(single_genome_size),
                format_bytes(per_batch_limit),
                format_bytes(max_index_bytes)
            );
            eprintln!("[batch] Including it as a single-genome batch anyway");

            // Add oversized genome as its own batch
            let mut oversized = GenomeBatch::new();
            oversized.add(genome);
            batches.push(oversized);
            continue;
        }

        // Best-fit: find the batch with least remaining space that can still fit this genome
        let mut best_fit_idx: Option<usize> = None;
        let mut best_fit_remaining = u64::MAX;

        for (idx, batch) in batches.iter().enumerate() {
            if !batch.would_exceed(&genome, per_batch_limit) {
                // Calculate remaining space after adding this genome
                let new_size = estimate_index_size(batch.total_bp + genome.total_bp);
                let remaining = per_batch_limit.saturating_sub(new_size);

                // Best-fit: prefer the batch with least remaining space (tightest fit)
                if remaining < best_fit_remaining {
                    best_fit_remaining = remaining;
                    best_fit_idx = Some(idx);
                }
            }
        }

        if let Some(idx) = best_fit_idx {
            batches[idx].add(genome);
        } else {
            // No batch has room, create a new one
            let mut new_batch = GenomeBatch::new();
            new_batch.add(genome);
            batches.push(new_batch);
        }
    }

    batches
}

/// Read sequences from a single source file that match the given prefixes
/// Stream-write genome batches from a single source file using Mutex-protected writers
/// This allows parallel processing of multiple source files
fn stream_write_batches_from_source_parallel(
    source_file: &Path,
    prefix_to_batches: &HashMap<String, Vec<usize>>,
    batch_writers: &[std::sync::Mutex<std::io::BufWriter<File>>],
) -> Result<()> {
    let file = File::open(source_file)
        .with_context(|| format!("Failed to open source FASTA: {}", source_file.display()))?;

    let source_str = source_file.to_string_lossy();
    let reader: Box<dyn BufRead> = if source_str.ends_with(".gz") || source_str.ends_with(".bgz") {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            noodles::bgzf::io::reader::Reader::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(256 * 1024, file))
    };

    let mut current_batch_indices: Vec<usize> = Vec::new();
    // Buffer to accumulate a complete sequence before writing
    let mut seq_buffer: Vec<u8> = Vec::with_capacity(64 * 1024);

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // Flush previous sequence to its batches
            if !seq_buffer.is_empty() && !current_batch_indices.is_empty() {
                for &batch_idx in &current_batch_indices {
                    let mut writer = batch_writers[batch_idx].lock().unwrap();
                    writer.write_all(&seq_buffer)?;
                }
                seq_buffer.clear();
            }

            // Extract prefix and find which batches need this sequence
            let name = line
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("");
            let prefix = extract_pansn_prefix(name);
            current_batch_indices = prefix_to_batches
                .get(&prefix)
                .cloned()
                .unwrap_or_default();
        }

        // Buffer the line
        if !current_batch_indices.is_empty() {
            seq_buffer.extend_from_slice(line.as_bytes());
            seq_buffer.push(b'\n');
        }
    }

    // Flush final sequence
    if !seq_buffer.is_empty() && !current_batch_indices.is_empty() {
        for &batch_idx in &current_batch_indices {
            let mut writer = batch_writers[batch_idx].lock().unwrap();
            writer.write_all(&seq_buffer)?;
        }
    }

    Ok(())
}

/// Write all genome batches by processing source files in parallel
/// Each source file is processed by a separate thread, writing to Mutex-protected batch files
pub fn write_all_genome_batches(
    batches: &[GenomeBatch],
    batch_dir: &Path,
) -> Result<Vec<PathBuf>> {
    use std::sync::Mutex;

    // Build prefix -> batch indices mapping
    let mut prefix_to_batches: HashMap<String, Vec<usize>> = HashMap::new();
    for (batch_idx, batch) in batches.iter().enumerate() {
        for genome in &batch.genomes {
            prefix_to_batches
                .entry(genome.prefix.clone())
                .or_default()
                .push(batch_idx);
        }
    }

    // Collect all unique source files
    let mut source_files: Vec<PathBuf> = batches
        .iter()
        .flat_map(|b| b.genomes.iter().map(|g| g.source_file.clone()))
        .collect();
    source_files.sort();
    source_files.dedup();

    // Create batch subdirectories and output files with Mutex protection
    let mut batch_paths: Vec<PathBuf> = Vec::with_capacity(batches.len());
    let mut batch_writers: Vec<Mutex<std::io::BufWriter<File>>> = Vec::with_capacity(batches.len());

    for i in 0..batches.len() {
        let batch_subdir = batch_dir.join(format!("batch_{}", i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("genomes.fa");
        let file = File::create(&batch_path)
            .with_context(|| format!("Failed to create batch FASTA: {}", batch_path.display()))?;
        batch_writers.push(Mutex::new(std::io::BufWriter::with_capacity(256 * 1024, file)));
        batch_paths.push(batch_path);
    }

    // Process source files in parallel
    let results: Vec<Result<()>> = source_files
        .par_iter()
        .map(|source_file| {
            stream_write_batches_from_source_parallel(source_file, &prefix_to_batches, &batch_writers)
        })
        .collect();

    // Check for errors
    for result in results {
        result?;
    }

    // Flush all writers
    for writer in &batch_writers {
        writer.lock().unwrap().flush()?;
    }

    Ok(batch_paths)
}
/// Report batch statistics and warn about oversized batches
pub fn report_batch_plan(batches: &[GenomeBatch], total_bp: u64, max_index_bytes: u64, quiet: bool) {
    if quiet {
        return;
    }

    let total_genomes: usize = batches.iter().map(|b| b.genomes.len()).sum();
    let total_pairs = batches.len() * batches.len();
    let self_pairs = batches.len();
    let cross_pairs = total_pairs - self_pairs;

    eprintln!(
        "[batch] Partitioned {} genomes ({}) into {} batches:",
        total_genomes,
        format_bytes(total_bp),
        batches.len()
    );

    // Find the two largest batches (for peak disk estimate during cross-batch alignment)
    let mut batch_sizes: Vec<u64> = batches.iter().map(|b| b.estimated_index_bytes).collect();
    batch_sizes.sort_by(|a, b| b.cmp(a)); // Sort descending

    let largest_batch = batch_sizes.first().copied().unwrap_or(0);
    let second_largest = batch_sizes.get(1).copied().unwrap_or(0);
    let peak_disk_estimate = largest_batch + second_largest;

    for (i, batch) in batches.iter().enumerate() {
        eprintln!(
            "[batch]   Batch {}: {} genomes, {}, est. index {}",
            i + 1,
            batch.genomes.len(),
            format_bytes(batch.total_bp),
            format_bytes(batch.estimated_index_bytes)
        );
    }

    eprintln!(
        "[batch] Will run {} batch alignments ({} self + {} cross-batch)",
        total_pairs, self_pairs, cross_pairs
    );

    // Warn if peak disk usage will exceed the limit
    if peak_disk_estimate > max_index_bytes {
        eprintln!(
            "[batch] WARNING: Peak disk usage (~{}) will exceed --batch-bytes limit ({})",
            format_bytes(peak_disk_estimate),
            format_bytes(max_index_bytes)
        );
        eprintln!(
            "[batch] NOTE: Minimum --batch-bytes needed for these genomes: {}",
            format_bytes(peak_disk_estimate + peak_disk_estimate / 10) // Add 10% margin
        );
        eprintln!(
            "[batch] TIP: Consider using a higher --batch-bytes limit, or split large genomes by chromosome"
        );
    }
}

/// Configuration for batch alignment
pub struct BatchAlignConfig {
    pub frequency: Option<usize>,
    pub threads: usize,
    pub min_alignment_length: u64,
    pub zstd_compress: bool,
    pub zstd_level: u32,
    pub keep_self: bool,
    pub quiet: bool,
}

/// Run batch alignment on a set of FASTA files
/// Returns a temp file containing aggregated PAF output
pub fn run_batch_alignment(
    fasta_files: &[String],
    max_index_bytes: u64,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::{self, FastGAIntegration};

    // partition_into_batches() handles dividing by 2 internally
    let per_batch_limit = max_index_bytes / 2;

    if !config.quiet {
        eprintln!(
            "[batch] Batch mode: max {} total index size ({} per batch)",
            format_bytes(max_index_bytes),
            format_bytes(per_batch_limit)
        );
    }

    // Determine temp directory
    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    // Create temp directory for batch FASTAs
    let batch_dir = temp_base.join(format!("sweepga_batch_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    // Parse genome sizes
    if !config.quiet {
        eprintln!("[batch] Scanning input files for genome sizes...");
    }
    let genomes = parse_genome_sizes(fasta_files)?;
    let total_bp: u64 = genomes.iter().map(|g| g.total_bp).sum();

    if genomes.is_empty() {
        anyhow::bail!("No genomes found in input files");
    }

    // Partition into batches
    let batches = partition_into_batches(genomes.clone(), max_index_bytes);

    // Check if any batch exceeds the per-batch limit - if so, use sequence-level splitting
    let any_oversized = batches
        .iter()
        .any(|b| b.estimated_index_bytes > per_batch_limit);

    if any_oversized {
        if !config.quiet {
            eprintln!(
                "[batch] Genome-level batching produces oversized batches, switching to sequence-level splitting..."
            );
        }
        // Clean up and use sequence-level batching instead
        let _ = std::fs::remove_dir_all(&batch_dir);
        return run_sequence_batch_alignment(fasta_files, max_index_bytes, config, tempdir);
    }

    if batches.len() == 1 {
        if !config.quiet {
            eprintln!(
                "[batch] All genomes ({}) fit in single batch (est. index {}), no batching needed",
                format_bytes(total_bp),
                format_bytes(batches[0].estimated_index_bytes)
            );
        }
        // Fall back to normal single-run alignment
        let _ = std::fs::remove_dir_all(&batch_dir);
        return run_single_batch_alignment(fasta_files, config, tempdir);
    }

    report_batch_plan(&batches, total_bp, max_index_bytes, config.quiet);

    // Write all batch FASTAs in a single pass through source files (streaming)
    if !config.quiet {
        eprintln!(
            "[batch] Writing {} batch FASTA files (streaming, single pass)...",
            batches.len()
        );
    }

    let batch_files = write_all_genome_batches(&batches, &batch_dir)?;

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
    );

    // Phase 1: Create GDBs for ALL batches upfront
    // ALNtoPAF needs GDB files for both query and target to resolve sequence names
    if !config.quiet {
        eprintln!("[batch] Creating GDB files for all batches...");
    }
    let mut gdb_bases: Vec<String> = Vec::new();
    for (i, batch_file) in batch_files.iter().enumerate() {
        if !config.quiet {
            eprintln!("[batch]   GDB for batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        gdb_bases.push(gdb_base);
    }

    // Create merged output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;
    let num_batches = batch_files.len();

    // Phase 2: For each target batch, create index, run alignments, cleanup index
    for i in 0..num_batches {
        // Build index for batch i (target) - GDB already exists
        if !config.quiet {
            eprintln!("[batch] Building index for batch {}...", i + 1);
        }
        fastga.create_index_only(&gdb_bases[i])?;

        // Track disk usage of index files
        let index_dir = std::path::Path::new(&gdb_bases[i])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::add_bytes(index_size);

        // Compress if requested
        if config.zstd_compress {
            fastga_integration::FastGAIntegration::compress_index(
                &gdb_bases[i],
                config.zstd_level,
            )?;
        }

        for j in 0..num_batches {
            // Skip self-alignment if not keeping self mappings and same batch
            if i == j && !config.keep_self && num_batches > 1 {
                // For self-batch, we still want inter-genome alignments within the batch
                // So we DO run self-alignment, filtering happens later
            }

            if !config.quiet {
                if i == j {
                    eprintln!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    eprintln!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            // Run alignment: query batch j against target batch i
            // Use direct PAF output to avoid ALNtoPAF conversion issues
            let paf_bytes = fastga.align_direct_paf(&batch_files[j], &batch_files[i])?;

            // Append to merged output
            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                eprintln!("[batch]   {} alignments", line_count);
            }

            // Clean up query index if FastGA auto-created one (j != i)
            // This keeps disk usage bounded to ~1 batch index at a time
            if j != i {
                let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[j]);
            }
        }

        // Clean up index (GIX only) for batch i to free disk space
        // Keep GDB files - they're needed for other batches as queries
        if !config.quiet {
            eprintln!("[batch] Cleaning up index for batch {}...", i + 1);
        }
        // Untrack disk usage before cleanup
        crate::disk_usage::remove_bytes(index_size);
        let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[i]);
    }

    // Phase 3: Clean up all GDB files and batch directory
    for gdb_base in &gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !config.quiet {
        eprintln!(
            "[batch] Completed batch alignment: {} total alignments",
            total_alignments
        );
    }

    Ok(merged_paf)
}

/// Run alignment without batching (single FastGA run)
/// Uses align_to_temp_paf_in_tempdir to ensure index files are created in tempdir
fn run_single_batch_alignment(
    fasta_files: &[String],
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::FastGAIntegration;

    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
    );

    // Use first file for self-alignment
    let input_path = std::fs::canonicalize(&fasta_files[0])
        .with_context(|| format!("Failed to resolve path: {}", fasta_files[0]))?;

    // Use _in_tempdir variant to ensure index files go to tempdir, not next to input
    fastga.align_to_temp_paf_in_tempdir(&input_path, &input_path)
}

/// Run batch alignment with sequence-level splitting
/// This is used when individual genomes are too large for genome-level batching
pub fn run_sequence_batch_alignment(
    fasta_files: &[String],
    max_index_bytes: u64,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::{self, FastGAIntegration};

    // partition_sequences_into_batches() handles dividing by 2 internally
    let _per_batch_limit = max_index_bytes / 2;

    // Determine temp directory
    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    // Create temp directory for batch FASTAs
    let batch_dir = temp_base.join(format!("sweepga_seqbatch_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    // Parse sequences
    if !config.quiet {
        eprintln!("[batch] Parsing individual sequences for sequence-level batching...");
    }
    let sequences = parse_sequences(fasta_files)?;
    let total_bp: u64 = sequences.iter().map(|s| s.length).sum();

    if sequences.is_empty() {
        anyhow::bail!("No sequences found in input files");
    }

    if !config.quiet {
        eprintln!(
            "[batch] Found {} sequences totaling {}",
            sequences.len(),
            format_bytes(total_bp)
        );
    }

    // Partition sequences into batches
    let batches = partition_sequences_into_batches(sequences, max_index_bytes);

    if batches.len() == 1 {
        if !config.quiet {
            eprintln!(
                "[batch] All sequences fit in single batch (est. index {})",
                format_bytes(batches[0].estimated_index_bytes)
            );
        }
        // Clean up and fall back to normal single-run alignment
        let _ = std::fs::remove_dir_all(&batch_dir);
        return run_single_batch_alignment(fasta_files, config, tempdir);
    }

    report_sequence_batch_plan(&batches, total_bp, max_index_bytes, config.quiet);

    // Write all batch FASTAs in a single pass through source files (streaming)
    if !config.quiet {
        eprintln!(
            "[batch] Writing {} batch FASTA files (streaming, single pass)...",
            batches.len()
        );
    }

    let batch_files = write_all_sequence_batches(&batches, &batch_dir)?;

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
    );

    // Phase 1: Create GDBs for ALL batches upfront
    if !config.quiet {
        eprintln!("[batch] Creating GDB files for all batches...");
    }
    let mut gdb_bases: Vec<String> = Vec::new();
    for (i, batch_file) in batch_files.iter().enumerate() {
        if !config.quiet {
            eprintln!("[batch]   GDB for batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        gdb_bases.push(gdb_base);
    }

    // Create merged output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;
    let num_batches = batch_files.len();

    // Phase 2: For each target batch, create index, run alignments, cleanup index
    for i in 0..num_batches {
        // Build index for batch i (target)
        if !config.quiet {
            eprintln!("[batch] Building index for batch {}...", i + 1);
        }
        fastga.create_index_only(&gdb_bases[i])?;

        // Track disk usage of index files
        let index_dir = std::path::Path::new(&gdb_bases[i])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::add_bytes(index_size);

        // Compress if requested
        if config.zstd_compress {
            fastga_integration::FastGAIntegration::compress_index(
                &gdb_bases[i],
                config.zstd_level,
            )?;
        }

        for j in 0..num_batches {
            if !config.quiet {
                if i == j {
                    eprintln!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    eprintln!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            // Run alignment: query batch j against target batch i
            let paf_bytes = fastga.align_direct_paf(&batch_files[j], &batch_files[i])?;

            // Append to merged output
            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                eprintln!("[batch]   {} alignments", line_count);
            }

            // Clean up query index if FastGA auto-created one (j != i)
            if j != i {
                let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[j]);
            }
        }

        // Clean up index for batch i to free disk space
        if !config.quiet {
            eprintln!("[batch] Cleaning up index for batch {}...", i + 1);
        }
        crate::disk_usage::remove_bytes(index_size);
        let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[i]);
    }

    // Phase 3: Clean up all GDB files and batch directory
    for gdb_base in &gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !config.quiet {
        eprintln!(
            "[batch] Completed sequence-level batch alignment: {} total alignments",
            total_alignments
        );
    }

    Ok(merged_paf)
}

/// Configuration for asymmetric batch alignment
pub struct AsymmetricBatchConfig {
    pub frequency: Option<usize>,
    pub threads: usize,
    pub min_alignment_length: u64,
    pub zstd_compress: bool,
    pub zstd_level: u32,
    pub quiet: bool,
    pub bidirectional: bool,
}

/// Run asymmetric batch alignment: queries → targets (optionally bidirectional)
/// Unlike all-vs-all mode, this only aligns queries to targets, not within each set.
/// This is more efficient when you have distinct query and target sets.
pub fn run_asymmetric_batch_alignment(
    query_files: &[String],
    target_files: &[String],
    max_index_bytes: u64,
    config: &AsymmetricBatchConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::{self, FastGAIntegration};

    // In asymmetric mode, query and target indexes don't coexist for the same batch
    // The target index exists, and FastGA may create a query index temporarily
    // So we still use max_index_bytes / 2 for safety
    let per_batch_limit = max_index_bytes / 2;

    if !config.quiet {
        eprintln!(
            "[batch] Asymmetric batch mode: max {} total index size ({} per batch)",
            format_bytes(max_index_bytes),
            format_bytes(per_batch_limit)
        );
        if config.bidirectional {
            eprintln!("[batch] Bidirectional mode: will align A→B and B→A");
        } else {
            eprintln!("[batch] Unidirectional mode: will align queries→targets only");
        }
    }

    // Determine temp directory
    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    // Create temp directory for batch FASTAs
    let batch_dir = temp_base.join(format!("sweepga_asymbatch_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    // Parse and partition queries
    if !config.quiet {
        eprintln!("[batch] Scanning query files for genome sizes...");
    }
    let query_genomes = parse_genome_sizes(query_files)?;
    if !config.quiet {
        eprintln!("[batch] Found {} query genomes, total {} bp",
            query_genomes.len(),
            query_genomes.iter().map(|g| g.total_bp).sum::<u64>());
        for g in &query_genomes {
            eprintln!("[batch]   {} -> {} bp (est. index {})",
                g.prefix, g.total_bp, format_bytes(estimate_index_size(g.total_bp)));
        }
    }
    let query_batches = partition_into_batches_with_limit(query_genomes, per_batch_limit);

    // Parse and partition targets
    if !config.quiet {
        eprintln!("[batch] Scanning target files for genome sizes...");
    }
    let target_genomes = parse_genome_sizes(target_files)?;
    if !config.quiet {
        eprintln!("[batch] Found {} target genomes, total {} bp",
            target_genomes.len(),
            target_genomes.iter().map(|g| g.total_bp).sum::<u64>());
        // Only show first 10 to avoid spam
        for (i, g) in target_genomes.iter().enumerate() {
            if i < 10 {
                eprintln!("[batch]   {} -> {} bp (est. index {})",
                    g.prefix, g.total_bp, format_bytes(estimate_index_size(g.total_bp)));
            } else if i == 10 {
                eprintln!("[batch]   ... and {} more genomes", target_genomes.len() - 10);
                break;
            }
        }
    }
    let target_batches = partition_into_batches_with_limit(target_genomes, per_batch_limit);

    // Check if any batch exceeds the limit - if so, switch to sequence-level splitting
    let query_oversized = query_batches.iter().any(|b| b.estimated_index_bytes > per_batch_limit);
    let target_oversized = target_batches.iter().any(|b| b.estimated_index_bytes > per_batch_limit);

    if query_oversized || target_oversized {
        if !config.quiet {
            eprintln!(
                "[batch] Genome-level batching produces oversized batches, switching to sequence-level splitting..."
            );
        }
        // Clean up and use sequence-level batching instead
        let _ = std::fs::remove_dir_all(&batch_dir);
        return run_asymmetric_sequence_batch_alignment(
            query_files,
            target_files,
            max_index_bytes,
            config,
            tempdir,
        );
    }

    let num_query_batches = query_batches.len();
    let num_target_batches = target_batches.len();
    let num_alignments = if config.bidirectional {
        num_query_batches * num_target_batches * 2
    } else {
        num_query_batches * num_target_batches
    };

    if !config.quiet {
        eprintln!(
            "[batch] Partitioned {} query genomes into {} batches",
            query_batches.iter().map(|b| b.genomes.len()).sum::<usize>(),
            num_query_batches
        );
        for (i, batch) in query_batches.iter().enumerate() {
            eprintln!("[batch]   Query batch {}: {} genomes, {} bp, est. index {}",
                i, batch.genomes.len(), batch.total_bp, format_bytes(batch.estimated_index_bytes));
        }
        eprintln!(
            "[batch] Partitioned {} target genomes into {} batches",
            target_batches.iter().map(|b| b.genomes.len()).sum::<usize>(),
            num_target_batches
        );
        for (i, batch) in target_batches.iter().enumerate() {
            eprintln!("[batch]   Target batch {}: {} genomes, {} bp, est. index {}",
                i, batch.genomes.len(), batch.total_bp, format_bytes(batch.estimated_index_bytes));
        }
        eprintln!("[batch] Per-batch limit: {}", format_bytes(per_batch_limit));
        eprintln!("[batch] Will run {} batch alignments", num_alignments);
    }

    // Create subdirectories for queries and targets
    let query_dir = batch_dir.join("queries");
    let target_dir = batch_dir.join("targets");
    std::fs::create_dir_all(&query_dir)?;
    std::fs::create_dir_all(&target_dir)?;

    // Write query batch FASTAs
    if !config.quiet {
        eprintln!(
            "[batch] Writing {} query batch FASTA files...",
            num_query_batches
        );
    }
    let query_batch_files = write_batches_to_dir(&query_batches, &query_dir, "query")?;

    // Write target batch FASTAs
    if !config.quiet {
        eprintln!(
            "[batch] Writing {} target batch FASTA files...",
            num_target_batches
        );
    }
    let target_batch_files = write_batches_to_dir(&target_batches, &target_dir, "target")?;

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
    );

    // Phase 1: Create GDBs for all batches
    if !config.quiet {
        eprintln!("[batch] Creating GDB files for all batches...");
    }
    let mut query_gdb_bases: Vec<String> = Vec::new();
    let mut target_gdb_bases: Vec<String> = Vec::new();

    for (i, batch_file) in query_batch_files.iter().enumerate() {
        if !config.quiet {
            eprintln!("[batch]   GDB for query batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        query_gdb_bases.push(gdb_base);
    }

    for (i, batch_file) in target_batch_files.iter().enumerate() {
        if !config.quiet {
            eprintln!("[batch]   GDB for target batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        target_gdb_bases.push(gdb_base);
    }

    // Create merged output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;

    // Phase 2: Run alignments (queries → targets)
    for t in 0..num_target_batches {
        // Build index for target batch t
        if !config.quiet {
            eprintln!("[batch] Building index for target batch {}...", t + 1);
        }
        fastga.create_index_only(&target_gdb_bases[t])?;

        // Track disk usage
        let index_dir = std::path::Path::new(&target_gdb_bases[t])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::add_bytes(index_size);

        // Compress if requested
        if config.zstd_compress {
            fastga_integration::FastGAIntegration::compress_index(
                &target_gdb_bases[t],
                config.zstd_level,
            )?;
        }

        // Align all query batches to this target batch
        for q in 0..num_query_batches {
            if !config.quiet {
                eprintln!(
                    "[batch] Aligning query batch {} → target batch {}...",
                    q + 1,
                    t + 1
                );
            }

            let paf_bytes =
                fastga.align_direct_paf(&query_batch_files[q], &target_batch_files[t])?;

            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                eprintln!("[batch]   {} alignments", line_count);
            }

            // Clean up query index if FastGA created one
            let _ = fastga_integration::FastGAIntegration::cleanup_index(&query_gdb_bases[q]);
        }

        // Clean up target index
        if !config.quiet {
            eprintln!("[batch] Cleaning up index for target batch {}...", t + 1);
        }
        crate::disk_usage::remove_bytes(index_size);
        let _ = fastga_integration::FastGAIntegration::cleanup_index(&target_gdb_bases[t]);
    }

    // Phase 2b: Run reverse alignments if bidirectional (targets → queries)
    if config.bidirectional {
        if !config.quiet {
            eprintln!("[batch] Starting reverse direction (targets → queries)...");
        }

        for q in 0..num_query_batches {
            // Build index for query batch q (now as target)
            if !config.quiet {
                eprintln!(
                    "[batch] Building index for query batch {} (as target)...",
                    q + 1
                );
            }
            fastga.create_index_only(&query_gdb_bases[q])?;

            let index_dir = std::path::Path::new(&query_gdb_bases[q])
                .parent()
                .unwrap_or(std::path::Path::new("."));
            let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
            crate::disk_usage::add_bytes(index_size);

            if config.zstd_compress {
                fastga_integration::FastGAIntegration::compress_index(
                    &query_gdb_bases[q],
                    config.zstd_level,
                )?;
            }

            // Align all target batches to this query batch
            for t in 0..num_target_batches {
                if !config.quiet {
                    eprintln!(
                        "[batch] Aligning target batch {} → query batch {}...",
                        t + 1,
                        q + 1
                    );
                }

                let paf_bytes =
                    fastga.align_direct_paf(&target_batch_files[t], &query_batch_files[q])?;

                let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
                total_alignments += line_count;

                merged_output.write_all(&paf_bytes)?;

                if !config.quiet {
                    eprintln!("[batch]   {} alignments", line_count);
                }

                let _ = fastga_integration::FastGAIntegration::cleanup_index(&target_gdb_bases[t]);
            }

            // Clean up query index
            if !config.quiet {
                eprintln!(
                    "[batch] Cleaning up index for query batch {} (as target)...",
                    q + 1
                );
            }
            crate::disk_usage::remove_bytes(index_size);
            let _ = fastga_integration::FastGAIntegration::cleanup_index(&query_gdb_bases[q]);
        }
    }

    // Phase 3: Clean up all GDB files and batch directory
    for gdb_base in &query_gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    for gdb_base in &target_gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !config.quiet {
        eprintln!(
            "[batch] Completed asymmetric batch alignment: {} total alignments",
            total_alignments
        );
    }

    Ok(merged_paf)
}

/// Run asymmetric batch alignment with sequence-level splitting
/// Used when individual genomes are too large for genome-level batching
///
/// Key optimization: For unidirectional alignment (queries → targets):
/// - Only partition and write out TARGET sequences (to control index size)
/// - Use original query files directly (no query partitioning needed)
/// - This saves disk space by avoiding duplicate query data
fn run_asymmetric_sequence_batch_alignment(
    query_files: &[String],
    target_files: &[String],
    max_index_bytes: u64,
    config: &AsymmetricBatchConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::{self, FastGAIntegration};

    let per_batch_limit = max_index_bytes / 2;

    if !config.quiet {
        eprintln!(
            "[batch] Asymmetric SEQUENCE-level batch mode: max {} total index size ({} per batch)",
            format_bytes(max_index_bytes),
            format_bytes(per_batch_limit)
        );
        if config.bidirectional {
            eprintln!("[batch] Bidirectional mode: will partition both queries and targets");
        } else {
            eprintln!("[batch] Unidirectional mode: only partitioning targets (queries used directly)");
        }
    }

    // Determine temp directory
    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    // Create temp directory for batch FASTAs
    let batch_dir = temp_base.join(format!("sweepga_asymseqbatch_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    // Parse sequences from target files (always needed)
    if !config.quiet {
        eprintln!("[batch] Parsing target sequences...");
    }
    let target_sequences = parse_sequences(target_files)?;
    let target_total_bp: u64 = target_sequences.iter().map(|s| s.length).sum();

    if !config.quiet {
        eprintln!(
            "[batch] Found {} target sequences totaling {}",
            target_sequences.len(),
            format_bytes(target_total_bp)
        );
    }

    // Partition target sequences into batches
    let target_batches = partition_sequences_into_batches(target_sequences, max_index_bytes);
    let num_target_batches = target_batches.len();

    // For bidirectional mode, we also need to partition queries (they become targets in reverse)
    let (query_batches, num_query_batches) = if config.bidirectional {
        if !config.quiet {
            eprintln!("[batch] Parsing query sequences (for bidirectional mode)...");
        }
        let query_sequences = parse_sequences(query_files)?;
        let query_total_bp: u64 = query_sequences.iter().map(|s| s.length).sum();

        if !config.quiet {
            eprintln!(
                "[batch] Found {} query sequences totaling {}",
                query_sequences.len(),
                format_bytes(query_total_bp)
            );
        }

        let batches = partition_sequences_into_batches(query_sequences, max_index_bytes);
        let num = batches.len();
        (Some(batches), num)
    } else {
        (None, 1) // Treat all queries as a single "batch" (the original files)
    };

    let num_alignments = if config.bidirectional {
        num_target_batches + num_query_batches // Each direction: all queries to each target batch
    } else {
        num_target_batches // One alignment per target batch
    };

    if !config.quiet {
        eprintln!(
            "[batch] Partitioned target sequences into {} batches:",
            num_target_batches
        );
        for (i, batch) in target_batches.iter().enumerate() {
            eprintln!(
                "[batch]   Target batch {}: {} sequences, {} bp, est. index {}",
                i,
                batch.sequences.len(),
                batch.total_bp,
                format_bytes(batch.estimated_index_bytes)
            );
        }
        if let Some(ref qb) = query_batches {
            eprintln!(
                "[batch] Partitioned query sequences into {} batches:",
                num_query_batches
            );
            for (i, batch) in qb.iter().enumerate() {
                eprintln!(
                    "[batch]   Query batch {}: {} sequences, {} bp, est. index {}",
                    i,
                    batch.sequences.len(),
                    batch.total_bp,
                    format_bytes(batch.estimated_index_bytes)
                );
            }
        } else {
            eprintln!("[batch] Query files will be used directly (no partitioning)");
        }
        eprintln!("[batch] Per-batch limit: {}", format_bytes(per_batch_limit));
        eprintln!("[batch] Will run {} batch alignments", num_alignments);
    }

    // Create target directory
    let target_dir = batch_dir.join("targets");
    std::fs::create_dir_all(&target_dir)?;

    // Write target batch FASTAs
    if !config.quiet {
        eprintln!(
            "[batch] Writing {} target batch FASTA files...",
            num_target_batches
        );
    }
    let target_batch_files = write_sequence_batches_to_dir(&target_batches, &target_dir, "target")?;

    // Write query batch FASTAs only if bidirectional
    let query_batch_files = if let Some(ref qb) = query_batches {
        let query_dir = batch_dir.join("queries");
        std::fs::create_dir_all(&query_dir)?;
        if !config.quiet {
            eprintln!(
                "[batch] Writing {} query batch FASTA files...",
                num_query_batches
            );
        }
        Some(write_sequence_batches_to_dir(qb, &query_dir, "query")?)
    } else {
        None
    };

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
    );

    // Phase 1: Create GDBs for target batches
    if !config.quiet {
        eprintln!("[batch] Creating GDB files for target batches...");
    }
    let mut target_gdb_bases: Vec<String> = Vec::new();

    for (i, batch_file) in target_batch_files.iter().enumerate() {
        if !config.quiet {
            eprintln!("[batch]   GDB for target batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        target_gdb_bases.push(gdb_base);
    }

    // For bidirectional mode, also create GDBs for query batches
    let query_gdb_bases: Option<Vec<String>> = if let Some(ref qbf) = query_batch_files {
        if !config.quiet {
            eprintln!("[batch] Creating GDB files for query batches...");
        }
        let mut bases = Vec::new();
        for (i, batch_file) in qbf.iter().enumerate() {
            if !config.quiet {
                eprintln!("[batch]   GDB for query batch {}...", i + 1);
            }
            let gdb_base = fastga.create_gdb_only(batch_file)?;
            bases.push(gdb_base);
        }
        Some(bases)
    } else {
        None
    };

    // Create merged output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;

    // Phase 2: Run alignments (queries → targets)
    // For unidirectional mode: use original query files directly
    // For bidirectional mode: use query batch files
    for t in 0..num_target_batches {
        // Build index for target batch t
        if !config.quiet {
            eprintln!("[batch] Building index for target batch {}...", t + 1);
        }
        fastga.create_index_only(&target_gdb_bases[t])?;

        // Track disk usage
        let index_dir = std::path::Path::new(&target_gdb_bases[t])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::add_bytes(index_size);

        if config.zstd_compress {
            fastga_integration::FastGAIntegration::compress_index(
                &target_gdb_bases[t],
                config.zstd_level,
            )?;
        }

        // Align queries to this target batch
        if let Some(ref qbf) = query_batch_files {
            // Bidirectional mode: align each query batch
            for (q, query_file) in qbf.iter().enumerate() {
                if !config.quiet {
                    eprintln!(
                        "[batch] Aligning query batch {} → target batch {}...",
                        q + 1,
                        t + 1
                    );
                }

                let paf_bytes = fastga.align_direct_paf(query_file, &target_batch_files[t])?;

                let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
                total_alignments += line_count;

                merged_output.write_all(&paf_bytes)?;

                if !config.quiet {
                    eprintln!("[batch]   {} alignments", line_count);
                }
            }
        } else {
            // Unidirectional mode: align original query files directly
            for (q, query_file) in query_files.iter().enumerate() {
                if !config.quiet {
                    if query_files.len() == 1 {
                        eprintln!(
                            "[batch] Aligning queries → target batch {}...",
                            t + 1
                        );
                    } else {
                        eprintln!(
                            "[batch] Aligning query file {} → target batch {}...",
                            q + 1,
                            t + 1
                        );
                    }
                }

                // Canonicalize the query file path
                let query_path = std::fs::canonicalize(query_file)
                    .with_context(|| format!("Failed to resolve path: {}", query_file))?;

                let paf_bytes = fastga.align_direct_paf(&query_path, &target_batch_files[t])?;

                let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
                total_alignments += line_count;

                merged_output.write_all(&paf_bytes)?;

                if !config.quiet {
                    eprintln!("[batch]   {} alignments", line_count);
                }
            }
        }

        // Clean up target index for this batch
        if !config.quiet {
            eprintln!("[batch] Cleaning up index for target batch {}...", t + 1);
        }
        crate::disk_usage::remove_bytes(index_size);
        let _ = fastga_integration::FastGAIntegration::cleanup_index(&target_gdb_bases[t]);
    }

    // Phase 2b: Run reverse alignments if bidirectional
    if config.bidirectional {
        let query_gdb_bases = query_gdb_bases.as_ref().expect("Query GDBs required for bidirectional");
        let query_batch_files = query_batch_files.as_ref().expect("Query batch files required for bidirectional");

        if !config.quiet {
            eprintln!("[batch] Starting reverse direction (targets → queries)...");
        }

        for (q, _) in query_batch_files.iter().enumerate() {
            // Build index for query batch q (now as target)
            if !config.quiet {
                eprintln!(
                    "[batch] Building index for query batch {} (as target)...",
                    q + 1
                );
            }
            fastga.create_index_only(&query_gdb_bases[q])?;

            let index_dir = std::path::Path::new(&query_gdb_bases[q])
                .parent()
                .unwrap_or(std::path::Path::new("."));
            let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
            crate::disk_usage::add_bytes(index_size);

            if config.zstd_compress {
                fastga_integration::FastGAIntegration::compress_index(
                    &query_gdb_bases[q],
                    config.zstd_level,
                )?;
            }

            // Align all target batches to this query batch
            for t in 0..num_target_batches {
                if !config.quiet {
                    eprintln!(
                        "[batch] Aligning target batch {} → query batch {}...",
                        t + 1,
                        q + 1
                    );
                }

                let paf_bytes =
                    fastga.align_direct_paf(&target_batch_files[t], &query_batch_files[q])?;

                let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
                total_alignments += line_count;

                merged_output.write_all(&paf_bytes)?;

                if !config.quiet {
                    eprintln!("[batch]   {} alignments", line_count);
                }
            }

            // Clean up query index
            if !config.quiet {
                eprintln!(
                    "[batch] Cleaning up index for query batch {} (as target)...",
                    q + 1
                );
            }
            crate::disk_usage::remove_bytes(index_size);
            let _ = fastga_integration::FastGAIntegration::cleanup_index(&query_gdb_bases[q]);
        }
    }

    // Phase 3: Clean up all GDB files and batch directory
    if let Some(ref qgdb) = query_gdb_bases {
        for gdb_base in qgdb {
            let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
        }
    }
    for gdb_base in &target_gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !config.quiet {
        eprintln!(
            "[batch] Completed asymmetric sequence-level batch alignment: {} total alignments",
            total_alignments
        );
    }

    Ok(merged_paf)
}

/// Write sequence batches to a directory, returning the paths to the FASTA files
fn write_sequence_batches_to_dir(
    batches: &[SequenceBatch],
    dir: &Path,
    prefix: &str,
) -> Result<Vec<PathBuf>> {
    use std::sync::Mutex;

    // Build sequence_name -> batch index mapping
    let mut seq_to_batch: HashMap<String, usize> = HashMap::new();
    for (batch_idx, batch) in batches.iter().enumerate() {
        for seq in &batch.sequences {
            seq_to_batch.insert(seq.name.clone(), batch_idx);
        }
    }

    // Collect all unique source files
    let mut source_files: Vec<PathBuf> = batches
        .iter()
        .flat_map(|b| b.sequences.iter().map(|s| s.source_file.clone()))
        .collect();
    source_files.sort();
    source_files.dedup();

    // Create batch files with Mutex protection
    let mut batch_paths: Vec<PathBuf> = Vec::with_capacity(batches.len());
    let mut batch_writers: Vec<Mutex<std::io::BufWriter<File>>> = Vec::with_capacity(batches.len());

    for i in 0..batches.len() {
        let batch_subdir = dir.join(format!("{}_{}", prefix, i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("sequences.fa");
        let file = File::create(&batch_path)
            .with_context(|| format!("Failed to create batch FASTA: {}", batch_path.display()))?;
        batch_writers.push(Mutex::new(std::io::BufWriter::with_capacity(256 * 1024, file)));
        batch_paths.push(batch_path);
    }

    // Process source files in parallel
    let results: Vec<Result<()>> = source_files
        .par_iter()
        .map(|source_file| {
            stream_write_sequences_from_source(source_file, &seq_to_batch, &batch_writers)
        })
        .collect();

    for result in results {
        result?;
    }

    // Flush all writers
    for writer in &batch_writers {
        writer.lock().unwrap().flush()?;
    }

    Ok(batch_paths)
}

/// Stream sequences from a source file to the appropriate batch files
fn stream_write_sequences_from_source(
    source_file: &Path,
    seq_to_batch: &HashMap<String, usize>,
    batch_writers: &[std::sync::Mutex<std::io::BufWriter<File>>],
) -> Result<()> {
    use std::io::Write;

    let file = File::open(source_file)
        .with_context(|| format!("Failed to open source FASTA: {}", source_file.display()))?;

    let reader: Box<dyn BufRead> = if source_file
        .to_string_lossy()
        .ends_with(".gz")
        || source_file.to_string_lossy().ends_with(".bgz")
    {
        Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut current_batch_idx: Option<usize> = None;
    let mut buffer = Vec::with_capacity(64 * 1024);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with('>') {
            // Flush previous sequence's buffer
            if let Some(batch_idx) = current_batch_idx {
                if !buffer.is_empty() {
                    let mut writer = batch_writers[batch_idx].lock().unwrap();
                    writer.write_all(&buffer)?;
                    buffer.clear();
                }
            }

            // Extract sequence name
            let name = trimmed
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("");

            // Look up which batch this sequence belongs to
            current_batch_idx = seq_to_batch.get(name).copied();

            if current_batch_idx.is_some() {
                buffer.extend_from_slice(line.as_bytes());
                buffer.push(b'\n');
            }
        } else if !trimmed.is_empty() {
            if current_batch_idx.is_some() {
                buffer.extend_from_slice(line.as_bytes());
                buffer.push(b'\n');
            }
        }
    }

    // Flush final buffer
    if let Some(batch_idx) = current_batch_idx {
        if !buffer.is_empty() {
            let mut writer = batch_writers[batch_idx].lock().unwrap();
            writer.write_all(&buffer)?;
        }
    }

    Ok(())
}

/// Partition genomes into batches with a specific per-batch limit (no /2 division)
fn partition_into_batches_with_limit(
    genomes: Vec<GenomeInfo>,
    per_batch_limit: u64,
) -> Vec<GenomeBatch> {
    let mut batches: Vec<GenomeBatch> = Vec::new();

    for genome in genomes {
        let single_genome_size = estimate_index_size(genome.total_bp);
        if single_genome_size > per_batch_limit {
            eprintln!(
                "[batch] WARNING: Genome {} ({}) has estimated index {} which exceeds per-batch limit {}",
                genome.prefix,
                format_bytes(genome.total_bp),
                format_bytes(single_genome_size),
                format_bytes(per_batch_limit)
            );
            eprintln!("[batch] Including it as a single-genome batch anyway");

            let mut oversized = GenomeBatch::new();
            oversized.add(genome);
            batches.push(oversized);
            continue;
        }

        let mut best_fit_idx: Option<usize> = None;
        let mut best_fit_remaining = u64::MAX;

        for (idx, batch) in batches.iter().enumerate() {
            if !batch.would_exceed(&genome, per_batch_limit) {
                let new_size = estimate_index_size(batch.total_bp + genome.total_bp);
                let remaining = per_batch_limit.saturating_sub(new_size);

                if remaining < best_fit_remaining {
                    best_fit_remaining = remaining;
                    best_fit_idx = Some(idx);
                }
            }
        }

        if let Some(idx) = best_fit_idx {
            batches[idx].add(genome);
        } else {
            let mut new_batch = GenomeBatch::new();
            new_batch.add(genome);
            batches.push(new_batch);
        }
    }

    batches
}

/// Write batches to a directory, returning the paths to the FASTA files
fn write_batches_to_dir(
    batches: &[GenomeBatch],
    dir: &Path,
    prefix: &str,
) -> Result<Vec<PathBuf>> {
    use std::sync::Mutex;

    // Build genome_prefix -> batch indices mapping
    let mut prefix_to_batches: HashMap<String, Vec<usize>> = HashMap::new();
    for (batch_idx, batch) in batches.iter().enumerate() {
        for genome in &batch.genomes {
            prefix_to_batches
                .entry(genome.prefix.clone())
                .or_default()
                .push(batch_idx);
        }
    }

    // Collect all unique source files
    let mut source_files: Vec<PathBuf> = batches
        .iter()
        .flat_map(|b| b.genomes.iter().map(|g| g.source_file.clone()))
        .collect();
    source_files.sort();
    source_files.dedup();

    // Create batch files with Mutex protection
    let mut batch_paths: Vec<PathBuf> = Vec::with_capacity(batches.len());
    let mut batch_writers: Vec<Mutex<std::io::BufWriter<File>>> = Vec::with_capacity(batches.len());

    for i in 0..batches.len() {
        let batch_subdir = dir.join(format!("{}_{}", prefix, i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("genomes.fa");
        let file = File::create(&batch_path)
            .with_context(|| format!("Failed to create batch FASTA: {}", batch_path.display()))?;
        batch_writers.push(Mutex::new(std::io::BufWriter::with_capacity(256 * 1024, file)));
        batch_paths.push(batch_path);
    }

    // Process source files in parallel
    let results: Vec<Result<()>> = source_files
        .par_iter()
        .map(|source_file| {
            stream_write_batches_from_source_parallel(source_file, &prefix_to_batches, &batch_writers)
        })
        .collect();

    for result in results {
        result?;
    }

    // Flush all writers
    for writer in &batch_writers {
        writer.lock().unwrap().flush()?;
    }

    Ok(batch_paths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_index_size() {
        // 100 Mbp -> ~1.3 GB
        let size_100m = estimate_index_size(100_000_000);
        assert!(size_100m > 1_200_000_000 && size_100m < 1_400_000_000);

        // 1 Gbp -> ~12.1 GB
        let size_1g = estimate_index_size(1_000_000_000);
        assert!(size_1g > 12_000_000_000 && size_1g < 12_200_000_000);
    }

    #[test]
    fn test_format_bytes() {
        assert_eq!(format_bytes(500), "500 bytes");
        assert_eq!(format_bytes(1500), "1.5 KB");
        assert_eq!(format_bytes(1_500_000), "1.5 MB");
        assert_eq!(format_bytes(1_500_000_000), "1.5 GB");
    }

    #[test]
    fn test_extract_pansn_prefix() {
        assert_eq!(extract_pansn_prefix("SGDref#1#chrI"), "SGDref#1#");
        assert_eq!(extract_pansn_prefix("genome#hap#seq"), "genome#hap#");
        assert_eq!(extract_pansn_prefix("simple"), "simple#");
    }
}
