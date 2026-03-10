//! Batch alignment for large genome sets
//!
//! When aligning many genomes, the aligner's resource usage can exceed available
//! disk space (fastga) or memory (wfmash). This module partitions genomes into
//! batches based on a cost model, runs the aligner on each batch pair, and
//! aggregates results.
//!
//! Cost models:
//!   FastGA (disk): GIX_size ≈ 0.1 GB + 12 bytes/bp
//!   Wfmash (memory): sketch ≈ 0.5 GB + 20 bytes/bp

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Estimated index overhead in bytes (fixed cost per index) — FastGA
const FASTGA_OVERHEAD_BYTES: u64 = 100_000_000; // 100 MB

/// Estimated bytes of index per basepair of genome — FastGA
const FASTGA_BYTES_PER_BP: f64 = 12.0;

/// Estimated memory overhead — wfmash
const WFMASH_OVERHEAD_BYTES: u64 = 500_000_000; // 500 MB

/// Estimated bytes of memory per basepair — wfmash
const WFMASH_BYTES_PER_BP: f64 = 20.0;

/// Default auto batch limit for FastGA (disk): ~50 GB
const FASTGA_AUTO_BYTES: u64 = 50_000_000_000;

/// Default auto batch limit for wfmash (memory): ~4 GB
const WFMASH_AUTO_BYTES: u64 = 4_000_000_000;

/// Resource limits for batch alignment.
#[derive(Debug, Clone)]
pub struct BatchLimits {
    /// Maximum disk bytes for index files per batch.
    pub max_disk_bytes: u64,
    /// Optional maximum memory bytes. When set, batches are also bounded by
    /// estimated peak RAM usage.
    pub max_memory_bytes: Option<u64>,
}

impl Default for BatchLimits {
    fn default() -> Self {
        Self {
            max_disk_bytes: FASTGA_AUTO_BYTES,
            max_memory_bytes: None,
        }
    }
}

/// Estimate peak memory usage for aligning a batch of genomes (wfmash model).
pub fn estimate_memory_usage(genome_bp: u64) -> u64 {
    WFMASH_OVERHEAD_BYTES + (genome_bp as f64 * WFMASH_BYTES_PER_BP) as u64
}

/// Trait for aligner-specific batch operations.
///
/// The batch loop calls these methods in order:
/// 1. `prepare_all()` — one-time setup for all batch files (e.g., create GDBs)
/// 2. For each target batch:
///    a. `prepare_target()` — build index for target
///    b. For each query batch: `align()` — run alignment
///    c. `cleanup_target()` — remove target index
/// 3. `cleanup_all()` — final cleanup (e.g., remove GDBs)
pub trait BatchAligner {
    /// One-time setup for all batch files (e.g., create GDB files).
    fn prepare_all(&self, batch_files: &[PathBuf], quiet: bool) -> Result<()>;

    /// Prepare a target batch for alignment (e.g., build index).
    fn prepare_target(&self, batch_idx: usize, quiet: bool) -> Result<()>;

    /// Run alignment of query batch against target batch. Returns PAF bytes.
    fn align(&self, query: &Path, target: &Path) -> Result<Vec<u8>>;

    /// Clean up target-specific resources (e.g., remove index files).
    fn cleanup_target(&self, batch_idx: usize, quiet: bool) -> Result<()>;

    /// Final cleanup of all resources (e.g., remove GDB files).
    fn cleanup_all(&self) -> Result<()>;

    /// Estimate resource cost in bytes for a batch of the given total basepairs.
    fn estimate_cost(&self, genome_bp: u64) -> u64;

    /// Label for the resource being limited (e.g., "disk" or "memory").
    fn resource_label(&self) -> &str;

    /// Default auto limit in bytes.
    fn auto_limit(&self) -> u64;

    /// Run a single-batch alignment (no partitioning needed).
    fn align_single(&self, fasta_files: &[String], tempdir: Option<&str>) -> Result<tempfile::NamedTempFile>;
}

/// FastGA batch aligner — limits disk usage from GIX index files.
pub struct FastGABatchAligner {
    fastga: crate::fastga_integration::FastGAIntegration,
    gdb_bases: std::cell::RefCell<Vec<String>>,
    zstd_compress: bool,
    zstd_level: u32,
}

impl FastGABatchAligner {
    pub fn new(
        frequency: Option<usize>,
        threads: usize,
        min_alignment_length: u64,
        tempdir: Option<String>,
        zstd_compress: bool,
        zstd_level: u32,
    ) -> Self {
        Self {
            fastga: crate::fastga_integration::FastGAIntegration::new(
                frequency,
                threads,
                min_alignment_length,
                tempdir,
            ),
            gdb_bases: std::cell::RefCell::new(Vec::new()),
            zstd_compress,
            zstd_level,
        }
    }
}

impl BatchAligner for FastGABatchAligner {
    fn prepare_all(&self, batch_files: &[PathBuf], quiet: bool) -> Result<()> {
        if !quiet {
            eprintln!("[batch] Creating GDB files for all batches...");
        }
        let mut gdb_bases = self.gdb_bases.borrow_mut();
        gdb_bases.clear();
        for (i, batch_file) in batch_files.iter().enumerate() {
            if !quiet {
                eprintln!("[batch]   GDB for batch {}...", i + 1);
            }
            let gdb_base = self.fastga.create_gdb_only(batch_file)?;
            gdb_bases.push(gdb_base);
        }
        Ok(())
    }

    fn prepare_target(&self, batch_idx: usize, quiet: bool) -> Result<()> {
        let gdb_bases = self.gdb_bases.borrow();
        if !quiet {
            eprintln!("[batch] Building index for batch {}...", batch_idx + 1);
        }
        self.fastga.create_index_only(&gdb_bases[batch_idx])?;

        // Track disk usage of index files
        let index_dir = std::path::Path::new(&gdb_bases[batch_idx])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::add_bytes(index_size);

        if self.zstd_compress {
            crate::fastga_integration::FastGAIntegration::compress_index(
                &gdb_bases[batch_idx],
                self.zstd_level,
            )?;
        }
        Ok(())
    }

    fn align(&self, query: &Path, target: &Path) -> Result<Vec<u8>> {
        self.fastga.align_direct_paf(query, target)
    }

    fn cleanup_target(&self, batch_idx: usize, quiet: bool) -> Result<()> {
        let gdb_bases = self.gdb_bases.borrow();
        if !quiet {
            eprintln!("[batch] Cleaning up index for batch {}...", batch_idx + 1);
        }
        // Estimate and untrack disk usage
        let index_dir = std::path::Path::new(&gdb_bases[batch_idx])
            .parent()
            .unwrap_or(std::path::Path::new("."));
        let index_size = crate::disk_usage::scan_fastga_index_files(index_dir).unwrap_or(0);
        crate::disk_usage::remove_bytes(index_size);
        let _ = crate::fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[batch_idx]);
        Ok(())
    }

    fn cleanup_all(&self) -> Result<()> {
        let gdb_bases = self.gdb_bases.borrow();
        for gdb_base in gdb_bases.iter() {
            let _ = crate::fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
        }
        Ok(())
    }

    fn estimate_cost(&self, genome_bp: u64) -> u64 {
        estimate_index_size(genome_bp)
    }

    fn resource_label(&self) -> &str {
        "disk"
    }

    fn auto_limit(&self) -> u64 {
        FASTGA_AUTO_BYTES
    }

    fn align_single(&self, fasta_files: &[String], _tempdir: Option<&str>) -> Result<tempfile::NamedTempFile> {
        let input_path = std::fs::canonicalize(&fasta_files[0])
            .with_context(|| format!("Failed to resolve path: {}", fasta_files[0]))?;
        self.fastga.align_to_temp_paf(&input_path, &input_path)
    }
}

/// Wfmash batch aligner — limits memory usage from sketching.
pub struct WfmashBatchAligner {
    num_threads: usize,
    min_alignment_length: Option<u64>,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
}

impl WfmashBatchAligner {
    pub fn new(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
    ) -> Self {
        Self {
            num_threads,
            min_alignment_length,
            map_pct_identity,
            temp_dir,
        }
    }

    /// Create a WfmashIntegration configured for a specific batch FASTA.
    fn create_wfmash_for(&self, fasta_path: &Path) -> Result<crate::wfmash_integration::WfmashIntegration> {
        let avg_seq = avg_seq_len_from_fai(fasta_path).ok();
        crate::wfmash_integration::WfmashIntegration::adaptive(
            self.num_threads,
            self.min_alignment_length,
            self.map_pct_identity.clone(),
            self.temp_dir.clone(),
            None, // segment_length: adaptive
            avg_seq,
            None, // sparsify
            None, // num_mappings
            None, // pairs_file
        )
    }
}

impl BatchAligner for WfmashBatchAligner {
    fn prepare_all(&self, _batch_files: &[PathBuf], _quiet: bool) -> Result<()> {
        Ok(()) // wfmash indexes on each call, no upfront prep needed
    }

    fn prepare_target(&self, _batch_idx: usize, _quiet: bool) -> Result<()> {
        Ok(()) // no-op
    }

    fn align(&self, query: &Path, target: &Path) -> Result<Vec<u8>> {
        use crate::aligner::Aligner;
        let wfmash = self.create_wfmash_for(target)?;
        wfmash.align_direct_paf(query, target)
    }

    fn cleanup_target(&self, _batch_idx: usize, _quiet: bool) -> Result<()> {
        Ok(()) // no-op
    }

    fn cleanup_all(&self) -> Result<()> {
        Ok(()) // no-op
    }

    fn estimate_cost(&self, genome_bp: u64) -> u64 {
        estimate_memory_usage(genome_bp)
    }

    fn resource_label(&self) -> &str {
        "memory"
    }

    fn auto_limit(&self) -> u64 {
        WFMASH_AUTO_BYTES
    }

    fn align_single(&self, fasta_files: &[String], _tempdir: Option<&str>) -> Result<tempfile::NamedTempFile> {
        use crate::aligner::Aligner;
        let input_path = std::fs::canonicalize(&fasta_files[0])
            .with_context(|| format!("Failed to resolve path: {}", fasta_files[0]))?;
        let wfmash = self.create_wfmash_for(&input_path)?;
        wfmash.align_to_temp_paf(&input_path, &input_path)
    }
}

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

/// A batch of genomes to align together
#[derive(Debug, Clone)]
pub struct GenomeBatch {
    /// Genomes in this batch
    pub genomes: Vec<GenomeInfo>,
    /// Total basepairs in batch
    pub total_bp: u64,
    /// Estimated resource cost in bytes
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
            estimated_index_bytes: FASTGA_OVERHEAD_BYTES,
        }
    }

    /// Add a genome to this batch, using the given cost function.
    pub fn add_with_cost(&mut self, genome: GenomeInfo, cost_fn: &dyn Fn(u64) -> u64) {
        self.total_bp += genome.total_bp;
        self.estimated_index_bytes = cost_fn(self.total_bp);
        self.genomes.push(genome);
    }

    /// Add a genome to this batch (FastGA cost model for backwards compat).
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

    /// Check if adding a genome would exceed the limit using a custom cost function.
    pub fn would_exceed_with_cost(&self, genome: &GenomeInfo, max_bytes: u64, cost_fn: &dyn Fn(u64) -> u64) -> bool {
        let new_total = self.total_bp + genome.total_bp;
        cost_fn(new_total) > max_bytes
    }

    /// Check if adding a genome would exceed BatchLimits (disk and/or memory).
    pub fn would_exceed_limits(&self, genome: &GenomeInfo, limits: &BatchLimits) -> bool {
        let new_total = self.total_bp + genome.total_bp;
        if estimate_index_size(new_total) > limits.max_disk_bytes {
            return true;
        }
        if let Some(max_mem) = limits.max_memory_bytes {
            if estimate_memory_usage(new_total) > max_mem {
                return true;
            }
        }
        false
    }
}

/// Estimate FastGA index size in bytes for a given genome size
pub fn estimate_index_size(genome_bp: u64) -> u64 {
    FASTGA_OVERHEAD_BYTES + (genome_bp as f64 * FASTGA_BYTES_PER_BP) as u64
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

/// Parse genome sizes from FASTA file(s)
/// Returns map of genome prefix -> GenomeInfo
pub fn parse_genome_sizes(fasta_files: &[String]) -> Result<Vec<GenomeInfo>> {
    let mut genomes: HashMap<String, GenomeInfo> = HashMap::new();

    for fasta_file in fasta_files {
        let path = Path::new(fasta_file);
        let file =
            File::open(path).with_context(|| format!("Failed to open FASTA: {}", fasta_file))?;

        let reader: Box<dyn BufRead> =
            if fasta_file.ends_with(".gz") || fasta_file.ends_with(".bgz") {
                Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
            } else {
                Box::new(BufReader::new(file))
            };

        let mut current_prefix: Option<String> = None;
        let mut current_bp: u64 = 0;

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            if trimmed.starts_with('>') {
                // Save previous genome's data
                if let Some(prefix) = current_prefix.take() {
                    let entry = genomes.entry(prefix.clone()).or_insert_with(|| GenomeInfo {
                        prefix: prefix.clone(),
                        total_bp: 0,
                        source_file: path.to_path_buf(),
                    });
                    entry.total_bp += current_bp;
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
            let entry = genomes.entry(prefix.clone()).or_insert_with(|| GenomeInfo {
                prefix: prefix.clone(),
                total_bp: 0,
                source_file: path.to_path_buf(),
            });
            entry.total_bp += current_bp;
        }
    }

    // Convert to sorted vector for deterministic ordering
    let mut genome_list: Vec<GenomeInfo> = genomes.into_values().collect();
    genome_list.sort_by(|a, b| a.prefix.cmp(&b.prefix));

    Ok(genome_list)
}

/// Compute average sequence length from a FASTA index (.fai) file.
fn avg_seq_len_from_fai(path: &Path) -> Result<u64> {
    let fai_path = path.with_extension(
        path.extension()
            .map(|e| format!("{}.fai", e.to_string_lossy()))
            .unwrap_or_else(|| "fai".to_string()),
    );
    let fai = std::fs::read_to_string(&fai_path)
        .with_context(|| format!("FASTA index not found: {}", fai_path.display()))?;
    let mut total_bases: u64 = 0;
    let mut num_seqs: u64 = 0;
    for line in fai.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            if let Ok(len) = fields[1].parse::<u64>() {
                total_bases += len;
                num_seqs += 1;
            }
        }
    }
    anyhow::ensure!(num_seqs > 0, "FASTA index {} contains no sequences", fai_path.display());
    Ok(total_bases / num_seqs)
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

/// Partition genomes into batches based on max index size (FastGA cost model).
pub fn partition_into_batches(genomes: Vec<GenomeInfo>, max_index_bytes: u64) -> Vec<GenomeBatch> {
    partition_into_batches_with_cost(genomes, max_index_bytes, &estimate_index_size)
}

/// Partition genomes into batches using a custom cost function.
pub fn partition_into_batches_with_cost(
    genomes: Vec<GenomeInfo>,
    max_bytes: u64,
    cost_fn: &dyn Fn(u64) -> u64,
) -> Vec<GenomeBatch> {
    let mut batches: Vec<GenomeBatch> = Vec::new();
    let mut current_batch = GenomeBatch::new();

    for genome in genomes {
        // Check if this genome alone exceeds the limit
        let single_genome_cost = cost_fn(genome.total_bp);
        if single_genome_cost > max_bytes {
            eprintln!(
                "[batch] WARNING: Genome {} ({}) has estimated cost {} which exceeds limit {}",
                genome.prefix,
                format_bytes(genome.total_bp),
                format_bytes(single_genome_cost),
                format_bytes(max_bytes)
            );
            eprintln!("[batch] Including it as a single-genome batch anyway");

            // Flush current batch if non-empty
            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
                current_batch = GenomeBatch::new();
            }

            // Add oversized genome as its own batch
            let mut oversized = GenomeBatch::new();
            oversized.add_with_cost(genome, cost_fn);
            batches.push(oversized);
            continue;
        }

        // Check if adding to current batch would exceed limit
        if current_batch.would_exceed_with_cost(&genome, max_bytes, cost_fn) {
            // Start new batch
            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
            }
            current_batch = GenomeBatch::new();
        }

        current_batch.add_with_cost(genome, cost_fn);
    }

    // Don't forget the last batch
    if !current_batch.genomes.is_empty() {
        batches.push(current_batch);
    }

    batches
}

/// Write a subset of genomes to a temporary FASTA file
pub fn write_batch_fasta(batch: &GenomeBatch, output_path: &Path) -> Result<()> {
    let mut output = File::create(output_path)
        .with_context(|| format!("Failed to create batch FASTA: {}", output_path.display()))?;

    // Group genomes by source file for efficient reading
    let mut by_source: HashMap<PathBuf, Vec<&str>> = HashMap::new();
    for genome in &batch.genomes {
        by_source
            .entry(genome.source_file.clone())
            .or_default()
            .push(&genome.prefix);
    }

    for (source_file, prefixes) in by_source {
        let file = File::open(&source_file)
            .with_context(|| format!("Failed to open source FASTA: {}", source_file.display()))?;

        let source_str = source_file.to_string_lossy();
        let reader: Box<dyn BufRead> =
            if source_str.ends_with(".gz") || source_str.ends_with(".bgz") {
                Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
            } else {
                Box::new(BufReader::new(file))
            };

        let mut include_current = false;

        for line in reader.lines() {
            let line = line?;

            if line.starts_with('>') {
                // Check if this sequence belongs to a genome we want
                let name = line
                    .trim_start_matches('>')
                    .split_whitespace()
                    .next()
                    .unwrap_or("");
                let prefix = extract_pansn_prefix(name);
                include_current = prefixes.iter().any(|p| *p == prefix);
            }

            if include_current {
                writeln!(output, "{}", line)?;
            }
        }
    }

    Ok(())
}

/// Report batch statistics, including the resource type being limited.
pub fn report_batch_plan(batches: &[GenomeBatch], total_bp: u64, resource_label: &str, quiet: bool) {
    if quiet {
        return;
    }

    let total_genomes: usize = batches.iter().map(|b| b.genomes.len()).sum();
    let total_pairs = batches.len() * batches.len();
    let self_pairs = batches.len();
    let cross_pairs = total_pairs - self_pairs;

    eprintln!(
        "[batch] Partitioned {} genomes ({}) into {} batches (limiting {}):",
        total_genomes,
        format_bytes(total_bp),
        batches.len(),
        resource_label,
    );

    for (i, batch) in batches.iter().enumerate() {
        eprintln!(
            "[batch]   Batch {}: {} genomes, {}, est. {} {}",
            i + 1,
            batch.genomes.len(),
            format_bytes(batch.total_bp),
            format_bytes(batch.estimated_index_bytes),
            resource_label,
        );
    }

    eprintln!(
        "[batch] Will run {} batch alignments ({} self + {} cross-batch)",
        total_pairs, self_pairs, cross_pairs
    );
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

/// Resolve "auto" batch bytes for a given aligner.
pub fn resolve_batch_bytes(value: &str, aligner: &dyn BatchAligner) -> Result<u64> {
    if value.eq_ignore_ascii_case("auto") {
        Ok(aligner.auto_limit())
    } else {
        parse_size_string(value)
    }
}

/// Parse a human-readable size string like "2G", "500M", "100k" into bytes.
pub fn parse_size_string(s: &str) -> Result<u64> {
    let s = s.trim();
    if s.is_empty() {
        anyhow::bail!("Empty size string");
    }

    let (num_str, multiplier) = match s.as_bytes().last() {
        Some(b'k' | b'K') => (&s[..s.len() - 1], 1_000u64),
        Some(b'm' | b'M') => (&s[..s.len() - 1], 1_000_000u64),
        Some(b'g' | b'G') => (&s[..s.len() - 1], 1_000_000_000u64),
        Some(b't' | b'T') => (&s[..s.len() - 1], 1_000_000_000_000u64),
        _ => (s, 1u64),
    };

    let value: f64 = num_str
        .parse()
        .with_context(|| format!("Invalid size value: {s}"))?;
    Ok((value * multiplier as f64) as u64)
}

/// Run batch alignment on a set of FASTA files using the given BatchAligner.
/// Returns a temp file containing aggregated PAF output.
pub fn run_batch_alignment_generic(
    fasta_files: &[String],
    max_bytes: u64,
    aligner: &dyn BatchAligner,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
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

    // Partition into batches using the aligner's cost model
    let cost_fn = |bp: u64| aligner.estimate_cost(bp);
    let batches = partition_into_batches_with_cost(genomes, max_bytes, &cost_fn);

    if batches.len() == 1 {
        if !config.quiet {
            eprintln!(
                "[batch] All genomes ({}) fit in single batch (est. {} {}), no batching needed",
                format_bytes(total_bp),
                format_bytes(batches[0].estimated_index_bytes),
                aligner.resource_label(),
            );
        }
        // Fall back to normal single-run alignment
        let _ = std::fs::remove_dir_all(&batch_dir);
        return aligner.align_single(fasta_files, tempdir);
    }

    report_batch_plan(&batches, total_bp, aligner.resource_label(), config.quiet);

    // Write batch FASTAs to separate subdirectories
    let mut batch_files: Vec<PathBuf> = Vec::new();
    for (i, batch) in batches.iter().enumerate() {
        let batch_subdir = batch_dir.join(format!("batch_{}", i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("genomes.fa");
        if !config.quiet {
            eprintln!(
                "[batch] Writing batch {} ({} genomes)...",
                i + 1,
                batch.genomes.len()
            );
        }
        write_batch_fasta(batch, &batch_path)?;
        batch_files.push(batch_path);
    }

    // Phase 1: Aligner-specific setup
    aligner.prepare_all(&batch_files, config.quiet)?;

    // Create merged output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;
    let num_batches = batch_files.len();

    // Phase 2: For each target batch, prepare, align all queries, cleanup
    for i in 0..num_batches {
        aligner.prepare_target(i, config.quiet)?;

        for j in 0..num_batches {
            if !config.quiet {
                if i == j {
                    eprintln!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    eprintln!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            let paf_bytes = aligner.align(&batch_files[j], &batch_files[i])?;

            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                eprintln!("[batch]   {} alignments", line_count);
            }
        }

        aligner.cleanup_target(i, config.quiet)?;
    }

    // Phase 3: Final cleanup
    aligner.cleanup_all()?;
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !config.quiet {
        eprintln!(
            "[batch] Completed batch alignment: {} total alignments",
            total_alignments
        );
    }

    Ok(merged_paf)
}

/// Run batch alignment on a set of FASTA files (FastGA-specific, backwards compatible).
/// Returns a temp file containing aggregated PAF output.
pub fn run_batch_alignment(
    fasta_files: &[String],
    max_index_bytes: u64,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    let aligner = FastGABatchAligner::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
        tempdir.map(String::from),
        config.zstd_compress,
        config.zstd_level,
    );
    run_batch_alignment_generic(fasta_files, max_index_bytes, &aligner, config, tempdir)
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
    fn test_estimate_memory_usage() {
        // 100 Mbp -> 0.5 GB + 2 GB = ~2.5 GB
        let mem = estimate_memory_usage(100_000_000);
        assert!(mem > 2_400_000_000 && mem < 2_600_000_000);
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

    #[test]
    fn test_parse_size_string() {
        assert_eq!(parse_size_string("500M").unwrap(), 500_000_000);
        assert_eq!(parse_size_string("2G").unwrap(), 2_000_000_000);
        assert_eq!(parse_size_string("1.5G").unwrap(), 1_500_000_000);
        assert_eq!(parse_size_string("100k").unwrap(), 100_000);
        assert_eq!(parse_size_string("1234").unwrap(), 1234);
        assert!(parse_size_string("").is_err());
        assert!(parse_size_string("abc").is_err());
    }

    #[test]
    fn test_partition_with_custom_cost() {
        // Use wfmash cost model
        let genomes = vec![
            GenomeInfo { prefix: "A#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "B#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "C#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
        ];
        // 100Mbp wfmash cost = 500MB + 2GB = 2.5GB; two genomes = ~4.5GB
        // With 4GB limit, each genome should be in its own batch
        let batches = partition_into_batches_with_cost(genomes, WFMASH_AUTO_BYTES, &estimate_memory_usage);
        assert!(batches.len() >= 2, "Expected at least 2 batches, got {}", batches.len());
    }
}
