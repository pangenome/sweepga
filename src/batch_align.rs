//! Batch alignment for large genome sets
//!
//! When aligning many genomes, the aligner's resource usage can exceed available
//! disk or memory. This module partitions genomes into batches based on total
//! sequence data size, runs the aligner on each batch pair, and aggregates results.
//!
//! ## GIXmake Index Size Limits
//!
//! FastGA's GIXmake indexer has practical size limits for k-mer index creation:
//! - **Safe range**: Batches with ≤40MB of sequence data typically succeed
//! - **Failure threshold**: Batches with ≥48MB often fail silently during index creation
//! - **Recommendation**: For yeast-sized genomes (~12MB each), keep batches to ≤3 genomes
//!
//! When GIXmake failures occur, `run_batch_alignment_with_budget()` automatically
//! reduces batch size and retries. Other batch functions provide error messages
//! suggesting manual --batch-bytes adjustment or switching to wfmash (no size limit).

use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Reasons why a batch alignment attempt might need to restart
#[derive(Debug)]
pub enum RestartReason {
    /// Disk budget exceeded during alignment
    BudgetExceeded,
    /// GIXmake index creation failed due to size limits
    IndexSizeLimitExceeded { batch_size_mb: u64 },
}

impl std::fmt::Display for RestartReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RestartReason::BudgetExceeded => write!(f, "disk budget exceeded"),
            RestartReason::IndexSizeLimitExceeded { batch_size_mb } => {
                write!(f, "GIXmake index size limit exceeded ({}MB)", batch_size_mb)
            }
        }
    }
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
        frequency: usize,
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
            log::info!("[batch] Creating GDB files for all batches...");
        }
        let mut gdb_bases = self.gdb_bases.borrow_mut();
        gdb_bases.clear();
        for (i, batch_file) in batch_files.iter().enumerate() {
            if !quiet {
                log::info!("[batch]   GDB for batch {}...", i + 1);
            }
            let gdb_base = self.fastga.create_gdb_only(batch_file)?;
            gdb_bases.push(gdb_base);
        }
        Ok(())
    }

    fn prepare_target(&self, batch_idx: usize, quiet: bool) -> Result<()> {
        let gdb_bases = self.gdb_bases.borrow();
        if !quiet {
            log::info!("[batch] Building index for batch {}...", batch_idx + 1);
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
            log::info!("[batch] Cleaning up index for batch {}...", batch_idx + 1);
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
    /// Mapping density (wfmash `-x`). Forwarded verbatim to every per-batch
    /// `WfmashIntegration`. `None` (or `Some(1.0)`) keeps all mappings.
    sparsify: Option<f64>,
}

impl WfmashBatchAligner {
    pub fn new(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
        sparsify: Option<f64>,
    ) -> Self {
        Self {
            num_threads,
            min_alignment_length,
            map_pct_identity,
            temp_dir,
            sparsify,
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
            None, // segment_length: adapt from avg_seq
            avg_seq,
            self.sparsify,
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
        }
    }

    pub fn add(&mut self, genome: GenomeInfo) {
        self.total_bp += genome.total_bp;
        self.genomes.push(genome);
    }
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

/// Partition genomes into batches so that each batch's total basepairs
/// does not exceed `max_bp`.  This is the user-facing limit: when
/// someone says `--batch-bytes 50m` they mean 50 MB of sequence data,
/// not some inflated cost estimate.
pub fn partition_into_batches_by_bp(genomes: Vec<GenomeInfo>, max_bp: u64) -> Vec<GenomeBatch> {
    let mut batches: Vec<GenomeBatch> = Vec::new();
    let mut current_batch = GenomeBatch::new();

    for genome in genomes {
        // A single genome bigger than the limit gets its own batch
        if genome.total_bp > max_bp {
            log::warn!(
                "[batch] Genome {} ({}) exceeds batch limit {}; including it as a single-genome batch anyway",
                genome.prefix,
                format_bytes(genome.total_bp),
                format_bytes(max_bp),
            );

            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
                current_batch = GenomeBatch::new();
            }

            let mut oversized = GenomeBatch::new();
            oversized.add(genome);
            batches.push(oversized);
            continue;
        }

        // Would adding this genome push us over the limit?
        if current_batch.total_bp + genome.total_bp > max_bp {
            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
            }
            current_batch = GenomeBatch::new();
        }

        current_batch.add(genome);
    }

    if !current_batch.genomes.is_empty() {
        batches.push(current_batch);
    }

    batches
}

/// Partition genomes into batches of at most `max_count` genomes each.
/// This is the `--batch-size` code path: the user directly controls how many
/// genomes go into each batch.
pub fn partition_into_batches_by_count(genomes: Vec<GenomeInfo>, max_count: usize) -> Vec<GenomeBatch> {
    let mut batches: Vec<GenomeBatch> = Vec::new();
    for chunk in genomes.chunks(max_count) {
        let mut batch = GenomeBatch::new();
        for g in chunk {
            batch.add(g.clone());
        }
        batches.push(batch);
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

/// Report batch statistics.
pub fn report_batch_plan(batches: &[GenomeBatch], total_bp: u64, quiet: bool) {
    if quiet {
        return;
    }

    let total_genomes: usize = batches.iter().map(|b| b.genomes.len()).sum();
    let total_pairs = batches.len() * batches.len();
    let self_pairs = batches.len();
    let cross_pairs = total_pairs - self_pairs;

    log::info!(
        "[batch] Partitioned {} genomes ({}) into {} batches:",
        total_genomes,
        format_bytes(total_bp),
        batches.len(),
    );

    for (i, batch) in batches.iter().enumerate() {
        log::info!(
            "[batch]   Batch {}: {} genomes, {}",
            i + 1,
            batch.genomes.len(),
            format_bytes(batch.total_bp),
        );
    }

    log::info!(
        "[batch] Will run {} batch alignments ({} self + {} cross-batch)",
        total_pairs, self_pairs, cross_pairs
    );
}

// ============================================================================
// Disk budget constants and computation
// ============================================================================

/// Empirical multiplier: GDB + BPS files are ~2× the input basepairs.
/// Calibrated on yeast genomes (~12 MB per genome).
const GDB_FACTOR: f64 = 2.0;

/// Empirical multiplier: each ktab file is ~1× the sequence size per thread.
/// Calibrated on yeast genomes with 8 threads.
const KTAB_PER_THREAD: f64 = 1.0;

/// Compute the maximum basepairs per batch given a disk budget.
///
/// The disk model is:
///   peak_disk = batch_fastas + all_gdbs + one_target_index + paf_reserve
///
/// The only term we control via batch size is `one_target_index`.
///
/// Returns `None` if the budget is too small for even a single genome's index.
pub fn compute_batch_bp_from_budget(
    total_bp: u64,
    genome_sizes: &[u64],
    n_threads: usize,
    zstd: bool,
    disk_budget: u64,
) -> Option<u64> {
    let zstd_factor: f64 = if zstd { 0.5 } else { 1.0 };
    let paf_reserve = total_bp / 10;

    let fixed_overhead = (total_bp as f64 * (1.0 + GDB_FACTOR)) as u64 + paf_reserve;
    let index_factor = n_threads as f64 * KTAB_PER_THREAD * zstd_factor;

    let largest_genome = genome_sizes.iter().copied().max().unwrap_or(0);

    // Check if we can fit at least the largest genome's index
    let min_index = (largest_genome as f64 * index_factor) as u64;
    if disk_budget < fixed_overhead.saturating_add(min_index) {
        return None;
    }

    let available_for_index = disk_budget.saturating_sub(fixed_overhead);
    let max_batch_bp = (available_for_index as f64 / index_factor) as u64;

    // Clamp to at least the largest single genome (can't split a genome)
    Some(max_batch_bp.max(largest_genome))
}

/// Estimate peak disk usage for a given total basepairs and thread count.
pub fn estimate_peak_disk(
    total_bp: u64,
    _n_genomes: usize,
    batch_bp: Option<u64>,
    n_threads: usize,
    zstd: bool,
) -> u64 {
    let zstd_factor: f64 = if zstd { 0.5 } else { 1.0 };
    let index_factor = n_threads as f64 * KTAB_PER_THREAD * zstd_factor;
    let paf_reserve = total_bp / 10;
    let fixed = (total_bp as f64 * (1.0 + GDB_FACTOR)) as u64 + paf_reserve;
    let index_bp = batch_bp.unwrap_or(total_bp);
    let index_cost = (index_bp as f64 * index_factor) as u64;
    fixed + index_cost
}

/// Resolve effective batch bytes from `--max-disk` and `--batch-bytes` flags.
///
/// Returns `Some(max_bp)` if batching should be used, `None` if the whole
/// input fits within the budget (or if neither flag was set).
pub fn resolve_batch_bytes_from_sizes(
    max_disk: Option<u64>,
    batch_bytes: Option<u64>,
    genome_sizes: &[u64],
    n_genomes: usize,
    n_threads: usize,
    zstd: bool,
    quiet: bool,
) -> Result<Option<u64>> {
    match (max_disk, batch_bytes) {
        (None, None) => Ok(None),
        (None, Some(bb)) => Ok(Some(bb)),
        (Some(_), Some(_)) => {
            anyhow::bail!("--max-disk conflicts with --batch-bytes; pass only one.")
        }
        (Some(budget), None) => {
            if genome_sizes.is_empty() {
                anyhow::bail!("No genomes found in input");
            }

            let total_bp: u64 = genome_sizes.iter().sum();
            let max_batch_bp = compute_batch_bp_from_budget(
                total_bp,
                genome_sizes,
                n_threads,
                zstd,
                budget,
            );

            match max_batch_bp {
                None => {
                    let largest = genome_sizes.iter().copied().max().unwrap_or(0);
                    let zstd_factor: f64 = if zstd { 0.5 } else { 1.0 };
                    let min_budget = (total_bp as f64 * (1.0 + GDB_FACTOR)) as u64
                        + total_bp / 10
                        + (largest as f64 * n_threads as f64 * KTAB_PER_THREAD * zstd_factor) as u64;
                    anyhow::bail!(
                        "Disk budget {} is too small for input data.\n\
                         Minimum budget: ~{}.\n\
                         Suggestions: increase --max-disk to at least {}, or use --zstd to halve index size.",
                        format_bytes(budget),
                        format_bytes(min_budget),
                        format_bytes(min_budget),
                    );
                }
                Some(bp) => {
                    let estimated_peak =
                        estimate_peak_disk(total_bp, n_genomes, Some(bp), n_threads, zstd);
                    let n_batches = if bp >= total_bp {
                        1
                    } else {
                        (total_bp + bp - 1) / bp
                    };

                    if !quiet {
                        log::info!("[budget] Disk budget: {}", format_bytes(budget));
                        log::info!(
                            "[budget] Estimated peak: {} ({:.0}% of budget)",
                            format_bytes(estimated_peak),
                            estimated_peak as f64 / budget as f64 * 100.0,
                        );
                        log::info!(
                            "[budget] Batch size: {} (~{} batches)",
                            format_bytes(bp),
                            n_batches,
                        );
                    }

                    // Pre-flight: check available disk space
                    if let Ok(available) = crate::disk_usage::available_disk_bytes("/tmp") {
                        if !quiet {
                            log::info!(
                                "[budget] Pre-flight: {} available on /tmp ({})",
                                format_bytes(available),
                                if available >= estimated_peak { "OK" } else { "WARNING: may be tight" },
                            );
                        }
                        if available < estimated_peak {
                            log::warn!(
                                "[budget] Available disk ({}) is less than estimated peak ({})",
                                format_bytes(available),
                                format_bytes(estimated_peak),
                            );
                        }
                    }

                    if bp >= total_bp {
                        if !quiet {
                            log::info!("[budget] All genomes fit within budget, no batching needed");
                        }
                        return Ok(None);
                    }

                    Ok(Some(bp))
                }
            }
        }
    }
}

/// Resolve effective batch bytes from `--max-disk` / `--batch-bytes` for FASTA input.
pub fn resolve_batch_bytes(
    max_disk: Option<u64>,
    batch_bytes: Option<u64>,
    fasta_files: &[String],
    n_threads: usize,
    zstd: bool,
    quiet: bool,
) -> Result<Option<u64>> {
    if max_disk.is_none() {
        // Fast path: nothing to compute if no budget was set.
        return resolve_batch_bytes_from_sizes(
            max_disk, batch_bytes, &[], 0, n_threads, zstd, quiet,
        );
    }
    let genomes = parse_genome_sizes(fasta_files)?;
    let genome_sizes: Vec<u64> = genomes.iter().map(|g| g.total_bp).collect();
    resolve_batch_bytes_from_sizes(
        max_disk,
        batch_bytes,
        &genome_sizes,
        genomes.len(),
        n_threads,
        zstd,
        quiet,
    )
}

/// Configuration for batch alignment. The aligner's own parameters
/// (threads, k-mer frequency, block length, zstd) live on the
/// `BatchAligner` instance passed alongside this config; this struct
/// only carries the orchestration-level knobs that the batch runner
/// itself consumes.
pub struct BatchAlignConfig {
    pub keep_self: bool,
    pub quiet: bool,
}

/// Maximum number of adaptive restarts before giving up.
const MAX_RESTARTS: usize = 5;

/// Budget threshold: trigger abort when tracked usage exceeds this fraction.
const BUDGET_THRESHOLD: f64 = 0.90;

/// Run batch alignment with disk budget enforcement and adaptive restart.
///
/// Monitors disk usage after each `prepare_target()` call. If the tracked usage
/// exceeds 90% of `disk_budget`, the current attempt is aborted, `max_batch_bp`
/// is halved, and the entire alignment restarts from scratch with smaller batches.
///
/// This is the v1 "simple restart" strategy: completed results from aborted
/// attempts are discarded to guarantee correctness without complex cross-batch
/// bookkeeping. A future v2 could preserve completed target batches.
pub fn run_batch_alignment_with_budget(
    fasta_files: &[String],
    disk_budget: u64,
    initial_batch_bp: u64,
    aligner: &dyn BatchAligner,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    // Parse genome sizes once (shared across restarts)
    if !config.quiet {
        log::info!("[batch] Scanning input files for genome sizes...");
    }
    let genomes = parse_genome_sizes(fasta_files)?;
    let total_bp: u64 = genomes.iter().map(|g| g.total_bp).sum();
    let genome_bp: Vec<u64> = genomes.iter().map(|g| g.total_bp).collect();

    if genomes.is_empty() {
        anyhow::bail!("No genomes found in input files");
    }

    let largest_genome = genome_bp.iter().copied().max().unwrap_or(0);
    let mut max_batch_bp = initial_batch_bp;
    let mut restarts = 0u32;
    let mut restart_reason = RestartReason::BudgetExceeded; // Default, updated on GIXmake failure

    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    loop {
        // Reset disk usage tracking for each attempt
        crate::disk_usage::reset();

        let batches = partition_into_batches_by_bp(genomes.clone(), max_batch_bp);

        if !config.quiet {
            log::info!(
                "[budget] Batch size: {} ({} batches)",
                format_bytes(max_batch_bp),
                batches.len(),
            );
        }

        // Single batch → no budget concerns from batching, use optimized path
        if batches.len() == 1 {
            if !config.quiet {
                log::info!(
                    "[batch] All genomes ({}) fit in single batch, no batching needed",
                    format_bytes(total_bp),
                );
            }
            return aligner.align_single(fasta_files, tempdir);
        }

        report_batch_plan(&batches, total_bp, config.quiet);

        // Create temp directory for this attempt
        let batch_dir = temp_base.join(format!("sweepga_batch_{}", std::process::id()));
        std::fs::create_dir_all(&batch_dir)?;

        // Write batch FASTAs
        let mut batch_files: Vec<PathBuf> = Vec::new();
        for (i, batch) in batches.iter().enumerate() {
            let batch_subdir = batch_dir.join(format!("batch_{}", i));
            std::fs::create_dir_all(&batch_subdir)?;
            let batch_path = batch_subdir.join("genomes.fa");
            if !config.quiet {
                log::info!(
                    "[batch] Writing batch {} ({} genomes)...",
                    i + 1,
                    batch.genomes.len(),
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

        let mut total_alignments = 0usize;
        let num_batches = batch_files.len();
        let mut budget_exceeded = false;

        // Phase 2: For each target batch, prepare, align all queries, cleanup
        for i in 0..num_batches {
            // Try to prepare target batch - this may fail due to GIXmake size limits
            match aligner.prepare_target(i, config.quiet) {
                Ok(()) => {
                    // Success, continue with batch processing
                },
                Err(e) => {
                    // Check if this is a GIXmake size limit error that should trigger restart
                    if let Some(error_chain) = e.chain().find_map(|err| {
                        err.downcast_ref::<crate::fastga_integration::IndexCreationError>()
                    }) {
                        if let crate::fastga_integration::IndexCreationError::SizeLimitExceeded { batch_size_mb, .. } = error_chain {
                            log::info!(
                                "[gixmake] Index creation failed for batch {} ({}MB) - batch too large for GIXmake",
                                i + 1, batch_size_mb
                            );

                            // Set restart reason and clean up
                            restart_reason = RestartReason::IndexSizeLimitExceeded {
                                batch_size_mb: *batch_size_mb
                            };
                            aligner.cleanup_all()?;
                            budget_exceeded = true; // Reuse existing restart mechanism
                            break;
                        }
                    }

                    // Not a GIXmake size limit error - propagate the original error
                    return Err(e);
                }
            }

            // Budget check after prepare_target — the index is the dominant cost
            let (exceeded, current, _) =
                crate::disk_usage::check_budget(disk_budget, BUDGET_THRESHOLD);
            if exceeded {
                log::info!(
                    "[budget] Disk usage {} exceeds {:.0}% of budget {}",
                    crate::disk_usage::format_bytes(current),
                    BUDGET_THRESHOLD * 100.0,
                    format_bytes(disk_budget),
                );
                aligner.cleanup_target(i, config.quiet)?;
                budget_exceeded = true;
                break;
            }

            if !config.quiet {
                log::info!(
                    "[batch] Target batch {} [disk: {} / {}]",
                    i + 1,
                    crate::disk_usage::format_bytes(current),
                    format_bytes(disk_budget),
                );
            }

            for j in 0..num_batches {
                if !config.quiet {
                    if i == j {
                        log::info!("[batch] Aligning batch {} to itself...", i + 1);
                    } else {
                        log::info!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                    }
                }

                let paf_bytes = aligner.align(&batch_files[j], &batch_files[i])?;
                let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
                total_alignments += line_count;
                merged_output.write_all(&paf_bytes)?;

                if !config.quiet {
                    log::info!("[batch]   {} alignments", line_count);
                }
            }

            aligner.cleanup_target(i, config.quiet)?;

            if !config.quiet {
                let current = crate::disk_usage::current_usage();
                log::info!(
                    "[batch] Completed target batch {} [disk: {} / {}]",
                    i + 1,
                    crate::disk_usage::format_bytes(current),
                    format_bytes(disk_budget),
                );
            }
        }

        // Phase 3: Final cleanup
        aligner.cleanup_all()?;
        let _ = std::fs::remove_dir_all(&batch_dir);

        if !budget_exceeded {
            // Success — report and verify
            if !config.quiet {
                log::info!(
                    "[batch] Completed batch alignment: {} total alignments",
                    total_alignments,
                );
                let peak = crate::disk_usage::peak_usage();
                log::info!(
                    "[budget] Peak disk usage: {} ({:.0}% of budget)",
                    crate::disk_usage::format_bytes(peak),
                    peak as f64 / disk_budget as f64 * 100.0,
                );
            }

            // Verify completeness
            drop(merged_output);
            let genome_prefixes: Vec<String> = batches
                .iter()
                .flat_map(|b| b.genomes.iter().map(|g| g.prefix.clone()))
                .collect();
            let verification = verify_batch_completeness(
                merged_paf.path(),
                &genome_prefixes,
                !config.keep_self,
            )?;
            if !config.quiet {
                log::info!(
                    "[batch] Verification: {}/{} genome pairs present",
                    verification.found_pairs, verification.expected_pairs,
                );
            }

            return Ok(merged_paf);
        }

        // Budget was exceeded — apply halving backoff
        restarts += 1;
        if restarts as usize > MAX_RESTARTS {
            match restart_reason {
                RestartReason::BudgetExceeded => {
                    anyhow::bail!(
                        "Exceeded max restarts ({}) while trying to fit within disk budget {}. \
                         Use --zstd to halve index size or increase --max-disk.",
                        MAX_RESTARTS,
                        format_bytes(disk_budget),
                    );
                }
                RestartReason::IndexSizeLimitExceeded { batch_size_mb } => {
                    anyhow::bail!(
                        "Exceeded max restarts ({}) due to GIXmake index size limits (last failed batch: {}MB). \
                         Try reducing --batch-bytes further or use a different aligner (wfmash) that doesn't have this limit.",
                        MAX_RESTARTS,
                        batch_size_mb
                    );
                }
            }
        }

        let old_batch_bp = max_batch_bp;
        max_batch_bp /= 2;

        // Can't go below the largest single genome
        if max_batch_bp < largest_genome {
            max_batch_bp = largest_genome;
            if old_batch_bp == largest_genome {
                anyhow::bail!(
                    "Cannot reduce batch size below largest genome ({}). \
                     Use --zstd or increase --max-disk.",
                    format_bytes(largest_genome),
                );
            }
        }

        match restart_reason {
            RestartReason::BudgetExceeded => {
                log::info!(
                    "[budget] Restart {}/{}: reducing batch from {} to {} (disk budget exceeded)",
                    restarts,
                    MAX_RESTARTS,
                    format_bytes(old_batch_bp),
                    format_bytes(max_batch_bp),
                );
            }
            RestartReason::IndexSizeLimitExceeded { batch_size_mb } => {
                log::info!(
                    "[gixmake] Restart {}/{}: reducing batch from {} to {} (GIXmake index size limit, {}MB batch failed)",
                    restarts,
                    MAX_RESTARTS,
                    format_bytes(old_batch_bp),
                    format_bytes(max_batch_bp),
                    batch_size_mb
                );
            }
        }

        // Reset restart reason for next iteration
        restart_reason = RestartReason::BudgetExceeded;
    }
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
        log::info!("[batch] Scanning input files for genome sizes...");
    }
    let genomes = parse_genome_sizes(fasta_files)?;
    let total_bp: u64 = genomes.iter().map(|g| g.total_bp).sum();

    if genomes.is_empty() {
        anyhow::bail!("No genomes found in input files");
    }

    // Partition into batches by total basepairs
    let batches = partition_into_batches_by_bp(genomes, max_bytes);

    if batches.len() == 1 {
        if !config.quiet {
            log::info!(
                "[batch] All genomes ({}) fit in single batch, no batching needed",
                format_bytes(total_bp),
            );
        }
        // Fall back to normal single-run alignment
        let _ = std::fs::remove_dir_all(&batch_dir);
        return aligner.align_single(fasta_files, tempdir);
    }

    report_batch_plan(&batches, total_bp, config.quiet);

    // Write batch FASTAs to separate subdirectories
    let mut batch_files: Vec<PathBuf> = Vec::new();
    for (i, batch) in batches.iter().enumerate() {
        let batch_subdir = batch_dir.join(format!("batch_{}", i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("genomes.fa");
        if !config.quiet {
            log::info!(
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
        // Try to prepare target batch - provide helpful error for GIXmake size limit failures
        if let Err(e) = aligner.prepare_target(i, config.quiet) {
            // Check if this is a GIXmake size limit error
            if let Some(error_chain) = e.chain().find_map(|err| {
                err.downcast_ref::<crate::fastga_integration::IndexCreationError>()
            }) {
                if let crate::fastga_integration::IndexCreationError::SizeLimitExceeded { batch_size_mb, suggested_limit, .. } = error_chain {
                    return Err(anyhow::anyhow!(
                        "GIXmake index creation failed for batch {} ({}MB). \
                         This batch size exceeds FastGA's index limit. \
                         Try --batch-bytes {}M or use run_batch_alignment_with_budget() for automatic retry.",
                        i + 1, batch_size_mb, suggested_limit
                    ));
                }
            }
            // Not a GIXmake error - propagate original error
            return Err(e);
        }

        for j in 0..num_batches {
            if !config.quiet {
                if i == j {
                    log::info!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    log::info!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            let paf_bytes = aligner.align(&batch_files[j], &batch_files[i])?;

            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                log::info!("[batch]   {} alignments", line_count);
            }
        }

        aligner.cleanup_target(i, config.quiet)?;
    }

    // Phase 3: Final cleanup
    aligner.cleanup_all()?;

    if !config.quiet {
        log::info!(
            "[batch] Completed batch alignment: {} total alignments",
            total_alignments
        );
    }

    // Phase 4: Verify completeness — check all expected genome pairs have alignments
    drop(merged_output); // flush before reading
    let genome_prefixes: Vec<String> = batches
        .iter()
        .flat_map(|b| b.genomes.iter().map(|g| g.prefix.clone()))
        .collect();
    let verification = verify_batch_completeness(
        merged_paf.path(),
        &genome_prefixes,
        !config.keep_self,
    )?;
    if !config.quiet {
        log::info!(
            "[batch] Verification: {}/{} genome pairs present",
            verification.found_pairs, verification.expected_pairs
        );
    }

    let _ = std::fs::remove_dir_all(&batch_dir);

    Ok(merged_paf)
}

/// Run batch alignment partitioned by genome count (`--batch-size`).
///
/// Each batch contains at most `max_count` genomes. The alignment loop
/// follows the same prepare/align/cleanup protocol as the byte-based path.
pub fn run_batch_alignment_by_count(
    fasta_files: &[String],
    max_count: usize,
    aligner: &dyn BatchAligner,
    config: &BatchAlignConfig,
    tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    let temp_base = if let Some(dir) = tempdir {
        PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        PathBuf::from(tmpdir)
    } else {
        PathBuf::from("/tmp")
    };

    let batch_dir = temp_base.join(format!("sweepga_batch_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    if !config.quiet {
        log::info!("[batch] Scanning input files for genome sizes...");
    }
    let genomes = parse_genome_sizes(fasta_files)?;
    let total_bp: u64 = genomes.iter().map(|g| g.total_bp).sum();

    if genomes.is_empty() {
        anyhow::bail!("No genomes found in input files");
    }

    let batches = partition_into_batches_by_count(genomes, max_count);

    if batches.len() == 1 {
        if !config.quiet {
            log::info!(
                "[batch] All genomes ({}) fit in single batch, no batching needed",
                format_bytes(total_bp),
            );
        }
        let _ = std::fs::remove_dir_all(&batch_dir);
        return aligner.align_single(fasta_files, tempdir);
    }

    report_batch_plan(&batches, total_bp, config.quiet);

    let mut batch_files: Vec<PathBuf> = Vec::new();
    for (i, batch) in batches.iter().enumerate() {
        let batch_subdir = batch_dir.join(format!("batch_{}", i));
        std::fs::create_dir_all(&batch_subdir)?;
        let batch_path = batch_subdir.join("genomes.fa");
        if !config.quiet {
            log::info!(
                "[batch] Writing batch {} ({} genomes)...",
                i + 1,
                batch.genomes.len()
            );
        }
        write_batch_fasta(batch, &batch_path)?;
        batch_files.push(batch_path);
    }

    aligner.prepare_all(&batch_files, config.quiet)?;

    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    let mut total_alignments = 0;
    let num_batches = batch_files.len();

    for i in 0..num_batches {
        // Try to prepare target batch - provide helpful error for GIXmake size limit failures
        if let Err(e) = aligner.prepare_target(i, config.quiet) {
            // Check if this is a GIXmake size limit error
            if let Some(error_chain) = e.chain().find_map(|err| {
                err.downcast_ref::<crate::fastga_integration::IndexCreationError>()
            }) {
                if let crate::fastga_integration::IndexCreationError::SizeLimitExceeded { batch_size_mb, suggested_limit, .. } = error_chain {
                    return Err(anyhow::anyhow!(
                        "GIXmake index creation failed for batch {} ({}MB). \
                         This batch size exceeds FastGA's index limit. \
                         Try --batch-bytes {}M or use the disk budget mode for automatic retry.",
                        i + 1, batch_size_mb, suggested_limit
                    ));
                }
            }
            // Not a GIXmake error - propagate original error
            return Err(e);
        }

        for j in 0..num_batches {
            if !config.quiet {
                if i == j {
                    log::info!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    log::info!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            let paf_bytes = aligner.align(&batch_files[j], &batch_files[i])?;
            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;
            merged_output.write_all(&paf_bytes)?;

            if !config.quiet {
                log::info!("[batch]   {} alignments", line_count);
            }
        }

        aligner.cleanup_target(i, config.quiet)?;
    }

    aligner.cleanup_all()?;

    if !config.quiet {
        log::info!(
            "[batch] Completed batch alignment: {} total alignments",
            total_alignments
        );
    }

    drop(merged_output);
    let genome_prefixes: Vec<String> = batches
        .iter()
        .flat_map(|b| b.genomes.iter().map(|g| g.prefix.clone()))
        .collect();
    let verification = verify_batch_completeness(
        merged_paf.path(),
        &genome_prefixes,
        !config.keep_self,
    )?;
    if !config.quiet {
        log::info!(
            "[batch] Verification: {}/{} genome pairs present",
            verification.found_pairs, verification.expected_pairs
        );
    }

    let _ = std::fs::remove_dir_all(&batch_dir);

    Ok(merged_paf)
}

/// Result of batch completeness verification.
///
/// The `missing` pair list is surfaced to the user via `log::warn!`
/// inside `verify_batch_completeness`; callers here only need the
/// `found_pairs` / `expected_pairs` counters.
#[derive(Debug, Clone)]
pub struct BatchVerification {
    /// Genome pairs found in the output.
    pub found_pairs: usize,
    /// Expected genome pairs.
    pub expected_pairs: usize,
}

/// Verify that a PAF file contains alignments for all expected genome pairs.
///
/// Scans the PAF output for query/target genome prefixes and checks that every
/// expected pair is represented. Missing pairs are logged as warnings.
///
/// `exclude_self` controls whether self-alignments (same genome prefix) are
/// expected. When `true`, pairs like A→A are not required.
pub fn verify_batch_completeness(
    paf_path: &Path,
    expected_genomes: &[String],
    exclude_self: bool,
) -> Result<BatchVerification> {
    // Build expected pairs
    let mut expected: HashSet<(String, String)> = HashSet::new();
    for q in expected_genomes {
        for t in expected_genomes {
            if exclude_self && q == t {
                continue;
            }
            expected.insert((q.clone(), t.clone()));
        }
    }

    // Scan PAF for found pairs
    let mut found: HashSet<(String, String)> = HashSet::new();

    let file = File::open(paf_path)
        .with_context(|| format!("Failed to open PAF for verification: {}", paf_path.display()))?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let query_prefix = extract_pansn_prefix(fields[0]);
        let target_prefix = extract_pansn_prefix(fields[5]);
        found.insert((query_prefix, target_prefix));
    }

    let missing: Vec<(String, String)> = expected
        .difference(&found)
        .cloned()
        .collect::<Vec<_>>();

    if !missing.is_empty() {
        log::warn!(
            "[batch] {}/{} genome pairs missing from output",
            missing.len(),
            expected.len()
        );
        for (q, t) in &missing {
            log::warn!("[batch]   missing: {} -> {}", q, t);
        }
    }

    Ok(BatchVerification {
        found_pairs: expected.intersection(&found).count(),
        expected_pairs: expected.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_verify_completeness_all_pairs_present() {
        // Create a minimal PAF file with all expected genome pairs
        let mut paf = tempfile::NamedTempFile::with_suffix(".paf").unwrap();
        // PAF: qname qlen qstart qend strand tname tlen tstart tend matches block_len mapq
        let pairs = [
            ("A#1#chr1", "B#1#chr1"),
            ("A#1#chr1", "C#1#chr1"),
            ("B#1#chr1", "A#1#chr1"),
            ("B#1#chr1", "C#1#chr1"),
            ("C#1#chr1", "A#1#chr1"),
            ("C#1#chr1", "B#1#chr1"),
        ];
        for (q, t) in &pairs {
            writeln!(
                paf,
                "{}\t1000\t0\t1000\t+\t{}\t1000\t0\t1000\t900\t1000\t60",
                q, t
            ).unwrap();
        }
        paf.flush().unwrap();

        let genomes = vec!["A#1#".to_string(), "B#1#".to_string(), "C#1#".to_string()];
        let result = verify_batch_completeness(paf.path(), &genomes, true).unwrap();
        assert_eq!(result.expected_pairs, 6); // 3 genomes, exclude self = 3*2
        assert_eq!(result.found_pairs, 6);
        assert_eq!(result.found_pairs, result.expected_pairs);
    }

    #[test]
    fn test_verify_completeness_missing_pair() {
        let mut paf = tempfile::NamedTempFile::with_suffix(".paf").unwrap();
        // Only A->B and B->A, missing C pairs
        for (q, t) in &[("A#1#chr1", "B#1#chr1"), ("B#1#chr1", "A#1#chr1")] {
            writeln!(
                paf,
                "{}\t1000\t0\t1000\t+\t{}\t1000\t0\t1000\t900\t1000\t60",
                q, t
            ).unwrap();
        }
        paf.flush().unwrap();

        let genomes = vec!["A#1#".to_string(), "B#1#".to_string(), "C#1#".to_string()];
        let result = verify_batch_completeness(paf.path(), &genomes, true).unwrap();
        assert_eq!(result.expected_pairs, 6);
        assert_eq!(result.found_pairs, 2);
        assert_eq!(result.expected_pairs - result.found_pairs, 4); // A->C, C->A, B->C, C->B
    }

    #[test]
    fn test_verify_completeness_with_self() {
        let mut paf = tempfile::NamedTempFile::with_suffix(".paf").unwrap();
        // Include self-alignments
        for (q, t) in &[
            ("A#1#chr1", "A#1#chr1"),
            ("A#1#chr1", "B#1#chr1"),
            ("B#1#chr1", "A#1#chr1"),
            ("B#1#chr1", "B#1#chr1"),
        ] {
            writeln!(
                paf,
                "{}\t1000\t0\t1000\t+\t{}\t1000\t0\t1000\t900\t1000\t60",
                q, t
            ).unwrap();
        }
        paf.flush().unwrap();

        let genomes = vec!["A#1#".to_string(), "B#1#".to_string()];
        // exclude_self = false means we expect self-pairs too
        let result = verify_batch_completeness(paf.path(), &genomes, false).unwrap();
        assert_eq!(result.expected_pairs, 4); // 2*2 = 4 (A->A, A->B, B->A, B->B)
        assert_eq!(result.found_pairs, 4);
        assert_eq!(result.found_pairs, result.expected_pairs);
    }

    #[test]
    fn test_partition_by_bp() {
        let genomes = vec![
            GenomeInfo { prefix: "A#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "B#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "C#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
        ];
        // 100Mbp each, 150Mbp limit -> 3 batches (one per genome, can't fit two)
        let batches = partition_into_batches_by_bp(genomes.clone(), 150_000_000);
        assert_eq!(batches.len(), 3, "Expected 3 batches, got {}", batches.len());

        // 200Mbp limit -> 2 batches (AB, C)
        let batches = partition_into_batches_by_bp(genomes.clone(), 200_000_000);
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].genomes.len(), 2);
        assert_eq!(batches[1].genomes.len(), 1);

        // 300Mbp limit -> 1 batch (all fit)
        let batches = partition_into_batches_by_bp(genomes, 300_000_000);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].genomes.len(), 3);
    }

    #[test]
    fn test_partition_by_count() {
        let genomes = vec![
            GenomeInfo { prefix: "A#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "B#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "C#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "D#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
            GenomeInfo { prefix: "E#1#".to_string(), total_bp: 100_000_000, source_file: PathBuf::from("a.fa") },
        ];

        // 2 genomes per batch -> 3 batches (AB, CD, E)
        let batches = partition_into_batches_by_count(genomes.clone(), 2);
        assert_eq!(batches.len(), 3);
        assert_eq!(batches[0].genomes.len(), 2);
        assert_eq!(batches[1].genomes.len(), 2);
        assert_eq!(batches[2].genomes.len(), 1);

        // 5 genomes per batch -> 1 batch
        let batches = partition_into_batches_by_count(genomes.clone(), 5);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].genomes.len(), 5);

        // 1 genome per batch -> 5 batches
        let batches = partition_into_batches_by_count(genomes, 1);
        assert_eq!(batches.len(), 5);
        for b in &batches {
            assert_eq!(b.genomes.len(), 1);
        }
    }

    #[test]
    fn test_batch_size_from_budget() {
        // 8 yeast genomes, ~12M bp each, 8 threads, no zstd
        let total_bp = 96_000_000u64;
        let genome_sizes = vec![12_000_000u64; 8];
        let n_threads = 8;

        // With 2 GB budget → should need no batching (all fits)
        let bp = compute_batch_bp_from_budget(total_bp, &genome_sizes, n_threads, false, 2_000_000_000);
        assert!(bp.is_some());
        assert!(bp.unwrap() >= total_bp, "Expected >= total_bp, got {}", bp.unwrap());

        // With 500 MB budget → should produce batches
        let bp = compute_batch_bp_from_budget(total_bp, &genome_sizes, n_threads, false, 500_000_000);
        assert!(bp.is_some());
        let bp_val = bp.unwrap();
        assert!(bp_val < total_bp, "Expected < total_bp, got {}", bp_val);
        assert!(bp_val >= 12_000_000, "Expected >= largest genome, got {}", bp_val);

        // With zstd → more room (halved index), so batch_bp should be larger
        let bp_no_zstd = compute_batch_bp_from_budget(total_bp, &genome_sizes, n_threads, false, 500_000_000).unwrap();
        let bp_zstd = compute_batch_bp_from_budget(total_bp, &genome_sizes, n_threads, true, 500_000_000).unwrap();
        assert!(bp_zstd >= bp_no_zstd, "zstd should allow larger batches");
    }

    #[test]
    fn test_batch_size_budget_too_small() {
        let total_bp = 96_000_000u64;
        let genome_sizes = vec![12_000_000u64; 8];
        // fixed_overhead ~= 96M * 3.0 + 9.6M ≈ 298M
        // min index = 12M * 8 threads = 96M
        // minimum total = ~394M
        // 10M is way too small
        let bp = compute_batch_bp_from_budget(total_bp, &genome_sizes, 8, false, 10_000_000);
        assert!(bp.is_none(), "Should return None for impossibly small budget");
    }

    #[test]
    fn test_backoff_halving() {
        // Simulate the halving backoff logic from run_batch_alignment_with_budget
        let largest_genome = 12_000_000u64;
        let mut max_batch_bp = 96_000_000u64;

        // Each halving should reduce by 2x
        for expected_restart in 1..=MAX_RESTARTS {
            let old = max_batch_bp;
            max_batch_bp /= 2;
            if max_batch_bp < largest_genome {
                max_batch_bp = largest_genome;
            }
            assert!(
                max_batch_bp <= old,
                "Restart {}: batch should not grow ({} > {})",
                expected_restart, max_batch_bp, old,
            );
        }

        // After 5 halvings: 96M → 48M → 24M → 12M (clamped) → 12M → 12M
        assert_eq!(max_batch_bp, largest_genome);
    }

    #[test]
    fn test_check_budget_threshold() {
        crate::disk_usage::reset();

        // Simulate 800 bytes of disk usage against a 1000-byte budget
        crate::disk_usage::add_bytes(800);
        let (exceeded, current, budget) = crate::disk_usage::check_budget(1000, 0.90);
        assert!(!exceeded, "800/1000 = 80% should not exceed 90% threshold");
        assert_eq!(current, 800);
        assert_eq!(budget, 1000);

        // Add more to cross the threshold
        crate::disk_usage::add_bytes(200); // now 1000 total
        let (exceeded, current, _) = crate::disk_usage::check_budget(1000, 0.90);
        assert!(exceeded, "1000/1000 = 100% should exceed 90% threshold");
        assert_eq!(current, 1000);

        crate::disk_usage::reset();
    }

    #[test]
    fn test_check_budget_just_below_threshold() {
        crate::disk_usage::reset();

        // 899 of 1000 = 89.9%, just below 90%
        crate::disk_usage::add_bytes(899);
        let (exceeded, _, _) = crate::disk_usage::check_budget(1000, 0.90);
        assert!(!exceeded, "89.9% should not exceed 90% threshold");

        // 901 of 1000 = 90.1%, just above
        crate::disk_usage::add_bytes(2); // now 901
        let (exceeded, _, _) = crate::disk_usage::check_budget(1000, 0.90);
        assert!(exceeded, "90.1% should exceed 90% threshold");

        crate::disk_usage::reset();
    }

    #[test]
    fn test_estimate_peak_disk() {
        let total_bp = 96_000_000u64;
        // No batching: index covers all genomes
        let peak_no_batch = estimate_peak_disk(total_bp, 8, None, 8, false);
        // With batching: index covers only one batch
        let peak_batched = estimate_peak_disk(total_bp, 8, Some(24_000_000), 8, false);

        assert!(peak_batched < peak_no_batch,
            "Batched peak {} should be less than unbatched peak {}",
            peak_batched, peak_no_batch);

        // Zstd should reduce peak
        let peak_zstd = estimate_peak_disk(total_bp, 8, Some(24_000_000), 8, true);
        assert!(peak_zstd < peak_batched,
            "Zstd peak {} should be less than non-zstd peak {}",
            peak_zstd, peak_batched);
    }
}
