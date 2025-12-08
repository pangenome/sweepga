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

/// Partition genomes into batches based on max index size
pub fn partition_into_batches(genomes: Vec<GenomeInfo>, max_index_bytes: u64) -> Vec<GenomeBatch> {
    let mut batches: Vec<GenomeBatch> = Vec::new();
    let mut current_batch = GenomeBatch::new();

    for genome in genomes {
        // Check if this genome alone exceeds the limit
        let single_genome_size = estimate_index_size(genome.total_bp);
        if single_genome_size > max_index_bytes {
            eprintln!(
                "[batch] WARNING: Genome {} ({}) has estimated index {} which exceeds limit {}",
                genome.prefix,
                format_bytes(genome.total_bp),
                format_bytes(single_genome_size),
                format_bytes(max_index_bytes)
            );
            eprintln!("[batch] Including it as a single-genome batch anyway");

            // Flush current batch if non-empty
            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
                current_batch = GenomeBatch::new();
            }

            // Add oversized genome as its own batch
            let mut oversized = GenomeBatch::new();
            oversized.add(genome);
            batches.push(oversized);
            continue;
        }

        // Check if adding to current batch would exceed limit
        if current_batch.would_exceed(&genome, max_index_bytes) {
            // Start new batch
            if !current_batch.genomes.is_empty() {
                batches.push(current_batch);
            }
            current_batch = GenomeBatch::new();
        }

        current_batch.add(genome);
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

/// Report batch statistics
pub fn report_batch_plan(batches: &[GenomeBatch], total_bp: u64, quiet: bool) {
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
    let batches = partition_into_batches(genomes, max_index_bytes);

    if batches.len() == 1 {
        if !config.quiet {
            eprintln!(
                "[batch] All genomes ({}) fit in single batch (est. index {}), no batching needed",
                format_bytes(total_bp),
                format_bytes(batches[0].estimated_index_bytes)
            );
        }
        // Fall back to normal single-run alignment
        return run_single_batch_alignment(fasta_files, config, tempdir);
    }

    report_batch_plan(&batches, total_bp, config.quiet);

    // Write batch FASTAs to separate subdirectories to avoid GDB conflicts
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

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
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
fn run_single_batch_alignment(
    fasta_files: &[String],
    config: &BatchAlignConfig,
    _tempdir: Option<&str>,
) -> Result<tempfile::NamedTempFile> {
    use crate::fastga_integration::FastGAIntegration;

    let fastga = FastGAIntegration::new(
        config.frequency,
        config.threads,
        config.min_alignment_length,
    );

    // Use first file for self-alignment
    let input_path = std::fs::canonicalize(&fasta_files[0])
        .with_context(|| format!("Failed to resolve path: {}", fasta_files[0]))?;

    fastga.align_to_temp_paf(&input_path, &input_path)
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
