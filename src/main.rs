mod agc;
mod aln_filter;
mod batch_align;
mod binary_paths;
mod compact_mapping;
mod disk_usage;
mod fastga_integration;
mod filter_types;
mod grouped_mappings;
mod knn_graph;
mod mapping;
mod mash;
mod paf;
mod paf_filter;
mod plane_sweep;
mod plane_sweep_core;
mod plane_sweep_exact;
mod plane_sweep_scaffold;
mod sequence_index;
mod tree_sparsify;
mod unified_filter;
mod union_find;

use anyhow::{Context, Result};
use clap::Parser;

use crate::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

/// Timing context for minimap2-style logging
struct TimingContext {
    start_time: Instant,
    start_cpu: f64,
}

impl TimingContext {
    fn new() -> Self {
        Self {
            start_time: Instant::now(),
            start_cpu: Self::cpu_time(),
        }
    }

    /// Get current CPU time (user + system) in seconds
    fn cpu_time() -> f64 {
        unsafe {
            let mut usage: libc::rusage = std::mem::zeroed();
            libc::getrusage(libc::RUSAGE_SELF, &mut usage);
            let user = usage.ru_utime.tv_sec as f64 + usage.ru_utime.tv_usec as f64 / 1_000_000.0;
            let system = usage.ru_stime.tv_sec as f64 + usage.ru_stime.tv_usec as f64 / 1_000_000.0;
            user + system
        }
    }

    /// Calculate elapsed wall time and CPU ratio
    fn stats(&self) -> (f64, f64) {
        let elapsed = self.start_time.elapsed().as_secs_f64();
        let cpu_used = Self::cpu_time() - self.start_cpu;
        let cpu_ratio = if elapsed > 0.0 {
            cpu_used / elapsed
        } else {
            0.0
        };
        (elapsed, cpu_ratio)
    }

    /// Log message with timing in minimap2 format (currently disabled)
    fn log(&self, _phase: &str, _message: &str) {
        // Logging disabled - eprintln statements commented out
    }
}

/// File type detected from content
#[derive(Debug, Clone, Copy, PartialEq)]
enum FileType {
    Fasta,
    Paf,
    Aln, // .1aln binary format
    Agc, // AGC compressed genome archive
}

/// Detect file type by reading first non-empty line or checking file extension
/// Handles .gz files automatically
fn detect_file_type(path: &str) -> Result<FileType> {
    // Check for .agc extension first (AGC archive)
    if path.ends_with(".agc") || path.ends_with(".AGC") {
        return Ok(FileType::Agc);
    }

    // Check for .1aln extension (binary format)
    if path.ends_with(".1aln") {
        return Ok(FileType::Aln);
    }

    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") || path.ends_with(".bgz") {
        // Use bgzf reader for .gz and .bgz files (handles both gzip and bgzip)
        Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Read first non-empty line
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            anyhow::bail!("Empty file: {path}");
        }
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            break;
        }
    }

    let trimmed = line.trim();

    // FASTA starts with >
    if trimmed.starts_with('>') {
        return Ok(FileType::Fasta);
    }

    // PAF is tab-delimited with 12+ fields
    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() >= 12 {
        // Check that numeric fields parse correctly (columns 1, 2, 3, 6, 7, 8, 9, 10, 11)
        if fields[1].parse::<u64>().is_ok()
            && fields[2].parse::<u64>().is_ok()
            && fields[3].parse::<u64>().is_ok()
            && fields[6].parse::<u64>().is_ok()
            && fields[9].parse::<u64>().is_ok()
            && fields[10].parse::<u64>().is_ok()
        {
            return Ok(FileType::Paf);
        }
    }

    anyhow::bail!("Could not detect file type for {path}: not FASTA (starts with >), PAF (12+ tab-delimited fields), or .1aln (binary)");
}

/// Method for calculating ANI from alignments
#[derive(Debug, Clone, Copy, PartialEq)]
enum AniMethod {
    All,                     // Use all alignments
    Orthogonal,              // Use best 1:1 mappings only
    NPercentile(f64, NSort), // Use alignments up to N-percentile (e.g., 50 for N50, 90 for N90)
}

/// Sorting method for N-percentile calculation
#[derive(Debug, Clone, Copy, PartialEq)]
enum NSort {
    Length,   // Sort by alignment length
    Identity, // Sort by identity
    Score,    // Sort by identity * log(length)
}

/// Parse a number that may have metric suffix (k/K=1000, m/M=1e6, g/G=1e9)
fn parse_metric_number(s: &str) -> Result<u64, String> {
    if s.is_empty() {
        return Err("Empty string".to_string());
    }

    let (num_part, suffix) = if s.ends_with(|c: char| c.is_ascii_alphabetic()) {
        let last_char = s.chars().last().unwrap();
        (&s[..s.len() - last_char.len_utf8()], Some(last_char))
    } else {
        (s, None)
    };

    let base: f64 = num_part
        .parse()
        .map_err(|e| format!("Invalid number: {e}"))?;

    let multiplier = match suffix {
        Some('k') | Some('K') => 1000.0,
        Some('m') | Some('M') => 1_000_000.0,
        Some('g') | Some('G') => 1_000_000_000.0,
        Some(c) => {
            return Err(format!(
                "Unknown suffix '{c}'. Use k/K (1000), m/M (1e6), or g/G (1e9)"
            ))
        }
        None => 1.0,
    };

    let result = base * multiplier;

    if result > u64::MAX as f64 {
        return Err(format!("Value {result} too large for u64"));
    }

    Ok(result as u64)
}

/// SweepGA - Fast genome alignment with plane sweep filtering
///
/// This tool wraps genome aligners (FastGA by default) and applies wfmash's filtering algorithms.
/// Can also process existing PAF files from any aligner.
/// File type (FASTA or PAF) is automatically detected from content.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    // ============================================================================
    // Input/Output
    // ============================================================================
    /// PAF (filter) or FASTA (align & filter). Multiple FASTAs for all-pairs alignment
    #[clap(value_name = "FILE", num_args = 0..,
           long_help = "Input files: FASTA (1+) or PAF (1 only), auto-detected\n\
                        \n  \
                        1 FASTA: align to self and filter\n  \
                        2+ FASTA: align all pairs and filter\n  \
                        1 PAF: filter alignments\n  \
                        stdin: auto-detect and process")]
    files: Vec<String>,

    /// Output file path (auto-detects format from extension: .paf or .1aln)
    #[clap(long = "output-file")]
    output_file: Option<String>,

    /// Output PAF format (default, kept for clarity and backwards compatibility)
    #[clap(long = "paf")]
    output_paf: bool,

    /// Output .1aln binary format instead of default PAF
    #[clap(long = "1aln")]
    output_1aln: bool,

    // ============================================================================
    // Alignment (FASTA input only)
    // ============================================================================
    /// Aligner to use for FASTA input
    #[clap(long = "aligner", default_value = "fastga", value_parser = ["fastga"],
           help_heading = "Alignment options")]
    aligner: String,

    /// FastGA k-mer frequency threshold (use k-mers occurring ≤ N times)
    #[clap(short = 'f', long = "frequency", help_heading = "Alignment options")]
    frequency: Option<usize>,

    /// Align all genome pairs separately (slower, uses more memory, but handles many genomes)
    #[clap(long = "all-pairs", help_heading = "Alignment options")]
    all_pairs: bool,

    /// Maximum index size per batch (e.g., "10G", "500M"). Partitions genomes to fit disk limits.
    /// Formula: index ≈ 0.1GB + 12 bytes/bp. Use when scratch space is limited.
    #[clap(long = "batch-bytes", value_parser = parse_metric_number, help_heading = "Alignment options")]
    batch_bytes: Option<u64>,

    /// Compress k-mer index with zstd for ~2x disk savings and faster I/O
    #[clap(long = "zstd", help_heading = "Alignment options")]
    zstd_compress: bool,

    /// Zstd compression level (1-19, higher = smaller files but slower). Default: 3
    #[clap(
        long = "zstd-level",
        default_value = "3",
        help_heading = "Alignment options"
    )]
    zstd_level: u32,

    // ============================================================================
    // Basic Filtering
    // ============================================================================
    /// Minimum block length
    #[clap(short = 'b', long = "block-length", default_value = "0", value_parser = parse_metric_number,
           help_heading = "Basic filtering")]
    block_length: u64,

    /// n:m-best mappings kept in query:target dimensions. 1:1 (orthogonal), use ∞/many for unbounded
    #[clap(
        short = 'n',
        long = "num-mappings",
        default_value = "1:1",
        help_heading = "Basic filtering"
    )]
    num_mappings: String,

    /// Maximum overlap ratio for plane sweep filtering
    #[clap(
        short = 'o',
        long = "overlap",
        default_value = "0.95",
        help_heading = "Basic filtering"
    )]
    overlap: f64,

    /// Scoring function for plane sweep
    #[clap(long = "scoring", default_value = "log-length-ani",
           value_parser = ["ani", "length", "length-ani", "log-length-ani", "matches"],
           help_heading = "Basic filtering")]
    scoring: String,

    /// Minimum identity threshold (0-1 fraction, 1-100%, or "aniN" for Nth percentile)
    #[clap(
        short = 'i',
        long = "min-identity",
        default_value = "0",
        help_heading = "Basic filtering"
    )]
    min_identity: String,

    /// Keep self-mappings (excluded by default)
    #[clap(long = "self", help_heading = "Basic filtering")]
    keep_self: bool,

    /// Disable all filtering
    #[clap(short = 'N', long = "no-filter", help_heading = "Basic filtering")]
    no_filter: bool,

    // ============================================================================
    // Scaffolding and Chaining
    // ============================================================================
    /// Scaffold jump (gap) distance. 0 = disable scaffolding, >0 = enable (accepts k/m/g suffix)
    #[clap(short = 'j', long = "scaffold-jump", default_value = "0", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    scaffold_jump: u64,

    /// Minimum scaffold chain length (accepts k/m/g suffix)
    #[clap(short = 's', long = "scaffold-mass", default_value = "10k", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    scaffold_mass: u64,

    /// Scaffold filter mode: "1:1" (best), "M:N" (M per query, N per target), "many" (unbounded)
    #[clap(
        short = 'm',
        long = "scaffold-filter",
        default_value = "many:many",
        help_heading = "Scaffolding and chaining"
    )]
    scaffold_filter: String,

    /// Scaffold chain overlap threshold for plane sweep filtering
    #[clap(
        short = 'O',
        long = "scaffold-overlap",
        default_value = "0.5",
        help_heading = "Scaffolding and chaining"
    )]
    scaffold_overlap: f64,

    /// Maximum Euclidean distance from scaffold anchor for rescue (0 = no rescue, accepts k/m/g suffix)
    #[clap(short = 'd', long = "scaffold-dist", default_value = "0", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    scaffold_dist: u64,

    /// Minimum scaffold identity threshold (0-1 fraction, 1-100%, "aniN", or defaults to -i)
    #[clap(
        short = 'Y',
        long = "min-scaffold-identity",
        default_value = "0",
        help_heading = "Scaffolding and chaining"
    )]
    min_scaffold_identity: String,

    /// Output scaffold chains only (for debugging)
    #[clap(long = "scaffolds-only", help_heading = "Scaffolding and chaining")]
    scaffolds_only: bool,

    // ============================================================================
    // Advanced Filtering
    // ============================================================================
    /// Keep this fraction or tree pattern (e.g., "0.5" or "tree:3" or "tree:3,2,0.1")
    #[clap(
        short = 'x',
        long = "sparsify",
        default_value = "1.0",
        help_heading = "Advanced filtering"
    )]
    sparsify: String,

    /// Method for calculating ANI: all, orthogonal, nX[-sort] (e.g. n50, n90-identity, n100-score)
    #[clap(
        long = "ani-method",
        default_value = "n100",
        help_heading = "Advanced filtering"
    )]
    ani_method: String,

    // ============================================================================
    // General Options
    // ============================================================================
    /// Number of threads for parallel processing
    #[clap(
        short = 't',
        long = "threads",
        default_value = "8",
        help_heading = "General options"
    )]
    threads: usize,

    /// Quiet mode (no progress output)
    #[clap(long = "quiet", help_heading = "General options")]
    quiet: bool,

    /// Temporary directory for intermediate files (defaults to TMPDIR env var, then /tmp)
    #[clap(long = "tempdir", help_heading = "General options")]
    tempdir: Option<String>,

    /// Check FastGA binary locations and exit (diagnostic tool)
    #[clap(long = "check-fastga", help_heading = "General options")]
    check_fastga: bool,

    /// Report disk usage statistics (current, peak, cumulative bytes written)
    #[clap(long = "disk-usage", help_heading = "General options")]
    disk_usage: bool,

    // ============================================================================
    // AGC Archive Options
    // ============================================================================
    /// Extract only samples matching this prefix (AGC input)
    #[clap(long = "agc-prefix", help_heading = "AGC archive options")]
    agc_prefix: Option<String>,

    /// Extract only these samples (comma-separated or @file with one per line)
    #[clap(long = "agc-samples", help_heading = "AGC archive options")]
    agc_samples: Option<String>,

    /// Query samples for asymmetric alignment (prefix or comma-separated/@file list)
    #[clap(long = "agc-queries", help_heading = "AGC archive options")]
    agc_queries: Option<String>,

    /// Target samples for asymmetric alignment (prefix or comma-separated/@file list)
    #[clap(long = "agc-targets", help_heading = "AGC archive options")]
    agc_targets: Option<String>,

    /// Temp directory for AGC extraction (defaults to --tempdir)
    #[clap(long = "agc-tempdir", help_heading = "AGC archive options")]
    agc_tempdir: Option<String>,

    // ============================================================================
    // Pair Selection (for incremental/targeted alignment)
    // ============================================================================
    /// File of sample pairs to align (TSV: query<tab>target, one pair per line)
    #[clap(long = "pairs", help_heading = "Pair selection")]
    pairs_file: Option<String>,

    /// File to append completed pairs (for resume capability)
    #[clap(long = "pairs-done", help_heading = "Pair selection")]
    pairs_done: Option<String>,

    /// File to write remaining pairs after run (for resume capability)
    #[clap(long = "pairs-remaining", help_heading = "Pair selection")]
    pairs_remaining: Option<String>,

    /// List all sample pairs and exit (for generating pairs files)
    #[clap(long = "list-pairs", help_heading = "Pair selection")]
    list_pairs: bool,

    /// Shuffle pair order (for load balancing across runs)
    #[clap(long = "shuffle-pairs", help_heading = "Pair selection")]
    shuffle_pairs: bool,

    /// Seed for pair shuffling (for reproducible ordering)
    #[clap(long = "shuffle-seed", help_heading = "Pair selection")]
    shuffle_seed: Option<u64>,

    /// Maximum number of pairs to process in this run (0 = unlimited)
    #[clap(long = "max-pairs", default_value = "0", help_heading = "Pair selection")]
    max_pairs: usize,

    /// Start index for pair range (0-based, use with shuffle-seed for stable ordering)
    #[clap(long = "pair-start", default_value = "0", help_heading = "Pair selection")]
    pair_start: usize,

    /// Sparsification strategy for pair selection (none, auto, random:<frac>, giant:<prob>, tree:<near>:<far>:<random>[:<kmer>])
    #[clap(long = "sparsify-pairs", default_value = "none", help_heading = "Pair selection")]
    sparsify_pairs: String,
}

fn parse_filter_mode(mode: &str, _filter_type: &str) -> (FilterMode, Option<usize>, Option<usize>) {
    // Handle various ways to specify "no filtering"
    let lower = mode.to_lowercase();
    match lower.as_str() {
        "1:1" => (FilterMode::OneToOne, Some(1), Some(1)),
        "1" | "1:∞" | "1:infinity" | "1:many" => (FilterMode::OneToMany, Some(1), None),
        "∞:1" | "infinity:1" | "many:1" => (FilterMode::ManyToMany, None, Some(1)),
        "many:many" | "∞:∞" | "infinity:infinity" | "many" | "∞" | "infinity" | "-1" | "-1:-1" => {
            (FilterMode::ManyToMany, None, None)
        }
        s if s.contains(':') => {
            // Parse custom like "10:5" or "2:3"
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() == 2 {
                // Check for infinity/many on each side
                let per_query = match parts[0] {
                    "∞" | "infinity" | "many" | "-1" => None,
                    n => n.parse::<usize>().ok().filter(|&x| x > 0), // Reject 0
                };
                let per_target = match parts[1] {
                    "∞" | "infinity" | "many" | "-1" => None,
                    n => n.parse::<usize>().ok().filter(|&x| x > 0), // Reject 0
                };

                // Determine filter mode based on limits
                let mode = match (per_query, per_target) {
                    (Some(1), Some(1)) => FilterMode::OneToOne,
                    (Some(1), None) => FilterMode::OneToMany,
                    _ => FilterMode::ManyToMany,
                };
                return (mode, per_query, per_target);
            }
            // eprintln!("Warning: Invalid {filter_type} filter '{s}', using default 1:1");
            (FilterMode::OneToOne, Some(1), Some(1))
        }
        _ => {
            // Try parsing as a single number
            if let Ok(n) = mode.parse::<usize>() {
                if n == 0 {
                    // eprintln!("Error: 0 is not a valid filter value. Use 1 for best mapping only.");
                    std::process::exit(1);
                }
                // Single number means filter on query axis only
                return (FilterMode::OneToMany, Some(n), None);
            }
            // eprintln!("Warning: Invalid {filter_type} filter '{mode}', using default 1:1");
            (FilterMode::OneToOne, Some(1), Some(1))
        }
    }
}

/// Parse ANI calculation method from string
fn parse_ani_method(method_str: &str) -> Option<AniMethod> {
    let lower = method_str.to_lowercase();

    match lower.as_str() {
        "all" => Some(AniMethod::All),
        "orthogonal" | "1:1" => Some(AniMethod::Orthogonal),
        _ => {
            // Try to parse N-percentile methods (e.g., n50, n90-identity, n100-score)
            if let Some(rest) = lower.strip_prefix('n') {
                // Split on hyphen to get percentile and optional sort method
                let parts: Vec<&str> = rest.split('-').collect();

                // Parse the percentile number
                if let Ok(percentile) = parts[0].parse::<f64>() {
                    if percentile > 0.0 && percentile <= 100.0 {
                        // Determine sort method
                        let sort_method = if parts.len() > 1 {
                            match parts[1] {
                                "length" => NSort::Length,
                                "identity" => NSort::Identity,
                                "score" => NSort::Score,
                                _ => return None,
                            }
                        } else {
                            // Default to identity for backwards compatibility
                            NSort::Identity
                        };

                        return Some(AniMethod::NPercentile(percentile, sort_method));
                    }
                }
            }
            None
        }
    }
}

/// Parse identity value from string (fraction, percentage, or aniN)
fn parse_identity_value(value: &str, ani_percentile: Option<f64>) -> Result<f64> {
    let lower = value.to_lowercase();

    if let Some(remainder) = lower.strip_prefix("ani") {
        // Parse aniN, aniN+X, or aniN-X
        if let Some(ani_value) = ani_percentile {
            if remainder.is_empty() {
                // Default "ani" means ani50
                return Ok(ani_value);
            }

            // Parse percentile number and optional offset
            // Find offset position if exists
            let (percentile_str, offset_part) = if let Some(plus_pos) = remainder.find('+') {
                (
                    &remainder[..plus_pos],
                    Some(('+', &remainder[plus_pos + 1..])),
                )
            } else if let Some(minus_pos) = remainder.find('-') {
                (
                    &remainder[..minus_pos],
                    Some(('-', &remainder[minus_pos + 1..])),
                )
            } else {
                (remainder, None)
            };

            // For now we only support ani50 (median), ignore the percentile number
            // Could extend to support ani25, ani75, etc. in future
            if !percentile_str.is_empty() && percentile_str != "50" {
                // eprintln!(
                //     "[sweepga] WARNING: Only ani50 (median) currently supported, using median"
                // );
            }

            // Apply offset if provided
            if let Some((sign, offset_str)) = offset_part {
                let offset: f64 = offset_str
                    .parse()
                    .map_err(|_| anyhow::anyhow!("Invalid ANI offset: {offset_str}"))?;
                match sign {
                    '+' => Ok((ani_value + offset / 100.0).min(1.0)),
                    '-' => Ok((ani_value - offset / 100.0).max(0.0)),
                    _ => Ok(ani_value),
                }
            } else {
                Ok(ani_value)
            }
        } else {
            anyhow::bail!("Cannot use ANI-based threshold without input alignments");
        }
    } else if let Ok(val) = value.parse::<f64>() {
        // Numeric value
        if val > 1.0 {
            Ok(val / 100.0) // Percentage to fraction
        } else {
            Ok(val) // Already fraction
        }
    } else {
        anyhow::bail!("Invalid identity value: {value}");
    }
}

/// Calculate ANI statistics between genome pairs using specified method
fn calculate_ani_stats(input_path: &str, method: AniMethod, quiet: bool) -> Result<f64> {
    use crate::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
    use tempfile::NamedTempFile;

    let final_input_path = match method {
        AniMethod::All => {
            // Use all alignments as-is
            input_path.to_string()
        }
        AniMethod::Orthogonal => {
            // First pass: Filter to get best 1:1 mappings only
            if !quiet {
                // eprintln!("[sweepga] Calculating ANI from best 1:1 mappings (min length: 1kb)");
            }
            let temp_filtered = NamedTempFile::new()?;
            let filtered_path = temp_filtered.path().to_str().unwrap().to_string();

            // Create config for 1:1 filtering (no scaffolding, just best mappings)
            let filter_config = FilterConfig {
                chain_gap: 2000,
                min_block_length: 1000, // Ignore very short alignments for ANI
                mapping_filter_mode: FilterMode::OneToOne,
                mapping_max_per_query: Some(1),
                mapping_max_per_target: Some(1),
                plane_sweep_secondaries: 0,
                scaffold_filter_mode: FilterMode::OneToOne,
                scaffold_max_per_query: Some(1),
                scaffold_max_per_target: Some(1),
                overlap_threshold: 0.95,
                sparsity: 1.0,
                no_merge: false,
                scaffold_gap: 10000,
                min_scaffold_length: 0, // No scaffolding for ANI calculation
                scaffold_overlap_threshold: 0.95,
                scaffold_max_deviation: 0,
                prefix_delimiter: '#',
                skip_prefix: false,
                scoring_function: ScoringFunction::Matches,
                min_identity: 0.0,
                min_scaffold_identity: 0.0,
            };

            let filter = PafFilter::new(filter_config);
            filter.filter_paf(input_path, &filtered_path)?;

            // Keep temp file alive until we're done
            Box::leak(Box::new(temp_filtered));
            filtered_path
        }
        AniMethod::NPercentile(_, _) => {
            // N-percentile methods are handled differently - we'll process all alignments but only use best ones
            input_path.to_string()
        }
    };

    // For N-percentile methods, we need to collect all alignments first
    if let AniMethod::NPercentile(percentile, sort_method) = method {
        return calculate_ani_n_percentile(input_path, percentile, sort_method, quiet);
    }

    // For All and Orthogonal methods, calculate directly
    let file = File::open(&final_input_path)?;
    let reader = BufReader::new(file);

    // Group alignments by genome pair (using PanSN prefixes)
    // Store (total_matches, total_block_length) for proper weighted calculation
    let mut genome_pairs: HashMap<(String, String), (f64, f64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Extract genome prefixes (up to and including last #)
        // This handles haplotypes correctly: HG002#1#chr1 -> HG002#1#
        let query_genome = if let Some(pos) = fields[0].rfind('#') {
            fields[0][..=pos].to_string()
        } else {
            fields[0].to_string()
        };
        let target_genome = if let Some(pos) = fields[5].rfind('#') {
            fields[5][..=pos].to_string()
        } else {
            fields[5].to_string()
        };

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        // Parse matches and block length
        let matches = fields[9].parse::<f64>().unwrap_or(0.0);
        let block_len = fields[10].parse::<f64>().unwrap_or(1.0);

        // Check for divergence tag (overrides matches-based identity)
        let mut final_matches = matches;
        for field in &fields[11..] {
            if let Some(div_str) = field.strip_prefix("dv:f:") {
                if let Ok(div) = div_str.parse::<f64>() {
                    // Convert divergence back to equivalent matches
                    final_matches = (1.0 - div) * block_len;
                    break;
                }
            }
        }

        let key = if query_genome < target_genome {
            (query_genome, target_genome)
        } else {
            (target_genome, query_genome)
        };

        // Accumulate total matches and total block length for weighted calculation
        let entry = genome_pairs.entry(key).or_insert((0.0, 0.0));
        entry.0 += final_matches;
        entry.1 += block_len;
    }

    if genome_pairs.is_empty() {
        // eprintln!("[sweepga] WARNING: No inter-genome alignments found for ANI calculation");
        return Ok(0.0); // Default to no filtering
    }

    // Calculate weighted average ANI for each genome pair
    // ANI = total_matches / total_block_length
    let mut ani_values: Vec<f64> = genome_pairs
        .iter()
        .map(|((_q, _t), (total_matches, total_length))| {
            if *total_length > 0.0 {
                total_matches / total_length
            } else {
                0.0
            }
        })
        .collect();

    ani_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Get 50th percentile (median)
    let median_idx = ani_values.len() / 2;
    let ani50 = if ani_values.len().is_multiple_of(2) && ani_values.len() > 1 {
        (ani_values[median_idx - 1] + ani_values[median_idx]) / 2.0
    } else {
        ani_values[median_idx]
    };

    // eprintln!(
    //     "[sweepga] ANI statistics from {} genome pairs:",
    //     genome_pairs.len()
    // );
    // eprintln!(
    //     "[sweepga]   Min: {:.1}%, Median: {:.1}%, Max: {:.1}%",
    //     ani_values.first().unwrap_or(&0.0) * 100.0,
    //     ani50 * 100.0,
    //     ani_values.last().unwrap_or(&0.0) * 100.0
    // );

    Ok(ani50)
}

/// Calculate ANI using N-percentile method - use best alignments covering N% of genome pairs
fn calculate_ani_n_percentile(
    input_path: &str,
    percentile: f64,
    sort_method: NSort,
    quiet: bool,
) -> Result<f64> {
    // Quiet mode or logging disabled

    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    // Collect all alignments with their lengths and identity
    #[allow(dead_code)]
    struct Alignment {
        query_genome: String,
        target_genome: String,
        matches: f64,
        block_length: f64,
        identity: f64,
        query_length: u64,
        target_length: u64,
    }

    let mut alignments = Vec::new();
    let mut genome_sizes: HashMap<String, u64> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Extract genome prefixes
        let query_genome = if let Some(pos) = fields[0].rfind('#') {
            fields[0][..=pos].to_string()
        } else {
            fields[0].to_string()
        };
        let target_genome = if let Some(pos) = fields[5].rfind('#') {
            fields[5][..=pos].to_string()
        } else {
            fields[5].to_string()
        };

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        // Parse query and target lengths
        let query_length = fields[1].parse::<u64>().unwrap_or(0);
        let target_length = fields[6].parse::<u64>().unwrap_or(0);

        // Track genome+chromosome sizes to avoid double counting
        let query_key = format!(
            "{}{}",
            query_genome,
            fields[0].rsplit('#').next().unwrap_or("")
        );
        let target_key = format!(
            "{}{}",
            target_genome,
            fields[5].rsplit('#').next().unwrap_or("")
        );
        genome_sizes.entry(query_key).or_insert(query_length);
        genome_sizes.entry(target_key).or_insert(target_length);

        let matches = fields[9].parse::<f64>().unwrap_or(0.0);
        let block_len = fields[10].parse::<f64>().unwrap_or(1.0);

        // Check for divergence tag
        let mut final_matches = matches;
        for field in &fields[11..] {
            if let Some(div_str) = field.strip_prefix("dv:f:") {
                if let Ok(div) = div_str.parse::<f64>() {
                    final_matches = (1.0 - div) * block_len;
                    break;
                }
            }
        }

        let identity = final_matches / block_len.max(1.0);

        alignments.push(Alignment {
            query_genome,
            target_genome,
            matches: final_matches,
            block_length: block_len,
            identity,
            query_length,
            target_length,
        });
    }

    if alignments.is_empty() {
        // eprintln!("[sweepga] WARNING: No inter-genome alignments found for ANI calculation");
        return Ok(0.0);
    }

    // Sort based on method
    match sort_method {
        NSort::Length => {
            // Sort by block length (descending)
            alignments.sort_by(|a, b| b.block_length.partial_cmp(&a.block_length).unwrap());
        }
        NSort::Identity => {
            // Sort by identity (descending)
            alignments.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap());
        }
        NSort::Score => {
            // Sort by identity * log(length) score (descending)
            alignments.sort_by(|a, b| {
                let score_a = a.identity * a.block_length.ln().max(1.0);
                let score_b = b.identity * b.block_length.ln().max(1.0);
                score_b.partial_cmp(&score_a).unwrap()
            });
        }
    }

    // Calculate total genome sizes (sum of all unique genome+chromosome lengths)
    let total_genome_size: f64 = genome_sizes.values().map(|&v| v as f64).sum();

    // Use genome size as denominator for N-percentile, not alignment length
    // This gives us true genome coverage percentage
    let n_threshold = total_genome_size * (percentile / 100.0);

    if !quiet {
        // eprintln!(
        //     "[sweepga] Total genome size: {:.1} Mb, N{} threshold: {:.1} Mb",
        //     total_genome_size / 1_000_000.0,
        //     percentile as i32,
        //     n_threshold / 1_000_000.0
        // );
    }

    // Take alignments until we reach N-percentile threshold
    let mut cumulative_length = 0.0;
    let mut genome_pairs: HashMap<(String, String), (f64, f64)> = HashMap::new();

    for alignment in alignments {
        cumulative_length += alignment.block_length;

        let key = if alignment.query_genome < alignment.target_genome {
            (alignment.query_genome, alignment.target_genome)
        } else {
            (alignment.target_genome, alignment.query_genome)
        };

        let entry = genome_pairs.entry(key).or_insert((0.0, 0.0));
        entry.0 += alignment.matches;
        entry.1 += alignment.block_length;

        if cumulative_length >= n_threshold {
            break;
        }
    }

    // Calculate weighted ANI for each genome pair
    let mut ani_values: Vec<f64> = genome_pairs
        .iter()
        .map(|((_q, _t), (total_matches, total_length))| {
            if *total_length > 0.0 {
                total_matches / total_length
            } else {
                0.0
            }
        })
        .collect();

    ani_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Get median
    let median_idx = ani_values.len() / 2;
    let ani50 = if ani_values.len().is_multiple_of(2) && ani_values.len() > 1 {
        (ani_values[median_idx - 1] + ani_values[median_idx]) / 2.0
    } else {
        ani_values[median_idx]
    };

    Ok(ani50)
}

/// Convert .1aln file to PAF using native reader (fast path)
fn aln_to_paf_native(aln_path: &str) -> Result<tempfile::NamedTempFile> {
    use fastga_rs::AlnReader;
    use std::io::Write;

    // Create temp file for PAF output
    let mut temp_paf = tempfile::NamedTempFile::with_suffix(".paf")?;

    // Open native reader
    let mut reader = AlnReader::open(aln_path)?;

    // Convert each record to PAF format
    while let Some(rec) = reader.read_record()? {
        let qname = reader.get_seq_name(rec.query_id, 0)?;
        let tname = reader.get_seq_name(rec.target_id, 1)?;

        let aln_len = (rec.query_end - rec.query_start) as usize;
        let matches = aln_len.saturating_sub(rec.diffs as usize);
        let identity = if aln_len > 0 {
            100.0 * (matches as f64) / (aln_len as f64)
        } else {
            0.0
        };

        // PAF format: qname qlen qstart qend strand tname tlen tstart tend matches alen mapq [tags...]
        writeln!(
            temp_paf,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t60\tid:f:{:.4}",
            qname,
            rec.query_len,
            rec.query_start,
            rec.query_end,
            if rec.reverse != 0 { '-' } else { '+' },
            tname,
            rec.target_len,
            rec.target_start,
            rec.target_end,
            matches,
            aln_len,
            identity
        )?;
    }

    temp_paf.flush()?;
    Ok(temp_paf)
}

/// Convert .1aln file to PAF using ALNtoPAF (fallback) or native reader (fast)
fn aln_to_paf(aln_path: &str, threads: usize) -> Result<tempfile::NamedTempFile> {
    // Try native reader first (2.3x faster)
    if let Ok(temp_paf) = aln_to_paf_native(aln_path) {
        return Ok(temp_paf);
    }

    // Fallback to ALNtoPAF (use embedded binary)
    let alnto_paf_bin = sweepga::binary_paths::get_embedded_binary_path("ALNtoPAF")?;

    // Create temp file for PAF output
    let temp_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let paf_path = temp_paf.path();

    // Run ALNtoPAF: ALNtoPAF [-mxsS] [-T<threads>] <alignment.1aln>
    // -x: output CIGAR with X's (required if we want to convert back to .1aln later)
    // It writes PAF to stdout
    let output = std::process::Command::new(&alnto_paf_bin)
        .arg("-x") // Generate CIGAR with X's
        .arg(format!("-T{threads}"))
        .arg(aln_path)
        .output()?;

    if !output.status.success() {
        anyhow::bail!(
            "ALNtoPAF conversion failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    // Write stdout to temp file
    std::fs::write(paf_path, &output.stdout)?;

    Ok(temp_paf)
}

// Note: Native .1aln writing is available via fastga-rs::AlnWriter
// For PAF → .1aln conversion in sweepga, we still use PAFtoALN for compatibility
// Future: implement direct PAF → .1aln conversion using AlnWriter

/// Extract PanSN genome prefix from a sequence name
/// Example: "SGDref#1#chrI" -> Some("SGDref#1#")
fn extract_genome_prefix(seq_name: &str) -> Option<String> {
    // Look for PanSN format: genome#haplotype#chromosome
    // We want to keep genome#haplotype# as the group
    let parts: Vec<&str> = seq_name.split('#').collect();
    if parts.len() >= 2 {
        // Keep first two parts plus trailing #
        Some(format!("{}#{}#", parts[0], parts[1]))
    } else {
        None
    }
}

/// Detect genome groups from a FASTA file by reading headers
fn detect_genome_groups(fasta_path: &Path) -> Result<Vec<String>> {
    use std::collections::BTreeSet;
    use std::io::{BufRead, BufReader};

    let file = File::open(fasta_path)?;
    let reader: Box<dyn BufRead> = if fasta_path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut groups = BTreeSet::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(stripped) = line.strip_prefix('>') {
            let seq_name = stripped.split_whitespace().next().unwrap_or("");
            if let Some(prefix) = extract_genome_prefix(seq_name) {
                groups.insert(prefix);
            }
        }
    }

    Ok(groups.into_iter().collect())
}

/// Write sequences belonging to a specific genome group to a new FASTA file
fn write_genome_fasta(input_path: &Path, output_path: &Path, genome_prefix: &str) -> Result<usize> {
    use std::io::{BufRead, BufReader, Write};

    let file = File::open(input_path)?;
    let reader: Box<dyn BufRead> = if input_path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(noodles::bgzf::io::reader::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut output = File::create(output_path)?;
    let mut writing = false;
    let mut seq_count = 0;

    for line in reader.lines() {
        let line = line?;
        if let Some(stripped) = line.strip_prefix('>') {
            let seq_name = stripped.split_whitespace().next().unwrap_or("");
            writing = if let Some(prefix) = extract_genome_prefix(seq_name) {
                prefix == genome_prefix
            } else {
                false
            };

            if writing {
                seq_count += 1;
                writeln!(output, "{line}")?;
            }
        } else if writing {
            writeln!(output, "{line}")?;
        }
    }

    output.flush()?;
    Ok(seq_count)
}

/// Perform all-pairs pairwise alignment for multiple FASTA files
/// If a single FASTA with multiple genomes, splits by PanSN prefix
#[allow(clippy::too_many_arguments)]
fn align_multiple_fastas(
    fasta_files: &[String],
    frequency: Option<usize>,
    all_pairs: bool,
    batch_bytes: Option<u64>,
    threads: usize,
    keep_self: bool,
    tempdir: Option<&str>,
    timing: &TimingContext,
    quiet: bool,
    min_alignment_length: u64,
    zstd_compress: bool,
    zstd_level: u32,
) -> Result<tempfile::NamedTempFile> {
    // Detect genome groups from input files
    let mut num_genomes = 0;
    for fasta_file in fasta_files {
        let path = Path::new(fasta_file);
        let groups = detect_genome_groups(path)?;
        num_genomes += groups.len();
        if !quiet {
            timing.log(
                "detect",
                &format!("{}: {} genome groups", fasta_file, groups.len()),
            );
        }
    }

    // Check for batch mode
    if let Some(max_index_bytes) = batch_bytes {
        if !quiet {
            timing.log(
                "batch",
                &format!(
                    "Batch mode enabled: max {} index size per batch",
                    batch_align::format_bytes(max_index_bytes)
                ),
            );
        }

        let config = batch_align::BatchAlignConfig {
            frequency,
            threads,
            min_alignment_length,
            zstd_compress,
            zstd_level,
            keep_self,
            quiet,
        };

        return batch_align::run_batch_alignment(fasta_files, max_index_bytes, &config, tempdir);
    }

    // Decide on alignment mode
    if all_pairs {
        // --all-pairs mode: split genomes and align each pair separately (bidirectional)
        return align_all_pairs_mode(
            fasta_files,
            frequency,
            threads,
            keep_self,
            tempdir,
            timing,
            quiet,
            min_alignment_length,
            zstd_compress,
            zstd_level,
        );
    }

    // Default mode: single FastGA run with adaptive frequency
    // Let fastga_integration.rs handle auto-detection via PanSN haplotype counting
    let effective_frequency = if let Some(f) = frequency {
        // User specified -f, use it as-is
        if !quiet {
            timing.log(
                "align",
                &format!("Using user-specified frequency threshold: {f}"),
            );
        }
        Some(f)
    } else {
        // No user-specified frequency: let fastga_integration.rs auto-detect
        // from PanSN naming convention
        None
    };

    if !quiet {
        timing.log(
            "align",
            &format!("Aligning {num_genomes} genomes in single FastGA run"),
        );
    }

    // Create FastGA integration with effective frequency
    let fastga = create_fastga_integration(
        effective_frequency,
        threads,
        min_alignment_length,
        tempdir.map(String::from),
    )?;

    // Run single FastGA alignment on the original input (all genomes together)
    // Use first file as both query and target for self-alignment
    let input_file = &fasta_files[0];
    // Canonicalize path to avoid empty parent directory issues in fastga-rs
    let input_path = std::fs::canonicalize(input_file)
        .with_context(|| format!("Failed to resolve path: {input_file}"))?;

    // If zstd compression is requested, pre-build and compress the index
    if zstd_compress {
        if !quiet {
            timing.log("index", "Building index for compression...");
        }
        let gdb_base = fastga.prepare_gdb(&input_path)?;
        // Strip .fa extension if present (fastga-rs returns path with .fa but index uses base name)
        let gdb_base_stripped = gdb_base
            .strip_suffix(".fa")
            .or_else(|| gdb_base.strip_suffix(".fna"))
            .or_else(|| gdb_base.strip_suffix(".fasta"))
            .unwrap_or(&gdb_base);
        fastga_integration::FastGAIntegration::compress_index(gdb_base_stripped, zstd_level)?;
        if !quiet {
            timing.log("index", "Index compressed with zstd");
        }
    }

    let temp_paf = fastga.align_to_temp_paf(&input_path, &input_path)?;

    if !quiet {
        let alignment_count = std::fs::read_to_string(temp_paf.path())?.lines().count();
        timing.log(
            "align",
            &format!("FastGA produced {alignment_count} alignments"),
        );
    }

    Ok(temp_paf)
}

/// Check if pair-mode processing is requested
fn is_pair_mode(args: &Args) -> bool {
    args.pairs_file.is_some()
        || args.list_pairs
        || args.shuffle_pairs
        || args.pair_start > 0
        || args.max_pairs > 0
        || args.pairs_done.is_some()
        || args.pairs_remaining.is_some()
        || args.sparsify_pairs != "none"
}

fn process_agc_archive(
    agc_path: &str,
    args: &Args,
    agc_tempdir: Option<&str>,
    timing: &TimingContext,
) -> Result<tempfile::NamedTempFile> {
    use agc::AgcSource;

    if !args.quiet {
        timing.log("agc", &format!("Opening AGC archive: {}", agc_path));
    }

    let mut agc = AgcSource::open(agc_path)?;

    // Determine temp directory for extraction
    let temp_base = if let Some(dir) = agc_tempdir {
        std::path::PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        std::path::PathBuf::from(tmpdir)
    } else {
        std::path::PathBuf::from("/tmp")
    };

    // Check for pair-mode processing (--list-pairs, --pairs, --shuffle-pairs, etc.)
    if is_pair_mode(args) || args.list_pairs {
        // Generate or read pairs
        let pairs = get_pairs_from_args(&mut agc, args, timing, Some(&temp_base))?;

        // Handle --list-pairs: output pairs and exit
        if args.list_pairs {
            // Apply shuffle and range before listing
            let pairs = apply_pair_filters(pairs, args, timing)?;
            for pair in &pairs {
                println!("{}", pair.to_tsv());
            }
            // Return empty PAF file (won't be used since we exit early)
            let output_paf = tempfile::Builder::new()
                .prefix("sweepga_pairs_")
                .suffix(".paf")
                .tempfile_in(&temp_base)?;
            return Ok(output_paf);
        }

        // Apply filters (shuffle, range, done filtering)
        let pairs = apply_pair_filters(pairs, args, timing)?;

        // Process pairs
        return process_agc_pairs(&mut agc, &pairs, args, &temp_base, timing);
    }

    // Standard mode: process all samples together

    // Determine which samples to process
    let samples = determine_agc_samples(&mut agc, args, timing)?;

    if samples.is_empty() {
        anyhow::bail!("No samples found in AGC archive matching the specified criteria");
    }

    if !args.quiet {
        timing.log("agc", &format!("Selected {} samples for alignment", samples.len()));
    }

    // Get sizes for batch planning
    let sizes = agc.get_sample_sizes_for(&samples)?;
    let total_bp: u64 = sizes.values().sum();

    if !args.quiet {
        timing.log(
            "agc",
            &format!(
                "Total sequence: {} bp ({} samples)",
                batch_align::format_bytes(total_bp),
                samples.len()
            ),
        );
    }

    // Check if we should use batch mode
    if let Some(max_index_bytes) = args.batch_bytes {
        // Batch mode: partition samples and process incrementally
        return process_agc_batched(
            &mut agc,
            &samples,
            &sizes,
            max_index_bytes,
            args,
            &temp_base,
            timing,
        );
    }

    // Non-batch mode: extract all samples to a single temp FASTA
    if !args.quiet {
        timing.log("agc", "Extracting all samples to temp FASTA...");
    }

    let temp_fasta = agc.extract_samples_to_temp(&samples, Some(&temp_base))?;
    let fasta_path = temp_fasta.path().to_path_buf();

    if !args.quiet {
        timing.log("agc", &format!("Extracted to: {}", fasta_path.display()));
    }

    // Run FastGA alignment using direct PAF output (bypasses ALNtoPAF segfaults)
    let fastga = create_fastga_integration(
        args.frequency,
        args.threads,
        args.block_length,
        args.tempdir.clone(),
    )?;

    if !args.quiet {
        timing.log("align", "Running FastGA alignment...");
    }

    let paf_bytes = fastga.align_direct_paf(&fasta_path, &fasta_path)?;

    // Write to temp file
    use std::io::Write;
    let mut temp_paf = tempfile::Builder::new()
        .prefix("sweepga_agc_")
        .suffix(".paf")
        .tempfile_in(&temp_base)?;
    temp_paf.write_all(&paf_bytes)?;
    temp_paf.flush()?;

    // temp_fasta will be cleaned up when it goes out of scope
    // But we need to keep it alive until alignment is done, which it is now

    Ok(temp_paf)
}

/// Determine which samples to extract from an AGC archive based on CLI args
fn determine_agc_samples(
    agc: &mut agc::AgcSource,
    args: &Args,
    timing: &TimingContext,
) -> Result<Vec<String>> {
    use agc::parse_sample_list;

    // Check for asymmetric alignment (queries vs targets)
    if args.agc_queries.is_some() || args.agc_targets.is_some() {
        let queries = get_agc_sample_set(agc, args.agc_queries.as_deref(), "queries", timing)?;
        let targets = get_agc_sample_set(agc, args.agc_targets.as_deref(), "targets", timing)?;

        if !args.quiet {
            timing.log(
                "agc",
                &format!("{} queries, {} targets", queries.len(), targets.len()),
            );
        }

        // For asymmetric alignment, we need both sets
        // Return the union for now (TODO: implement true asymmetric alignment)
        let mut all_samples: std::collections::HashSet<String> = queries.into_iter().collect();
        all_samples.extend(targets);
        return Ok(all_samples.into_iter().collect());
    }

    // Check for explicit sample list
    if let Some(ref sample_spec) = args.agc_samples {
        let samples = parse_sample_list(sample_spec)?;
        if !args.quiet {
            timing.log("agc", &format!("Using {} samples from list", samples.len()));
        }
        return Ok(samples);
    }

    // Check for prefix filter
    if let Some(ref prefix) = args.agc_prefix {
        let samples = agc.list_samples_with_prefix(prefix);
        if !args.quiet {
            timing.log(
                "agc",
                &format!("{} samples matching prefix '{}'", samples.len(), prefix),
            );
        }
        return Ok(samples);
    }

    // Default: all samples
    let samples = agc.list_samples();
    if !args.quiet {
        timing.log("agc", &format!("Using all {} samples", samples.len()));
    }
    Ok(samples)
}

/// Get a set of samples from AGC based on a spec (prefix or list)
fn get_agc_sample_set(
    agc: &mut agc::AgcSource,
    spec: Option<&str>,
    _name: &str,
    _timing: &TimingContext,
) -> Result<Vec<String>> {
    match spec {
        None => {
            // No spec means all samples
            Ok(agc.list_samples())
        }
        Some(s) if s.starts_with('@') || s.contains(',') => {
            // Explicit list
            agc::parse_sample_list(s)
        }
        Some(prefix) => {
            // Prefix match
            Ok(agc.list_samples_with_prefix(prefix))
        }
    }
}

/// A pair of samples to align (query, target)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SamplePair {
    pub query: String,
    pub target: String,
}

impl SamplePair {
    pub fn new(query: String, target: String) -> Self {
        Self { query, target }
    }

    /// Parse from TSV line (query\ttarget)
    pub fn from_tsv(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() >= 2 {
            Some(Self::new(parts[0].to_string(), parts[1].to_string()))
        } else {
            None
        }
    }

    /// Format as TSV line
    pub fn to_tsv(&self) -> String {
        format!("{}\t{}", self.query, self.target)
    }
}

/// Read pairs from a TSV file
fn read_pairs_file(path: &str) -> Result<Vec<SamplePair>> {
    let file = File::open(path).context(format!("Failed to open pairs file: {}", path))?;
    let reader = BufReader::new(file);
    let mut pairs = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(pair) = SamplePair::from_tsv(trimmed) {
            pairs.push(pair);
        }
    }

    Ok(pairs)
}

/// Generate all pairs from a list of samples (symmetric, no self-pairs)
fn generate_all_pairs(samples: &[String]) -> Vec<SamplePair> {
    let mut pairs = Vec::new();
    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
            pairs.push(SamplePair::new(samples[i].clone(), samples[j].clone()));
        }
    }
    pairs
}

/// Generate cartesian product of queries x targets
fn generate_cartesian_pairs(queries: &[String], targets: &[String]) -> Vec<SamplePair> {
    let mut pairs = Vec::new();
    for q in queries {
        for t in targets {
            if q != t {
                pairs.push(SamplePair::new(q.clone(), t.clone()));
            }
        }
    }
    pairs
}

/// Filter out pairs that are already in the done file
fn filter_done_pairs(pairs: Vec<SamplePair>, done_path: &str) -> Result<Vec<SamplePair>> {
    let done_pairs: std::collections::HashSet<SamplePair> = if Path::new(done_path).exists() {
        read_pairs_file(done_path)?.into_iter().collect()
    } else {
        std::collections::HashSet::new()
    };

    Ok(pairs.into_iter().filter(|p| !done_pairs.contains(p)).collect())
}

/// Shuffle pairs with optional seed for reproducibility
fn shuffle_pairs(mut pairs: Vec<SamplePair>, seed: Option<u64>) -> Vec<SamplePair> {
    use rand::seq::SliceRandom;
    use rand::SeedableRng;

    if let Some(seed) = seed {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        pairs.shuffle(&mut rng);
    } else {
        let mut rng = rand::thread_rng();
        pairs.shuffle(&mut rng);
    }
    pairs
}

/// Write pairs to a TSV file
fn write_pairs_file(path: &str, pairs: &[SamplePair]) -> Result<()> {
    use std::io::Write;
    let mut file = File::create(path).context(format!("Failed to create pairs file: {}", path))?;
    for pair in pairs {
        writeln!(file, "{}", pair.to_tsv())?;
    }
    Ok(())
}

/// Append a pair to a file (for done tracking)
fn append_pair_to_file(path: &str, pair: &SamplePair) -> Result<()> {
    use std::io::Write;
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .context(format!("Failed to open pairs file for append: {}", path))?;
    writeln!(file, "{}", pair.to_tsv())?;
    Ok(())
}

/// Generate or read pairs based on CLI arguments
fn get_pairs_from_args(
    agc: &mut agc::AgcSource,
    args: &Args,
    timing: &TimingContext,
    temp_base: Option<&std::path::Path>,
) -> Result<Vec<SamplePair>> {
    // Option 1: Read pairs from file
    if let Some(ref pairs_file) = args.pairs_file {
        let pairs = read_pairs_file(pairs_file)?;
        if !args.quiet {
            timing.log("pairs", &format!("Read {} pairs from file", pairs.len()));
        }
        return Ok(pairs);
    }

    // Option 2: Cartesian product of queries x targets
    if args.agc_queries.is_some() && args.agc_targets.is_some() {
        let queries = get_agc_sample_set(agc, args.agc_queries.as_deref(), "queries", timing)?;
        let targets = get_agc_sample_set(agc, args.agc_targets.as_deref(), "targets", timing)?;

        let pairs = generate_cartesian_pairs(&queries, &targets);
        if !args.quiet {
            timing.log(
                "pairs",
                &format!(
                    "Generated {} pairs ({} queries x {} targets)",
                    pairs.len(),
                    queries.len(),
                    targets.len()
                ),
            );
        }
        return Ok(pairs);
    }

    // Get samples list
    let samples = if let Some(ref prefix) = args.agc_prefix {
        agc.list_samples_with_prefix(prefix)
    } else if let Some(ref sample_spec) = args.agc_samples {
        agc::parse_sample_list(sample_spec)?
    } else {
        agc.list_samples()
    };

    // Parse sparsification strategy
    let strategy: knn_graph::SparsificationStrategy = args
        .sparsify_pairs
        .parse()
        .map_err(|e| anyhow::anyhow!("Invalid sparsification strategy: {}", e))?;

    // Option 3: Sparsified pairs using minhash/knn
    if strategy != knn_graph::SparsificationStrategy::None {
        if !args.quiet {
            timing.log("pairs", &format!("Sparsification: {}", strategy));
        }

        // For sparsification, we need to compute mash sketches
        // Stream one sample at a time to avoid materializing all sequences in memory

        // Determine k-mer size from strategy
        let kmer_size = match &strategy {
            knn_graph::SparsificationStrategy::TreeSampling(_, _, _, Some(k)) => *k,
            _ => mash::DEFAULT_KMER_SIZE,
        };

        if !args.quiet {
            timing.log(
                "mash",
                &format!(
                    "Computing sketches for {} samples (k={})...",
                    samples.len(),
                    kmer_size
                ),
            );
        }

        // Stream: extract each sample directly to memory, compute sketch, discard sequence
        // Only sketches are kept in memory (~8KB each vs ~12MB per yeast genome)
        // No temp files needed - RAGC extracts directly to Vec<u8>
        let mut sketches: Vec<mash::KmerSketch> = Vec::with_capacity(samples.len());
        for (i, sample) in samples.iter().enumerate() {
            // Extract directly to memory (no temp file)
            let seq = agc.extract_sample_to_bytes(sample)?;

            // Compute sketch and discard sequence immediately
            let sketch = mash::KmerSketch::from_sequence(&seq, kmer_size, mash::DEFAULT_SKETCH_SIZE);
            sketches.push(sketch);
            drop(seq); // Explicitly drop to free memory

            // Progress indicator for large sets
            if !args.quiet && samples.len() > 100 && (i + 1) % 50 == 0 {
                timing.log("mash", &format!("  sketched {}/{} samples", i + 1, samples.len()));
            }
        }

        // Select pairs using pre-computed sketches (no sequence data needed)
        let pair_indices = knn_graph::select_pairs_from_sketches(&sketches, &strategy);

        // Convert indices to SamplePair
        let pairs: Vec<SamplePair> = pair_indices
            .into_iter()
            .map(|(i, j)| SamplePair::new(samples[i].clone(), samples[j].clone()))
            .collect();

        if !args.quiet {
            let total_possible = samples.len() * (samples.len() - 1) / 2;
            let percentage = (pairs.len() as f64 / total_possible as f64) * 100.0;
            timing.log(
                "pairs",
                &format!(
                    "Selected {} pairs from {} samples ({:.1}% of {} possible)",
                    pairs.len(),
                    samples.len(),
                    percentage,
                    total_possible
                ),
            );
        }

        return Ok(pairs);
    }

    // Option 4: All pairs from samples (no sparsification)
    let pairs = generate_all_pairs(&samples);
    if !args.quiet {
        timing.log(
            "pairs",
            &format!(
                "Generated {} pairs from {} samples",
                pairs.len(),
                samples.len()
            ),
        );
    }
    Ok(pairs)
}

/// Apply filtering to pairs: shuffle, range selection, done filtering
fn apply_pair_filters(
    pairs: Vec<SamplePair>,
    args: &Args,
    timing: &TimingContext,
) -> Result<Vec<SamplePair>> {
    let mut pairs = pairs;

    // Apply shuffle if requested
    if args.shuffle_pairs {
        pairs = shuffle_pairs(pairs, args.shuffle_seed);
        if !args.quiet {
            let seed_str = args
                .shuffle_seed
                .map(|s| format!(" (seed={})", s))
                .unwrap_or_default();
            timing.log("pairs", &format!("Shuffled pairs{}", seed_str));
        }
    }

    // Apply range selection (start + max)
    if args.pair_start > 0 || args.max_pairs > 0 {
        let total = pairs.len();
        let start = args.pair_start.min(total);
        let end = if args.max_pairs > 0 {
            (start + args.max_pairs).min(total)
        } else {
            total
        };
        pairs = pairs[start..end].to_vec();
        if !args.quiet {
            timing.log(
                "pairs",
                &format!("Selected pairs {}..{} of {}", start, end, total),
            );
        }
    }

    // Filter out done pairs if specified
    if let Some(ref done_path) = args.pairs_done {
        let before = pairs.len();
        pairs = filter_done_pairs(pairs, done_path)?;
        let filtered = before - pairs.len();
        if !args.quiet && filtered > 0 {
            timing.log(
                "pairs",
                &format!("Filtered {} done pairs, {} remaining", filtered, pairs.len()),
            );
        }
    }

    // Write remaining pairs to file if requested
    if let Some(ref remaining_path) = args.pairs_remaining {
        write_pairs_file(remaining_path, &pairs)?;
        if !args.quiet {
            timing.log(
                "pairs",
                &format!("Wrote {} remaining pairs to {}", pairs.len(), remaining_path),
            );
        }
    }

    Ok(pairs)
}

/// Process AGC archive pair-by-pair with checkpointing
#[allow(clippy::too_many_arguments)]
fn process_agc_pairs(
    agc: &mut agc::AgcSource,
    pairs: &[SamplePair],
    args: &Args,
    temp_base: &std::path::Path,
    timing: &TimingContext,
) -> Result<tempfile::NamedTempFile> {
    use std::io::Write;

    if pairs.is_empty() {
        if !args.quiet {
            timing.log("pairs", "No pairs to process");
        }
        // Return empty PAF file
        let output_paf = tempfile::Builder::new()
            .prefix("sweepga_pairs_")
            .suffix(".paf")
            .tempfile_in(temp_base)?;
        return Ok(output_paf);
    }

    if !args.quiet {
        timing.log("pairs", &format!("Processing {} pairs", pairs.len()));
    }

    // Create FastGA integration
    let fastga = create_fastga_integration(
        args.frequency,
        args.threads,
        args.block_length,
        Some(temp_base.to_string_lossy().to_string()),
    )?;

    // Create work directory for extractions
    let work_dir = temp_base.join(format!("sweepga_pairs_{}", std::process::id()));
    std::fs::create_dir_all(&work_dir)?;

    // Create output PAF file
    let mut output_paf = tempfile::Builder::new()
        .prefix("sweepga_pairs_")
        .suffix(".paf")
        .tempfile_in(temp_base)?;

    let mut total_alignments = 0usize;
    let mut processed = 0;

    for pair in pairs {
        processed += 1;
        if !args.quiet {
            eprintln!(
                "[pairs] ({}/{}) {} -> {}",
                processed,
                pairs.len(),
                pair.query,
                pair.target
            );
        }
        let prev_total = total_alignments;

        // Extract query and target to temp FASTA files
        let query_fasta = work_dir.join(format!("{}.fa", pair.query.replace('#', "_")));
        let target_fasta = work_dir.join(format!("{}.fa", pair.target.replace('#', "_")));

        // Extract samples (skip if same and file exists)
        if pair.query == pair.target {
            agc.extract_sample_to_fasta(&pair.query, &query_fasta)?;
            // Run self-alignment
            let paf_bytes = fastga.align_direct_paf(&query_fasta, &query_fasta)?;
            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;
            output_paf.as_file_mut().write_all(&paf_bytes)?;

            // Cleanup
            let _ = std::fs::remove_file(&query_fasta);
        } else {
            // Extract both samples
            agc.extract_sample_to_fasta(&pair.query, &query_fasta)?;
            agc.extract_sample_to_fasta(&pair.target, &target_fasta)?;

            // Create index for target and align query to it
            let target_gdb = fastga.create_gdb_only(&target_fasta)?;
            fastga.create_index_only(&target_gdb)?;

            if args.zstd_compress {
                fastga_integration::FastGAIntegration::compress_index(&target_gdb, args.zstd_level)?;
            }

            let paf_bytes = fastga.align_direct_paf(&query_fasta, &target_fasta)?;
            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;
            output_paf.as_file_mut().write_all(&paf_bytes)?;

            // Cleanup
            let _ = fastga_integration::FastGAIntegration::cleanup_all(&target_gdb);
            let _ = std::fs::remove_file(&query_fasta);
            let _ = std::fs::remove_file(&target_fasta);
        }

        // Record completion if done file specified
        if let Some(ref done_path) = args.pairs_done {
            append_pair_to_file(done_path, pair)?;
        }

        if !args.quiet {
            let pair_alignments = total_alignments - prev_total;
            eprintln!(
                "[pairs]   {} alignments (total: {})",
                pair_alignments,
                total_alignments
            );
        }
    }

    // Cleanup work directory
    let _ = std::fs::remove_dir_all(&work_dir);

    if !args.quiet {
        timing.log(
            "pairs",
            &format!(
                "Completed {} pairs, {} total alignments",
                pairs.len(),
                total_alignments
            ),
        );
    }

    Ok(output_paf)
}

/// Process AGC archive in batches to limit disk usage
fn process_agc_batched(
    agc: &mut agc::AgcSource,
    samples: &[String],
    sizes: &std::collections::HashMap<String, u64>,
    max_index_bytes: u64,
    args: &Args,
    temp_base: &std::path::Path,
    _timing: &TimingContext,
) -> Result<tempfile::NamedTempFile> {
    use std::io::Write;

    if !args.quiet {
        eprintln!(
            "[batch] Batch mode: max {} index size",
            batch_align::format_bytes(max_index_bytes)
        );
    }

    // Partition samples into batches based on estimated index size
    // Index estimate: 0.1GB + 12 bytes per bp
    let mut batches: Vec<Vec<String>> = Vec::new();
    let mut current_batch: Vec<String> = Vec::new();
    let mut current_size: u64 = 0;

    for sample in samples {
        let sample_size = sizes.get(sample).copied().unwrap_or(0);
        let estimated_index = 100_000_000 + sample_size * 12; // 0.1GB + 12 bytes/bp

        if current_size + estimated_index > max_index_bytes && !current_batch.is_empty() {
            batches.push(std::mem::take(&mut current_batch));
            current_size = 0;
        }

        current_batch.push(sample.clone());
        current_size += estimated_index;
    }

    if !current_batch.is_empty() {
        batches.push(current_batch);
    }

    let num_batches = batches.len();
    if !args.quiet {
        eprintln!("[batch] Partitioned {} samples into {} batches", samples.len(), num_batches);
    }

    // Phase 1: Extract all batches to temp FASTA files
    let batch_dir = temp_base.join(format!("sweepga_agc_batches_{}", std::process::id()));
    std::fs::create_dir_all(&batch_dir)?;

    let mut batch_files: Vec<std::path::PathBuf> = Vec::new();
    for (i, batch) in batches.iter().enumerate() {
        let batch_file = batch_dir.join(format!("batch_{}.fa", i));
        agc.extract_samples_to_fasta(batch, &batch_file)?;
        if !args.quiet {
            eprintln!("[batch] Extracted batch {} ({} samples)", i + 1, batch.len());
        }
        batch_files.push(batch_file);
    }

    // Create FastGA integration
    // Use temp_base (not batch_dir) so TMPDIR isn't set to a directory we'll delete
    let fastga = create_fastga_integration(
        args.frequency,
        args.threads,
        args.block_length,
        Some(temp_base.to_string_lossy().to_string()),
    )?;

    // Phase 2: Create GDB for all batches
    if !args.quiet {
        eprintln!("[batch] Creating GDB files...");
    }
    let mut gdb_bases: Vec<String> = Vec::new();
    for (i, batch_file) in batch_files.iter().enumerate() {
        if !args.quiet {
            eprintln!("[batch]   GDB for batch {}...", i + 1);
        }
        let gdb_base = fastga.create_gdb_only(batch_file)?;
        gdb_bases.push(gdb_base);
    }

    // Create merged output file
    let mut output_paf = tempfile::Builder::new()
        .prefix("sweepga_agc_")
        .suffix(".paf")
        .tempfile_in(temp_base)?;

    let mut total_alignments = 0;

    // Phase 3: All-pairs alignment between batches
    // For each target batch i, create index, run alignments, cleanup index
    for i in 0..num_batches {
        if !args.quiet {
            eprintln!("[batch] Building index for batch {}...", i + 1);
        }
        fastga.create_index_only(&gdb_bases[i])?;

        // Compress if requested
        if args.zstd_compress {
            fastga_integration::FastGAIntegration::compress_index(&gdb_bases[i], args.zstd_level)?;
        }

        for j in 0..num_batches {
            if !args.quiet {
                if i == j {
                    eprintln!("[batch] Aligning batch {} to itself...", i + 1);
                } else {
                    eprintln!("[batch] Aligning batch {} vs batch {}...", j + 1, i + 1);
                }
            }

            // Run alignment: query batch j against target batch i
            // Use direct PAF output to avoid ALNtoPAF conversion issues
            let paf_bytes = fastga.align_direct_paf(&batch_files[j], &batch_files[i])?;

            let line_count = paf_bytes.iter().filter(|&&b| b == b'\n').count();
            total_alignments += line_count;

            output_paf.as_file_mut().write_all(&paf_bytes)?;

            if !args.quiet {
                eprintln!("[batch]   {} alignments", line_count);
            }

            // Clean up query index if FastGA auto-created one (j != i)
            if j != i {
                let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[j]);
            }
        }

        // Clean up index for batch i to free disk space
        if !args.quiet {
            eprintln!("[batch] Cleaning up index for batch {}...", i + 1);
        }
        let _ = fastga_integration::FastGAIntegration::cleanup_index(&gdb_bases[i]);
    }

    // Phase 4: Clean up all GDB files and batch directory
    for gdb_base in &gdb_bases {
        let _ = fastga_integration::FastGAIntegration::cleanup_all(gdb_base);
    }
    let _ = std::fs::remove_dir_all(&batch_dir);

    if !args.quiet {
        eprintln!("[batch] Completed batch alignment: {} total alignments", total_alignments);
    }

    Ok(output_paf)
}

/// Align all genome pairs separately in both directions (for --all-pairs mode)
#[allow(clippy::too_many_arguments)]
fn align_all_pairs_mode(
    fasta_files: &[String],
    frequency: Option<usize>,
    threads: usize,
    keep_self: bool,
    tempdir: Option<&str>,
    timing: &TimingContext,
    quiet: bool,
    min_alignment_length: u64,
    zstd_compress: bool,
    zstd_level: u32,
) -> Result<tempfile::NamedTempFile> {
    use std::io::Write;

    // Determine temp directory
    let temp_base = if let Some(dir) = tempdir {
        std::path::PathBuf::from(dir)
    } else if let Ok(tmpdir) = std::env::var("TMPDIR") {
        std::path::PathBuf::from(tmpdir)
    } else {
        std::path::PathBuf::from("/tmp")
    };

    // Create temp directory for genome FASTAs
    let temp_genome_dir = temp_base.join(format!("sweepga_genomes_{}", std::process::id()));
    std::fs::create_dir_all(&temp_genome_dir)?;

    if !quiet {
        timing.log(
            "split",
            &format!("Using temp directory: {}", temp_genome_dir.display()),
        );
    }

    // Collect all genome groups from all input files
    let mut all_groups: Vec<(String, String)> = Vec::new(); // (genome_prefix, source_fasta)

    for fasta_file in fasta_files {
        let path = Path::new(fasta_file);
        let groups = detect_genome_groups(path)?;

        for group in groups {
            all_groups.push((group, fasta_file.clone()));
        }
    }

    // Write per-genome FASTA files
    let mut genome_files: HashMap<String, std::path::PathBuf> = HashMap::new();

    for (genome_prefix, source_file) in &all_groups {
        if genome_files.contains_key(genome_prefix) {
            continue; // Already processed
        }

        let genome_name = genome_prefix.trim_end_matches('#').replace('#', "_");
        let genome_path = temp_genome_dir.join(format!("{genome_name}.fa"));

        let seq_count = write_genome_fasta(Path::new(source_file), &genome_path, genome_prefix)?;

        if !quiet {
            timing.log(
                "split",
                &format!(
                    "Wrote {} sequences for {} to {}",
                    seq_count,
                    genome_prefix,
                    genome_path.display()
                ),
            );
        }

        genome_files.insert(genome_prefix.clone(), genome_path);
    }

    // Get unique genome prefixes in deterministic order
    let mut genome_prefixes: Vec<String> = genome_files.keys().cloned().collect();
    genome_prefixes.sort();

    if !quiet {
        timing.log(
            "index",
            &format!("Building indices for {} genomes", genome_prefixes.len()),
        );
    }

    // Create FastGA integration
    let fastga = create_fastga_integration(
        frequency,
        threads,
        min_alignment_length,
        tempdir.map(String::from),
    )?;

    // Step 1: Build GDB and GIX indices for all genomes (once each)
    for genome_prefix in &genome_prefixes {
        let fasta_path = &genome_files[genome_prefix];

        if !quiet {
            timing.log(
                "index",
                &format!("Building index for {}", genome_prefix.trim_end_matches('#')),
            );
        }

        // Create .gdb and .gix files
        let gdb_base = fastga.prepare_gdb(fasta_path)?;

        // Compress index if requested
        if zstd_compress {
            // Strip .fa extension if present (fastga-rs returns path with .fa but index uses base name)
            let gdb_base_stripped = gdb_base
                .strip_suffix(".fa")
                .or_else(|| gdb_base.strip_suffix(".fna"))
                .or_else(|| gdb_base.strip_suffix(".fasta"))
                .unwrap_or(&gdb_base);
            fastga_integration::FastGAIntegration::compress_index(gdb_base_stripped, zstd_level)?;
        }
    }

    if !quiet {
        timing.log(
            "align",
            &format!("Aligning {} genomes pairwise", genome_prefixes.len()),
        );
    }

    // Create merged PAF output file
    let merged_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    let mut merged_output = File::create(merged_paf.path())?;

    // Align all pairs in both directions (complete matrix)
    let mut total_pairs = 0;
    let mut total_alignments = 0;

    for i in 0..genome_prefixes.len() {
        for j in 0..genome_prefixes.len() {
            // Skip self-alignments
            if i == j {
                continue;
            }

            let genome_i = &genome_prefixes[i];
            let genome_j = &genome_prefixes[j];

            let fasta_i = &genome_files[genome_i];
            let fasta_j = &genome_files[genome_j];

            if !quiet {
                timing.log(
                    "align",
                    &format!(
                        "Aligning {} vs {}",
                        genome_i.trim_end_matches('#'),
                        genome_j.trim_end_matches('#')
                    ),
                );
            }

            // Align using pre-built GDB/GIX indices (FastGA auto-detects them)
            let temp_paf = fastga.align_to_temp_paf(fasta_i, fasta_j)?;

            // Append to merged output
            let paf_content = std::fs::read_to_string(temp_paf.path())?;
            merged_output.write_all(paf_content.as_bytes())?;
            merged_output.flush()?; // Flush after each alignment to prevent corruption

            let alignment_count = paf_content.lines().count();
            total_alignments += alignment_count;
            total_pairs += 1;
        }

        // Handle self-alignments if requested
        if keep_self {
            let genome_i = &genome_prefixes[i];
            let fasta_i = &genome_files[genome_i];

            if !quiet {
                timing.log(
                    "align",
                    &format!("Self-aligning {}", genome_i.trim_end_matches('#')),
                );
            }

            let temp_paf = fastga.align_to_temp_paf(fasta_i, fasta_i)?;

            let paf_content = std::fs::read_to_string(temp_paf.path())?;
            merged_output.write_all(paf_content.as_bytes())?;
            merged_output.flush()?; // Flush after each alignment to prevent corruption

            let alignment_count = paf_content.lines().count();
            total_alignments += alignment_count;
            total_pairs += 1;
        }
    }

    // Final flush (redundant but safe)
    merged_output.flush()?;

    if !quiet {
        timing.log(
            "align",
            &format!(
                "Completed {total_pairs} pairwise alignments, {total_alignments} total alignments"
            ),
        );
    }

    // Cleanup temp genome files and indices
    for genome_path in genome_files.values() {
        // Remove FASTA file
        let _ = std::fs::remove_file(genome_path);

        // Remove GIX files based on FASTA path (GDB is now embedded in .1aln)
        let gdb_base = genome_path.with_extension("");
        let _ = std::fs::remove_file(format!("{}.gix", gdb_base.display()));
        let _ = std::fs::remove_file(format!("{}.ktab", gdb_base.display()));
    }

    let _ = std::fs::remove_dir(&temp_genome_dir);

    Ok(merged_paf)
}

/// Create FastGA integration with optional frequency parameter
fn create_fastga_integration(
    frequency: Option<usize>,
    num_threads: usize,
    min_alignment_length: u64,
    temp_dir: Option<String>,
) -> Result<fastga_integration::FastGAIntegration> {
    use crate::fastga_integration::FastGAIntegration;
    Ok(FastGAIntegration::new(
        frequency,
        num_threads,
        min_alignment_length,
        temp_dir,
    ))
}

fn main() -> Result<()> {
    // Set OUT_DIR to help fastga-rs find its binaries
    // This ensures FAtoGDB, GIXmake, and GIXrm are found when FastGA needs them
    if let Some(cargo_home) = std::env::var("CARGO_HOME")
        .ok()
        .or_else(|| std::env::var("HOME").ok().map(|h| format!("{h}/.cargo")))
    {
        let lib_dir = format!("{cargo_home}/lib/sweepga");
        if std::path::Path::new(&lib_dir).exists() {
            std::env::set_var("OUT_DIR", &lib_dir);
        }
    }

    let args = Args::parse();

    // Handle --check-fastga diagnostic flag
    if args.check_fastga {
        println!("=== FastGA Binary Locations ===\n");

        let bins = ["FastGA", "ALNtoPAF", "PAFtoALN"];
        let mut all_found = true;

        for bin in &bins {
            match binary_paths::get_embedded_binary_path(bin) {
                Ok(path) => {
                    println!("✓ {}: {}", bin, path.display());
                    if path.exists() {
                        let metadata = std::fs::metadata(&path)?;
                        println!("  Size: {} bytes", metadata.len());
                        #[cfg(unix)]
                        {
                            use std::os::unix::fs::PermissionsExt;
                            println!(
                                "  Executable: {}",
                                metadata.permissions().mode() & 0o111 != 0
                            );
                        }

                        let path_str = path.to_string_lossy();
                        if path_str.contains("/build/fastga-rs-") {
                            println!("  Location: Embedded (development build)");
                        } else if path_str.contains("/.cargo/lib/") {
                            println!("  Location: Installed (cargo install)");
                        } else {
                            println!("  Location: System PATH");
                        }
                    }
                }
                Err(e) => {
                    println!("✗ {bin}: NOT FOUND");
                    println!("  Error: {e}");
                    all_found = false;
                }
            }
            println!();
        }

        if all_found {
            println!("All FastGA binaries found and ready to use!");
            std::process::exit(0);
        } else {
            println!("Some FastGA binaries are missing. Please run:");
            println!("  cargo build --release");
            println!("  ./install.sh");
            std::process::exit(1);
        }
    }

    // Show help if no arguments provided and stdin is a terminal (interactive mode)
    if args.files.is_empty() {
        use std::io::IsTerminal;
        if std::io::stdin().is_terminal() {
            Args::parse_from(["sweepga", "-h"]);
        }
    }

    let timing = TimingContext::new();

    // Print startup banner
    if !args.quiet {
        use chrono::Local;
        let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");
        let cmd_line: Vec<String> = std::env::args().collect();
        timing.log("start", &format!("{} | {}", timestamp, cmd_line.join(" ")));
    }

    // Set TMPDIR if --tempdir is specified (FastGA uses this via Rust's tempfile crate)
    if let Some(ref tempdir) = args.tempdir {
        std::env::set_var("TMPDIR", tempdir);
        if !args.quiet {
            timing.log("tempdir", &format!("Using temp directory: {tempdir}"));
        }
    }

    // Track alignment time separately
    let mut alignment_time: Option<f64> = None;

    // Detect file types early so we can use them for output conversion
    let input_file_types: Vec<FileType> = if !args.files.is_empty() {
        let file_types: Result<Vec<FileType>> =
            args.files.iter().map(|f| detect_file_type(f)).collect();
        file_types?
    } else {
        vec![]
    };

    // Enable .1aln workflow when:
    // - User explicitly requested .1aln output with --1aln flag
    // - Input is .1aln or FASTA (not PAF)
    // - OR output file extension is .1aln
    // GDB is now preserved using AlnWriter::create_with_gdb (via open_write_from)
    // When using stdin (no files), treat as PAF input
    let input_is_paf = args.files.is_empty()
        || (!input_file_types.is_empty() && input_file_types.iter().all(|ft| *ft == FileType::Paf));
    let want_1aln_output = args.output_1aln
        || args
            .output_file
            .as_ref()
            .is_some_and(|f| f.ends_with(".1aln"));
    let use_1aln_workflow = !input_is_paf && want_1aln_output;

    if use_1aln_workflow {
        // PURE .1ALN WORKFLOW - FastGA produces .1aln, filter as .1aln, output .1aln
        if !args.quiet {
            timing.log("detect", "Using .1aln workflow (FastGA native format)");
        }

        // Step 1: Get .1aln input (either from file or from FastGA alignment)
        let (_temp_1aln, aln_input_path) =
            if !input_file_types.is_empty() && input_file_types[0] == FileType::Aln {
                // Input is already .1aln
                if !args.quiet {
                    timing.log("detect", &format!("Input: {} (.1aln)", args.files[0]));
                }
                (None, args.files[0].clone())
            } else if !input_file_types.is_empty() && input_file_types[0] == FileType::Fasta {
                // FASTA input - use FastGA to produce .1aln
                if args.files.len() > 1 {
                    // eprintln!(
                    //     "[sweepga] ERROR: Multiple FASTA files not supported in .1aln workflow yet"
                    // );
                    // eprintln!("[sweepga] Use --paf flag to use PAF workflow for multiple files");
                    std::process::exit(1);
                }

                let alignment_start = Instant::now();
                // Canonicalize path to avoid empty parent directory issues in fastga-rs
                let path = std::fs::canonicalize(&args.files[0])
                    .with_context(|| format!("Failed to resolve path: {}", args.files[0]))?;
                let fastga = create_fastga_integration(
                    args.frequency,
                    args.threads,
                    args.block_length,
                    args.tempdir.clone(),
                )?;

                if !args.quiet {
                    timing.log("align", &format!("Running FastGA on {}", args.files[0]));
                }

                let temp_1aln = fastga.align_to_temp_1aln(&path, &path)?;
                alignment_time = Some(alignment_start.elapsed().as_secs_f64());

                if !args.quiet {
                    timing.log(
                        "align",
                        &format!("FastGA complete ({:.1}s)", alignment_time.unwrap()),
                    );
                }

                let aln_path = temp_1aln.path().to_string_lossy().into_owned();

                // Track disk usage of alignment file
                disk_usage::track_file_created(&aln_path);
                // Also track .1gdb if it exists
                let gdb_path = aln_path.replace(".1aln", ".1gdb");
                if std::path::Path::new(&gdb_path).exists() {
                    disk_usage::track_file_created(&gdb_path);
                }

                (Some(temp_1aln), aln_path)
            } else {
                // eprintln!("[sweepga] ERROR: No valid input provided");
                std::process::exit(1);
            };

        // Step 2: Parse filter config
        let (plane_sweep_mode, plane_sweep_query_limit, plane_sweep_target_limit) =
            parse_filter_mode(&args.num_mappings, "plane sweep");
        let (scaffold_filter_mode, scaffold_max_per_query, scaffold_max_per_target) =
            parse_filter_mode(&args.scaffold_filter, "scaffold");
        let scoring_function = match args.scoring.as_str() {
            "ani" | "identity" => ScoringFunction::Identity,
            "length" => ScoringFunction::Length,
            "length-ani" | "length-identity" => ScoringFunction::LengthIdentity,
            "matches" => ScoringFunction::Matches,
            "log-length-ani" | "log-length-identity" => ScoringFunction::LogLengthIdentity,
            _ => ScoringFunction::LogLengthIdentity,
        };

        // Parse sparsification strategy
        use tree_sparsify::SparsificationStrategy;
        let sparsification_strategy = args.sparsify.parse::<SparsificationStrategy>()?;
        let sparsity_fraction = match sparsification_strategy {
            SparsificationStrategy::Fraction(f) => f,
            SparsificationStrategy::Tree(_, _, _) => 1.0, // Tree filtering handled separately
        };

        let filter_config = FilterConfig {
            chain_gap: args.scaffold_jump,
            min_block_length: args.block_length,
            mapping_filter_mode: plane_sweep_mode,
            mapping_max_per_query: plane_sweep_query_limit,
            mapping_max_per_target: plane_sweep_target_limit,
            plane_sweep_secondaries: 0,
            scaffold_filter_mode,
            scaffold_max_per_query,
            scaffold_max_per_target,
            overlap_threshold: args.overlap,
            sparsity: sparsity_fraction,
            no_merge: true,
            scaffold_gap: args.scaffold_jump,
            min_scaffold_length: args.scaffold_mass,
            scaffold_overlap_threshold: args.scaffold_overlap,
            scaffold_max_deviation: args.scaffold_dist,
            prefix_delimiter: '#',
            skip_prefix: false,
            scoring_function,
            min_identity: 0.0,
            min_scaffold_identity: 0.0,
        };

        // Step 2.6: Apply tree filtering if requested (natively on .1aln format)
        let tree_filtered_input =
            if let SparsificationStrategy::Tree(k_nearest, k_farthest, random_fraction) =
                sparsification_strategy
            {
                // Apply tree filtering directly to .1aln
                let temp_tree_filtered = tempfile::NamedTempFile::with_suffix(".1aln")?;

                use tree_sparsify::apply_tree_filter_to_1aln;
                apply_tree_filter_to_1aln(
                    &aln_input_path,
                    temp_tree_filtered
                        .path()
                        .to_str()
                        .context("Invalid temp path")?,
                    k_nearest,
                    k_farthest,
                    random_fraction,
                    args.quiet,
                )?;

                // Use tree-filtered output as input for plane sweep filtering
                Some(temp_tree_filtered)
            } else {
                None
            };

        let final_filter_input = if let Some(ref tf) = tree_filtered_input {
            tf.path().to_path_buf()
        } else {
            std::path::PathBuf::from(&aln_input_path)
        };

        // Step 3: Filter .1aln directly using unified_filter (format-preserving)
        use crate::unified_filter::filter_file;
        if let Some(ref output_file) = args.output_file {
            filter_file(
                final_filter_input,
                output_file,
                &filter_config,
                false,
                args.keep_self,
            )?;
        } else {
            // Write to temp file then copy to stdout
            let temp_output = tempfile::NamedTempFile::with_suffix(".1aln")?;
            filter_file(
                final_filter_input,
                temp_output.path(),
                &filter_config,
                false,
                args.keep_self,
            )?;

            // Copy binary to stdout
            use std::io::copy;
            let mut file = std::fs::File::open(temp_output.path())?;
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();
            copy(&mut file, &mut handle)?;
        }

        if !args.quiet {
            let (total_elapsed, _) = timing.stats();
            if let Some(align_time) = alignment_time {
                let filtering_time = total_elapsed - align_time;
                timing.log("done", &format!("Alignment: {align_time:.1}s, Filtering: {filtering_time:.1}s, Total: {total_elapsed:.1}s"));
            } else {
                timing.log("done", &format!("Total: {total_elapsed:.1}s"));
            }
        }
        return Ok(());
    }

    // Detect file types and route accordingly
    let (_temp_paf, input_path) = if !args.files.is_empty() {
        let file_types = &input_file_types;

        if !args.quiet {
            for (i, (file, ftype)) in args.files.iter().zip(file_types.iter()).enumerate() {
                timing.log(
                    "detect",
                    &format!("Input {}: {} ({:?})", i + 1, file, ftype),
                );
            }
        }

        // Check if all files are FASTA
        let all_fasta = file_types.iter().all(|ft| *ft == FileType::Fasta);
        let fasta_count = file_types
            .iter()
            .filter(|ft| **ft == FileType::Fasta)
            .count();

        match (args.files.len(), all_fasta) {
            (1, true) => {
                // Single FASTA - check if it has multiple genomes
                // Canonicalize path to avoid empty parent directory issues in fastga-rs
                let path = std::fs::canonicalize(&args.files[0])
                    .with_context(|| format!("Failed to resolve path: {}", args.files[0]))?;
                let groups = detect_genome_groups(&path)?;

                if groups.len() > 1 {
                    // Multiple genomes in single FASTA - use pairwise alignment
                    if !args.quiet {
                        timing.log(
                            "detect",
                            &format!(
                                "Detected {} genomes in {}, using pairwise alignment",
                                groups.len(),
                                args.files[0]
                            ),
                        );
                    }

                    let alignment_start = Instant::now();
                    let temp_paf = align_multiple_fastas(
                        &args.files,
                        args.frequency,
                        args.all_pairs,
                        args.batch_bytes,
                        args.threads,
                        args.keep_self,
                        args.tempdir.as_deref(),
                        &timing,
                        args.quiet,
                        args.block_length,
                        args.zstd_compress,
                        args.zstd_level,
                    )?;

                    alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "All pairwise alignments complete ({:.1}s)",
                                alignment_time.unwrap()
                            ),
                        );
                    }

                    let paf_path = temp_paf.path().to_string_lossy().into_owned();
                    (Some(temp_paf), paf_path)
                } else {
                    // Single genome - self-alignment
                    let alignment_start = Instant::now();

                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "Running {} self-alignment on {}",
                                args.aligner, args.files[0]
                            ),
                        );
                    }

                    let fastga = create_fastga_integration(
                        args.frequency,
                        args.threads,
                        args.block_length,
                        args.tempdir.clone(),
                    )?;
                    let temp_paf = fastga.align_to_temp_paf(&path, &path)?;

                    alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "{} alignment complete ({:.1}s)",
                                args.aligner,
                                alignment_time.unwrap()
                            ),
                        );
                    }

                    let paf_path = temp_paf.path().to_string_lossy().into_owned();
                    (Some(temp_paf), paf_path)
                }
            }
            (2, true) => {
                // Two FASTAs - check if they're the same or different
                let groups1 = detect_genome_groups(Path::new(&args.files[0]))?;
                let groups2 = detect_genome_groups(Path::new(&args.files[1]))?;

                if groups1.len() > 1 || groups2.len() > 1 {
                    // Multiple genomes across files - use pairwise alignment
                    if !args.quiet {
                        timing.log(
                            "detect",
                            &format!(
                                "Detected {} + {} genomes, using pairwise alignment",
                                groups1.len(),
                                groups2.len()
                            ),
                        );
                    }

                    let alignment_start = Instant::now();
                    let temp_paf = align_multiple_fastas(
                        &args.files,
                        args.frequency,
                        args.all_pairs,
                        args.batch_bytes,
                        args.threads,
                        args.keep_self,
                        args.tempdir.as_deref(),
                        &timing,
                        args.quiet,
                        args.block_length,
                        args.zstd_compress,
                        args.zstd_level,
                    )?;

                    alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "All pairwise alignments complete ({:.1}s)",
                                alignment_time.unwrap()
                            ),
                        );
                    }

                    let paf_path = temp_paf.path().to_string_lossy().into_owned();
                    (Some(temp_paf), paf_path)
                } else {
                    // Simple pairwise alignment (legacy behavior)
                    let alignment_start = Instant::now();

                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "Running {} alignment: {} -> {}",
                                args.aligner, args.files[0], args.files[1]
                            ),
                        );
                    }

                    // Canonicalize paths to avoid empty parent directory issues in fastga-rs
                    let target = std::fs::canonicalize(&args.files[0])
                        .with_context(|| format!("Failed to resolve path: {}", args.files[0]))?;
                    let query = std::fs::canonicalize(&args.files[1])
                        .with_context(|| format!("Failed to resolve path: {}", args.files[1]))?;
                    let fastga = create_fastga_integration(
                        args.frequency,
                        args.threads,
                        args.block_length,
                        args.tempdir.clone(),
                    )?;
                    let temp_paf = fastga.align_to_temp_paf(&target, &query)?;

                    alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                    if !args.quiet {
                        timing.log(
                            "align",
                            &format!(
                                "{} alignment complete ({:.1}s)",
                                args.aligner,
                                alignment_time.unwrap()
                            ),
                        );
                    }

                    let paf_path = temp_paf.path().to_string_lossy().into_owned();
                    (Some(temp_paf), paf_path)
                }
            }
            (n, true) if n > 2 => {
                // Multiple FASTAs - always use pairwise alignment
                if !args.quiet {
                    timing.log(
                        "detect",
                        &format!(
                            "{fasta_count} FASTA files provided, using all-pairs pairwise alignment"
                        ),
                    );
                }

                let alignment_start = Instant::now();
                let temp_paf = align_multiple_fastas(
                    &args.files,
                    args.frequency,
                    args.all_pairs,
                    args.batch_bytes,
                    args.threads,
                    args.keep_self,
                    args.tempdir.as_deref(),
                    &timing,
                    args.quiet,
                    args.block_length,
                    args.zstd_compress,
                    args.zstd_level,
                )?;

                alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                if !args.quiet {
                    timing.log(
                        "align",
                        &format!(
                            "All pairwise alignments complete ({:.1}s)",
                            alignment_time.unwrap()
                        ),
                    );
                }

                let paf_path = temp_paf.path().to_string_lossy().into_owned();
                (Some(temp_paf), paf_path)
            }
            (1, false) if file_types[0] == FileType::Paf => {
                // Filter existing PAF
                (None, args.files[0].clone())
            }
            (1, false) if file_types[0] == FileType::Aln => {
                // Convert .1aln to PAF for filtering
                if !args.quiet {
                    timing.log(
                        "convert",
                        &format!("Converting .1aln to PAF: {}", args.files[0]),
                    );
                }

                let temp_paf = aln_to_paf(&args.files[0], args.threads)?;
                let paf_path = temp_paf.path().to_str().unwrap().to_string();

                if !args.quiet {
                    timing.log("convert", "Conversion complete");
                }

                (Some(temp_paf), paf_path)
            }
            (1, false) if file_types[0] == FileType::Agc => {
                // AGC archive - extract samples and align
                let alignment_start = Instant::now();

                let agc_tempdir = args.agc_tempdir.as_deref().or(args.tempdir.as_deref());

                let temp_paf = process_agc_archive(
                    &args.files[0],
                    &args,
                    agc_tempdir,
                    &timing,
                )?;

                alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                if !args.quiet {
                    timing.log(
                        "align",
                        &format!(
                            "AGC alignment complete ({:.1}s)",
                            alignment_time.unwrap()
                        ),
                    );
                }

                let paf_path = temp_paf.path().to_string_lossy().into_owned();
                (Some(temp_paf), paf_path)
            }
            _ => {
                anyhow::bail!("Invalid file combination: expected FASTA file(s) for alignment, 1 PAF for filtering, 1 .1aln for filtering, or 1 .agc archive");
            }
        }
    } else {
        // Read from stdin - save to temp file for auto-detection and two-pass processing
        let temp = tempfile::NamedTempFile::new()?;
        let temp_path = temp.path().to_str().unwrap().to_string();

        {
            use std::io::{self, BufRead, Write};
            let stdin = io::stdin();
            let mut temp_file = std::fs::File::create(&temp_path)?;

            for line in stdin.lock().lines() {
                writeln!(temp_file, "{}", line?)?;
            }
            temp_file.flush()?;
        }

        // Check if stdin had any content
        let file_size = std::fs::metadata(&temp_path)?.len();
        if file_size == 0 && !args.no_filter {
            use clap::CommandFactory;
            Args::command().print_help()?;
            std::process::exit(0);
        }

        // Detect stdin type
        let file_type = detect_file_type(&temp_path)?;
        if !args.quiet {
            timing.log("detect", &format!("stdin: {file_type:?}"));
        }

        match file_type {
            FileType::Fasta => {
                // Self-align stdin FASTA
                let alignment_start = Instant::now();

                if !args.quiet {
                    timing.log(
                        "align",
                        &format!("Running {} self-alignment on stdin", args.aligner),
                    );
                }

                // Canonicalize path to avoid empty parent directory issues in fastga-rs
                let path = std::fs::canonicalize(&temp_path)
                    .with_context(|| format!("Failed to resolve path: {temp_path}"))?;
                let fastga = create_fastga_integration(
                    args.frequency,
                    args.threads,
                    args.block_length,
                    args.tempdir.clone(),
                )?;
                let temp_paf = fastga.align_to_temp_paf(&path, &path)?;

                alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                if !args.quiet {
                    timing.log(
                        "align",
                        &format!(
                            "{} alignment complete ({:.1}s)",
                            args.aligner,
                            alignment_time.unwrap()
                        ),
                    );
                }

                let paf_path = temp_paf.path().to_string_lossy().into_owned();
                // Keep both temp files alive
                Box::leak(Box::new(temp));
                (Some(temp_paf), paf_path)
            }
            FileType::Paf => {
                // Filter stdin PAF
                (Some(temp), temp_path)
            }
            FileType::Aln => {
                // Convert stdin .1aln to PAF for filtering
                if !args.quiet {
                    timing.log("convert", "Converting .1aln from stdin to PAF");
                }

                let temp_paf = aln_to_paf(&temp_path, args.threads)?;
                let paf_path = temp_paf.path().to_str().unwrap().to_string();

                if !args.quiet {
                    timing.log("convert", "Conversion complete");
                }

                // Keep both temp files alive
                Box::leak(Box::new(temp));
                (Some(temp_paf), paf_path)
            }
            FileType::Agc => {
                anyhow::bail!("AGC archives cannot be read from stdin - please provide a file path");
            }
        }
    };

    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // Handle no-filter mode - just copy input to stdout
    if args.no_filter {
        use std::io::{BufRead, Write};

        let input = crate::paf::open_paf_input(&input_path)?;
        let stdout = std::io::stdout();
        let mut output = stdout.lock();

        for line in input.lines() {
            writeln!(output, "{}", line?)?;
        }

        return Ok(());
    }

    // Parse -n parameter for plane sweep
    let (plane_sweep_mode, plane_sweep_query_limit, plane_sweep_target_limit) =
        parse_filter_mode(&args.num_mappings, "plane sweep");

    // Parse scaffold filter mode
    let (scaffold_filter_mode, scaffold_max_per_query, scaffold_max_per_target) =
        parse_filter_mode(&args.scaffold_filter, "scaffold");

    // Parse scoring function
    let scoring_function = match args.scoring.as_str() {
        "ani" | "identity" => ScoringFunction::Identity,
        "length" => ScoringFunction::Length,
        "length-ani" | "length-identity" => ScoringFunction::LengthIdentity,
        "matches" => ScoringFunction::Matches,
        "log-length-ani" | "log-length-identity" => ScoringFunction::LogLengthIdentity,
        _ => ScoringFunction::LogLengthIdentity,
    };

    // Parse sparsification strategy (PAF workflow)
    use tree_sparsify::SparsificationStrategy;
    let sparsification_strategy = args.sparsify.parse::<SparsificationStrategy>()?;
    let sparsity_fraction = match sparsification_strategy {
        SparsificationStrategy::Fraction(f) => f,
        SparsificationStrategy::Tree(_, _, _) => 1.0, // Tree filtering will be applied to PAF
    };

    // Placeholder for identity values - will be calculated after we have input path
    let temp_config = FilterConfig {
        chain_gap: args.scaffold_jump,
        min_block_length: args.block_length,
        mapping_filter_mode: plane_sweep_mode,
        mapping_max_per_query: plane_sweep_query_limit,
        mapping_max_per_target: plane_sweep_target_limit,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode,
        scaffold_max_per_query,
        scaffold_max_per_target,
        overlap_threshold: args.overlap,
        sparsity: sparsity_fraction,
        no_merge: true,
        scaffold_gap: args.scaffold_jump,
        min_scaffold_length: args.scaffold_mass,
        scaffold_overlap_threshold: args.scaffold_overlap,
        scaffold_max_deviation: args.scaffold_dist,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function,
        min_identity: 0.0,          // Will be set later
        min_scaffold_identity: 0.0, // Will be set later
    };

    // Parse ANI calculation method
    let ani_method = parse_ani_method(&args.ani_method);
    if ani_method.is_none() {
        // eprintln!(
        //     "[sweepga] WARNING: Unknown ANI method '{}', using 'n50-identity'",
        //     args.ani_method
        // );
    }
    let ani_method = ani_method.unwrap_or(AniMethod::NPercentile(50.0, NSort::Identity));

    // Now calculate ANI if needed for identity thresholds
    let ani_percentile = if args.min_identity.to_lowercase().contains("ani")
        || args.min_scaffold_identity.to_lowercase().contains("ani")
    {
        Some(calculate_ani_stats(&input_path, ani_method, args.quiet)?)
    } else {
        None
    };

    // Parse identity thresholds
    let min_identity = parse_identity_value(&args.min_identity, ani_percentile)?;
    let min_scaffold_identity = if args.min_scaffold_identity.is_empty() {
        min_identity // If empty string, use min_identity
    } else {
        parse_identity_value(&args.min_scaffold_identity, ani_percentile)?
    };

    // Only report thresholds if they're non-zero
    if !args.quiet && (min_identity > 0.0 || min_scaffold_identity > 0.0) {
        if min_identity > 0.0 {
            timing.log(
                "config",
                &format!("Mapping identity threshold: {:.1}%", min_identity * 100.0),
            );
        }
        if min_scaffold_identity > 0.0 && min_scaffold_identity != min_identity {
            timing.log(
                "config",
                &format!(
                    "Scaffold identity threshold: {:.1}%",
                    min_scaffold_identity * 100.0
                ),
            );
        }
    }

    // Create final config with calculated identity values
    let mut config = temp_config;
    config.min_identity = min_identity;
    config.min_scaffold_identity = min_scaffold_identity;

    // Determine output format based on flags and file extensions
    // Default: PAF output (unless --1aln flag or .1aln extension)
    let output_1aln = if let Some(ref outfile) = args.output_file {
        // Auto-detect from output file extension (.1aln → binary, .paf → text)
        outfile.ends_with(".1aln")
    } else {
        // Use --1aln flag (default is false = PAF output)
        args.output_1aln
    };

    // Apply filtering - always to temp file, then copy to stdout
    if !args.quiet {
        timing.log("parse", &format!("Parsing input PAF: {input_path}"));
    }

    let output_temp = tempfile::NamedTempFile::with_suffix(".paf")?;
    let (output_file, output_path_buf) = output_temp.keep()?;
    let output_path = output_path_buf.to_str().unwrap().to_string();
    drop(output_file); // Close the file handle so filter_paf can open it

    // Apply tree-based sparsification if requested
    let tree_filtered_path = if let SparsificationStrategy::Tree(k_nearest, k_farthest, rand_frac) =
        sparsification_strategy
    {
        if !args.quiet {
            timing.log(
                "tree-filter",
                &format!(
                    "Applying tree sparsification: tree:{}{}{}",
                    k_nearest,
                    if k_farthest > 0 {
                        format!(",{k_farthest}")
                    } else {
                        String::new()
                    },
                    if rand_frac > 0.0 {
                        format!(",{rand_frac}")
                    } else {
                        String::new()
                    }
                ),
            );
        }

        // Create temporary file for tree-filtered PAF
        let tree_temp = tempfile::NamedTempFile::with_suffix(".tree.paf")?;
        let (tree_file, tree_path_buf) = tree_temp.keep()?;
        let tree_path = tree_path_buf.to_str().unwrap().to_string();
        drop(tree_file); // Close so apply_tree_filter_to_paf can open it

        // Apply tree filtering
        tree_sparsify::apply_tree_filter_to_paf(
            &input_path,
            &tree_path,
            k_nearest,
            k_farthest,
            rand_frac,
        )?;

        Some(tree_path)
    } else {
        None
    };

    // Use tree-filtered PAF if available, otherwise use original input
    let filter_input_path = tree_filtered_path.as_ref().unwrap_or(&input_path);

    // Note: -f (no_filter) implies --self (keep self-mappings)
    let filter = PafFilter::new(config)
        .with_keep_self(args.keep_self || args.no_filter)
        .with_scaffolds_only(args.scaffolds_only);
    filter.filter_paf(filter_input_path, &output_path)?;

    // Convert output format if requested
    // Note: PAFtoALN automatically appends .paf to the input filename, so we need to strip it
    let paftoaln_input_path = if output_1aln {
        output_path
            .strip_suffix(".paf")
            .unwrap_or(&output_path)
            .to_string()
    } else {
        output_path.clone()
    };

    let (final_output_path, _aln_temp) = if output_1aln {
        // Convert PAF to 1aln using PAFtoALN
        // Requires the original FASTA or .1aln file(s) for sequence metadata
        if input_file_types.is_empty() {
            if !args.quiet {
                // eprintln!(
                //     "[sweepga] WARNING: .1aln output requires FASTA or .1aln input (not PAF), falling back to PAF output"
                // );
            }
            (output_path.clone(), None::<tempfile::NamedTempFile>)
        } else {
            let first_type = input_file_types[0];
            if first_type != FileType::Fasta && first_type != FileType::Aln {
                if !args.quiet {
                    // eprintln!(
                    //     "[sweepga] WARNING: .1aln output requires FASTA or .1aln input (not PAF), falling back to PAF output"
                    // );
                }
                (output_path.clone(), None::<tempfile::NamedTempFile>)
            } else {
                // Convert filtered PAF to .1aln using PAFtoALN
                // Find PAFtoALN binary (embedded)
                let paftoaln_bin = sweepga::binary_paths::get_embedded_binary_path("PAFtoALN")
                    .map_err(|_| anyhow::anyhow!("PAFtoALN tool not found. Cannot convert PAF to .1aln. Use --paf flag to output PAF format instead."))?;

                if !args.quiet {
                    timing.log("convert", "Converting filtered PAF to .1aln");
                }

                // PAFtoALN usage: PAFtoALN [-T<int>] <alignments> <source1.fa|source1.1gdb> [<source2.fa|source2.1gdb>]
                // Note: PAFtoALN automatically appends .paf to the alignments filename for input
                // and creates <alignments>.1aln for output
                let mut cmd = std::process::Command::new(&paftoaln_bin);
                cmd.arg(format!("-T{}", args.threads))
                    .arg(&paftoaln_input_path);

                // Pass FASTA or .1gdb file(s) for sequence metadata
                // If input is .1aln, pass matching .1gdb; if FASTA, pass FASTA
                for input_file in &args.files {
                    if input_file.ends_with(".1aln") {
                        // Input is .1aln - pass matching .1gdb
                        let gdb_path = input_file
                            .strip_suffix(".1aln")
                            .unwrap_or(input_file)
                            .to_string()
                            + ".1gdb";
                        cmd.arg(&gdb_path);
                    } else {
                        // Input is FASTA - pass it directly
                        cmd.arg(input_file);
                    }
                }

                // Suppress PAFtoALN's normal output unless in verbose mode
                if args.quiet {
                    cmd.stdout(std::process::Stdio::null())
                        .stderr(std::process::Stdio::null());
                }

                let status = cmd.status()?;
                if !status.success() {
                    anyhow::bail!("PAFtoALN conversion failed");
                }

                // PAFtoALN creates the output file as <input_path>.1aln
                let aln_path = format!("{paftoaln_input_path}.1aln");

                if !args.quiet {
                    timing.log("convert", "PAF to .1aln conversion complete");
                }

                (aln_path, None::<tempfile::NamedTempFile>)
            }
        }
    } else {
        (output_path.clone(), None::<tempfile::NamedTempFile>)
    };

    // Handle output based on format and destination
    use std::io::{BufRead, BufReader, Write};

    // Write output (PAF or .1aln) to file or stdout
    if let Some(output_file) = &args.output_file {
        std::fs::copy(&final_output_path, output_file)?;
    } else if output_1aln {
        // .1aln is binary - copy bytes directly to stdout
        use std::io::copy;
        let mut file = std::fs::File::open(&final_output_path)?;
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        copy(&mut file, &mut handle)?;
    } else {
        // PAF is text - write line by line to stdout
        let file = std::fs::File::open(&final_output_path)?;
        let reader = BufReader::new(file);
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();

        for line in reader.lines() {
            writeln!(handle, "{}", line?)?;
        }
    }

    // Clean up temp files
    let _ = std::fs::remove_file(&output_path);
    if output_1aln {
        let _ = std::fs::remove_file(&final_output_path);
    }

    if !args.quiet {
        let (total_elapsed, _) = timing.stats();
        let filtering_time = if let Some(align_time) = alignment_time {
            total_elapsed - align_time
        } else {
            total_elapsed
        };

        // Format the complete message
        if let Some(align_time) = alignment_time {
            timing.log("done", &format!("Alignment: {align_time:.1}s, Filtering: {filtering_time:.1}s, Total: {total_elapsed:.1}s"));
        } else {
            timing.log("done", &format!("Total: {total_elapsed:.1}s"));
        }
    }

    // Report disk usage if requested
    if args.disk_usage {
        disk_usage::log_summary();
    }

    Ok(())
}
