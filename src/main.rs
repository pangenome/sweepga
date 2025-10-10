mod aln_filter;
mod compact_mapping;
mod fastga_integration;
mod grouped_mappings;
mod mapping;
mod paf;
mod paf_filter;
mod plane_sweep;
mod plane_sweep_core;
mod plane_sweep_exact;
mod sequence_index;
mod unified_filter;
mod union_find;

use anyhow::Result;
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

    /// Log message with timing in minimap2 format
    fn log(&self, phase: &str, message: &str) {
        let (elapsed, cpu_ratio) = self.stats();
        eprintln!("[sweepga::{phase}::{elapsed:.3}*{cpu_ratio:.2}] {message}");
    }
}

/// File type detected from content
#[derive(Debug, Clone, Copy, PartialEq)]
enum FileType {
    Fasta,
    Paf,
    Aln, // .1aln binary format
}

/// Detect file type by reading first non-empty line or checking file extension
/// Handles .gz files automatically
fn detect_file_type(path: &str) -> Result<FileType> {
    // Check for .1aln extension first (binary format)
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
            anyhow::bail!("Empty file: {}", path);
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

    anyhow::bail!("Could not detect file type for {}: not FASTA (starts with >), PAF (12+ tab-delimited fields), or .1aln (binary)", path);
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
    /// PAF (filter) or FASTA (align & filter). Multiple FASTAs for all-pairs alignment
    #[clap(value_name = "FILE", num_args = 0..,
           long_help = "Input files: FASTA (1+) or PAF (1 only), auto-detected\n\
                        \n  \
                        1 FASTA: align to self and filter\n  \
                        2+ FASTA: align all pairs and filter\n  \
                        1 PAF: filter alignments\n  \
                        stdin: auto-detect and process")]
    files: Vec<String>,

    /// Aligner to use for FASTA input
    #[clap(long = "aligner", default_value = "fastga", value_parser = ["fastga"])]
    aligner: String,

    /// FastGA k-mer frequency threshold (use k-mers occurring ≤ N times)
    #[clap(short = 'f', long = "frequency")]
    frequency: Option<usize>,

    /// Align all genome pairs separately (slower, uses more memory, but handles many genomes)
    #[clap(long = "all-pairs")]
    all_pairs: bool,

    /// Minimum block length
    #[clap(short = 'b', long = "block-length", default_value = "0", value_parser = parse_metric_number)]
    block_length: u64,

    /// Maximum overlap ratio for plane sweep filtering
    #[clap(short = 'o', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Keep this fraction of mappings
    #[clap(short = 'x', long = "sparsify", default_value = "1.0")]
    sparsify: f64,

    /// n:m-best mappings kept in query:target dimensions. 1:1 (orthogonal), use ∞/many for unbounded
    #[clap(short = 'n', long = "num-mappings", default_value = "1:1")]
    num_mappings: String,

    /// Scaffold filter: "1:1" (best), "M:N" (top M per query, N per target; ∞/many for unbounded)
    #[clap(
        short = 'm',
        long = "scaffold-filter",
        default_value = "many:many",
        hide = true
    )]
    scaffold_filter: String,

    /// Scaffold jump (gap) distance. 0 = disable scaffolding (plane sweep only), >0 = enable scaffolding
    #[clap(short = 'j', long = "scaffold-jump", default_value = "0", value_parser = parse_metric_number, hide = true)]
    scaffold_jump: u64,

    /// Minimum scaffold length when scaffolding is enabled
    #[clap(short = 's', long = "scaffold-mass", default_value = "10k", value_parser = parse_metric_number, hide = true)]
    scaffold_mass: u64,

    /// Scaffold chain overlap threshold
    #[clap(
        short = 'O',
        long = "scaffold-overlap",
        default_value = "0.5",
        hide = true
    )]
    scaffold_overlap: f64,

    /// Maximum distance from scaffold anchor (0 = no rescue, only keep scaffold members)
    #[clap(short = 'd', long = "scaffold-dist", default_value = "20k", value_parser = parse_metric_number, hide = true)]
    scaffold_dist: u64,

    /// Scoring function for plane sweep
    #[clap(long = "scoring", default_value = "log-length-ani",
           value_parser = ["ani", "length", "length-ani", "log-length-ani", "matches"])]
    scoring: String,

    /// Method for calculating ANI: all, orthogonal, nX[-sort] (e.g. n50, n90-identity, n100-score)
    #[clap(long = "ani-method", default_value = "n100")]
    ani_method: String,

    /// Minimum identity threshold (0-1 fraction, 1-100%, or "aniN" for Nth percentile)
    #[clap(short = 'i', long = "min-identity", default_value = "0")]
    min_identity: String,

    /// Minimum scaffold identity threshold (0-1 fraction, 1-100%, "aniN", or defaults to -i)
    #[clap(
        short = 'Y',
        long = "min-scaffold-identity",
        default_value = "0",
        hide = true
    )]
    min_scaffold_identity: String,

    /// Disable all filtering
    #[clap(short = 'N', long = "no-filter")]
    no_filter: bool,

    /// Keep self-mappings (excluded by default)
    #[clap(long = "self")]
    keep_self: bool,

    /// Output scaffold chains only (for debugging)
    #[clap(long = "scaffolds-only", hide = true)]
    scaffolds_only: bool,

    /// Quiet mode (no progress output)
    #[clap(long = "quiet")]
    quiet: bool,

    /// Number of threads for parallel processing
    #[clap(short = 't', long = "threads", default_value = "8")]
    threads: usize,

    /// Output PAF format instead of default .1aln (text instead of binary)
    #[clap(long = "paf")]
    output_paf: bool,

    /// Output file path (auto-detects format from extension: .paf or .1aln)
    #[clap(long = "output-file")]
    output_file: Option<String>,

    /// Temporary directory for intermediate files (defaults to TMPDIR env var, then /tmp)
    #[clap(long = "tempdir")]
    tempdir: Option<String>,
}

fn parse_filter_mode(mode: &str, filter_type: &str) -> (FilterMode, Option<usize>, Option<usize>) {
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
            eprintln!("Warning: Invalid {filter_type} filter '{s}', using default 1:1");
            (FilterMode::OneToOne, Some(1), Some(1))
        }
        _ => {
            // Try parsing as a single number
            if let Ok(n) = mode.parse::<usize>() {
                if n == 0 {
                    eprintln!("Error: 0 is not a valid filter value. Use 1 for best mapping only.");
                    std::process::exit(1);
                }
                // Single number means filter on query axis only
                return (FilterMode::OneToMany, Some(n), None);
            }
            eprintln!("Warning: Invalid {filter_type} filter '{mode}', using default 1:1");
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
                eprintln!(
                    "[sweepga] WARNING: Only ani50 (median) currently supported, using median"
                );
            }

            // Apply offset if provided
            if let Some((sign, offset_str)) = offset_part {
                let offset: f64 = offset_str
                    .parse()
                    .map_err(|_| anyhow::anyhow!("Invalid ANI offset: {}", offset_str))?;
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
        anyhow::bail!("Invalid identity value: {}", value);
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
                eprintln!("[sweepga] Calculating ANI from best 1:1 mappings (min length: 1kb)");
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
        eprintln!("[sweepga] WARNING: No inter-genome alignments found for ANI calculation");
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

    eprintln!(
        "[sweepga] ANI statistics from {} genome pairs:",
        genome_pairs.len()
    );
    eprintln!(
        "[sweepga]   Min: {:.1}%, Median: {:.1}%, Max: {:.1}%",
        ani_values.first().unwrap_or(&0.0) * 100.0,
        ani50 * 100.0,
        ani_values.last().unwrap_or(&0.0) * 100.0
    );

    Ok(ani50)
}

/// Calculate ANI using N-percentile method - use best alignments covering N% of genome pairs
fn calculate_ani_n_percentile(
    input_path: &str,
    percentile: f64,
    sort_method: NSort,
    quiet: bool,
) -> Result<f64> {
    if !quiet {
        let method_name = match sort_method {
            NSort::Length => format!("longest alignments (N{} by length)", percentile as i32),
            NSort::Identity => format!(
                "highest identity alignments (N{} by identity)",
                percentile as i32
            ),
            NSort::Score => format!(
                "best scoring alignments (N{} by identity × log(length))",
                percentile as i32
            ),
        };
        eprintln!("[sweepga] Calculating ANI from {method_name}");
    }

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
        eprintln!("[sweepga] WARNING: No inter-genome alignments found for ANI calculation");
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
        eprintln!(
            "[sweepga] Total genome size: {:.1} Mb, N{} threshold: {:.1} Mb",
            total_genome_size / 1_000_000.0,
            percentile as i32,
            n_threshold / 1_000_000.0
        );
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

    // Calculate total alignment length included
    let total_included: f64 = genome_pairs.values().map(|(_, len)| len).sum();
    let coverage_pct = (total_included / total_genome_size) * 100.0;

    let method_str = match sort_method {
        NSort::Length => format!("N{} by length", percentile as i32),
        NSort::Identity => format!("N{} by identity", percentile as i32),
        NSort::Score => format!("N{} by score", percentile as i32),
    };
    eprintln!(
        "[sweepga] ANI statistics from {} genome pairs ({}):",
        genome_pairs.len(),
        method_str
    );
    eprintln!(
        "[sweepga]   Coverage: {:.1}% of genomes ({:.1} Mb aligned / {:.1} Mb sum of genome sizes)",
        coverage_pct,
        total_included / 1_000_000.0,
        total_genome_size / 1_000_000.0
    );
    eprintln!(
        "[sweepga]   Min: {:.1}%, Median: {:.1}%, Max: {:.1}%",
        ani_values.first().unwrap_or(&0.0) * 100.0,
        ani50 * 100.0,
        ani_values.last().unwrap_or(&0.0) * 100.0
    );

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

    // Fallback to ALNtoPAF
    let alnto_paf_bin = std::env::current_exe()
        .ok()
        .and_then(|exe| exe.parent().map(|p| p.join("ALNtoPAF")))
        .filter(|p| p.exists())
        .or_else(|| {
            let local = std::path::PathBuf::from("./ALNtoPAF");
            if local.exists() {
                Some(local)
            } else {
                None
            }
        })
        .or_else(|| {
            let deps = std::path::PathBuf::from("deps/fastga/ALNtoPAF");
            if deps.exists() {
                Some(deps)
            } else {
                None
            }
        })
        .ok_or_else(|| {
            anyhow::anyhow!(
                "ALNtoPAF tool not found and native reader failed. Cannot convert .1aln to PAF."
            )
        })?;

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
    threads: usize,
    keep_self: bool,
    tempdir: Option<&str>,
    timing: &TimingContext,
    quiet: bool,
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
        );
    }

    // Default mode: single FastGA run with adaptive frequency
    let effective_frequency = if let Some(f) = frequency {
        // User specified -f, use it as-is
        if !quiet {
            timing.log(
                "align",
                &format!("Using user-specified frequency threshold: {f}"),
            );
        }
        Some(f)
    } else if num_genomes > 1 {
        // Auto-set frequency to at least num_genomes to avoid filtering out valid k-mers
        let auto_freq = num_genomes;
        if !quiet {
            timing.log(
                "align",
                &format!("Setting frequency threshold to {auto_freq} (number of genome groups)"),
            );
        }
        Some(auto_freq)
    } else {
        // Single genome, use FastGA default
        frequency
    };

    if !quiet {
        timing.log(
            "align",
            &format!("Aligning {num_genomes} genomes in single FastGA run"),
        );
    }

    // Create FastGA integration with effective frequency
    let fastga = create_fastga_integration(effective_frequency, threads)?;

    // Run single FastGA alignment on the original input (all genomes together)
    // Use first file as both query and target for self-alignment
    let input_file = &fasta_files[0];
    let temp_paf = fastga.align_to_temp_paf(Path::new(input_file), Path::new(input_file))?;

    if !quiet {
        let alignment_count = std::fs::read_to_string(temp_paf.path())?.lines().count();
        timing.log(
            "align",
            &format!("FastGA produced {alignment_count} alignments"),
        );
    }

    Ok(temp_paf)
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
    let fastga = create_fastga_integration(frequency, threads)?;

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
        // We don't need to store the result - FastGA will auto-detect them by the FASTA path
        let _gdb_base = fastga.prepare_gdb(fasta_path)?;
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
) -> Result<fastga_integration::FastGAIntegration> {
    use crate::fastga_integration::FastGAIntegration;
    Ok(FastGAIntegration::new(frequency, num_threads))
}

fn main() -> Result<()> {
    let args = Args::parse();
    let timing = TimingContext::new();

    // Print startup banner
    if !args.quiet {
        use chrono::Local;
        let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");
        let cmd_line: Vec<String> = std::env::args().collect();
        timing.log("start", &format!("{} | {}", timestamp, cmd_line.join(" ")));
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
    // - Input is .1aln or FASTA (not PAF)
    // - User hasn't explicitly requested PAF output with --paf flag
    let input_is_paf =
        !input_file_types.is_empty() && input_file_types.iter().all(|ft| *ft == FileType::Paf);
    let want_paf_output = args.output_paf
        || args
            .output_file
            .as_ref()
            .map_or(false, |f| f.ends_with(".paf"));
    let use_1aln_workflow = !input_is_paf && !want_paf_output;

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
                    eprintln!(
                        "[sweepga] ERROR: Multiple FASTA files not supported in .1aln workflow yet"
                    );
                    eprintln!("[sweepga] Use --paf flag to use PAF workflow for multiple files");
                    std::process::exit(1);
                }

                let alignment_start = Instant::now();
                let path = Path::new(&args.files[0]);
                let fastga = create_fastga_integration(args.frequency, args.threads)?;

                if !args.quiet {
                    timing.log("align", &format!("Running FastGA on {}", args.files[0]));
                }

                let temp_1aln = fastga.align_to_temp_1aln(path, path)?;
                alignment_time = Some(alignment_start.elapsed().as_secs_f64());

                if !args.quiet {
                    timing.log(
                        "align",
                        &format!("FastGA complete ({:.1}s)", alignment_time.unwrap()),
                    );
                }

                let aln_path = temp_1aln.path().to_string_lossy().into_owned();
                (Some(temp_1aln), aln_path)
            } else {
                eprintln!("[sweepga] ERROR: No valid input provided");
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
            sparsity: args.sparsify,
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

        // Step 3: Filter .1aln directly using unified_filter (format-preserving)
        use crate::unified_filter::filter_file;
        if let Some(ref output_file) = args.output_file {
            filter_file(&aln_input_path, output_file, &filter_config, false)?;
        } else {
            // Write to temp file then copy to stdout
            let temp_output = tempfile::NamedTempFile::with_suffix(".1aln")?;
            filter_file(&aln_input_path, temp_output.path(), &filter_config, false)?;

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
                let path = Path::new(&args.files[0]);
                let groups = detect_genome_groups(path)?;

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
                        args.threads,
                        args.keep_self,
                        args.tempdir.as_deref(),
                        &timing,
                        args.quiet,
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

                    let fastga = create_fastga_integration(args.frequency, args.threads)?;
                    let temp_paf = fastga.align_to_temp_paf(path, path)?;

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
                        args.threads,
                        args.keep_self,
                        args.tempdir.as_deref(),
                        &timing,
                        args.quiet,
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

                    let target = Path::new(&args.files[0]);
                    let query = Path::new(&args.files[1]);
                    let fastga = create_fastga_integration(args.frequency, args.threads)?;
                    let temp_paf = fastga.align_to_temp_paf(target, query)?;

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
                    args.threads,
                    args.keep_self,
                    args.tempdir.as_deref(),
                    &timing,
                    args.quiet,
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
            _ => {
                anyhow::bail!("Invalid file combination: expected FASTA file(s) for alignment, 1 PAF for filtering, or 1 .1aln for filtering");
            }
        }
    } else {
        // Check if stdin is available
        use std::io::IsTerminal;
        let stdin_available = !std::io::stdin().is_terminal();

        if !stdin_available && !args.no_filter {
            use clap::CommandFactory;
            Args::command().print_help()?;
            std::process::exit(0);
        }

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

                let path = Path::new(&temp_path);
                let fastga = create_fastga_integration(args.frequency, args.threads)?;
                let temp_paf = fastga.align_to_temp_paf(path, path)?;

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
        sparsity: args.sparsify,
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
        eprintln!(
            "[sweepga] WARNING: Unknown ANI method '{}', using 'n50-identity'",
            args.ani_method
        );
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
    // Default: .1aln for FASTA/Aln input, PAF for PAF input
    let output_1aln = if args.output_paf {
        // User explicitly requested PAF output
        false
    } else if let Some(ref outfile) = args.output_file {
        // Auto-detect from output file extension
        !outfile.ends_with(".paf")
    } else if !input_file_types.is_empty() && input_file_types[0] == FileType::Paf {
        // PAF input → PAF output only
        false
    } else if !input_file_types.is_empty()
        && (input_file_types[0] == FileType::Fasta || input_file_types[0] == FileType::Aln)
    {
        // FASTA or .1aln input → default to .1aln output
        true
    } else {
        // Fallback: PAF output
        false
    };

    // Apply filtering - always to temp file, then copy to stdout
    if !args.quiet {
        timing.log("parse", &format!("Parsing input PAF: {input_path}"));
    }

    let output_temp = tempfile::NamedTempFile::with_suffix(".paf")?;
    let (output_file, output_path_buf) = output_temp.keep()?;
    let output_path = output_path_buf.to_str().unwrap().to_string();
    drop(output_file); // Close the file handle so filter_paf can open it

    // Note: -f (no_filter) implies --self (keep self-mappings)
    let filter = PafFilter::new(config)
        .with_keep_self(args.keep_self || args.no_filter)
        .with_scaffolds_only(args.scaffolds_only);
    filter.filter_paf(&input_path, &output_path)?;

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
                eprintln!(
                    "[sweepga] WARNING: .1aln output requires FASTA or .1aln input (not PAF), falling back to PAF output"
                );
            }
            (output_path.clone(), None::<tempfile::NamedTempFile>)
        } else {
            let first_type = input_file_types[0];
            if first_type != FileType::Fasta && first_type != FileType::Aln {
                if !args.quiet {
                    eprintln!(
                        "[sweepga] WARNING: .1aln output requires FASTA or .1aln input (not PAF), falling back to PAF output"
                    );
                }
                (output_path.clone(), None::<tempfile::NamedTempFile>)
            } else {
                // Convert filtered PAF to .1aln using PAFtoALN
                // Find PAFtoALN binary (bundled with fastga-rs)
                let paftoaln_bin = fastga_rs::embedded::get_binary_path("PAFtoALN")
                    .ok()
                    .and_then(|p| if p.exists() { Some(p) } else { None })
                    .ok_or_else(|| anyhow::anyhow!("PAFtoALN tool not found. Cannot convert PAF to .1aln. Use --paf flag to output PAF format instead."))?;

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

    Ok(())
}
