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
mod union_find;

use anyhow::Result;
use clap::Parser;

use crate::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

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

/// SweepGA - Fast genome alignment with sophisticated filtering
///
/// This tool wraps genome aligners (FastGA by default) and applies wfmash's filtering algorithms.
/// Can also process existing PAF files from any aligner.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// FASTA files: either one (self-alignment) or two (target then query)
    /// If no FASTA provided, reads PAF from stdin or -i flag
    #[clap(value_name = "FASTA", num_args = 0..=2)]
    fasta_files: Vec<String>,

    /// Input PAF file (alternative to FASTA files)
    #[clap(short = 'i', long = "input", conflicts_with = "fasta_files")]
    input: Option<String>,

    /// Output PAF file (stdout if not specified)
    #[clap(short = 'o', long = "output")]
    output: Option<String>,

    /// Aligner to use for FASTA input
    #[clap(long = "aligner", default_value = "fastga", value_parser = ["fastga"])]
    aligner: String,

    /// Minimum block length
    #[clap(short = 'b', long = "block-length", default_value = "0", value_parser = parse_metric_number)]
    block_length: u64,

    /// Maximum overlap ratio for plane sweep filtering
    #[clap(short = 'p', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Keep this fraction of mappings
    #[clap(short = 'x', long = "sparsify", default_value = "1.0")]
    sparsify: f64,

    /// Mapping filter: "1:1" (best), "M:N" (top M per query, N per target; ∞/many for unbounded)
    #[clap(short = 'n', long = "num-mappings", default_value = "many:many")]
    num_mappings: String,

    /// Scaffold filter: "1:1" (best), "M:N" (top M per query, N per target; ∞/many for unbounded)
    #[clap(short = 'm', long = "scaffold-filter", default_value = "1:1")]
    scaffold_filter: String,

    /// Scaffold jump (gap) distance. 0 = disable scaffolding (plane sweep only), >0 = enable scaffolding
    #[clap(short = 'j', long = "scaffold-jump", default_value = "10k", value_parser = parse_metric_number)]
    scaffold_jump: u64,

    /// Minimum scaffold length when scaffolding is enabled
    #[clap(short = 's', long = "scaffold-mass", default_value = "10k", value_parser = parse_metric_number)]
    scaffold_mass: u64,

    /// Scaffold chain overlap threshold
    #[clap(short = 'O', long = "scaffold-overlap", default_value = "0.5")]
    scaffold_overlap: f64,

    /// Maximum distance from scaffold anchor (0 = no rescue, only keep scaffold members)
    #[clap(short = 'd', long = "scaffold-dist", default_value = "20k", value_parser = parse_metric_number)]
    scaffold_dist: u64,

    /// Scoring function for plane sweep
    #[clap(long = "scoring", default_value = "log-length-identity",
           value_parser = ["identity", "length", "length-identity", "log-length-identity", "matches"])]
    scoring: String,

    /// Minimum identity threshold (0-1 fraction, 1-100%, or "aniN" for Nth percentile)
    #[clap(short = 'y', long = "min-identity", default_value = "ani50")]
    min_identity: String,

    /// Method for calculating ANI: all, orthogonal, nX[-sort] (e.g. n50, n90-identity, n100-score)
    #[clap(long = "ani-method", default_value = "n50")]
    ani_method: String,

    /// Minimum scaffold identity threshold (0-1 fraction or 1-100%, defaults to -y)
    #[clap(short = 'Y', long = "min-scaffold-identity")]
    min_scaffold_identity: Option<f64>,

    /// Disable all filtering
    #[clap(short = 'f', long = "no-filter")]
    no_filter: bool,

    /// Keep self-mappings (excluded by default)
    #[clap(long = "self")]
    keep_self: bool,

    /// Output scaffold chains only (for debugging)
    #[clap(long = "scaffolds-only")]
    scaffolds_only: bool,

    /// Quiet mode (no progress output)
    #[clap(long = "quiet")]
    quiet: bool,

    /// Number of threads for parallel processing
    #[clap(short = 't', long = "threads", default_value = "8")]
    threads: usize,
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
            eprintln!(
                "Warning: Invalid {filter_type} filter '{s}', using default 1:1"
            );
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
            eprintln!(
                "Warning: Invalid {filter_type} filter '{mode}', using default 1:1"
            );
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

    if lower.starts_with("ani") {
        // Parse aniN, aniN+X, or aniN-X
        if let Some(ani_value) = ani_percentile {
            let remainder = &lower[3..];
            if remainder.is_empty() {
                // Default "ani" means ani50
                return Ok(ani_value);
            }

            // Parse percentile number and optional offset
            // Find offset position if exists
            let (percentile_str, offset_part) = if let Some(plus_pos) = remainder.find('+') {
                (&remainder[..plus_pos], Some(('+', &remainder[plus_pos+1..])))
            } else if let Some(minus_pos) = remainder.find('-') {
                (&remainder[..minus_pos], Some(('-', &remainder[minus_pos+1..])))
            } else {
                (remainder, None)
            };

            // For now we only support ani50 (median), ignore the percentile number
            // Could extend to support ani25, ani75, etc. in future
            if !percentile_str.is_empty() && percentile_str != "50" {
                eprintln!("[sweepga] WARNING: Only ani50 (median) currently supported, using median");
            }

            // Apply offset if provided
            if let Some((sign, offset_str)) = offset_part {
                let offset: f64 = offset_str.parse()
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
            Ok(val / 100.0)  // Percentage to fraction
        } else {
            Ok(val)  // Already fraction
        }
    } else {
        anyhow::bail!("Invalid identity value: {}", value);
    }
}

/// Calculate ANI statistics between genome pairs using specified method
fn calculate_ani_stats(input_path: &str, method: AniMethod, quiet: bool) -> Result<f64> {
    use crate::paf_filter::{PafFilter, FilterConfig, FilterMode, ScoringFunction};
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
            let mut temp_filtered = NamedTempFile::new()?;
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
        return Ok(0.0);  // Default to no filtering
    }

    // Calculate weighted average ANI for each genome pair
    // ANI = total_matches / total_block_length
    let mut ani_values: Vec<f64> = genome_pairs.iter().map(|((_q, _t), (total_matches, total_length))| {
        if *total_length > 0.0 {
            total_matches / total_length
        } else {
            0.0
        }
    }).collect();

    ani_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Get 50th percentile (median)
    let median_idx = ani_values.len() / 2;
    let ani50 = if ani_values.len() % 2 == 0 && ani_values.len() > 1 {
        (ani_values[median_idx - 1] + ani_values[median_idx]) / 2.0
    } else {
        ani_values[median_idx]
    };

    eprintln!("[sweepga] ANI statistics from {} genome pairs:", genome_pairs.len());
    eprintln!("[sweepga]   Min: {:.1}%, Median: {:.1}%, Max: {:.1}%",
             ani_values.first().unwrap_or(&0.0) * 100.0,
             ani50 * 100.0,
             ani_values.last().unwrap_or(&0.0) * 100.0);

    Ok(ani50)
}

/// Calculate ANI using N-percentile method - use best alignments covering N% of genome pairs
fn calculate_ani_n_percentile(input_path: &str, percentile: f64, sort_method: NSort, quiet: bool) -> Result<f64> {
    if !quiet {
        let method_name = match sort_method {
            NSort::Length => format!("longest alignments (N{} by length)", percentile as i32),
            NSort::Identity => format!("highest identity alignments (N{} by identity)", percentile as i32),
            NSort::Score => format!("best scoring alignments (N{} by identity × log(length))", percentile as i32),
        };
        eprintln!("[sweepga] Calculating ANI from {}", method_name);
    }

    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    // Collect all alignments with their lengths and identity
    struct Alignment {
        query_genome: String,
        target_genome: String,
        matches: f64,
        block_length: f64,
        identity: f64,
    }

    let mut alignments = Vec::new();

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

    // Calculate total length and threshold for N-percentile
    let total_length: f64 = alignments.iter().map(|a| a.block_length).sum();
    let n_threshold = total_length * (percentile / 100.0);

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
    let mut ani_values: Vec<f64> = genome_pairs.iter().map(|((_q, _t), (total_matches, total_length))| {
        if *total_length > 0.0 {
            total_matches / total_length
        } else {
            0.0
        }
    }).collect();

    ani_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Get median
    let median_idx = ani_values.len() / 2;
    let ani50 = if ani_values.len() % 2 == 0 && ani_values.len() > 1 {
        (ani_values[median_idx - 1] + ani_values[median_idx]) / 2.0
    } else {
        ani_values[median_idx]
    };

    let method_str = match sort_method {
        NSort::Length => format!("N{} by length", percentile as i32),
        NSort::Identity => format!("N{} by identity", percentile as i32),
        NSort::Score => format!("N{} by score", percentile as i32),
    };
    eprintln!("[sweepga] ANI statistics from {} genome pairs ({}):", genome_pairs.len(), method_str);
    eprintln!("[sweepga]   Min: {:.1}%, Median: {:.1}%, Max: {:.1}%",
             ani_values.first().unwrap_or(&0.0) * 100.0,
             ani50 * 100.0,
             ani_values.last().unwrap_or(&0.0) * 100.0);

    Ok(ani50)
}

fn main() -> Result<()> {
    let start_time = Instant::now();
    let mut args = Args::parse();

    // Track alignment time separately
    let mut alignment_time: Option<f64> = None;

    // Handle FASTA input mode (alignment generation)
    let _temp_paf = if !args.fasta_files.is_empty() {
        let alignment_start = Instant::now();
        match args.aligner.as_str() {
            "fastga" => {
                use crate::fastga_integration::FastGAIntegration;
                use std::path::Path;

                let (targets, queries) = match args.fasta_files.len() {
                    1 => {
                        // Self-alignment
                        let path = Path::new(&args.fasta_files[0]);
                        (path, path)
                    }
                    2 => {
                        // Two files: first is target, second is query
                        (Path::new(&args.fasta_files[0]), Path::new(&args.fasta_files[1]))
                    }
                    _ => {
                        anyhow::bail!("Expected 1 or 2 FASTA files, got {}", args.fasta_files.len());
                    }
                };

                if !args.quiet {
                    if args.fasta_files.len() == 1 {
                        eprintln!("[sweepga] Running {} self-alignment on {}...", args.aligner, args.fasta_files[0]);
                    } else {
                        eprintln!("[sweepga] Running {} alignment: {} -> {}...",
                                 args.aligner, args.fasta_files[0], args.fasta_files[1]);
                    }
                }

                // Configure FastGA with minimal parameters (just threads)
                // FastGA will use its own defaults for identity and alignment length
                let fastga = FastGAIntegration::new(args.threads);

                // Run alignment and get temp PAF file
                let temp_paf = fastga.align_to_temp_paf(targets, queries)?;

                alignment_time = Some(alignment_start.elapsed().as_secs_f64());
                if !args.quiet {
                    eprintln!("[sweepga] {} alignment complete ({:.1}s). Applying filtering...",
                             args.aligner, alignment_time.unwrap());
                }

                // Update args to use the generated PAF file
                args.input = Some(temp_paf.path().to_string_lossy().into_owned());

                Some(temp_paf)
            }
            _ => {
                anyhow::bail!("Unknown aligner: {}", args.aligner);
            }
        }
    } else {
        None
    };

    // Check if stdin is available when no input specified
    let stdin_available = if args.input.is_none() {
        use std::io::IsTerminal;
        !std::io::stdin().is_terminal()
    } else {
        false
    };

    // If no input specified and no stdin, print help
    if args.input.is_none() && !stdin_available && !args.no_filter && args.fasta_files.is_empty() {
        use clap::CommandFactory;
        Args::command().print_help()?;
        std::process::exit(0);
    }

    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // Handle no-filter mode - just copy input to output
    if args.no_filter {
        use std::io::{self, BufRead, BufReader, Write};

        let input: Box<dyn BufRead> = if let Some(ref path) = args.input {
            Box::new(BufReader::new(std::fs::File::open(path)?))
        } else {
            Box::new(BufReader::new(io::stdin()))
        };

        let mut output: Box<dyn Write> = if let Some(ref path) = args.output {
            Box::new(std::fs::File::create(path)?)
        } else {
            Box::new(io::stdout())
        };

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
        "identity" => ScoringFunction::Identity,
        "length" => ScoringFunction::Length,
        "length-identity" => ScoringFunction::LengthIdentity,
        "matches" => ScoringFunction::Matches,
        "log-length-identity" | _ => ScoringFunction::LogLengthIdentity,
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
        min_identity: 0.0,  // Will be set later
        min_scaffold_identity: 0.0,  // Will be set later
    };

    // Handle input: if stdin, save to temp file for two-pass processing
    let (_input_temp, input_path) = if let Some(ref path) = args.input {
        (None, path.clone())
    } else {
        // Read from stdin into a temp file for two-pass processing
        let temp = tempfile::NamedTempFile::new()?;
        let temp_path = temp.path().to_str().unwrap().to_string();

        // Copy stdin to temp file
        {
            use std::io::{self, BufRead, Write};
            let stdin = io::stdin();
            let mut temp_file = std::fs::File::create(&temp_path)?;

            for line in stdin.lock().lines() {
                writeln!(temp_file, "{}", line?)?;
            }
        }

        (Some(temp), temp_path)
    };

    // Parse ANI calculation method
    let ani_method = parse_ani_method(&args.ani_method);
    if ani_method.is_none() {
        eprintln!("[sweepga] WARNING: Unknown ANI method '{}', using 'n50-identity'", args.ani_method);
    }
    let ani_method = ani_method.unwrap_or(AniMethod::NPercentile(50.0, NSort::Identity));

    // Now calculate ANI if needed for identity thresholds
    let ani_percentile = if args.min_identity.to_lowercase().contains("ani") {
        Some(calculate_ani_stats(&input_path, ani_method, args.quiet)?)
    } else {
        None
    };

    // Parse identity thresholds
    let min_identity = parse_identity_value(&args.min_identity, ani_percentile)?;
    if !args.quiet {
        eprintln!("[sweepga] Using minimum identity threshold: {:.1}%", min_identity * 100.0);
    }

    let min_scaffold_identity = if let Some(scaffold_val) = args.min_scaffold_identity {
        if scaffold_val > 1.0 {
            scaffold_val / 100.0
        } else {
            scaffold_val
        }
    } else {
        min_identity
    };

    // Create final config with calculated identity values
    let mut config = temp_config;
    config.min_identity = min_identity;
    config.min_scaffold_identity = min_scaffold_identity;

    let (_output_temp, output_path) = if let Some(ref path) = args.output {
        (None, path.clone())
    } else {
        // Use temp file then copy to stdout
        let temp = tempfile::NamedTempFile::new()?;
        let path = temp.path().to_str().unwrap().to_string();
        (Some(temp), path)
    };

    // Apply filtering
    // Note: -f (no_filter) implies --self (keep self-mappings)
    let filter = PafFilter::new(config)
        .with_keep_self(args.keep_self || args.no_filter)
        .with_scaffolds_only(args.scaffolds_only);
    filter.filter_paf(&input_path, &output_path)?;

    // If output was to stdout, copy temp file to stdout
    if args.output.is_none() {
        use std::io::{self, BufRead, BufReader, Write};
        let file = std::fs::File::open(&output_path)?;
        let reader = BufReader::new(file);
        let stdout = io::stdout();
        let mut handle = stdout.lock();

        for line in reader.lines() {
            writeln!(handle, "{}", line?)?;
        }
    }

    if !args.quiet {
        let total_elapsed = start_time.elapsed().as_secs_f64();
        let filtering_time = if let Some(align_time) = alignment_time {
            total_elapsed - align_time
        } else {
            total_elapsed
        };

        // Format the complete message
        if let Some(align_time) = alignment_time {
            eprintln!("[sweepga] Complete. Alignment: {:.1}s, Filtering: {:.1}s, Total: {:.1}s",
                       align_time, filtering_time, total_elapsed);
        } else {
            eprintln!("[sweepga] Filtering complete. Total: {:.1}s", total_elapsed);
        }
    }

    Ok(())
}
