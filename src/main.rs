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
use indicatif::{ProgressBar, ProgressStyle};

use crate::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

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

/// Calculate ANI statistics between genome pairs
fn calculate_ani_stats(input_path: &str) -> Result<f64> {
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    // Group alignments by genome pair (using PanSN prefixes)
    let mut genome_pairs: HashMap<(String, String), Vec<f64>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Extract genome prefixes (before first #)
        let query_genome = fields[0].split('#').next().unwrap_or(fields[0]).to_string();
        let target_genome = fields[5].split('#').next().unwrap_or(fields[5]).to_string();

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        // Parse matches and block length
        let matches = fields[9].parse::<f64>().unwrap_or(0.0);
        let block_len = fields[10].parse::<f64>().unwrap_or(1.0);
        let identity = matches / block_len.max(1.0);

        // Check for divergence tag
        let mut final_identity = identity;
        for field in &fields[11..] {
            if let Some(div_str) = field.strip_prefix("dv:f:") {
                if let Ok(div) = div_str.parse::<f64>() {
                    final_identity = 1.0 - div;
                    break;
                }
            }
        }

        let key = if query_genome < target_genome {
            (query_genome, target_genome)
        } else {
            (target_genome, query_genome)
        };

        genome_pairs.entry(key).or_default().push(final_identity);
    }

    if genome_pairs.is_empty() {
        eprintln!("[sweepga] WARNING: No inter-genome alignments found for ANI calculation");
        return Ok(0.0);  // Default to no filtering
    }

    // Calculate weighted average ANI for each genome pair
    let mut ani_values: Vec<f64> = genome_pairs.iter().map(|((_q, _t), identities)| {
        if identities.is_empty() {
            0.0
        } else {
            identities.iter().sum::<f64>() / identities.len() as f64
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

fn main() -> Result<()> {
    let mut args = Args::parse();

    // Handle FASTA input mode (alignment generation)
    let _temp_paf = if !args.fasta_files.is_empty() {
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

                if !args.quiet {
                    eprintln!("[sweepga] {} alignment complete. Applying filtering...", args.aligner);
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

    // Progress indicator
    let progress = if !args.quiet {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        pb.set_message("[sweepga] Filtering PAF records...");
        Some(pb)
    } else {
        None
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

    // Now calculate ANI if needed for identity thresholds
    let ani_percentile = if args.min_identity.to_lowercase().contains("ani") {
        Some(calculate_ani_stats(&input_path)?)
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
    let filter = PafFilter::new(config)
        .with_keep_self(args.keep_self)
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

    if let Some(pb) = progress {
        pb.finish_and_clear();
        if !args.quiet {
            eprintln!("[sweepga] Filtering complete");
        }
    }

    Ok(())
}
