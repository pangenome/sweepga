mod fastga_integration;
mod mapping;
mod paf;
mod paf_filter;
mod plane_sweep;
mod plane_sweep_exact;
mod union_find;

use anyhow::Result;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};

use crate::paf_filter::{FilterConfig, FilterMode, PafFilter};

/// Parse a number that may have metric suffix (k/K=1000, m/M=1e6, g/G=1e9)
fn parse_metric_number(s: &str) -> Result<u32, String> {
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

    if result > u32::MAX as f64 {
        return Err(format!("Value {result} too large for u32"));
    }

    Ok(result as u32)
}

/// SweepGA - Fast genome alignment with sophisticated filtering
///
/// This tool applies wfmash's filtering algorithms to PAF input from FASTGA or other aligners
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

    /// Minimum block length
    #[clap(short = 'b', long = "block-length", default_value = "0", value_parser = parse_metric_number)]
    block_length: u32,

    /// Maximum overlap ratio for plane sweep filtering
    #[clap(short = 'p', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Keep this fraction of mappings
    #[clap(short = 'x', long = "sparsify", default_value = "1.0")]
    sparsify: f64,

    /// Mapping filter: "1:1" (best), "M:N" (top M per query, N per target; ∞/many for unbounded), "N" (no filter)
    #[clap(short = 'n', long = "num-mappings", default_value = "1:1")]
    num_mappings: String,

    /// Scaffold filter: "1:1" (best), "M:N" (top M per query, N per target; ∞/many for unbounded)
    #[clap(long = "scaffold-filter", default_value = "1:1")]
    scaffold_filter: String,

    /// Scaffold jump (gap) distance. 0 = disable scaffolding (plane sweep only), >0 = enable scaffolding
    #[clap(short = 'j', long = "scaffold-jump", default_value = "10k", value_parser = parse_metric_number)]
    scaffold_jump: u32,

    /// Minimum scaffold length when scaffolding is enabled
    #[clap(short = 's', long = "scaffold-mass", default_value = "10k", value_parser = parse_metric_number)]
    scaffold_mass: u32,

    /// Scaffold chain overlap threshold
    #[clap(short = 'O', long = "scaffold-overlap", default_value = "0.5")]
    scaffold_overlap: f64,

    /// Maximum distance from scaffold anchor
    #[clap(short = 'd', long = "scaffold-dist", default_value = "100k", value_parser = parse_metric_number)]
    scaffold_dist: u32,

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

fn main() -> Result<()> {
    let mut args = Args::parse();

    // Handle FASTA input mode (FastGA alignment)
    let _temp_paf = if !args.fasta_files.is_empty() {
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
                eprintln!("Running FastGA self-alignment on {}...", args.fasta_files[0]);
            } else {
                eprintln!("Running FastGA alignment: {} -> {}...",
                         args.fasta_files[0], args.fasta_files[1]);
            }
        }

        // Configure FastGA with minimal parameters (just threads)
        // FastGA will use its own defaults for identity and alignment length
        let fastga = FastGAIntegration::new(args.threads);

        // Run alignment and get temp PAF file
        let temp_paf = fastga.align_to_temp_paf(targets, queries)?;

        if !args.quiet {
            eprintln!("FastGA complete. Processing alignments from temp file: {:?}", temp_paf.path());
        }

        // Update args to use the generated PAF file
        args.input = Some(temp_paf.path().to_string_lossy().into_owned());

        Some(temp_paf)
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

    // Set up filter configuration
    let config = FilterConfig {
        chain_gap: args.scaffold_jump, // Use scaffold_jump for merging into chains
        min_block_length: args.block_length,
        mapping_filter_mode: plane_sweep_mode,
        mapping_max_per_query: plane_sweep_query_limit,
        mapping_max_per_target: plane_sweep_target_limit,
        plane_sweep_secondaries: 0, // Not used with new format
        scaffold_filter_mode,
        scaffold_max_per_query,
        scaffold_max_per_target,
        overlap_threshold: args.overlap,
        sparsity: args.sparsify,
        no_merge: true, // Never merge primary mappings (CIGAR strings become invalid)
        scaffold_gap: args.scaffold_jump,
        min_scaffold_length: args.scaffold_mass,
        scaffold_overlap_threshold: args.scaffold_overlap,
        scaffold_max_deviation: args.scaffold_dist,
        prefix_delimiter: '#', // Default PanSN delimiter
        skip_prefix: true,     // true = don't group by prefix (filter per sequence, not per genome)
    };

    // Progress indicator
    let progress = if !args.quiet {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        pb.set_message("Filtering PAF records...");
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
        pb.finish_with_message("Filtering complete");
    }

    Ok(())
}
