mod mapping;
mod paf;
mod paf_filter;
mod union_find;
mod plane_sweep;
mod plane_sweep_exact;

use anyhow::Result;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};

use crate::paf_filter::{PafFilter, FilterConfig, FilterMode};

/// SweepGA - Fast genome alignment with sophisticated filtering
///
/// This tool applies wfmash's filtering algorithms to PAF input from FASTGA or other aligners
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Input PAF file
    #[clap(short = 'i', long = "input")]
    input: Option<String>,

    /// Output PAF file (stdout if not specified)
    #[clap(short = 'o', long = "output")]
    output: Option<String>,

    /// Minimum block length [0]
    #[clap(short = 'b', long = "block-length", default_value = "0")]
    block_length: u32,

    /// Maximum overlap ratio for plane sweep filtering [0.95]
    #[clap(short = 'p', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Keep this fraction of mappings [1.0]
    #[clap(short = 'x', long = "sparsify", default_value = "1.0")]
    sparsify: f64,

    /// Primary plane sweep filter (default: N = no filtering)
    /// Format: "1:1", "1" (=1:∞), "N" (=N:N, no filtering), or "M:N"
    #[clap(short = 'm', long = "mapping-filter", default_value = "N")]
    mapping_filter: String,

    /// Number of mappings to keep per position in plane sweep (default: 0 = best only)
    /// Use -n 1 to keep best + 1 secondary, -n -1 to keep all non-overlapping
    #[clap(short = 'n', long = "num-mappings", default_value = "0")]
    num_mappings: i32,

    /// Scaffold filter when -s > 0 (default: 1:1)
    /// Format: "1:1", "1" (=1:∞), "N" (=N:N, no filtering), or "M:N"
    #[clap(long = "scaffold-filter", default_value = "1:1")]
    scaffold_filter: String,

    /// Scaffold jump (gap) distance [100000]
    #[clap(short = 'j', long = "scaffold-jump", default_value = "100000")]
    scaffold_jump: u32,

    /// Minimum scaffold length [10000]
    #[clap(short = 's', long = "scaffold-mass", default_value = "10000")]
    scaffold_mass: u32,

    /// Scaffold chain overlap threshold [0.5]
    #[clap(short = 'O', long = "scaffold-overlap", default_value = "0.5")]
    scaffold_overlap: f64,

    /// Maximum distance from scaffold anchor [100000]
    #[clap(short = 'd', long = "scaffold-dist", default_value = "100000")]
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
    match mode.to_lowercase().as_str() {
        "1:1" => (FilterMode::OneToOne, Some(1), Some(1)),
        "1" | "1:n" | "1:∞" => (FilterMode::OneToMany, Some(1), None),
        "n:1" | "∞:1" => (FilterMode::ManyToMany, None, Some(1)),  // N:1 is ManyToMany with target limit
        "n" | "n:n" | "∞" | "∞:∞" => (FilterMode::ManyToMany, None, None),
        s if s.contains(':') => {
            // Parse custom like "10:5"
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() == 2 {
                let per_query = parts[0].parse::<usize>().ok();
                let per_target = parts[1].parse::<usize>().ok();
                if per_query.is_some() || per_target.is_some() {
                    let mode = match (per_query, per_target) {
                        (Some(1), Some(1)) => FilterMode::OneToOne,
                        (Some(1), _) => FilterMode::OneToMany,
                        _ => FilterMode::ManyToMany,  // Any other combination uses ManyToMany
                    };
                    return (mode, per_query, per_target);
                }
            }
            eprintln!("Warning: Invalid {} filter '{}', using default N:N", filter_type, s);
            (FilterMode::ManyToMany, None, None)
        }
        _ => {
            eprintln!("Warning: Invalid {} filter '{}', using default N:N", filter_type, mode);
            (FilterMode::ManyToMany, None, None)
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    // If no input specified, print help and exit
    if args.input.is_none() && !args.no_filter {
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
    let plane_sweep_secondaries = if args.num_mappings < 0 {
        usize::MAX  // Keep all non-overlapping
    } else {
        args.num_mappings as usize  // Keep best + this many secondaries
    };

    // Parse mapping filter mode (for backward compatibility)
    let (mapping_filter_mode, mapping_max_per_query, mapping_max_per_target) =
        if args.mapping_filter == "N" {
            // Use plane sweep as default with -n parameter
            (FilterMode::OneToMany, Some(1 + plane_sweep_secondaries), None)
        } else {
            parse_filter_mode(&args.mapping_filter, "mapping")
        };

    // Parse scaffold filter mode
    let (scaffold_filter_mode, scaffold_max_per_query, scaffold_max_per_target) = parse_filter_mode(&args.scaffold_filter, "scaffold");

    // Set up filter configuration
    let config = FilterConfig {
        chain_gap: args.scaffold_jump, // Use scaffold_jump for merging into chains
        min_block_length: args.block_length,
        mapping_filter_mode,
        mapping_max_per_query,
        mapping_max_per_target,
        plane_sweep_secondaries,  // Add plane sweep parameter
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
        prefix_delimiter: '#',  // Default PanSN delimiter
        skip_prefix: true,      // true = group by prefix (default behavior for PanSN)
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

    // Input path is required (unless --no-filter)
    let input_path = if let Some(ref path) = args.input {
        path.clone()
    } else {
        // This shouldn't happen due to earlier check, but handle it gracefully
        eprintln!("Error: Input file is required");
        std::process::exit(1);
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