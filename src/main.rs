mod mapping;
mod paf;
mod paf_filter;
mod union_find;

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
    /// Input PAF file (stdin if not specified)
    #[clap(short = 'i', long = "input")]
    input: Option<String>,

    /// Output PAF file (stdout if not specified)
    #[clap(short = 'o', long = "output")]
    output: Option<String>,

    /// Chain jump (gap) distance [2000]
    #[clap(short = 'c', long = "chain-jump", default_value = "2000")]
    chain_jump: u32,

    /// Minimum block length [0]
    #[clap(short = 'l', long = "block-length", default_value = "0")]
    block_length: u32,

    /// Maximum mappings per segment (controls "many" in filtering modes)
    #[clap(short = 'n', long = "mappings")]
    mappings: Option<usize>,

    /// Maximum overlap ratio [0.95]
    #[clap(short = 'O', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Keep this fraction of mappings [1.0]
    #[clap(short = 'x', long = "sparsify", default_value = "1.0")]
    sparsify: f64,

    /// Filtering mode: "1:1", "1:∞"/"1"/"map" (plane sweep), or "N:N"/"N" (no filtering)
    #[clap(short = 'm', long = "mode", default_value = "1")]
    mode: String,

    /// Keep all fragment mappings (no merge)
    #[clap(short = 'M', long = "no-merge")]
    no_merge: bool,

    /// Scaffold jump (gap) distance [100000]
    #[clap(short = 'j', long = "scaffold-jump", default_value = "100000")]
    scaffold_jump: u32,

    /// Minimum scaffold length [10000]
    #[clap(short = 'S', long = "scaffold-mass", default_value = "10000")]
    scaffold_mass: u32,

    /// Scaffold chain overlap threshold [0.5]
    #[clap(long = "scaffold-overlap", default_value = "0.5")]
    scaffold_overlap: f64,

    /// Maximum distance from scaffold anchor [100000]
    #[clap(short = 'D', long = "scaffold-dist", default_value = "100000")]
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

fn main() -> Result<()> {
    let args = Args::parse();

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

    // Parse filtering mode
    let (filter_mode, max_per_query, max_per_target) = match args.mode.as_str() {
        "1:1" => (FilterMode::OneToOne, Some(1), Some(1)),
        "1" | "1:∞" | "1:inf" | "1:N" | "map" => (FilterMode::OneToMany, Some(1), args.mappings),
        "N" | "N:N" => (FilterMode::ManyToMany, args.mappings, args.mappings),
        _ => {
            eprintln!("Invalid mode: {}. Use '1:1', '1' (or '1:∞', '1:inf', 'map'), or 'N' (or 'N:N')", args.mode);
            std::process::exit(1);
        }
    };

    // Set up filter configuration
    let config = FilterConfig {
        chain_gap: args.chain_jump,
        min_block_length: args.block_length,
        mapping_filter_mode: filter_mode,
        mapping_max_per_query: max_per_query,
        mapping_max_per_target: max_per_target,
        scaffold_filter_mode: FilterMode::OneToMany, // Default scaffold filtering (like wfmash)
        scaffold_max_per_query: None,  // No limit - keep all non-overlapping scaffolds
        scaffold_max_per_target: None,  // No limit on targets
        overlap_threshold: args.overlap,
        sparsity: args.sparsify,
        no_merge: args.no_merge,
        scaffold_gap: args.scaffold_jump,
        min_scaffold_length: args.scaffold_mass,
        scaffold_overlap_threshold: args.scaffold_overlap,
        scaffold_max_deviation: args.scaffold_dist,
        prefix_delimiter: '#',  // Default PanSN delimiter
        skip_prefix: true,      // Group by prefix by default
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

    // Handle stdin/stdout if paths not specified
    // Keep temp files in scope so they don't get deleted
    let (_input_temp, input_path) = if let Some(ref path) = args.input {
        (None, path.clone())
    } else {
        // Write stdin to temp file
        use std::io::{self, BufRead, Write};
        let temp = tempfile::NamedTempFile::new()?;
        let stdin = io::stdin();
        {
            let mut writer = std::io::BufWriter::new(&temp);
            for line in stdin.lock().lines() {
                writeln!(writer, "{}", line?)?;
            }
            writer.flush()?;
        }
        let path = temp.path().to_str().unwrap().to_string();
        (Some(temp), path)
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