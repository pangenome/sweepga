mod compact_filter;

use anyhow::Result;
use clap::Parser;
use compact_filter::CompactPafFilter;

/// SweepGA Compact - Memory-efficient PAF filtering
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Input PAF file
    #[clap(short = 'i', long = "input")]
    input: String,

    /// Output PAF file
    #[clap(short = 'o', long = "output")]
    output: String,

    /// Minimum block length [0]
    #[clap(short = 'b', long = "block-length", default_value = "0")]
    block_length: u32,

    /// Filter specification: "N" (no filter), "1:1", "1", etc.
    #[clap(short = 'm', long = "filter", default_value = "N")]
    filter: String,

    /// Overlap threshold [0.95]
    #[clap(short = 'p', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Prefix delimiter for grouping [#]
    #[clap(short = 'Y', long = "delimiter", default_value = "#")]
    delimiter: char,

    /// Keep self-mappings
    #[clap(long = "self")]
    keep_self: bool,

    /// Disable prefix grouping
    #[clap(long = "no-prefix")]
    no_prefix: bool,

    /// Debug mode - just read and write without filtering
    #[clap(long = "debug")]
    debug: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Create compact filter
    let mut filter = CompactPafFilter::new(
        args.block_length,
        args.keep_self,
        args.delimiter,
        !args.no_prefix,  // group_by_prefix
    );

    // Read PAF and build indices
    filter.read_paf(&args.input)?;

    // Get offsets to keep
    let kept_offsets = if args.debug {
        // Debug mode: keep all mappings (no filtering)
        eprintln!("DEBUG MODE: Writing all records without filtering");
        filter.get_all_offsets()
    } else {
        // Apply plane sweep filtering
        filter.apply_plane_sweep(&args.filter, args.overlap)?
    };

    // Write filtered output
    filter.write_filtered_output(&args.input, &args.output, kept_offsets)?;

    eprintln!("Done!");
    Ok(())
}