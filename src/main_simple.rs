mod mapping;
mod paf;
mod plane_sweep;

use anyhow::Result;
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::paf::PafReader;
use crate::plane_sweep::apply_plane_sweep;

/// SweepGA - Fast genome alignment filtering with plane sweep
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Input PAF file (stdin if not specified)
    #[clap(short = 'i', long = "input")]
    input: Option<String>,

    /// Output PAF file (stdout if not specified)
    #[clap(short = 'o', long = "output")]
    output: Option<String>,

    /// Minimum block length [0]
    #[clap(short = 'b', long = "block-length", default_value = "0")]
    block_length: u32,

    /// Plane sweep filter: "N" (no filter), "1:1", "1" (1:∞), etc.
    #[clap(short = 'm', long = "mapping-filter", default_value = "N")]
    filter: String,

    /// Overlap threshold for plane sweep [0.95]
    #[clap(short = 'p', long = "overlap", default_value = "0.95")]
    overlap: f64,

    /// Prefix delimiter for grouping [#]
    #[clap(short = 'Y', long = "group-prefix", default_value = "#")]
    delimiter: char,

    /// Keep self-mappings
    #[clap(long = "self")]
    keep_self: bool,

    /// Disable prefix grouping
    #[clap(long = "no-prefix")]
    no_prefix: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Read PAF input
    eprintln!("Reading PAF input...");
    let records = if let Some(ref path) = args.input {
        let file = File::open(path)?;
        let mut reader = PafReader::new(file);
        reader.read_all()?
    } else {
        let stdin = std::io::stdin();
        let mut reader = PafReader::new(stdin);
        reader.read_all()?
    };
    eprintln!("Read {} records", records.len());

    // Separate into vectors and filter
    let mut all_records = Vec::new();
    let mut mappings = Vec::new();
    let mut aux_data = Vec::new();
    let mut indices = Vec::new();

    for (idx, (paf, mapping, aux)) in records.into_iter().enumerate() {
        // Skip if below min block length
        if mapping.block_length < args.block_length {
            continue;
        }

        // Skip self-mappings unless requested
        if !args.keep_self && paf.query_name == paf.ref_name {
            continue;
        }

        all_records.push(paf);
        mappings.push(mapping);
        aux_data.push(aux);
        indices.push(idx);
    }

    eprintln!("After initial filters: {} mappings", mappings.len());

    // Check identity distribution
    let mut id_buckets = vec![0; 11];
    for mapping in &mappings {
        let identity = mapping.identity();
        let bucket = (identity / 10.0).min(10.0) as usize;
        id_buckets[bucket] += 1;
    }
    eprintln!("Identity distribution:");
    for (i, count) in id_buckets.iter().enumerate() {
        if *count > 0 {
            eprintln!("  {}0-{}0%: {}", i, i+1, count);
        }
    }

    // Apply plane sweep filtering
    if args.filter != "N" && args.filter != "N:N" {
        eprintln!("Applying plane sweep filter: {}", args.filter);
        apply_plane_sweep(
            &mut mappings,
            &mut aux_data,
            &args.filter,
            args.overlap,
            args.delimiter,
            !args.no_prefix,  // skip_prefix = true means group by prefix
        );
        eprintln!("After plane sweep: {} mappings", mappings.len());
    }

    // Track which records to keep after filtering
    let kept_indices: std::collections::HashSet<usize> = if mappings.len() < all_records.len() {
        // Plane sweep filtered some records, need to track which ones remain
        (0..mappings.len()).collect()
    } else {
        // No filtering applied or all records kept
        (0..all_records.len()).collect()
    };

    // Write output
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.output {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };

    // Write kept records
    for idx in kept_indices {
        if idx < all_records.len() {
            writeln!(writer, "{}", all_records[idx])?;
        }
    }

    Ok(())
}