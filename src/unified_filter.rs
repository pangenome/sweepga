//! Unified filtering module - works with both .1aln and PAF formats
//!
//! This module provides format-agnostic filtering by:
//! 1. Reading alignments into RecordMeta (same structure used by PAF filter)
//! 2. Using the SAME apply_filters() logic from PafFilter
//! 3. Writing passing records in the original format using ranks
//!
//! Key features:
//! - No PAF conversion for .1aln files (unless --paf flag is used)
//! - Preserves exact input format
//! - Uses SAME filtering logic as PAF filter (no duplication!)

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;

use crate::mapping::ChainStatus;
use crate::paf_filter::{FilterConfig, PafFilter, RecordMeta};

/// Extract RecordMeta from .1aln file (analogous to PAF extract_metadata)
pub fn extract_1aln_metadata<P: AsRef<Path>>(
    path: P,
) -> Result<(Vec<RecordMeta>, HashMap<String, i64>)> {
    let path_str = path.as_ref().to_str().context("Invalid path")?;

    // Open .1aln reader
    let mut reader = fastga_rs::AlnReader::open(path_str).context("Failed to open .1aln file")?;

    // Get all sequence names upfront (efficient bulk lookup)
    let id_to_name = reader.get_all_seq_names();

    // Build reverse mapping for writing
    let mut name_to_id = HashMap::new();
    for (id, name) in &id_to_name {
        name_to_id.insert(name.clone(), *id);
    }

    eprintln!(
        "[unified_filter] Loaded {} sequence names from .1aln",
        id_to_name.len()
    );

    // Read all alignments and create RecordMeta
    let mut metadata = Vec::new();
    let mut rank = 0;

    // Note: We need to read X records manually since AlnReader doesn't expose them
    // For now, use the mismatches field which comes from the 'D' record
    // TODO: Sum X values if they differ from D

    while let Some(aln) = reader.read_alignment()? {
        // Get actual sequence names (not numeric IDs)
        let query_id_num: i64 = aln.query_name.parse().unwrap_or(-1);
        let target_id_num: i64 = aln.target_name.parse().unwrap_or(-1);

        let query_name_full = id_to_name
            .get(&query_id_num)
            .cloned()
            .unwrap_or_else(|| aln.query_name.clone());
        let target_name_full = id_to_name
            .get(&target_id_num)
            .cloned()
            .unwrap_or_else(|| aln.target_name.clone());

        // Truncate to first word (before first space) to match PAF format
        // FASTA headers can contain descriptions after the ID, but PAF only uses the ID
        let query_name = query_name_full
            .split_whitespace()
            .next()
            .unwrap_or(&query_name_full)
            .to_string();
        let target_name = target_name_full
            .split_whitespace()
            .next()
            .unwrap_or(&target_name_full)
            .to_string();

        // Calculate identity and matches from .1aln data
        // The .1aln format (via fastga-rs AlnReader) provides:
        // - aln.matches: Number of matching bases (calculated from X records)
        //   matches = identity * query_span, where identity uses ALNtoPAF formula
        // - aln.mismatches: Total diffs from X records (or D as fallback)
        // - aln.block_len: Query span (query_end - query_start)
        //
        // Note: fastga-rs now correctly reads X records and calculates matches
        // using the same formula as ALNtoPAF: (sum(X) - del) / query_span / 2.0
        //
        // IMPORTANT: Do NOT use aln.identity() method! It recalculates identity
        // as matches/(matches+mismatches+gaps) which is NOT compatible with
        // the ALNtoPAF formula. Instead, derive identity from matches field.

        let query_span = (aln.query_end - aln.query_start) as u64;
        let target_span = (aln.target_end - aln.target_start) as u64;

        // PAF format column 10: alignment block length (query_span + target_span)
        // This matches ALNtoPAF output format
        let block_length = query_span + target_span;

        // Use matches from fastga-rs (calculated using X records)
        let matches = aln.matches as u64;

        // Derive identity from matches: identity = matches / query_span
        // This preserves the ALNtoPAF-compatible identity calculation
        let identity = if query_span > 0 {
            matches as f64 / query_span as f64
        } else {
            0.0
        };

        // Create RecordMeta (same structure as PAF filter!)
        metadata.push(RecordMeta {
            rank,
            query_name,
            target_name,
            query_start: aln.query_start as u64,
            query_end: aln.query_end as u64,
            target_start: aln.target_start as u64,
            target_end: aln.target_end as u64,
            block_length,
            identity,
            matches,
            alignment_length: block_length, // Total alignment length including gaps
            strand: aln.strand,
            chain_id: None,
            chain_status: ChainStatus::Unassigned,
            discard: false,
            overlapped: false,
        });

        rank += 1;
    }

    eprintln!(
        "[unified_filter] Read {} alignments from .1aln",
        metadata.len()
    );

    Ok((metadata, name_to_id))
}

/// Write filtered .1aln using passing ranks
/// IMPORTANT: Also copies the .1gdb file from input to output to preserve sequence names
pub fn write_1aln_filtered<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    passing_ranks: &HashMap<usize, RecordMeta>,
    _name_to_id: &HashMap<String, i64>,
) -> Result<()> {
    let path_str = input_path.as_ref().to_str().context("Invalid path")?;

    // Re-open input for reading
    let mut reader = fastga_rs::AlnReader::open(path_str)?;

    // Create output writer with GDB copied from input
    // This preserves sequence names when filtering
    let mut writer = fastga_rs::AlnWriter::create_with_gdb(
        output_path.as_ref(),
        input_path.as_ref(),
        true, // binary
    )?;

    let mut rank = 0;
    let mut written = 0;

    // Use unsafe to access the internal OneFile for low-level reading
    unsafe {
        let input_file = &mut reader.file;

        // Read through input file, looking for 'A' records
        loop {
            let line_type = input_file.read_line();

            if line_type == '\0' {
                break; // EOF
            }

            if line_type == 'A' {
                // Found an alignment record
                if passing_ranks.contains_key(&rank) {
                    // This record passed filtering - copy it directly with trace data
                    writer.copy_alignment_record_from_file(input_file)?;
                    written += 1;
                } else {
                    // Skip this alignment - read past all its associated records
                    loop {
                        let next_type = input_file.read_line();
                        if next_type == '\0' {
                            break; // EOF
                        }
                        if next_type == 'T' {
                            // Read past X record too
                            input_file.read_line();
                            break;
                        }
                        if next_type == 'A' {
                            // Hit next alignment without trace data
                            break;
                        }
                    }
                }
                rank += 1;
            }
        }
    }

    // Explicitly finalize the output file
    writer.finalize();

    eprintln!(
        "[unified_filter] Wrote {written} alignments to .1aln (GDB and trace data preserved)"
    );
    Ok(())
}

/// Main unified filtering function - works for both .1aln and PAF
pub fn filter_file<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    config: &FilterConfig,
    force_paf_output: bool,
    keep_self: bool,
) -> Result<()> {
    let input_str = input_path.as_ref().to_str().context("Invalid input path")?;

    // Determine input format by checking file content
    // .1aln files are binary OneCode format starting with specific magic bytes
    let is_1aln = {
        use std::io::Read;
        let mut file = std::fs::File::open(&input_path)?;
        let mut header = [0u8; 16];
        match file.read(&mut header) {
            Ok(n) if n >= 2 => {
                // OneCode format starts with '1' followed by version (ASCII '3')
                // or has specific binary signatures
                header[0] == b'1' && header[1] == b' '
            }
            _ => {
                // If we can't read header, fall back to extension check
                input_str.ends_with(".1aln")
            }
        }
    };

    if is_1aln {
        // .1aln input workflow
        eprintln!("[unified_filter] Reading .1aln metadata...");
        let (metadata, name_to_id) = extract_1aln_metadata(&input_path)?;

        // Use SAME filtering logic as PAF!
        eprintln!("[unified_filter] Applying filters...");
        let filter = PafFilter::new(config.clone()).with_keep_self(keep_self);
        let passing_ranks = filter.apply_filters(metadata)?;

        eprintln!(
            "[unified_filter] {} records passed filtering",
            passing_ranks.len()
        );

        // Determine output format
        let output_str = output_path
            .as_ref()
            .to_str()
            .context("Invalid output path")?;
        let output_1aln = !force_paf_output && !output_str.ends_with(".paf");

        if output_1aln {
            // Write .1aln output
            write_1aln_filtered(&input_path, &output_path, &passing_ranks, &name_to_id)?;
        } else {
            // Write PAF output - need to convert .1aln → PAF first
            // For now, use existing ALNtoPAF tool or implement direct conversion
            anyhow::bail!("1aln → PAF output not yet implemented in unified filter");
        }
    } else {
        // PAF input workflow - use existing PAF filter directly
        let filter = PafFilter::new(config.clone()).with_keep_self(keep_self);
        let input_str = input_path.as_ref();
        let output_str = output_path.as_ref();
        filter.filter_paf(input_str, output_str)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf_filter::{FilterMode, ScoringFunction};

    #[test]
    fn test_unified_1aln_filtering() {
        // Skip test if test_output.1aln doesn't exist
        if !std::path::Path::new("test_output.1aln").exists() {
            eprintln!("Skipping test - test_output.1aln not found");
            return;
        }

        let config = FilterConfig {
            chain_gap: 0,
            min_block_length: 0,
            mapping_filter_mode: FilterMode::OneToOne,
            mapping_max_per_query: Some(1),
            mapping_max_per_target: Some(1),
            plane_sweep_secondaries: 0,
            scaffold_filter_mode: FilterMode::ManyToMany,
            scaffold_max_per_query: None,
            scaffold_max_per_target: None,
            overlap_threshold: 0.95,
            sparsity: 1.0,
            no_merge: true,
            scaffold_gap: 0,
            min_scaffold_length: 0,
            scaffold_overlap_threshold: 0.95,
            scaffold_max_deviation: 0,
            prefix_delimiter: '#',
            skip_prefix: false,
            scoring_function: ScoringFunction::LogLengthIdentity,
            min_identity: 0.0,
            min_scaffold_identity: 0.0,
        };

        filter_file(
            "test_output.1aln",
            "test_filtered_result.1aln",
            &config,
            false,
            false, // keep_self
        )
        .unwrap();

        // Verify output exists
        assert!(std::path::Path::new("test_filtered_result.1aln").exists());

        eprintln!("✓ Test passed - filtered output created successfully");

        // Clean up
        let _ = std::fs::remove_file("test_filtered_result.1aln");
    }
}
