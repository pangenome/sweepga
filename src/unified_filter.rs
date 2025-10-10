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

use crate::paf_filter::{FilterConfig, PafFilter, RecordMeta};
use crate::mapping::ChainStatus;

/// Extract RecordMeta from .1aln file (analogous to PAF extract_metadata)
pub fn extract_1aln_metadata<P: AsRef<Path>>(path: P) -> Result<(Vec<RecordMeta>, HashMap<String, i64>)> {
    let path_str = path.as_ref().to_str().context("Invalid path")?;

    // Open .1aln reader
    let mut reader = fastga_rs::AlnReader::open(path_str)
        .context("Failed to open .1aln file")?;

    // Get all sequence names upfront (efficient bulk lookup)
    let id_to_name = reader.get_all_seq_names();

    // Build reverse mapping for writing
    let mut name_to_id = HashMap::new();
    for (id, name) in &id_to_name {
        name_to_id.insert(name.clone(), *id);
    }

    eprintln!("[unified_filter] Loaded {} sequence names from .1aln", id_to_name.len());

    // Read all alignments and create RecordMeta
    let mut metadata = Vec::new();
    let mut rank = 0;

    while let Some(aln) = reader.read_alignment()? {
        // Get actual sequence names (not numeric IDs)
        let query_id_num: i64 = aln.query_name.parse().unwrap_or(-1);
        let target_id_num: i64 = aln.target_name.parse().unwrap_or(-1);

        let query_name = id_to_name.get(&query_id_num)
            .cloned()
            .unwrap_or_else(|| aln.query_name.clone());
        let target_name = id_to_name.get(&target_id_num)
            .cloned()
            .unwrap_or_else(|| aln.target_name.clone());

        // Calculate identity
        let block_len = aln.block_len as u64;
        let matches = aln.matches as u64;
        let identity = if block_len > 0 {
            matches as f64 / block_len as f64
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
            block_length: block_len,
            identity,
            matches,
            alignment_length: block_len, // For .1aln, block_len = alignment_length
            strand: aln.strand,
            chain_id: None,
            chain_status: ChainStatus::Unassigned,
            discard: false,
            overlapped: false,
        });

        rank += 1;
    }

    eprintln!("[unified_filter] Read {} alignments from .1aln", metadata.len());

    Ok((metadata, name_to_id))
}

/// Write filtered .1aln using passing ranks
pub fn write_1aln_filtered<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    passing_ranks: &HashMap<usize, RecordMeta>,
    _name_to_id: &HashMap<String, i64>,
) -> Result<()> {
    let path_str = input_path.as_ref().to_str().context("Invalid path")?;

    // Re-open input for reading
    let mut reader = fastga_rs::AlnReader::open(path_str)?;

    // Create output writer
    let mut writer = fastga_rs::AlnWriter::create(output_path.as_ref(), true)?;

    let mut rank = 0;
    let mut written = 0;

    // Read through input, write passing records
    while let Some(aln) = reader.read_alignment()? {
        if passing_ranks.contains_key(&rank) {
            // This record passed filtering - write it
            // Note: sequence names are already numeric IDs in .1aln, so we write as-is
            writer.write_alignment(&aln)?;
            written += 1;
        }
        rank += 1;
    }

    eprintln!("[unified_filter] Wrote {} alignments to .1aln", written);
    Ok(())
}

/// Main unified filtering function - works for both .1aln and PAF
pub fn filter_file<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    config: &FilterConfig,
    force_paf_output: bool,
) -> Result<()> {
    let input_str = input_path.as_ref().to_str().context("Invalid input path")?;

    // Determine input format
    let is_1aln = input_str.ends_with(".1aln");

    if is_1aln {
        // .1aln input workflow
        eprintln!("[unified_filter] Reading .1aln metadata...");
        let (metadata, name_to_id) = extract_1aln_metadata(&input_path)?;

        // Use SAME filtering logic as PAF!
        eprintln!("[unified_filter] Applying filters...");
        let filter = PafFilter::new(config.clone());
        let passing_ranks = filter.apply_filters(metadata)?;

        eprintln!("[unified_filter] {} records passed filtering", passing_ranks.len());

        // Determine output format
        let output_str = output_path.as_ref().to_str().context("Invalid output path")?;
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
        let filter = PafFilter::new(config.clone());
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
    #[ignore] // Requires test data
    fn test_unified_1aln_filtering() {
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

        filter_file("test.1aln", "test_output.1aln", &config, false).unwrap();
    }
}
