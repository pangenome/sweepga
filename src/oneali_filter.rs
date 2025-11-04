use crate::alignment_filter_core::{
    build_identity_matrix, extract_genome_prefix, filter_by_genome_pairs,
    select_giant_component_pairs, select_tree_pairs, AlignmentInfo,
};
/// .1aln format filtering using core filtering logic
///
/// This module provides filtering for FastGA's native .1aln format,
/// avoiding the need to convert to PAF.
use anyhow::Result;
use std::collections::HashSet;

/// Apply tree-based filtering to .1aln file
///
/// # Arguments
/// * `input_path` - Path to input .1aln file
/// * `output_path` - Path to output filtered .1aln file
/// * `k_nearest` - Number of nearest neighbors
/// * `k_farthest` - Number of farthest neighbors
/// * `random_fraction` - Fraction of random pairs
///
/// # Returns
/// Number of alignments kept
pub fn apply_tree_filter(
    input_path: &str,
    output_path: &str,
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> Result<usize> {
    // Read all alignments and extract info
    let mut reader = fastga_rs::AlnReader::open(input_path)?;
    let mut alignments_info = Vec::new();
    let mut total_count = 0;

    while let Some(aln) = reader.read_alignment()? {
        total_count += 1;

        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        let alignment_length = (aln.query_end - aln.query_start) as usize;
        let identity = if alignment_length > 0 {
            aln.matches as f64 / alignment_length as f64
        } else {
            0.0
        };

        alignments_info.push(AlignmentInfo {
            query_genome,
            target_genome,
            identity,
            alignment_length,
        });
    }

    eprintln!(
        "[sweepga] Tree filtering .1aln: read {} alignments",
        total_count
    );

    // Build identity matrix
    let identity_matrix = build_identity_matrix(&alignments_info);

    // Select tree pairs
    let selected_pairs =
        select_tree_pairs(&identity_matrix, k_nearest, k_farthest, random_fraction);

    // Filter alignments
    let keep_indices: HashSet<usize> = filter_by_genome_pairs(&alignments_info, &selected_pairs)
        .into_iter()
        .collect();

    eprintln!(
        "[sweepga] Tree filtering .1aln: keeping {} alignments (tree:{}{}{})",
        keep_indices.len(),
        k_nearest,
        if k_farthest > 0 {
            format!(",{}", k_farthest)
        } else {
            String::new()
        },
        if random_fraction > 0.0 {
            format!(",{}", random_fraction)
        } else {
            String::new()
        }
    );

    // Write filtered alignments
    let mut reader = fastga_rs::AlnReader::open(input_path)?;
    let mut writer = fastga_rs::AlnWriter::create_with_gdb(
        std::path::Path::new(output_path),
        std::path::Path::new(input_path),
        true, // binary
    )?;

    let mut rank = 0;
    let mut written = 0;

    while let Some(aln) = reader.read_alignment()? {
        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Always keep self-comparisons
        let should_keep = if query_genome == target_genome {
            false // Actually skip self-comparisons
        } else {
            keep_indices.contains(&rank)
        };

        if should_keep {
            writer.write_alignment(&aln)?;
            written += 1;
        }

        if query_genome != target_genome {
            rank += 1; // Only increment rank for non-self alignments
        }
    }

    writer.finalize();

    eprintln!(
        "[sweepga] Wrote {} tree-filtered alignments to .1aln",
        written
    );

    Ok(written)
}

/// Apply giant component filtering to .1aln file
///
/// # Arguments
/// * `input_path` - Path to input .1aln file
/// * `output_path` - Path to output filtered .1aln file
/// * `connectivity_prob` - Desired connectivity probability
///
/// # Returns
/// Number of alignments kept
pub fn apply_giant_component_filter(
    input_path: &str,
    output_path: &str,
    connectivity_prob: f64,
) -> Result<usize> {
    // Read all alignments and extract info
    let mut reader = fastga_rs::AlnReader::open(input_path)?;
    let mut alignments_info = Vec::new();
    let mut total_count = 0;

    while let Some(aln) = reader.read_alignment()? {
        total_count += 1;

        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        let alignment_length = (aln.query_end - aln.query_start) as usize;
        let identity = if alignment_length > 0 {
            aln.matches as f64 / alignment_length as f64
        } else {
            0.0
        };

        alignments_info.push(AlignmentInfo {
            query_genome,
            target_genome,
            identity,
            alignment_length,
        });
    }

    eprintln!(
        "[sweepga] Giant component filtering .1aln: read {} alignments",
        total_count
    );

    // Build identity matrix
    let identity_matrix = build_identity_matrix(&alignments_info);

    // Count unique genomes
    let mut genomes = HashSet::new();
    for (g1, g2) in identity_matrix.keys() {
        genomes.insert(g1);
        genomes.insert(g2);
    }

    // Select pairs using giant component strategy
    let selected_pairs = select_giant_component_pairs(&identity_matrix, connectivity_prob);

    eprintln!(
        "[sweepga] Giant component filtering .1aln: {} genomes, connectivity_prob={:.3}",
        genomes.len(),
        connectivity_prob
    );

    // Filter alignments
    let keep_indices: HashSet<usize> = filter_by_genome_pairs(&alignments_info, &selected_pairs)
        .into_iter()
        .collect();

    eprintln!(
        "[sweepga] Giant component filtering .1aln: keeping {} alignments ({:.2}%)",
        keep_indices.len(),
        (keep_indices.len() as f64 / alignments_info.len() as f64) * 100.0
    );

    // Write filtered alignments
    let mut reader = fastga_rs::AlnReader::open(input_path)?;
    let mut writer = fastga_rs::AlnWriter::create_with_gdb(
        std::path::Path::new(output_path),
        std::path::Path::new(input_path),
        true, // binary
    )?;

    let mut rank = 0;
    let mut written = 0;

    while let Some(aln) = reader.read_alignment()? {
        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Skip self-comparisons
        let should_keep = if query_genome == target_genome {
            false
        } else {
            keep_indices.contains(&rank)
        };

        if should_keep {
            writer.write_alignment(&aln)?;
            written += 1;
        }

        if query_genome != target_genome {
            rank += 1;
        }
    }

    writer.finalize();

    eprintln!(
        "[sweepga] Wrote {} giant-component-filtered alignments to .1aln",
        written
    );

    Ok(written)
}
