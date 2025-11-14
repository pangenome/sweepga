/// Giant component (connectivity) based filtering using Erdős-Rényi random graph theory
///
/// This module implements sparsification that maintains graph connectivity by randomly
/// sampling edges with a probability calculated to ensure a single giant component.
use anyhow::Result;
use rand::Rng;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Compute edge probability for Erdős-Rényi random graph giant component
///
/// For a random graph G(n,p) to maintain a single giant component with probability
/// `connectivity_prob`, we need: p = (log n - log(-log(connectivity_prob)))/n
///
/// This is based on the sharp threshold result:
/// If p = (log n + c)/n, then P_connected → e^(-e^(-c)) as n → ∞
///
/// # Arguments
/// * `n` - Number of nodes (sequences)
/// * `connectivity_prob` - Desired probability that graph has giant component (0 < x < 1)
///
/// # Returns
/// Edge probability p for the random graph
pub fn compute_connectivity_probability(n: usize, connectivity_prob: f64) -> f64 {
    if n <= 1 {
        return 1.0;
    }

    // Clamp connectivity probability to reasonable range
    let x = connectivity_prob.clamp(0.001, 0.999);

    // For very small n, use simpler heuristics
    if n <= 10 {
        return match n {
            2 => 1.0, // 2 nodes: need the single edge
            3 => 0.8, // 3 nodes: need high probability
            4 => 0.7, // 4 nodes
            5 => 0.6, // 5 nodes
            _ => 0.5, // 6-10 nodes
        };
    }

    let n_f = n as f64;
    let log_n = n_f.ln();

    // c = -log(-log(x))
    let c = -(-x.ln()).ln();

    // p = (log n + c) / n
    let p = (log_n + c) / n_f;

    // Ensure p is reasonable (not too small or too large)
    p.clamp(0.001, 1.0)
}

/// Extract genome prefix from PanSN format sequence name
/// Format: SAMPLE#HAPLOTYPE#CONTIG or GENOME#...
fn extract_genome_prefix(name: &str) -> &str {
    name.split('#').next().unwrap_or(name)
}

/// Apply giant component (connectivity) filtering to PAF file
///
/// Randomly samples alignments based on Erdős-Rényi random graph theory to maintain
/// a single giant component with the specified probability.
///
/// # Arguments
/// * `input_path` - Path to input PAF file
/// * `output_path` - Path to output filtered PAF file
/// * `connectivity_prob` - Desired probability of maintaining giant component (0 < x < 1)
///
/// # Returns
/// Result with number of alignments kept
pub fn apply_giant_component_filter_to_paf(
    input_path: &str,
    output_path: &str,
    connectivity_prob: f64,
) -> Result<usize> {
    // Read PAF to count unique genomes
    let input_file = File::open(input_path)?;
    let reader = BufReader::new(input_file);

    let mut genomes = std::collections::HashSet::new();
    let mut total_alignments = 0;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        total_alignments += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let query_genome = extract_genome_prefix(fields[0]);
        let target_genome = extract_genome_prefix(fields[5]);

        genomes.insert(query_genome.to_string());
        genomes.insert(target_genome.to_string());
    }

    let n_genomes = genomes.len();
    let edge_prob = compute_connectivity_probability(n_genomes, connectivity_prob);

    // eprintln!(
    //     "[sweepga] Giant component filtering: {} genomes, connectivity_prob={:.3}, edge_prob={:.4}",
    //     n_genomes, connectivity_prob, edge_prob
    // );
    // eprintln!(
    //     "[sweepga] Giant component filtering: {} total alignments",
    //     total_alignments
    // );

    // Re-read and filter with random sampling
    let input_file = File::open(input_path)?;
    let reader = BufReader::new(input_file);
    let mut output_file = File::create(output_path)?;

    let mut rng = rand::thread_rng();
    let mut kept = 0;

    for line in reader.lines() {
        let line = line?;

        // Keep comments and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            writeln!(output_file, "{}", line)?;
            continue;
        }

        // Randomly sample with edge probability
        if rng.gen::<f64>() < edge_prob {
            writeln!(output_file, "{}", line)?;
            kept += 1;
        }
    }

    // eprintln!(
    //     "[sweepga] Giant component filtering: keeping {} alignments ({:.2}%)",
    //     kept,
    //     (kept as f64 / total_alignments as f64) * 100.0
    // );

    Ok(kept)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_connectivity_probability() {
        // Test small graphs
        assert_eq!(compute_connectivity_probability(2, 0.95), 1.0);
        assert_eq!(compute_connectivity_probability(3, 0.95), 0.8);

        // Test larger graphs - should decrease with size
        let p_100 = compute_connectivity_probability(100, 0.95);
        let p_1000 = compute_connectivity_probability(1000, 0.95);
        assert!(p_100 > p_1000);

        // Test higher connectivity requires higher edge probability
        let p_95 = compute_connectivity_probability(100, 0.95);
        let p_99 = compute_connectivity_probability(100, 0.99);
        assert!(p_99 > p_95);
    }

    #[test]
    fn test_extract_genome_prefix() {
        assert_eq!(extract_genome_prefix("HG002#1#chr1"), "HG002");
        assert_eq!(extract_genome_prefix("GRCh38#0#chr22"), "GRCh38");
        assert_eq!(extract_genome_prefix("simple"), "simple");
    }
}
