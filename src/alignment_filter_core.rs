use rand::Rng;
/// Core alignment filtering logic - format agnostic
///
/// This module contains the core filtering algorithms that work on abstract
/// alignment representations, independent of file format (.1aln or PAF).
use std::collections::{HashMap, HashSet};

/// Abstract representation of an alignment for filtering purposes
#[derive(Debug, Clone)]
pub struct AlignmentInfo {
    pub query_genome: String,
    pub target_genome: String,
    pub identity: f64,
    pub alignment_length: usize,
}

/// Extract genome prefix from PanSN format sequence name
/// Format: SAMPLE#HAPLOTYPE#CONTIG -> SAMPLE
pub fn extract_genome_prefix(name: &str) -> String {
    name.split('#').next().unwrap_or(name).to_string()
}

/// Build genome-pair identity matrix from alignments
///
/// Returns map of (genome1, genome2) -> weighted average identity
pub fn build_identity_matrix(alignments: &[AlignmentInfo]) -> HashMap<(String, String), f64> {
    let mut genome_pairs: HashMap<(String, String), (f64, f64)> = HashMap::new();

    for aln in alignments {
        // Skip self-comparisons
        if aln.query_genome == aln.target_genome {
            continue;
        }

        // Canonical ordering
        let key = if aln.query_genome < aln.target_genome {
            (aln.query_genome.clone(), aln.target_genome.clone())
        } else {
            (aln.target_genome.clone(), aln.query_genome.clone())
        };

        let entry = genome_pairs.entry(key).or_insert((0.0, 0.0));
        entry.0 += aln.identity * aln.alignment_length as f64; // weighted matches
        entry.1 += aln.alignment_length as f64; // total length
    }

    // Compute weighted average identity
    genome_pairs
        .into_iter()
        .map(|(key, (weighted_matches, total_length))| {
            let avg_identity = if total_length > 0.0 {
                weighted_matches / total_length
            } else {
                0.0
            };
            (key, avg_identity)
        })
        .collect()
}

/// Select genome pairs using tree-based k-nearest neighbor strategy
pub fn select_tree_pairs(
    identity_matrix: &HashMap<(String, String), f64>,
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> HashSet<(String, String)> {
    // Collect all unique genomes
    let mut genomes = HashSet::new();
    for (g1, g2) in identity_matrix.keys() {
        genomes.insert(g1.clone());
        genomes.insert(g2.clone());
    }

    let mut selected = HashSet::new();

    // For each genome, select k nearest and k farthest neighbors
    for genome in &genomes {
        // Get all pairs involving this genome
        let mut pairs: Vec<_> = identity_matrix
            .iter()
            .filter(|((g1, g2), _)| g1 == genome || g2 == genome)
            .map(|((g1, g2), &identity)| {
                let other = if g1 == genome { g2 } else { g1 };
                (other.clone(), identity, (g1.clone(), g2.clone()))
            })
            .collect();

        // Sort by identity (descending for nearest, ascending for farthest)
        pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        // Select k nearest (highest identity)
        for (_, _, pair) in pairs.iter().take(k_nearest) {
            selected.insert(pair.clone());
        }

        // Select k farthest (lowest identity)
        if k_farthest > 0 {
            for (_, _, pair) in pairs.iter().rev().take(k_farthest) {
                selected.insert(pair.clone());
            }
        }
    }

    // Add random fraction of remaining pairs
    if random_fraction > 0.0 {
        let mut rng = rand::thread_rng();
        for (pair, _) in identity_matrix.iter() {
            if !selected.contains(pair) && rng.gen::<f64>() < random_fraction {
                selected.insert(pair.clone());
            }
        }
    }

    selected
}

/// Compute edge probability for Erdős-Rényi giant component
pub fn compute_giant_component_probability(n_genomes: usize, connectivity_prob: f64) -> f64 {
    if n_genomes <= 1 {
        return 1.0;
    }

    let x = connectivity_prob.clamp(0.001, 0.999);

    // Small graph heuristics
    if n_genomes <= 10 {
        return match n_genomes {
            2 => 1.0,
            3 => 0.8,
            4 => 0.7,
            5 => 0.6,
            _ => 0.5,
        };
    }

    let n_f = n_genomes as f64;
    let log_n = n_f.ln();
    let c = -(-x.ln()).ln();
    let p = (log_n + c) / n_f;

    p.clamp(0.001, 1.0)
}

/// Select genome pairs using giant component (connectivity) strategy
pub fn select_giant_component_pairs(
    identity_matrix: &HashMap<(String, String), f64>,
    connectivity_prob: f64,
) -> HashSet<(String, String)> {
    // Count unique genomes
    let mut genomes = HashSet::new();
    for (g1, g2) in identity_matrix.keys() {
        genomes.insert(g1);
        genomes.insert(g2);
    }

    let n_genomes = genomes.len();
    let edge_prob = compute_giant_component_probability(n_genomes, connectivity_prob);

    // Random sampling
    let mut rng = rand::thread_rng();
    identity_matrix
        .keys()
        .filter(|_| rng.gen::<f64>() < edge_prob)
        .cloned()
        .collect()
}

/// Filter alignments based on selected genome pairs
///
/// Returns indices of alignments to keep
pub fn filter_by_genome_pairs(
    alignments: &[AlignmentInfo],
    selected_pairs: &HashSet<(String, String)>,
) -> Vec<usize> {
    alignments
        .iter()
        .enumerate()
        .filter_map(|(idx, aln)| {
            // Skip self-comparisons
            if aln.query_genome == aln.target_genome {
                return None;
            }

            // Check if this pair is selected
            let pair = if aln.query_genome < aln.target_genome {
                (aln.query_genome.clone(), aln.target_genome.clone())
            } else {
                (aln.target_genome.clone(), aln.query_genome.clone())
            };

            if selected_pairs.contains(&pair) {
                Some(idx)
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_genome_prefix() {
        assert_eq!(extract_genome_prefix("HG002#1#chr1"), "HG002");
        assert_eq!(extract_genome_prefix("GRCh38#0#chr22"), "GRCh38");
        assert_eq!(extract_genome_prefix("simple"), "simple");
    }

    #[test]
    fn test_build_identity_matrix() {
        let alignments = vec![
            AlignmentInfo {
                query_genome: "A".to_string(),
                target_genome: "B".to_string(),
                identity: 0.95,
                alignment_length: 1000,
            },
            AlignmentInfo {
                query_genome: "A".to_string(),
                target_genome: "B".to_string(),
                identity: 0.90,
                alignment_length: 1000,
            },
        ];

        let matrix = build_identity_matrix(&alignments);
        assert_eq!(matrix.len(), 1);

        let key = ("A".to_string(), "B".to_string());
        assert!(matrix.contains_key(&key));
        assert!((matrix[&key] - 0.925).abs() < 0.001); // weighted average
    }
}
