// Tree-based sparsification for PAF alignments
// Port of allwave's knn_graph module for filtering existing alignments

use anyhow::{anyhow, Context, Result};
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

/// Sparsification strategy for PAF alignment filtering
#[derive(Debug, Clone, PartialEq)]
pub enum SparsificationStrategy {
    /// Simple fraction-based sparsification (keep X% of alignments)
    Fraction(f64),
    /// Tree-based selection: (k_nearest, k_farthest, random_fraction)
    /// - k_nearest: number of most similar neighbors (required)
    /// - k_farthest: number of most dissimilar neighbors (optional, default 0)
    /// - random_fraction: fraction of random pairs (optional, default 0.0)
    Tree(usize, usize, f64),
}

impl FromStr for SparsificationStrategy {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        // Try to parse as float first
        if let Ok(frac) = s.parse::<f64>() {
            if !(0.0..=1.0).contains(&frac) {
                return Err(anyhow!("Fraction must be between 0.0 and 1.0, got {}", frac));
            }
            return Ok(SparsificationStrategy::Fraction(frac));
        }

        // Parse tree:neighbor[,stranger[,random]]
        if let Some(value) = s.strip_prefix("tree:") {
            let parts: Vec<&str> = value.split(',').collect();
            if parts.is_empty() {
                return Err(anyhow!("tree: requires at least neighbor count"));
            }

            let k_nearest = parts[0]
                .parse::<usize>()
                .context("tree: neighbor must be a positive integer")?;

            let k_farthest = if parts.len() > 1 {
                parts[1]
                    .parse::<usize>()
                    .context("tree: stranger must be a positive integer")?
            } else {
                0
            };

            let rand_frac = if parts.len() > 2 {
                let f = parts[2]
                    .parse::<f64>()
                    .context("tree: random must be a number between 0.0 and 1.0")?;
                if !(0.0..=1.0).contains(&f) {
                    return Err(anyhow!("tree: random must be between 0.0 and 1.0"));
                }
                f
            } else {
                0.0
            };

            return Ok(SparsificationStrategy::Tree(k_nearest, k_farthest, rand_frac));
        }

        Err(anyhow!("Invalid sparsification pattern '{}'. Use a fraction (0.0-1.0) or tree:neighbor[,stranger[,random]]", s))
    }
}

/// Extract PanSN genome prefix from sequence name
/// Example: "HG002#1#chr1" -> "HG002#1#"
fn extract_genome_prefix(seq_name: &str) -> String {
    // Look for PanSN format: genome#haplotype#chromosome
    let parts: Vec<&str> = seq_name.split('#').collect();
    if parts.len() >= 2 {
        format!("{}#{}#", parts[0], parts[1])
    } else {
        // Fallback: use whole name as genome
        seq_name.to_string()
    }
}

/// PAF alignment record (minimal fields needed for filtering)
#[derive(Debug, Clone)]
pub struct PafAlignment {
    pub query_name: String,
    pub target_name: String,
    pub matches: u64,
    pub block_length: u64,
    pub identity: f64,
}

/// Build genome-pair identity matrix from PAF alignments
/// Returns map of (genome1, genome2) -> weighted average identity
pub fn build_identity_matrix(alignments: &[PafAlignment]) -> HashMap<(String, String), f64> {
    let mut genome_pairs: HashMap<(String, String), (f64, f64)> = HashMap::new();

    for aln in alignments {
        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        // Canonical key (sorted order)
        let key = if query_genome < target_genome {
            (query_genome, target_genome)
        } else {
            (target_genome, query_genome)
        };

        // Accumulate total matches and total block length for weighted calculation
        let entry = genome_pairs.entry(key).or_insert((0.0, 0.0));
        entry.0 += aln.matches as f64;
        entry.1 += aln.block_length as f64;
    }

    // Calculate weighted average ANI for each genome pair
    genome_pairs
        .into_iter()
        .map(|(key, (total_matches, total_length))| {
            let identity = if total_length > 0.0 {
                total_matches / total_length
            } else {
                0.0
            };
            (key, identity)
        })
        .collect()
}

/// Select k-nearest and k-farthest neighbors for each genome
/// Returns set of (genome1, genome2) pairs to keep (canonical order)
pub fn select_tree_pairs(
    identity_matrix: &HashMap<(String, String), f64>,
    k_nearest: usize,
    k_farthest: usize,
) -> HashSet<(String, String)> {
    // Collect all unique genomes
    let mut genomes = HashSet::new();
    for (g1, g2) in identity_matrix.keys() {
        genomes.insert(g1.clone());
        genomes.insert(g2.clone());
    }

    let mut selected_pairs = HashSet::new();

    // For each genome, select k-nearest and k-farthest neighbors
    for genome in &genomes {
        // Get all pairs involving this genome
        let mut neighbors: Vec<(&String, f64)> = identity_matrix
            .iter()
            .filter_map(|((g1, g2), &identity)| {
                if g1 == genome {
                    Some((g2, identity))
                } else if g2 == genome {
                    Some((g1, identity))
                } else {
                    None
                }
            })
            .collect();

        // Sort by identity (descending for nearest, ascending for farthest)
        neighbors.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        // Select k-nearest (highest identity)
        for (neighbor, _) in neighbors.iter().take(k_nearest) {
            let pair = if genome < *neighbor {
                (genome.clone(), (*neighbor).clone())
            } else {
                ((*neighbor).clone(), genome.clone())
            };
            selected_pairs.insert(pair);
        }

        // Select k-farthest (lowest identity)
        if k_farthest > 0 {
            neighbors.reverse();
            for (neighbor, _) in neighbors.iter().take(k_farthest) {
                let pair = if genome < *neighbor {
                    (genome.clone(), (*neighbor).clone())
                } else {
                    ((*neighbor).clone(), genome.clone())
                };
                selected_pairs.insert(pair);
            }
        }
    }

    selected_pairs
}

/// Filter PAF alignments using tree-based sparsification
/// Returns indices of alignments to keep
pub fn filter_tree_based(
    alignments: &[PafAlignment],
    k_nearest: usize,
    k_farthest: usize,
    _random_fraction: f64, // TODO: implement random pairs
) -> Vec<usize> {
    // Build identity matrix from alignments
    let identity_matrix = build_identity_matrix(alignments);

    // Select genome pairs to keep
    let selected_pairs = select_tree_pairs(&identity_matrix, k_nearest, k_farthest);

    // Filter alignments based on selected pairs
    let mut keep_indices = Vec::new();
    for (i, aln) in alignments.iter().enumerate() {
        let query_genome = extract_genome_prefix(&aln.query_name);
        let target_genome = extract_genome_prefix(&aln.target_name);

        // Skip self-comparisons
        if query_genome == target_genome {
            continue;
        }

        // Check if this pair is selected
        let pair = if query_genome < target_genome {
            (query_genome, target_genome)
        } else {
            (target_genome, query_genome)
        };

        if selected_pairs.contains(&pair) {
            keep_indices.push(i);
        }
    }

    keep_indices
}

/// Apply tree-based sparsification to a PAF file
/// Reads PAF, filters based on tree strategy, writes to output file
pub fn apply_tree_filter_to_paf(
    input_path: &str,
    output_path: &str,
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> Result<()> {
    use std::fs::File;
    use std::io::{BufRead, BufReader, Write};

    // Parse PAF file into alignments
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    let mut alignments = Vec::new();
    let mut paf_lines = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Parse PAF fields
        let query_name = fields[0].to_string();
        let target_name = fields[5].to_string();
        let matches = fields[9].parse::<u64>().unwrap_or(0);
        let block_length = fields[10].parse::<u64>().unwrap_or(1);
        let identity = matches as f64 / block_length as f64;

        alignments.push(PafAlignment {
            query_name,
            target_name,
            matches,
            block_length,
            identity,
        });
        paf_lines.push(line);
    }

    eprintln!("[sweepga] Tree filtering: {} total alignments", alignments.len());

    // Apply tree filtering
    let keep_indices = filter_tree_based(&alignments, k_nearest, k_farthest, random_fraction);

    eprintln!("[sweepga] Tree filtering: keeping {} alignments (tree:{}{}{})",
        keep_indices.len(),
        k_nearest,
        if k_farthest > 0 { format!(",{}", k_farthest) } else { String::new() },
        if random_fraction > 0.0 { format!(",{}", random_fraction) } else { String::new() }
    );

    // Write filtered PAF
    let mut output = File::create(output_path)?;
    for &idx in &keep_indices {
        writeln!(output, "{}", paf_lines[idx])?;
    }
    output.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fraction() {
        let s = SparsificationStrategy::from_str("0.5").unwrap();
        assert_eq!(s, SparsificationStrategy::Fraction(0.5));
    }

    #[test]
    fn test_parse_tree() {
        let s = SparsificationStrategy::from_str("tree:3").unwrap();
        assert_eq!(s, SparsificationStrategy::Tree(3, 0, 0.0));

        let s = SparsificationStrategy::from_str("tree:3,2").unwrap();
        assert_eq!(s, SparsificationStrategy::Tree(3, 2, 0.0));

        let s = SparsificationStrategy::from_str("tree:3,2,0.1").unwrap();
        assert_eq!(s, SparsificationStrategy::Tree(3, 2, 0.1));
    }

    #[test]
    fn test_extract_genome_prefix() {
        assert_eq!(extract_genome_prefix("HG002#1#chr1"), "HG002#1#");
        assert_eq!(extract_genome_prefix("HG002#2#chr2"), "HG002#2#");
        assert_eq!(extract_genome_prefix("NA12878#1#chrX"), "NA12878#1#");
        assert_eq!(extract_genome_prefix("simple"), "simple");
    }
}
