//! K-nearest neighbor graph construction for intelligent sequence alignment sampling
//!
//! This module implements k-nearest neighbor graphs based on mash distances between sequences.
//! Each sequence is connected to its k most similar sequences, creating a structured graph
//! that maintains good connectivity while being much sparser than all-vs-all comparisons.
//!
//! Ported from allwave.

#![allow(dead_code)]

use crate::mash;
use std::collections::HashSet;

/// Mash distance parameters shared by all strategies that use sketching.
#[derive(Debug, Clone, PartialEq)]
pub struct MashParams {
    pub kmer_size: usize,
    pub sketch_size: usize,
}

impl Default for MashParams {
    fn default() -> Self {
        MashParams {
            kmer_size: mash::DEFAULT_KMER_SIZE,
            sketch_size: mash::DEFAULT_SKETCH_SIZE,
        }
    }
}

/// Sparsification strategy for pair selection
#[derive(Debug, Clone, PartialEq)]
pub enum SparsificationStrategy {
    /// No sparsification - all pairs
    None,
    /// Automatic selection based on sample count
    Auto,
    /// Random fraction of pairs
    Random(f64),
    /// Giant component connectivity guarantee
    Connectivity(f64),
    /// Tree-based: k_nearest, k_farthest, random_fraction
    TreeSampling(usize, usize, f64),
    /// Wfmash mapping density (-x flag). None = auto (ln(n)/n*10), Some = explicit fraction.
    /// Only valid with wfmash backend; errors if used with fastga.
    WfmashDensity(Option<f64>),
}

impl std::str::FromStr for SparsificationStrategy {
    type Err = String;

    /// Parse a `--sparsify` CLI value. Every variant uses `:` as separator.
    ///
    /// - `none` / `all`
    /// - `auto`
    /// - `<frac>` or `random:<frac>`                       (0 < frac ≤ 1)
    /// - `giant:<prob>` / `connectivity:<prob>`            (0 < prob < 1)
    /// - `tree:<near>[:<far>[:<random>]]` / `knn:…`
    /// - `wfmash:auto` / `wfmash:<frac>`
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Bare float shorthand for `random:<frac>`.
        if let Ok(frac) = s.parse::<f64>() {
            if frac <= 0.0 || frac > 1.0 {
                return Err("Bare-float sparsify must be between 0 and 1".to_string());
            }
            return Ok(SparsificationStrategy::Random(frac));
        }

        match s {
            "none" | "all" => Ok(SparsificationStrategy::None),
            "auto" => Ok(SparsificationStrategy::Auto),
            s if s.starts_with("random:") => {
                let fraction: f64 = s[7..]
                    .parse()
                    .map_err(|_| "Invalid random fraction".to_string())?;
                if fraction <= 0.0 || fraction > 1.0 {
                    return Err("Random fraction must be between 0 and 1".to_string());
                }
                Ok(SparsificationStrategy::Random(fraction))
            }
            s if s.starts_with("giant:") || s.starts_with("connectivity:") => {
                let colon_pos = s.find(':').unwrap();
                let prob: f64 = s[colon_pos + 1..]
                    .parse()
                    .map_err(|_| "Invalid giant component probability".to_string())?;
                if prob <= 0.0 || prob >= 1.0 {
                    return Err("Giant component probability must be between 0 and 1".to_string());
                }
                Ok(SparsificationStrategy::Connectivity(prob))
            }
            s if s.starts_with("tree:") || s.starts_with("knn:") => {
                let colon_pos = s.find(':').unwrap();
                let body = &s[colon_pos + 1..];
                let parts: Vec<&str> = body.split(':').collect();
                if parts.is_empty() || parts.len() > 3 {
                    return Err(
                        "Invalid tree format. Use: tree:<k_nearest>[:<k_farthest>[:<random_fraction>]]"
                            .to_string(),
                    );
                }

                let k_nearest: usize = parts[0]
                    .parse()
                    .map_err(|_| "Invalid k nearest count".to_string())?;
                let k_farthest: usize = if parts.len() > 1 {
                    parts[1]
                        .parse()
                        .map_err(|_| "Invalid k farthest count".to_string())?
                } else {
                    0
                };
                let random_frac: f64 = if parts.len() > 2 {
                    parts[2]
                        .parse()
                        .map_err(|_| "Invalid random fraction".to_string())?
                } else {
                    0.0
                };

                if k_nearest == 0 && k_farthest == 0 {
                    return Err(
                        "At least one of k_nearest or k_farthest must be greater than 0"
                            .to_string(),
                    );
                }
                if !(0.0..=1.0).contains(&random_frac) {
                    return Err("Random fraction must be between 0 and 1".to_string());
                }

                Ok(SparsificationStrategy::TreeSampling(
                    k_nearest,
                    k_farthest,
                    random_frac,
                ))
            }
            s if s.starts_with("wfmash:") => {
                let val = &s[7..];
                if val == "auto" {
                    Ok(SparsificationStrategy::WfmashDensity(None))
                } else {
                    let frac: f64 = val
                        .parse()
                        .map_err(|_| "Invalid wfmash density fraction".to_string())?;
                    if frac <= 0.0 || frac > 1.0 {
                        return Err(
                            "Wfmash density fraction must be between 0 and 1".to_string()
                        );
                    }
                    Ok(SparsificationStrategy::WfmashDensity(Some(frac)))
                }
            }
            _ => Err(format!(
                "Invalid sparsification strategy '{}'. Use: none, all, auto, giant:<probability>, connectivity:<probability>, random:<fraction>, tree:<near>:<far>:<random>, knn:<near>:<far>:<random>, wfmash:auto, or wfmash:<fraction>",
                s
            )),
        }
    }
}

impl std::fmt::Display for SparsificationStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SparsificationStrategy::None => write!(f, "none"),
            SparsificationStrategy::Auto => write!(f, "auto"),
            SparsificationStrategy::Random(frac) => write!(f, "random:{}", frac),
            SparsificationStrategy::Connectivity(prob) => write!(f, "giant:{}", prob),
            SparsificationStrategy::TreeSampling(near, far, rand) => {
                write!(f, "tree:{}:{}:{}", near, far, rand)
            }
            SparsificationStrategy::WfmashDensity(frac) => match frac {
                Some(f_val) => write!(f, "wfmash:{}", f_val),
                None => write!(f, "wfmash:auto"),
            },
        }
    }
}

impl SparsificationStrategy {
    /// Get a human-readable description of the strategy for logging.
    pub fn description(&self) -> String {
        match self {
            SparsificationStrategy::None => "all pairs (no sparsification)".to_string(),
            SparsificationStrategy::Auto => "auto (context-dependent)".to_string(),
            SparsificationStrategy::Random(f) => format!("random {:.1}%", f * 100.0),
            SparsificationStrategy::Connectivity(p) => {
                format!("giant component (p={:.2})", p)
            }
            SparsificationStrategy::TreeSampling(near, far, rand) => {
                format!(
                    "tree sampling (k_near={}, k_far={}, rand={:.1}%)",
                    near,
                    far,
                    rand * 100.0,
                )
            }
            SparsificationStrategy::WfmashDensity(frac) => match frac {
                Some(f) => format!("wfmash density ({:.1}%)", f * 100.0),
                None => "wfmash density (auto: ln(n)/n*10)".to_string(),
            },
        }
    }

    /// Compute the wfmash auto-density fraction for a given number of haplotypes.
    /// Uses the pggb heuristic: ln(n)/n * 10.
    /// Returns None if the fraction >= 1.0 (no sparsification needed).
    pub fn wfmash_auto_density(n_haps: usize) -> Option<f64> {
        if n_haps <= 1 {
            return None;
        }
        let n = n_haps as f64;
        let frac = n.ln() / n * 10.0;
        if frac >= 1.0 {
            None
        } else {
            Some(frac)
        }
    }
}

/// Extract sequence pairs using k-nearest neighbor graph with optional stranger-joining
/// Returns pairs as (index_i, index_j) tuples
pub fn extract_tree_pairs(
    sequences: &[Vec<u8>],
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
    mash_params: &MashParams,
) -> Vec<(usize, usize)> {
    if sequences.len() < 2 {
        return Vec::new();
    }

    // Compute distance matrix
    let distance_matrix = mash::compute_distance_matrix_with_params(
        sequences,
        mash_params.kmer_size,
        mash_params.sketch_size,
    );

    extract_tree_pairs_from_matrix(&distance_matrix, k_nearest, k_farthest, random_fraction)
}

/// Extract tree pairs from pre-computed distance matrix
pub fn extract_tree_pairs_from_matrix(
    distance_matrix: &[Vec<f64>],
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> Vec<(usize, usize)> {
    let n = distance_matrix.len();
    if n < 2 {
        return Vec::new();
    }

    let mut all_pairs = Vec::new();

    if k_nearest > 0 {
        let nearest_pairs = build_knn_graph(distance_matrix, k_nearest, false);
        all_pairs.extend(nearest_pairs);
    }

    if k_farthest > 0 {
        let farthest_pairs = build_knn_graph(distance_matrix, k_farthest, true);
        all_pairs.extend(farthest_pairs);
    }

    if random_fraction > 0.0 {
        let random_pairs = generate_random_pairs(n, random_fraction);
        all_pairs.extend(random_pairs);
    }

    // Normalize to canonical order and remove duplicates
    let mut canonical: Vec<(usize, usize)> = all_pairs
        .into_iter()
        .map(|(i, j)| if i < j { (i, j) } else { (j, i) })
        .collect();
    canonical.sort_unstable();
    canonical.dedup();

    canonical
}

/// Extract tree pairs and random pairs separately (for iterative alignment with early stopping)
/// Returns (tree_pairs, random_pairs) where tree_pairs should be processed first to guarantee connectivity
#[allow(clippy::type_complexity)]
pub fn extract_tree_pairs_separated(
    sequences: &[Vec<u8>],
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
    mash_params: &MashParams,
) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
    if sequences.len() < 2 {
        return (Vec::new(), Vec::new());
    }

    let distance_matrix = mash::compute_distance_matrix_with_params(
        sequences,
        mash_params.kmer_size,
        mash_params.sketch_size,
    );

    let mut tree_pairs = Vec::new();

    if k_nearest > 0 {
        let nearest_pairs = build_knn_graph(&distance_matrix, k_nearest, false);
        tree_pairs.extend(nearest_pairs);
    }

    if k_farthest > 0 {
        let farthest_pairs = build_knn_graph(&distance_matrix, k_farthest, true);
        tree_pairs.extend(farthest_pairs);
    }

    // Normalize and dedupe tree pairs
    let tree_set: HashSet<(usize, usize)> = tree_pairs
        .into_iter()
        .map(|(i, j)| if i < j { (i, j) } else { (j, i) })
        .collect();

    let tree_pairs: Vec<(usize, usize)> = tree_set.iter().copied().collect();

    // Generate random pairs
    let mut random_pairs = if random_fraction > 0.0 {
        generate_random_pairs(sequences.len(), random_fraction)
    } else {
        Vec::new()
    };

    // Remove any random pairs that are already in tree pairs
    random_pairs.retain(|pair| {
        let canonical = if pair.0 < pair.1 {
            *pair
        } else {
            (pair.1, pair.0)
        };
        !tree_set.contains(&canonical)
    });

    (tree_pairs, random_pairs)
}

/// Build k-nearest or k-farthest neighbor graph from distance matrix
fn build_knn_graph(
    distance_matrix: &[Vec<f64>],
    k_neighbors: usize,
    farthest: bool,
) -> Vec<(usize, usize)> {
    let n = distance_matrix.len();
    let mut pairs = Vec::new();

    for (i, row) in distance_matrix.iter().enumerate() {
        let mut neighbors: Vec<(f64, usize)> =
            (0..n).filter(|&j| i != j).map(|j| (row[j], j)).collect();

        if farthest {
            neighbors.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        } else {
            neighbors.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        }

        let k_actual = k_neighbors.min(neighbors.len());
        for &(_, j) in neighbors.iter().take(k_actual) {
            pairs.push((i, j));
        }
    }

    pairs
}

/// Generate random pairs with deterministic hashing based on indices
pub fn generate_random_pairs(n: usize, fraction: f64) -> Vec<(usize, usize)> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::Hasher;

    let mut pairs = Vec::new();
    let threshold = (fraction * u64::MAX as f64) as u64;

    for i in 0..n {
        for j in i + 1..n {
            let mut hasher = DefaultHasher::new();
            hasher.write_usize(i);
            hasher.write_usize(j);
            let hash = hasher.finish();

            if hash <= threshold {
                pairs.push((i, j));
            }
        }
    }

    pairs
}

/// Estimate the number of pairs in a tree-based sampling strategy
pub fn estimate_tree_pair_count(
    n: usize,
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> usize {
    let nearest_pairs = n * k_nearest.min(n.saturating_sub(1));
    let farthest_pairs = n * k_farthest.min(n.saturating_sub(1));
    let total_possible = n * (n - 1) / 2; // Undirected pairs
    let random_pairs = (total_possible as f64 * random_fraction).round() as usize;
    (nearest_pairs + farthest_pairs + random_pairs).min(total_possible)
}

/// Select pairs based on sparsification strategy
/// Returns pairs as (sample_index_i, sample_index_j) tuples
pub fn select_pairs(
    sample_count: usize,
    sequences: Option<&[Vec<u8>]>,
    strategy: &SparsificationStrategy,
    mash_params: &MashParams,
) -> Vec<(usize, usize)> {
    match strategy {
        SparsificationStrategy::None => {
            // All pairs
            let mut pairs = Vec::new();
            for i in 0..sample_count {
                for j in i + 1..sample_count {
                    pairs.push((i, j));
                }
            }
            pairs
        }
        SparsificationStrategy::Auto => {
            if sample_count <= 10 {
                // Very small: all pairs
                let mut pairs = Vec::new();
                for i in 0..sample_count {
                    for j in i + 1..sample_count {
                        pairs.push((i, j));
                    }
                }
                pairs
            } else if sample_count <= 50 {
                // Medium: giant-component pair selection
                select_pairs(
                    sample_count,
                    sequences,
                    &SparsificationStrategy::Connectivity(0.99),
                    mash_params,
                )
            } else {
                // Large: tree-based k-NN sampling
                if let Some(seqs) = sequences {
                    extract_tree_pairs(seqs, 5, 2, 0.05, mash_params)
                } else {
                    generate_random_pairs(sample_count, 0.1)
                }
            }
        }
        SparsificationStrategy::Random(fraction) => generate_random_pairs(sample_count, *fraction),
        SparsificationStrategy::Connectivity(prob) => {
            // Giant component: use enough pairs to guarantee connectivity with given probability
            // For n nodes, we need roughly n*ln(n) / 2 edges for connectivity
            // Adjust by probability factor
            let target_edges = ((sample_count as f64) * (sample_count as f64).ln() / 2.0
                * (-prob.ln()))
            .ceil() as usize;
            let total_possible = sample_count * (sample_count - 1) / 2;
            let fraction = (target_edges as f64 / total_possible as f64).min(1.0);

            if let Some(seqs) = sequences {
                // Use tree sampling with enough nearest neighbors
                let k_nearest = ((fraction * sample_count as f64).ceil() as usize).max(2);
                extract_tree_pairs(seqs, k_nearest, 1, 0.01, mash_params)
            } else {
                generate_random_pairs(sample_count, fraction)
            }
        }
        SparsificationStrategy::TreeSampling(k_nearest, k_farthest, random_frac) => {
            if let Some(seqs) = sequences {
                extract_tree_pairs(seqs, *k_nearest, *k_farthest, *random_frac, mash_params)
            } else {
                // Fallback without sequences - just use random
                generate_random_pairs(sample_count, *random_frac)
            }
        }
        SparsificationStrategy::WfmashDensity(_) => {
            // WfmashDensity is a mapping-level strategy, not pair selection.
            // When used for pair selection, it means "all pairs" (density is applied
            // at the aligner level, not at pair selection).
            let mut pairs = Vec::new();
            for i in 0..sample_count {
                for j in i + 1..sample_count {
                    pairs.push((i, j));
                }
            }
            pairs
        }
    }
}

/// Select pairs based on sparsification strategy using pre-computed sketches
/// This avoids materializing all sequences in memory at once
pub fn select_pairs_from_sketches(
    sketches: &[mash::KmerSketch],
    strategy: &SparsificationStrategy,
) -> Vec<(usize, usize)> {
    let sample_count = sketches.len();

    match strategy {
        SparsificationStrategy::None => {
            let mut pairs = Vec::new();
            for i in 0..sample_count {
                for j in i + 1..sample_count {
                    pairs.push((i, j));
                }
            }
            pairs
        }
        SparsificationStrategy::Auto => {
            if sample_count <= 10 {
                let mut pairs = Vec::new();
                for i in 0..sample_count {
                    for j in i + 1..sample_count {
                        pairs.push((i, j));
                    }
                }
                pairs
            } else if sample_count <= 50 {
                select_pairs_from_sketches(
                    sketches,
                    &SparsificationStrategy::Connectivity(0.99),
                )
            } else {
                let distance_matrix = mash::distance_matrix_from_sketches(sketches);
                extract_tree_pairs_from_matrix(&distance_matrix, 5, 2, 0.05)
            }
        }
        SparsificationStrategy::Random(fraction) => generate_random_pairs(sample_count, *fraction),
        SparsificationStrategy::Connectivity(prob) => {
            let target_edges = ((sample_count as f64) * (sample_count as f64).ln() / 2.0
                * (-prob.ln()))
            .ceil() as usize;
            let total_possible = sample_count * (sample_count - 1) / 2;
            let fraction = (target_edges as f64 / total_possible as f64).min(1.0);

            let k_nearest = ((fraction * sample_count as f64).ceil() as usize).max(2);
            let distance_matrix = mash::distance_matrix_from_sketches(sketches);
            extract_tree_pairs_from_matrix(&distance_matrix, k_nearest, 1, 0.01)
        }
        SparsificationStrategy::TreeSampling(k_nearest, k_farthest, random_frac) => {
            let distance_matrix = mash::distance_matrix_from_sketches(sketches);
            extract_tree_pairs_from_matrix(&distance_matrix, *k_nearest, *k_farthest, *random_frac)
        }
        SparsificationStrategy::WfmashDensity(_) => {
            // WfmashDensity = all pairs (density applied at aligner level)
            let mut pairs = Vec::new();
            for i in 0..sample_count {
                for j in i + 1..sample_count {
                    pairs.push((i, j));
                }
            }
            pairs
        }
    }
}

/// Merge MinHash sketches by unioning their minimizer vectors and
/// re-truncating to `sketch_size`. MinHash bottom-k has the mergeability
/// property that the bottom-k of the union equals the bottom-k of the
/// merged bottom-k sketches — so this produces exactly the sketch we
/// would have gotten had we computed it over the concatenated sequence,
/// without re-hashing.
pub fn merge_sketches(
    parts: &[&mash::KmerSketch],
    sketch_size: usize,
) -> mash::KmerSketch {
    let k = parts.first().map(|s| s.k).unwrap_or(mash::DEFAULT_KMER_SIZE);
    let mut minimizers: Vec<u64> = parts
        .iter()
        .flat_map(|s| s.minimizers.iter().copied())
        .collect();
    let length: usize = parts.iter().map(|s| s.length).sum();
    minimizers.sort_unstable();
    minimizers.dedup();
    minimizers.truncate(sketch_size);
    mash::KmerSketch { minimizers, k, length }
}

/// Expand PanSN haplotype-level pairs into contig-level pairs:
/// the cross-product of contigs across selected haplotype pairs, plus every
/// intra-haplotype contig pair (same-sample contigs — always kept together).
///
/// Output pairs are normalized `(min, max)` and deduplicated.
pub fn expand_haplotype_pairs(
    hap_pairs: &[(usize, usize)],
    hap_groups: &[Vec<usize>],
) -> Vec<(usize, usize)> {
    use std::collections::HashSet;
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    for &(hi, hj) in hap_pairs {
        for &ci in &hap_groups[hi] {
            for &cj in &hap_groups[hj] {
                if ci == cj {
                    continue;
                }
                let pair = if ci < cj { (ci, cj) } else { (cj, ci) };
                seen.insert(pair);
            }
        }
    }
    for contigs in hap_groups {
        for a in 0..contigs.len() {
            for b in a + 1..contigs.len() {
                let (ci, cj) = (contigs[a], contigs[b]);
                let pair = if ci < cj { (ci, cj) } else { (cj, ci) };
                seen.insert(pair);
            }
        }
    }
    let mut out: Vec<(usize, usize)> = seen.into_iter().collect();
    out.sort_unstable();
    out
}

/// Sparsify at PanSN haplotype granularity, then project back to contig
/// pairs.
///
/// For pangenome inputs where "sample#haplotype#contig" names pack multiple
/// contigs per haplotype, running `select_pairs_from_sketches` directly at
/// the contig level loses biological structure — a random fraction of
/// contig pairs gives an uneven sample across haplotypes. This helper
/// groups contigs by haplotype, merges their MinHash sketches into one
/// per-haplotype sketch, runs the strategy on the haplotype graph, and
/// then expands each selected haplotype pair to the full cross-product of
/// contig pairs (plus all intra-haplotype contig pairs, which are
/// free-to-keep).
///
/// Falls back to plain contig-level `select_pairs_from_sketches` when the
/// input isn't PanSN-named (detected as: the number of distinct haplotype
/// groups equals the number of input contigs).
pub fn select_pairs_haplotype_aware(
    names: &[&str],
    contig_sketches: &[mash::KmerSketch],
    strategy: &SparsificationStrategy,
    sketch_size: usize,
) -> Vec<(usize, usize)> {
    use crate::pansn::{group_indices_by_pansn, PanSnLevel};
    let hap_groups = group_indices_by_pansn(names.iter().copied(), PanSnLevel::Haplotype);
    if hap_groups.len() == contig_sketches.len() {
        return select_pairs_from_sketches(contig_sketches, strategy);
    }
    let hap_sketches: Vec<mash::KmerSketch> = hap_groups
        .iter()
        .map(|idxs| {
            let parts: Vec<&mash::KmerSketch> =
                idxs.iter().map(|&i| &contig_sketches[i]).collect();
            merge_sketches(&parts, sketch_size)
        })
        .collect();
    let hap_pairs = select_pairs_from_sketches(&hap_sketches, strategy);
    expand_haplotype_pairs(&hap_pairs, &hap_groups)
}

/// Sketch-free variant of [`select_pairs_haplotype_aware`] for the
/// `None` / `Random(_)` / `WfmashDensity(_)` strategies that don't
/// need per-node distances. The strategy is applied at the haplotype level
/// — so `Random(0.05)` samples 5% of *haplotype pairs* (and then expands
/// each to contig pairs), which matches the biologically-meaningful unit
/// rather than the per-contig pair count that can be wildly uneven.
pub fn select_pairs_haplotype_aware_no_sketch(
    names: &[&str],
    strategy: &SparsificationStrategy,
    mash_params: &MashParams,
) -> Vec<(usize, usize)> {
    use crate::pansn::{group_indices_by_pansn, PanSnLevel};
    let n = names.len();
    let hap_groups = group_indices_by_pansn(names.iter().copied(), PanSnLevel::Haplotype);
    if hap_groups.len() == n {
        return select_pairs(n, None, strategy, mash_params);
    }
    let hap_pairs = select_pairs(hap_groups.len(), None, strategy, mash_params);
    expand_haplotype_pairs(&hap_pairs, &hap_groups)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sparsification_none() {
        let s: SparsificationStrategy = "none".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::None);
    }

    #[test]
    fn test_parse_sparsification_random() {
        let s: SparsificationStrategy = "random:0.5".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::Random(0.5));
    }

    #[test]
    fn test_parse_sparsification_giant() {
        let s: SparsificationStrategy = "giant:0.99".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::Connectivity(0.99));
    }

    #[test]
    fn test_parse_sparsification_tree() {
        let s: SparsificationStrategy = "tree:5:2:0.1".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::TreeSampling(5, 2, 0.1));
    }

    #[test]
    fn test_build_knn_graph() {
        let distances = vec![
            vec![0.0, 0.1, 0.9],
            vec![0.1, 0.0, 0.8],
            vec![0.9, 0.8, 0.0],
        ];

        let pairs = build_knn_graph(&distances, 1, false);
        assert_eq!(pairs.len(), 3);
        assert!(pairs.contains(&(0, 1))); // seq0's nearest is seq1
        assert!(pairs.contains(&(1, 0))); // seq1's nearest is seq0
    }

    #[test]
    fn test_estimate_pair_count() {
        assert_eq!(estimate_tree_pair_count(4, 1, 0, 0.0), 4);
        // With 4 samples, max possible undirected pairs = 4 choose 2 = 6
        // k_nearest=1 + k_farthest=1 would give 8 directed edges, but clamped to 6
        assert_eq!(estimate_tree_pair_count(4, 1, 1, 0.0), 6);
    }

    #[test]
    fn test_select_pairs_none() {
        let pairs = select_pairs(4, None, &SparsificationStrategy::None, &MashParams::default());
        assert_eq!(pairs.len(), 6); // 4 choose 2 = 6
    }

    #[test]
    fn test_merge_sketches_is_bottom_k_union() {
        let a = mash::KmerSketch::from_sequence(b"ACGTACGTACGTACGT", 5, 100);
        let b = mash::KmerSketch::from_sequence(b"TTTTGGGGAAAACCCC", 5, 100);
        let merged = merge_sketches(&[&a, &b], 100);
        assert_eq!(merged.k, 5);
        assert_eq!(merged.length, a.length + b.length);
        // Every merged minimizer must come from one of the parts
        let union: std::collections::HashSet<u64> =
            a.minimizers.iter().chain(b.minimizers.iter()).copied().collect();
        for m in &merged.minimizers {
            assert!(union.contains(m));
        }
        // Sketch is sorted & deduped
        for w in merged.minimizers.windows(2) {
            assert!(w[0] < w[1]);
        }
    }

    #[test]
    fn test_expand_haplotype_pairs_cross_and_intra() {
        // hap_0 = indices [0,1], hap_1 = [2,3,4]
        let groups = vec![vec![0, 1], vec![2, 3, 4]];
        let pairs = expand_haplotype_pairs(&[(0, 1)], &groups);
        // 2*3 cross + 1 intra(hap0) + 3 intra(hap1) = 10 pairs
        assert_eq!(pairs.len(), 10);
        assert!(pairs.contains(&(0, 1))); // intra hap_0
        assert!(pairs.contains(&(2, 3))); // intra hap_1
        assert!(pairs.contains(&(0, 2))); // cross
        assert!(pairs.contains(&(1, 4))); // cross
    }

    #[test]
    fn test_haplotype_aware_none_strategy_matches_all_pairs() {
        // 4 contigs across 2 haplotypes; None strategy at hap level expands
        // to all 6 contig pairs.
        let names = ["S1#1#a", "S1#1#b", "S2#1#a", "S2#1#b"];
        let pairs = select_pairs_haplotype_aware_no_sketch(
            &names,
            &SparsificationStrategy::None,
            &MashParams::default(),
        );
        assert_eq!(pairs.len(), 6);
    }

    #[test]
    fn test_haplotype_aware_non_pansn_falls_back_to_contig() {
        // Non-PanSN names: every name is its own haplotype, so behavior must
        // be identical to direct contig-level selection.
        let names = ["chr1", "chr2", "chr3"];
        let pairs = select_pairs_haplotype_aware_no_sketch(
            &names,
            &SparsificationStrategy::None,
            &MashParams::default(),
        );
        assert_eq!(pairs.len(), 3); // 3 choose 2
    }

    #[test]
    fn test_select_pairs_random() {
        // Random fraction should produce approximately that fraction of pairs
        let pairs = select_pairs(10, None, &SparsificationStrategy::Random(0.5), &MashParams::default());
        let total_possible = 10 * 9 / 2; // 45 pairs
                                         // Allow some variance due to randomness
        assert!(pairs.len() > 0);
        assert!(pairs.len() <= total_possible);
    }

    #[test]
    fn test_generate_random_pairs() {
        let pairs = generate_random_pairs(5, 0.4);
        let total_possible = 5 * 4 / 2; // 10 pairs
                                        // Should get roughly 40% = 4 pairs (allow variance)
        assert!(pairs.len() > 0);
        assert!(pairs.len() <= total_possible);

        // All pairs should be valid (i < j)
        for (i, j) in &pairs {
            assert!(i < j);
        }
    }

    #[test]
    fn test_build_knn_graph_k2() {
        // Test with k=2 nearest neighbors
        let distances = vec![
            vec![0.0, 0.1, 0.2, 0.9],
            vec![0.1, 0.0, 0.15, 0.8],
            vec![0.2, 0.15, 0.0, 0.7],
            vec![0.9, 0.8, 0.7, 0.0],
        ];

        let pairs = build_knn_graph(&distances, 2, false);
        // Each node selects 2 nearest, but some overlap
        assert!(pairs.len() >= 4); // At least k edges per node (with overlap)

        // Node 0's nearest are 1 and 2
        assert!(pairs.contains(&(0, 1)) || pairs.contains(&(1, 0)));
        assert!(pairs.contains(&(0, 2)) || pairs.contains(&(2, 0)));
    }

    #[test]
    fn test_build_knn_graph_farthest() {
        // Test farthest neighbor selection
        let distances = vec![
            vec![0.0, 0.1, 0.9],
            vec![0.1, 0.0, 0.8],
            vec![0.9, 0.8, 0.0],
        ];

        let pairs = build_knn_graph(&distances, 1, true); // farthest=true
                                                          // Node 0's farthest is 2, node 1's farthest is 2, node 2's farthest is 0
        assert!(pairs.contains(&(0, 2)) || pairs.contains(&(2, 0)));
    }

    #[test]
    fn test_extract_tree_pairs_from_matrix() {
        let distances = vec![
            vec![0.0, 0.1, 0.5, 0.9],
            vec![0.1, 0.0, 0.4, 0.8],
            vec![0.5, 0.4, 0.0, 0.3],
            vec![0.9, 0.8, 0.3, 0.0],
        ];

        // k_nearest=1, k_farthest=0, random=0
        let pairs = extract_tree_pairs_from_matrix(&distances, 1, 0, 0.0);
        assert!(!pairs.is_empty());

        // All pairs should be valid
        for (i, j) in &pairs {
            assert_ne!(i, j);
        }
    }

    #[test]
    fn test_parse_sparsification_aliases() {
        // "all" is an alias for None
        let s: SparsificationStrategy = "all".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::None);

        // "connectivity:" is an alias for "giant:"
        let s: SparsificationStrategy = "connectivity:0.95".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::Connectivity(0.95));

        // "knn:" is an alias for "tree:"
        let s: SparsificationStrategy = "knn:3:1:0.1".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::TreeSampling(3, 1, 0.1));
    }

    #[test]
    fn test_parse_sparsification_wfmash() {
        let s: SparsificationStrategy = "wfmash:auto".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::WfmashDensity(None));

        let s: SparsificationStrategy = "wfmash:0.3".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::WfmashDensity(Some(0.3)));

        // Invalid: out of range
        let result: Result<SparsificationStrategy, _> = "wfmash:0.0".parse();
        assert!(result.is_err());

        let result: Result<SparsificationStrategy, _> = "wfmash:1.5".parse();
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_invalid_sparsification() {
        let result: Result<SparsificationStrategy, _> = "invalid".parse();
        assert!(result.is_err());

        let result: Result<SparsificationStrategy, _> = "random:1.5".parse(); // > 1.0
        assert!(result.is_err());

        let result: Result<SparsificationStrategy, _> = "random:-0.1".parse(); // < 0
        assert!(result.is_err());
    }

    #[test]
    fn test_description() {
        assert_eq!(
            SparsificationStrategy::None.description(),
            "all pairs (no sparsification)"
        );
        assert_eq!(
            SparsificationStrategy::Auto.description(),
            "auto (context-dependent)"
        );
        assert!(SparsificationStrategy::WfmashDensity(None)
            .description()
            .contains("auto"));
        assert!(SparsificationStrategy::WfmashDensity(Some(0.3))
            .description()
            .contains("30.0%"));
    }

    #[test]
    fn test_wfmash_auto_density() {
        // n=1 -> None
        assert_eq!(SparsificationStrategy::wfmash_auto_density(1), None);
        // n=2 -> ln(2)/2*10 ~ 3.47 >= 1.0 -> None
        assert_eq!(SparsificationStrategy::wfmash_auto_density(2), None);
        // n=100 -> ln(100)/100*10 ~ 0.46 -> Some(0.46)
        let frac = SparsificationStrategy::wfmash_auto_density(100).unwrap();
        assert!(frac > 0.4 && frac < 0.5);
    }

    #[test]
    fn test_display_roundtrip() {
        let strategies = vec![
            SparsificationStrategy::None,
            SparsificationStrategy::Auto,
            SparsificationStrategy::Random(0.5),
            SparsificationStrategy::Connectivity(0.99),
            SparsificationStrategy::TreeSampling(5, 2, 0.1),
            SparsificationStrategy::WfmashDensity(None),
            SparsificationStrategy::WfmashDensity(Some(0.3)),
        ];
        for strategy in strategies {
            let s = strategy.to_string();
            let parsed: SparsificationStrategy = s.parse().unwrap();
            assert_eq!(strategy, parsed, "Roundtrip failed for {}", s);
        }
    }

    #[test]
    fn test_estimate_pair_count_larger() {
        // With 10 samples, 10 choose 2 = 45 possible pairs
        assert_eq!(estimate_tree_pair_count(10, 0, 0, 0.0), 0);
        assert_eq!(estimate_tree_pair_count(10, 1, 0, 0.0), 10);
        assert_eq!(estimate_tree_pair_count(10, 2, 0, 0.0), 20);

        // With random fraction
        assert_eq!(estimate_tree_pair_count(10, 0, 0, 0.5), 23); // 45 * 0.5 rounded

        // Clamped to max pairs
        assert_eq!(estimate_tree_pair_count(10, 10, 10, 1.0), 45); // Can't exceed n choose 2
    }

    #[test]
    fn test_select_pairs_with_sequences() {
        // Test that tree sampling works with actual sequences
        let sequences: Vec<Vec<u8>> = vec![
            b"ATCGATCGATCGATCGATCGATCG".to_vec(),
            b"ATCGATCGATCGATCGATCGATCG".to_vec(), // Same as 0
            b"GCTAGCTAGCTAGCTAGCTAGCTA".to_vec(), // Different
            b"GCTAGCTAGCTAGCTAGCTAGCTA".to_vec(), // Same as 2
        ];

        let pairs = select_pairs(
            4,
            Some(&sequences),
            &SparsificationStrategy::TreeSampling(1, 0, 0.0),
            &MashParams::default(),
        );

        // Should connect similar sequences
        assert!(!pairs.is_empty());

        // 0 and 1 are identical, should be connected
        // 2 and 3 are identical, should be connected
        let has_01 = pairs.contains(&(0, 1)) || pairs.contains(&(1, 0));
        let has_23 = pairs.contains(&(2, 3)) || pairs.contains(&(3, 2));
        assert!(has_01 || has_23); // At least one identical pair should be connected
    }
}
