//! Mash distance computation using MinHash k-mer sketching
//!
//! This module provides functionality to compute pairwise mash distances between sequences
//! using k-mer sketching with MinHash. Ported from allwave.

#![allow(dead_code)]

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

/// Default k-mer size for mash distance computation
pub const DEFAULT_KMER_SIZE: usize = 15;

/// Default sketch size (number of minimizers to keep)
pub const DEFAULT_SKETCH_SIZE: usize = 1000;

/// A k-mer sketch using MinHash
#[derive(Debug, Clone)]
pub struct KmerSketch {
    /// Minimizers (smallest hash values) for this sequence
    pub minimizers: Vec<u64>,
    /// K-mer size used
    pub k: usize,
    /// Original sequence length
    pub length: usize,
}

impl KmerSketch {
    /// Create a new sketch from a sequence
    pub fn from_sequence(sequence: &[u8], k: usize, sketch_size: usize) -> Self {
        let minimizers = sketch_sequence(sequence, k, sketch_size);
        Self {
            minimizers,
            k,
            length: sequence.len(),
        }
    }

    /// Compute Jaccard index with another sketch
    pub fn jaccard(&self, other: &KmerSketch) -> f64 {
        if self.k != other.k {
            return 0.0;
        }

        let set1: HashSet<_> = self.minimizers.iter().collect();
        let set2: HashSet<_> = other.minimizers.iter().collect();

        let intersection_size = set1.intersection(&set2).count();
        let union_size = set1.union(&set2).count();

        if union_size == 0 {
            0.0
        } else {
            intersection_size as f64 / union_size as f64
        }
    }

    /// Compute mash distance from Jaccard index
    pub fn mash_distance(&self, other: &KmerSketch) -> f64 {
        let jaccard = self.jaccard(other);
        if jaccard <= 0.0 {
            1.0 // Maximum distance
        } else {
            // Mash distance: -1/k * ln(2*J/(1+J))
            // where J is Jaccard index
            let k = self.k as f64;
            let ratio = (2.0 * jaccard) / (1.0 + jaccard);
            if ratio <= 0.0 {
                1.0
            } else {
                (-1.0 / k) * ratio.ln()
            }
        }
    }
}

/// Create MinHash sketch from sequence
fn sketch_sequence(sequence: &[u8], k: usize, sketch_size: usize) -> Vec<u64> {
    if sequence.len() < k {
        return Vec::new();
    }

    let mut hashes = Vec::new();

    // Extract all k-mers and hash them
    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];

        // Skip k-mers containing non-ACGT characters
        if kmer.iter().any(|&b| !is_dna_base(b)) {
            continue;
        }

        // Hash the k-mer (considering both orientations)
        let hash_fwd = hash_kmer(kmer);
        let hash_rev = hash_kmer(&reverse_complement_kmer(kmer));

        // Use canonical k-mer (lexicographically smaller)
        let canonical_hash = hash_fwd.min(hash_rev);
        hashes.push(canonical_hash);
    }

    // Sort and take the smallest hashes (MinHash)
    hashes.sort_unstable();
    hashes.truncate(sketch_size);
    hashes
}

/// Hash a k-mer using a simple polynomial rolling hash
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Check if byte represents a DNA base
fn is_dna_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Compute reverse complement of a k-mer
fn reverse_complement_kmer(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b, // Keep non-DNA characters as-is
        })
        .collect()
}

/// Compute all-vs-all mash distances between sequences
pub fn compute_distance_matrix(sequences: &[Vec<u8>]) -> Vec<Vec<f64>> {
    compute_distance_matrix_with_params(sequences, DEFAULT_KMER_SIZE, DEFAULT_SKETCH_SIZE)
}

/// Compute all-vs-all mash distances with custom parameters
pub fn compute_distance_matrix_with_params(
    sequences: &[Vec<u8>],
    k: usize,
    sketch_size: usize,
) -> Vec<Vec<f64>> {
    let n = sequences.len();
    let mut matrix = vec![vec![0.0; n]; n];

    // Create sketches for all sequences
    let sketches: Vec<_> = sequences
        .iter()
        .map(|seq| KmerSketch::from_sequence(seq, k, sketch_size))
        .collect();

    // Compute pairwise distances
    for i in 0..n {
        for j in i + 1..n {
            let distance = sketches[i].mash_distance(&sketches[j]);
            matrix[i][j] = distance;
            matrix[j][i] = distance; // Symmetric
        }
    }

    matrix
}

/// Compute sketches for all sequences in parallel
pub fn compute_sketches_parallel(
    sequences: &[Vec<u8>],
    k: usize,
    sketch_size: usize,
) -> Vec<KmerSketch> {
    use rayon::prelude::*;

    sequences
        .par_iter()
        .map(|seq| KmerSketch::from_sequence(seq, k, sketch_size))
        .collect()
}

/// Compute distance matrix from pre-computed sketches
pub fn distance_matrix_from_sketches(sketches: &[KmerSketch]) -> Vec<Vec<f64>> {
    let n = sketches.len();
    let mut matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in i + 1..n {
            let distance = sketches[i].mash_distance(&sketches[j]);
            matrix[i][j] = distance;
            matrix[j][i] = distance;
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_sketch() {
        let seq = b"ATCGATCGATCG";
        let sketch = KmerSketch::from_sequence(seq, 4, 10);
        assert!(!sketch.minimizers.is_empty());
        assert_eq!(sketch.k, 4);
        assert_eq!(sketch.length, seq.len());
    }

    #[test]
    fn test_jaccard_identical() {
        let seq = b"ATCGATCGATCG";
        let sketch1 = KmerSketch::from_sequence(seq, 4, 10);
        let sketch2 = KmerSketch::from_sequence(seq, 4, 10);

        let jaccard = sketch1.jaccard(&sketch2);
        assert!((jaccard - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mash_distance_identical() {
        let seq = b"ATCGATCGATCG";
        let sketch1 = KmerSketch::from_sequence(seq, 4, 10);
        let sketch2 = KmerSketch::from_sequence(seq, 4, 10);

        let distance = sketch1.mash_distance(&sketch2);
        assert!(distance < 1e-10); // Should be very close to 0
    }

    #[test]
    fn test_reverse_complement() {
        let kmer = b"ATCG";
        let rc = reverse_complement_kmer(kmer);
        assert_eq!(rc, b"CGAT");
    }

    #[test]
    fn test_distance_matrix() {
        let sequences = vec![
            b"ATCGATCGATCGATCG".to_vec(),
            b"ATCGATCGATCGATCG".to_vec(),
            b"GGGGGGGGGGGGGGGG".to_vec(),
        ];

        let matrix = compute_distance_matrix(&sequences);
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);

        // Diagonal should be 0
        assert!(matrix[0][0] < 1e-6);
        assert!(matrix[1][1] < 1e-6);
        assert!(matrix[2][2] < 1e-6);

        // Identical sequences should have distance ~0
        assert!(matrix[0][1] < 1e-6);
        assert!(matrix[1][0] < 1e-6);

        // Different sequences should have non-zero distance
        assert!(matrix[0][2] > 0.0);
        assert!(matrix[2][0] > 0.0);
    }

    #[test]
    fn test_short_sequence() {
        // Sequence shorter than k-mer size should produce empty sketch
        let seq = b"AT";
        let sketch = KmerSketch::from_sequence(seq, 4, 10);
        assert!(sketch.minimizers.is_empty());
    }

    #[test]
    fn test_similar_sequences() {
        // Two similar sequences should have small distance
        let seq1 = b"ATCGATCGATCGATCGATCGATCG";
        let seq2 = b"ATCGATCGATCGATCGATCGATCG"; // Same
        let seq3 = b"ATCGATCGATCGTTCGATCGATCG"; // One mutation

        let sketch1 = KmerSketch::from_sequence(seq1, 6, 100);
        let sketch2 = KmerSketch::from_sequence(seq2, 6, 100);
        let sketch3 = KmerSketch::from_sequence(seq3, 6, 100);

        // Identical should have distance 0
        assert!(sketch1.mash_distance(&sketch2) < 1e-10);

        // One mutation should have small but non-zero distance
        let dist = sketch1.mash_distance(&sketch3);
        assert!(dist > 0.0);
        assert!(dist < 0.5); // Should still be fairly similar
    }

    #[test]
    fn test_is_dna_base() {
        assert!(is_dna_base(b'A'));
        assert!(is_dna_base(b'C'));
        assert!(is_dna_base(b'G'));
        assert!(is_dna_base(b'T'));
        assert!(is_dna_base(b'a'));
        assert!(is_dna_base(b't'));
        assert!(!is_dna_base(b'N'));
        assert!(!is_dna_base(b'X'));
    }

    #[test]
    fn test_sketch_with_ns() {
        // Sequences with N characters should skip those k-mers
        let seq = b"ATCGNNNNATCG";
        let sketch = KmerSketch::from_sequence(seq, 4, 10);
        // Should still produce some minimizers from valid k-mers
        assert!(!sketch.minimizers.is_empty());
    }

    #[test]
    fn test_jaccard_different_k() {
        let seq1 = b"ATCGATCGATCG";
        let seq2 = b"ATCGATCGATCG";
        let sketch1 = KmerSketch::from_sequence(seq1, 4, 10);
        let sketch2 = KmerSketch::from_sequence(seq2, 5, 10);

        // Different k-mer sizes should return 0 Jaccard
        assert_eq!(sketch1.jaccard(&sketch2), 0.0);
    }

    #[test]
    fn test_distance_matrix_from_sketches() {
        let sequences = vec![
            b"ATCGATCGATCGATCG".to_vec(),
            b"ATCGATCGATCGATCG".to_vec(),
            b"GCTAGCTAGCTAGCTA".to_vec(),
        ];

        let sketches: Vec<_> = sequences
            .iter()
            .map(|s| KmerSketch::from_sequence(s, 4, 100))
            .collect();

        let matrix = distance_matrix_from_sketches(&sketches);

        // Same tests as compute_distance_matrix
        assert_eq!(matrix.len(), 3);
        assert!(matrix[0][1] < 1e-6); // Identical
        assert!(matrix[0][2] > 0.0); // Different
    }
}
