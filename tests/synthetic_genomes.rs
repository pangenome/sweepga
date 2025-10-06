use rand::rngs::StdRng;
use rand::seq::{IteratorRandom, SliceRandom};
use rand::{Rng, SeedableRng};
/// Generate synthetic genomes with controlled mutations for testing
/// This ensures FastGA can find homologous regions to align
use std::fs;
use std::path::{Path, PathBuf};

/// Generate a stable random DNA sequence with a fixed seed
pub fn generate_base_sequence(length: usize, seed: u64) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    let bases = ['A', 'C', 'G', 'T'];

    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}

/// Mutation types
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum MutationType {
    Substitution,
    Insertion,
    Deletion,
}

/// Apply controlled mutations to create a derived sequence
#[allow(dead_code)]
pub fn mutate_sequence(base: &str, num_mutations: usize, seed: u64) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut sequence: Vec<char> = base.chars().collect();
    let bases = ['A', 'C', 'G', 'T'];

    // Generate random positions for mutations (avoiding duplicates)
    let mut positions: Vec<usize> = (0..sequence.len()).collect();
    positions.shuffle(&mut rng);
    positions.truncate(num_mutations);
    positions.sort_unstable();

    // Apply mutations from end to beginning to maintain positions
    for &pos in positions.iter().rev() {
        if pos >= sequence.len() {
            continue; // Skip if out of bounds
        }

        let mutation_type = match rng.gen_range(0..3) {
            0 => MutationType::Substitution,
            1 => MutationType::Insertion,
            _ => MutationType::Deletion,
        };

        match mutation_type {
            MutationType::Substitution => {
                // Change to a different base
                let current = sequence[pos];
                let new_base = bases
                    .iter()
                    .filter(|&&b| b != current)
                    .choose(&mut rng)
                    .copied()
                    .unwrap_or('N');
                sequence[pos] = new_base;
            }
            MutationType::Insertion => {
                // Insert 1-5 random bases
                let insert_len = rng.gen_range(1..=5);
                for _ in 0..insert_len {
                    let new_base = bases[rng.gen_range(0..4)];
                    sequence.insert(pos, new_base);
                }
            }
            MutationType::Deletion => {
                // Delete 1-5 bases (but don't exceed bounds)
                let delete_len = rng.gen_range(1..=5).min(sequence.len() - pos);
                for _ in 0..delete_len {
                    if pos < sequence.len() {
                        sequence.remove(pos);
                    }
                }
            }
        }
    }

    sequence.iter().collect()
}

/// Generate a pair of related sequences for testing alignment
#[allow(dead_code)]
pub fn generate_test_pair(length: usize, divergence: f64) -> (String, String) {
    // Use fixed seeds for reproducibility
    let base = generate_base_sequence(length, 42);

    // Calculate number of mutations based on divergence
    let num_mutations = (length as f64 * divergence) as usize;

    // Create mutated version
    let mutated = mutate_sequence(&base, num_mutations, 43);

    (base, mutated)
}

/// Create test FASTA files with various sequence lengths
#[allow(dead_code)]
pub fn create_scaled_test_files(dir: &Path) -> Vec<(PathBuf, PathBuf, usize)> {
    let sizes = vec![1000, 10000, 50000, 100000];
    let mut results = Vec::new();

    for size in sizes {
        let file1 = dir.join(format!("test_{size}bp_1.fa"));
        let file2 = dir.join(format!("test_{size}bp_2.fa"));

        // Generate sequences with 5% divergence
        let (seq1, seq2) = generate_test_pair(size, 0.05);

        // Write FASTA files
        fs::write(&file1, format!(">seq_{size}_1\n{seq1}\n")).unwrap();
        fs::write(&file2, format!(">seq_{size}_2\n{seq2}\n")).unwrap();

        results.push((file1, file2, size));
    }

    results
}

/// Create multi-chromosome test genome
#[allow(dead_code)]
pub fn create_multichrom_genome(path: &Path, size_per_chr: usize, num_chromosomes: usize) {
    let mut content = String::new();

    for i in 1..=num_chromosomes {
        // Each chromosome has a different seed for variety
        let seq = generate_base_sequence(size_per_chr, (i * 100) as u64);
        content.push_str(&format!(">chr{i}\n{seq}\n"));
    }

    fs::write(path, content).unwrap();
}

/// Create a genome with repeats for testing repeat detection
#[allow(dead_code)]
pub fn create_genome_with_repeats(path: &Path, base_size: usize) {
    let mut sequence = generate_base_sequence(base_size, 99);

    // Add tandem repeats
    let repeat_unit = "ACGTACGTACGT";
    for _ in 0..10 {
        sequence.push_str(repeat_unit);
    }

    // Add interspersed repeats
    let transposon = "AAAATTTTGGGGCCCC";
    let positions = [base_size / 4, base_size / 2, 3 * base_size / 4];

    let mut chars: Vec<char> = sequence.chars().collect();
    for pos in positions.iter().rev() {
        for ch in transposon.chars() {
            chars.insert(*pos, ch);
        }
    }

    let final_seq: String = chars.iter().collect();
    fs::write(path, format!(">genome_with_repeats\n{final_seq}\n")).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    #[allow(dead_code)]
    fn test_stable_sequence_generation() {
        // Same seed should produce same sequence
        let seq1 = generate_base_sequence(100, 42);
        let seq2 = generate_base_sequence(100, 42);
        assert_eq!(seq1, seq2, "Same seed should produce same sequence");

        // Different seeds should produce different sequences
        let seq3 = generate_base_sequence(100, 43);
        assert_ne!(
            seq1, seq3,
            "Different seeds should produce different sequences"
        );
    }

    #[test]
    #[allow(dead_code)]
    fn test_mutations() {
        let base = "ACGTACGTACGTACGT";
        let mutated = mutate_sequence(base, 3, 42);

        // Should have some differences
        assert_ne!(base, mutated, "Mutated sequence should differ from base");

        // Check that length doesn't change too drastically
        // With 3 mutations on a 16-char string, could have up to 3 insertions
        // or 3 deletions, so allow +/- 3 characters
        let len_diff = (mutated.len() as i32 - base.len() as i32).abs();
        assert!(
            len_diff <= 10, // Allow reasonable change for indels
            "Length change unexpected: {} -> {}",
            base.len(),
            mutated.len()
        );
    }

    #[test]
    #[allow(dead_code)]
    fn test_scaled_generation() {
        let temp_dir = TempDir::new().unwrap();
        let files = create_scaled_test_files(temp_dir.path());

        assert_eq!(files.len(), 4, "Should create 4 size variants");

        for (file1, file2, size) in files {
            assert!(file1.exists(), "File 1 should exist");
            assert!(file2.exists(), "File 2 should exist");

            let content1 = fs::read_to_string(&file1).unwrap();
            let content2 = fs::read_to_string(&file2).unwrap();

            // Check approximate size (accounting for header, newlines, and mutations)
            // Content should be at least 90% of expected size (mutations can delete)
            assert!(
                content1.len() > size * 9 / 10,
                "File 1 too small: {} bytes for {} bp sequence",
                content1.len(),
                size
            );
            assert!(
                content2.len() > size * 9 / 10,
                "File 2 too small: {} bytes for {} bp sequence",
                content2.len(),
                size
            );
        }
    }
}
