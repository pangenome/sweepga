/// Utility functions for testing
use std::fs;
use std::path::{Path, PathBuf};

/// Generate a random DNA sequence of given length
pub fn generate_dna_sequence(length: usize) -> String {
    use rand::{thread_rng, Rng};
    let bases = ['A', 'C', 'G', 'T'];
    let mut rng = thread_rng();

    (0..length)
        .map(|_| bases[rng.gen_range(0..4)])
        .collect()
}

/// Generate a FASTA file with specified sequences
pub fn create_fasta_file(path: &Path, sequences: &[(&str, &str)]) {
    let mut content = String::new();
    for (name, seq) in sequences {
        content.push_str(&format!(">{name}\n{seq}\n"));
    }
    fs::write(path, content).expect("Failed to write FASTA file");
}

/// Generate test genome with repeats and variations
pub fn generate_test_genome(path: &Path, size: usize, num_chromosomes: usize) {
    let mut sequences = Vec::new();

    for i in 1..=num_chromosomes {
        let chr_name = format!("chr{i}");
        let mut seq = generate_dna_sequence(size / num_chromosomes);

        // Add some repeats within the sequence
        if i == 1 {
            // Add tandem repeat
            let repeat = "ACGTACGTACGTACGT";
            seq.push_str(repeat);
            seq.push_str(repeat);
        }

        // Add some variation in chromosome 2
        if i == 2 {
            // Introduce some mutations
            let mut chars: Vec<char> = seq.chars().collect();
            for j in (0..chars.len()).step_by(100) {
                if j < chars.len() {
                    chars[j] = if chars[j] == 'A' { 'G' } else { 'A' };
                }
            }
            seq = chars.iter().collect();
        }

        sequences.push((chr_name, seq));
    }

    let content: String = sequences
        .iter()
        .map(|(name, seq)| format!(">{name}\n{seq}\n"))
        .collect();

    fs::write(path, content).expect("Failed to write test genome");
}

/// Create genomes with known homology for testing alignments
pub fn create_homologous_genomes(dir: &Path) -> (PathBuf, PathBuf) {
    let genome1 = dir.join("genome1.fa");
    let genome2 = dir.join("genome2.fa");

    // Create base sequence
    let base_seq = generate_dna_sequence(10000);

    // Genome 1: original
    create_fasta_file(&genome1, &[
        ("chr1", &base_seq),
        ("chr2", &generate_dna_sequence(8000)),
    ]);

    // Genome 2: with variations
    let variant_seq = base_seq.clone();

    // Introduce SNPs every 100 bases
    let mut chars: Vec<char> = variant_seq.chars().collect();
    for i in (0..chars.len()).step_by(100) {
        if i < chars.len() {
            chars[i] = match chars[i] {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => 'N',
            };
        }
    }

    // Add insertion at position 5000
    let insertion = "AAAAAAAAAA";
    chars.splice(5000..5000, insertion.chars());

    // Add deletion at position 7000 (remove 10 bases)
    if chars.len() > 7010 {
        chars.drain(7000..7010);
    }

    let variant_seq: String = chars.iter().collect();

    // Add an inverted segment
    let inverted = generate_dna_sequence(1000)
        .chars()
        .rev()
        .collect::<String>();

    create_fasta_file(&genome2, &[
        ("chr1_variant", &variant_seq),
        ("chr2_inverted", &inverted),
        ("chr3_novel", &generate_dna_sequence(5000)),
    ]);

    (genome1, genome2)
}

/// Parse PAF file and return basic statistics
pub struct PafStats {
    pub num_alignments: usize,
    pub num_queries: usize,
    pub num_targets: usize,
    pub has_extended_cigar: bool,
    pub has_scaffolds: bool,
    pub total_aligned_bases: u64,
    pub average_identity: f64,
}

pub fn analyze_paf(path: &Path) -> PafStats {
    let content = fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<&str> = content.lines().collect();

    let mut queries = std::collections::HashSet::new();
    let mut targets = std::collections::HashSet::new();
    let mut total_bases = 0u64;
    let mut total_identity = 0.0;
    let mut has_cigar = false;
    let mut has_scaffolds = false;

    for line in &lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 12 {
            queries.insert(fields[0]);
            targets.insert(fields[5]);

            if let Ok(bases) = fields[10].parse::<u64>() {
                total_bases += bases;
            }

            if let (Ok(matches), Ok(aligned)) = (fields[9].parse::<f64>(), fields[10].parse::<f64>()) {
                if aligned > 0.0 {
                    total_identity += matches / aligned;
                }
            }

            if line.contains("cg:Z:") && line.contains("=") {
                has_cigar = true;
            }

            if line.contains("st:Z:scaffold") {
                has_scaffolds = true;
            }
        }
    }

    PafStats {
        num_alignments: lines.len(),
        num_queries: queries.len(),
        num_targets: targets.len(),
        has_extended_cigar: has_cigar,
        has_scaffolds,
        total_aligned_bases: total_bases,
        average_identity: if lines.is_empty() { 0.0 } else { total_identity / lines.len() as f64 },
    }
}

/// Create test configuration for sweepga
pub struct TestConfig {
    pub threads: usize,
    pub block_length: Option<String>,
    pub scaffold_jump: Option<String>,
    pub scaffold_mass: Option<String>,
    pub overlap: Option<f64>,
}

impl Default for TestConfig {
    fn default() -> Self {
        TestConfig {
            threads: 2,
            block_length: None,
            scaffold_jump: None,
            scaffold_mass: None,
            overlap: None,
        }
    }
}

impl TestConfig {
    pub fn to_args(&self) -> Vec<String> {
        let mut args = vec!["-t".to_string(), self.threads.to_string()];

        if let Some(ref bl) = self.block_length {
            args.push("-b".to_string());
            args.push(bl.clone());
        }

        if let Some(ref sj) = self.scaffold_jump {
            args.push("-j".to_string());
            args.push(sj.clone());
        }

        if let Some(ref sm) = self.scaffold_mass {
            args.push("-s".to_string());
            args.push(sm.clone());
        }

        if let Some(overlap) = self.overlap {
            args.push("-p".to_string());
            args.push(overlap.to_string());
        }

        args
    }
}