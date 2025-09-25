/// Comprehensive integration tests for FastGA-RS
/// These tests verify that FastGA correctly generates alignments from FASTA inputs

use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Helper to count lines in a file
fn count_lines(path: &Path) -> usize {
    fs::read_to_string(path)
        .expect("Failed to read file")
        .lines()
        .count()
}

/// Helper to check if PAF output contains extended CIGAR
fn has_extended_cigar(path: &Path) -> bool {
    fs::read_to_string(path)
        .expect("Failed to read PAF")
        .lines()
        .any(|line| line.contains("cg:Z:") && line.contains("="))
}

/// Helper to run sweepga command
fn run_sweepga(args: &[&str]) -> Result<String, String> {
    let output = Command::new("cargo")
        .args(&["run", "--bin", "sweepga", "--release", "--quiet", "--"])
        .args(args)
        .output()
        .expect("Failed to run sweepga");

    if output.status.success() {
        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    } else {
        Err(String::from_utf8_lossy(&output.stderr).to_string())
    }
}

#[test]
fn test_fastga_self_alignment() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("self_align.paf");

    // Run self-alignment on yeast data (minimal filtering to keep more alignments)
    let result = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-t", "2",
        "--self",  // Include self-mappings
        "-n", "N",  // No mapping filter
        "-j", "0",  // No scaffolding
        "-o", output_path.to_str().unwrap(),
    ]);

    assert!(result.is_ok(), "FastGA self-alignment failed: {:?}", result);
    assert!(output_path.exists(), "Output PAF not created");

    let line_count = count_lines(&output_path);
    assert!(line_count > 0, "No alignments produced");
    assert!(line_count > 1000, "Too few alignments for yeast self-alignment");

    assert!(has_extended_cigar(&output_path), "Missing extended CIGAR format");
}

#[test]
fn test_fastga_pairwise_alignment() {
    // Import from parent module (tests/)
    #[path = "synthetic_genomes.rs"]
    mod synthetic_genomes;
    use synthetic_genomes::generate_test_pair;

    let temp_dir = TempDir::new().unwrap();

    // Create two test FASTA files with related sequences
    let fasta1 = temp_dir.path().join("seq1.fa");
    let fasta2 = temp_dir.path().join("seq2.fa");
    let output = temp_dir.path().join("pairwise.paf");

    // Generate related sequences (10kb with 5% divergence)
    let (seq1, seq2) = generate_test_pair(10000, 0.05);
    fs::write(&fasta1, format!(">seq1\n{}\n", seq1)).unwrap();
    fs::write(&fasta2, format!(">seq2\n{}\n", seq2)).unwrap();

    // Run pairwise alignment
    let result = run_sweepga(&[
        fasta1.to_str().unwrap(),
        fasta2.to_str().unwrap(),
        "-t", "1",
        "-o", output.to_str().unwrap(),
    ]);

    assert!(result.is_ok(), "Pairwise alignment failed");
    assert!(output.exists(), "Output not created");

    // Should produce at least one alignment between similar sequences
    let content = fs::read_to_string(&output).unwrap();
    assert!(!content.is_empty(), "No alignments produced");
    assert!(content.contains("seq1"), "Missing query sequence name");
    assert!(content.contains("seq2"), "Missing target sequence name");
}

#[test]
fn test_thread_parameter() {
    let temp_dir = TempDir::new().unwrap();
    let output1 = temp_dir.path().join("t1.paf");
    let output4 = temp_dir.path().join("t4.paf");

    // Create small test file
    let test_fa = temp_dir.path().join("test.fa");
    fs::write(&test_fa, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n").unwrap();

    // Run with 1 thread
    let result1 = run_sweepga(&[
        test_fa.to_str().unwrap(),
        "-t", "1",
        "--self",  // Include self-mappings for consistency
        "-o", output1.to_str().unwrap(),
    ]);

    // Run with 4 threads
    let result4 = run_sweepga(&[
        test_fa.to_str().unwrap(),
        "-t", "4",
        "--self",  // Include self-mappings for consistency
        "-o", output4.to_str().unwrap(),
    ]);

    assert!(result1.is_ok() && result4.is_ok(), "Thread tests failed");

    // Both should produce output
    assert!(output1.exists() && output4.exists(), "Outputs not created");

    // Results should be deterministic (same number of alignments)
    let lines1 = count_lines(&output1);
    let lines4 = count_lines(&output4);
    assert_eq!(lines1, lines4, "Thread count affected output ({} vs {})", lines1, lines4);
}

#[test]
fn test_filtering_with_fastga() {
    let temp_dir = TempDir::new().unwrap();
    let loose = temp_dir.path().join("loose.paf");
    let strict = temp_dir.path().join("strict.paf");

    // Test with different block length filters
    let result_loose = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-b", "0",  // No minimum block length
        "-t", "2",
        "--self",
        "-o", loose.to_str().unwrap(),
    ]);

    let result_strict = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-b", "10k",  // 10kb minimum
        "-t", "2",
        "--self",
        "-o", strict.to_str().unwrap(),
    ]);

    assert!(result_loose.is_ok() && result_strict.is_ok(), "Filtering tests failed");

    let loose_count = count_lines(&loose);
    let strict_count = count_lines(&strict);

    assert!(strict_count < loose_count,
            "Strict filter should produce fewer alignments ({} >= {})",
            strict_count, loose_count);
    assert!(strict_count > 0, "Strict filter removed all alignments");
}

#[test]
fn test_scaffold_filtering() {
    let temp_dir = TempDir::new().unwrap();
    let with_scaffold = temp_dir.path().join("scaffold.paf");
    let no_scaffold = temp_dir.path().join("no_scaffold.paf");

    // With scaffolding
    let result1 = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-j", "100k",  // Enable scaffolding with 100kb jump
        "-s", "50k",   // Minimum scaffold mass
        "-t", "2",
        "--self",
        "-o", with_scaffold.to_str().unwrap(),
    ]);

    // Without scaffolding
    let result2 = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-j", "0",  // Disable scaffolding
        "-t", "2",
        "--self",
        "-o", no_scaffold.to_str().unwrap(),
    ]);

    assert!(result1.is_ok() && result2.is_ok(), "Scaffold tests failed");

    // Check for scaffold annotations
    let scaffold_content = fs::read_to_string(&with_scaffold).unwrap();
    assert!(scaffold_content.contains("st:Z:scaffold"),
            "Missing scaffold annotations");
    assert!(scaffold_content.contains("ch:Z:chain"),
            "Missing chain annotations");
}

#[test]
fn test_empty_input_handling() {
    let temp_dir = TempDir::new().unwrap();
    let empty_fa = temp_dir.path().join("empty.fa");
    let output = temp_dir.path().join("empty_out.paf");

    // Create empty FASTA
    fs::write(&empty_fa, "").unwrap();

    // This should handle gracefully
    let result = run_sweepga(&[
        empty_fa.to_str().unwrap(),
        "-t", "1",
        "-o", output.to_str().unwrap(),
    ]);

    // Should either succeed with no output or fail gracefully
    // Not asserting success as behavior may vary, but shouldn't panic
    if result.is_ok() && output.exists() {
        let content = fs::read_to_string(&output).unwrap();
        assert!(content.is_empty() || content.lines().count() == 0,
                "Empty input produced alignments");
    }
}

#[test]
fn test_large_sequence_handling() {
    #[path = "synthetic_genomes.rs"]
    mod synthetic_genomes;
    use synthetic_genomes::generate_base_sequence;

    let temp_dir = TempDir::new().unwrap();
    let large_fa = temp_dir.path().join("large.fa");
    let output = temp_dir.path().join("large_out.paf");

    // Create a large sequence with internal repeats (100kb)
    let base = generate_base_sequence(50000, 999);
    // Duplicate part of it to create self-similarity
    let seq = format!("{}{}", &base[..25000], &base);
    fs::write(&large_fa, format!(">large_seq\n{}\n", seq)).unwrap();

    // Should handle without hanging or crashing
    let result = run_sweepga(&[
        large_fa.to_str().unwrap(),
        "-t", "2",
        "--self",  // Include self-mappings
        "-o", output.to_str().unwrap(),
    ]);

    assert!(result.is_ok(), "Failed on large sequence");
    if output.exists() {
        // Should find self-alignments in the duplicated regions
        let content = fs::read_to_string(&output).unwrap();
        assert!(!content.is_empty(), "Should find alignments in duplicated regions");
        assert!(has_extended_cigar(&output), "Invalid output for large sequence");
    }
}

#[test]
fn test_multisequence_fasta() {
    #[path = "synthetic_genomes.rs"]
    mod synthetic_genomes;
    use synthetic_genomes::{generate_base_sequence, mutate_sequence};

    let temp_dir = TempDir::new().unwrap();
    let multi_fa = temp_dir.path().join("multi.fa");
    let output = temp_dir.path().join("multi_out.paf");

    // Create multi-sequence FASTA with related sequences
    let base = generate_base_sequence(5000, 123);
    let seq1 = base.clone();
    let seq2 = mutate_sequence(&base, 50, 124);  // 1% divergence
    let seq3 = mutate_sequence(&base, 100, 125); // 2% divergence

    fs::write(&multi_fa, format!(
        ">seq1\n{}\n>seq2\n{}\n>seq3\n{}\n",
        seq1, seq2, seq3
    )).unwrap();

    let result = run_sweepga(&[
        multi_fa.to_str().unwrap(),
        "-t", "1",
        "--self",  // Include self-mappings
        "-o", output.to_str().unwrap(),
    ]);

    assert!(result.is_ok(), "Failed on multi-sequence FASTA");

    if output.exists() {
        let content = fs::read_to_string(&output).unwrap();
        // Should have alignments between the related sequences
        assert!(!content.is_empty(), "Should produce alignments");

        // Check that we found alignments between sequences
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines.len() > 0, "Should have at least one alignment");
    }
}

#[test]
#[ignore]  // This test requires significant time
fn test_performance_regression() {
    use std::time::Instant;

    let temp_dir = TempDir::new().unwrap();
    let output = temp_dir.path().join("perf.paf");

    let start = Instant::now();

    let result = run_sweepga(&[
        "data/scerevisiae8.fa",
        "-t", "4",
        "--self",
        "-o", output.to_str().unwrap(),
    ]);

    let duration = start.elapsed();

    assert!(result.is_ok(), "Performance test failed");

    // Yeast self-alignment should complete in reasonable time
    assert!(duration.as_secs() < 60,
            "Alignment took too long: {:?}", duration);

    // Should produce expected number of alignments
    let line_count = count_lines(&output);
    assert!(line_count > 1000 && line_count < 5000,
            "Unexpected alignment count: {}", line_count);
}