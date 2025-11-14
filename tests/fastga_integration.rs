#![allow(clippy::uninlined_format_args)]
/// Comprehensive integration tests for FastGA-RS
/// These tests verify that FastGA correctly generates alignments from FASTA inputs
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

#[path = "synthetic_genomes.rs"]
mod synthetic_genomes;

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
    // Use the compiled binary directly instead of cargo run
    // This assumes tests are run with `cargo test --release`
    let sweepga_path = if cfg!(debug_assertions) {
        "target/debug/sweepga"
    } else {
        "target/release/sweepga"
    };

    let output = Command::new(sweepga_path)
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
// FastGA binaries are now available via binary_paths
fn test_fastga_self_alignment() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("self_align.paf");

    // Run self-alignment on yeast data (minimal filtering to keep more alignments)
    let result = run_sweepga(&[
        "data/scerevisiae8.fa.gz",
        "-t",
        "2",
        "--self", // Include self-mappings
        "-n",
        "N", // No mapping filter
        "-j",
        "0",     // No scaffolding
        "--paf", // Output PAF format (with extended CIGAR)
    ]);

    assert!(result.is_ok(), "FastGA self-alignment failed: {result:?}");

    // Write stdout to output file
    fs::write(&output_path, result.unwrap()).unwrap();
    assert!(output_path.exists(), "Output PAF not created");

    let line_count = count_lines(&output_path);
    assert!(line_count > 0, "No alignments produced");
    // With 8 yeast genomes, we expect at least the self-mappings for each chromosome
    // There are ~17 chromosomes per genome, so 8 * 17 = ~136 minimum
    assert!(
        line_count > 100,
        "Too few alignments for yeast self-alignment (got {line_count})"
    );

    assert!(
        has_extended_cigar(&output_path),
        "Missing extended CIGAR format"
    );
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_fastga_pairwise_alignment() {
    // Import from synthetic_genomes module
    use self::synthetic_genomes::generate_test_pair;

    let temp_dir = TempDir::new().unwrap();

    // Create two test FASTA files with related sequences
    let fasta1 = temp_dir.path().join("seq1.fa");
    let fasta2 = temp_dir.path().join("seq2.fa");
    let output = temp_dir.path().join("pairwise.paf");

    // Generate related sequences (20kb with 1% divergence for better alignment)
    let (seq1, seq2) = generate_test_pair(20000, 0.01);
    fs::write(&fasta1, format!(">seq1\n{seq1}\n")).unwrap();
    fs::write(&fasta2, format!(">seq2\n{seq2}\n")).unwrap();

    // Debug: Check if files were created correctly
    assert!(
        fasta1.exists() && fasta2.exists(),
        "Input files not created"
    );
    let seq1_size = fs::metadata(&fasta1).unwrap().len();
    let seq2_size = fs::metadata(&fasta2).unwrap().len();
    // Sequences should be around 10kb (with header, at least 9.5kb)
    assert!(
        seq1_size > 9500 && seq2_size > 9500,
        "Input files too small: {seq1_size} and {seq2_size}"
    );

    // Run pairwise alignment
    eprintln!(
        "Running alignment: {} vs {}",
        fasta1.display(),
        fasta2.display()
    );
    let result = run_sweepga(&[
        fasta1.to_str().unwrap(),
        fasta2.to_str().unwrap(),
        "-t",
        "1",
        "-i",
        "0",     // Disable identity threshold to avoid filtering test alignment
        "--paf", // Output PAF format
    ]);

    if let Err(ref e) = result {
        eprintln!("Alignment error: {e}");
        // Try to preserve the files for debugging
        let _ = fs::copy(&fasta1, "/tmp/debug_seq1.fa");
        let _ = fs::copy(&fasta2, "/tmp/debug_seq2.fa");
    }

    assert!(result.is_ok(), "Pairwise alignment failed: {result:?}");

    // Write stdout to output file
    let content = result.unwrap();
    fs::write(&output, &content).unwrap();
    assert!(output.exists(), "Output not created");

    // Should produce at least one alignment between similar sequences
    assert!(
        !content.is_empty(),
        "No alignments produced (file has {} bytes)",
        content.len()
    );
    assert!(content.contains("seq1"), "Missing query sequence name");
    assert!(content.contains("seq2"), "Missing target sequence name");
}

#[test]
// FastGA binaries are now available via binary_paths
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
        "-t",
        "1",
        "--self", // Include self-mappings for consistency
        "--paf",  // Request PAF output for test comparison
    ]);

    // Run with 4 threads
    let result4 = run_sweepga(&[
        test_fa.to_str().unwrap(),
        "-t",
        "4",
        "--self", // Include self-mappings for consistency
        "--paf",  // Request PAF output for test comparison
    ]);

    if let Err(ref e) = result1 {
        eprintln!("result1 error: {e}");
    }
    if let Err(ref e) = result4 {
        eprintln!("result4 error: {e}");
    }

    assert!(
        result1.is_ok() && result4.is_ok(),
        "Thread tests failed: result1={result1:?}, result4={result4:?}"
    );

    // Write outputs to files
    fs::write(&output1, result1.unwrap()).unwrap();
    fs::write(&output4, result4.unwrap()).unwrap();

    // Both should produce output
    assert!(output1.exists() && output4.exists(), "Outputs not created");

    // Results should be deterministic (same number of alignments)
    let lines1 = count_lines(&output1);
    let lines4 = count_lines(&output4);
    assert_eq!(
        lines1, lines4,
        "Thread count affected output ({lines1} vs {lines4})"
    );
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_filtering_with_fastga() {
    let temp_dir = TempDir::new().unwrap();
    let loose = temp_dir.path().join("loose.paf");
    let strict = temp_dir.path().join("strict.paf");

    // Test with different block length filters
    let result_loose = run_sweepga(&[
        "data/scerevisiae8.fa.gz",
        "-b",
        "0", // No minimum block length
        "-t",
        "2",
        "--self",
        "--paf",
    ]);

    let result_strict = run_sweepga(&[
        "data/scerevisiae8.fa.gz",
        "-b",
        "10k", // 10kb minimum
        "-t",
        "2",
        "--self",
        "--paf",
    ]);

    assert!(
        result_loose.is_ok() && result_strict.is_ok(),
        "Filtering tests failed"
    );

    // Write outputs to files
    fs::write(&loose, result_loose.unwrap()).unwrap();
    fs::write(&strict, result_strict.unwrap()).unwrap();

    let loose_count = count_lines(&loose);
    let strict_count = count_lines(&strict);

    assert!(
        strict_count < loose_count,
        "Strict filter should produce fewer alignments ({strict_count} >= {loose_count})"
    );
    assert!(strict_count > 0, "Strict filter removed all alignments");
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_scaffold_filtering() {
    let temp_dir = TempDir::new().unwrap();
    let with_scaffold = temp_dir.path().join("scaffold.paf");
    let no_scaffold = temp_dir.path().join("no_scaffold.paf");

    // With scaffolding
    let result1 = run_sweepga(&[
        "data/scerevisiae8.fa.gz",
        "-j",
        "100k", // Enable scaffolding with 100kb jump
        "-s",
        "50k", // Minimum scaffold mass
        "-t",
        "2",
        "--self",
        "--paf",
    ]);

    // Without scaffolding
    let result2 = run_sweepga(&[
        "data/scerevisiae8.fa.gz",
        "-j",
        "0", // Disable scaffolding
        "-t",
        "2",
        "--self",
        "--paf",
    ]);

    assert!(result1.is_ok() && result2.is_ok(), "Scaffold tests failed");

    // Write outputs to files
    fs::write(&with_scaffold, result1.unwrap()).unwrap();
    fs::write(&no_scaffold, result2.unwrap()).unwrap();

    // Check for scaffold annotations
    let scaffold_content = fs::read_to_string(&with_scaffold).unwrap();
    assert!(
        scaffold_content.contains("st:Z:scaffold"),
        "Missing scaffold annotations"
    );
    assert!(
        scaffold_content.contains("ch:Z:chain"),
        "Missing chain annotations"
    );
}

#[test]
fn test_empty_input_handling() {
    let temp_dir = TempDir::new().unwrap();
    let empty_fa = temp_dir.path().join("empty.fa");
    let output = temp_dir.path().join("empty_out.paf");

    // Create empty FASTA
    fs::write(&empty_fa, "").unwrap();

    // This should handle gracefully
    let result = run_sweepga(&[empty_fa.to_str().unwrap(), "-t", "1", "--paf"]);

    // Should either succeed with no output or fail gracefully
    // Not asserting success as behavior may vary, but shouldn't panic
    if let Ok(content) = result {
        fs::write(&output, content).unwrap();
        if output.exists() {
            let content = fs::read_to_string(&output).unwrap();
            assert!(
                content.is_empty() || content.lines().count() == 0,
                "Empty input produced alignments"
            );
        }
    }
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_large_sequence_handling() {
    use self::synthetic_genomes::generate_base_sequence;

    let temp_dir = TempDir::new().unwrap();
    let large_fa = temp_dir.path().join("large.fa");
    let output = temp_dir.path().join("large_out.paf");

    // Create a large sequence with internal repeats (100kb)
    let base = generate_base_sequence(50000, 999);
    // Duplicate part of it to create self-similarity
    let seq = format!("{}{}", &base[..25000], &base);
    fs::write(&large_fa, format!(">large_seq\n{seq}\n")).unwrap();

    // Should handle without hanging or crashing
    let result = run_sweepga(&[
        large_fa.to_str().unwrap(),
        "-t",
        "2",
        "--self", // Include self-mappings
        "--paf",  // Request PAF output for test
    ]);

    if let Err(ref e) = result {
        eprintln!("Large sequence test error: {e}");
    }

    assert!(result.is_ok(), "Failed on large sequence: {result:?}");

    // Write stdout to output file
    fs::write(&output, result.unwrap()).unwrap();
    if output.exists() {
        // Should find self-alignments in the duplicated regions
        let content = fs::read_to_string(&output).unwrap();
        assert!(
            !content.is_empty(),
            "Should find alignments in duplicated regions"
        );
        assert!(
            has_extended_cigar(&output),
            "Invalid output for large sequence"
        );
    }
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_multisequence_fasta() {
    use self::synthetic_genomes::{generate_base_sequence, mutate_sequence};

    let temp_dir = TempDir::new().unwrap();
    let multi_fa = temp_dir.path().join("multi.fa");
    let output = temp_dir.path().join("multi_out.paf");

    // Create multi-sequence FASTA with related sequences
    // Use 20kb sequences based on what works for pairwise alignment
    let base = generate_base_sequence(20000, 123);
    let seq1 = base.clone();
    let seq2 = mutate_sequence(&base, 200, 124); // 1% divergence (200/20000)
    let seq3 = mutate_sequence(&base, 400, 125); // 2% divergence (400/20000)

    fs::write(
        &multi_fa,
        format!(">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n"),
    )
    .unwrap();

    let result = run_sweepga(&[
        multi_fa.to_str().unwrap(),
        "-t",
        "1",
        "--self", // Include self-mappings
        "--paf",  // Request PAF output for test
    ]);

    if let Err(ref e) = result {
        eprintln!("Multi-sequence test error: {e}");
    }

    assert!(result.is_ok(), "Failed on multi-sequence FASTA: {result:?}");

    // Write stdout to output file
    fs::write(&output, result.unwrap()).unwrap();
    if output.exists() {
        let content = fs::read_to_string(&output).unwrap();
        // Should have alignments between the related sequences
        assert!(!content.is_empty(), "Should produce alignments");

        // Check that we found alignments between sequences
        let lines: Vec<&str> = content.lines().collect();
        assert!(!lines.is_empty(), "Should have at least one alignment");
    }
}

#[test]
// FastGA binaries are now available via binary_paths
fn test_performance_regression() {
    use std::time::Instant;

    let temp_dir = TempDir::new().unwrap();
    let output = temp_dir.path().join("perf.paf");

    let start = Instant::now();

    let result = run_sweepga(&["data/scerevisiae8.fa.gz", "-t", "4", "--self", "--paf"]);

    let duration = start.elapsed();

    assert!(result.is_ok(), "Performance test failed");

    // Write stdout to output file
    fs::write(&output, result.unwrap()).unwrap();

    // Yeast self-alignment should complete in reasonable time
    assert!(
        duration.as_secs() < 60,
        "Alignment took too long: {duration:?}"
    );

    // Should produce expected number of alignments (8 yeast genomes self-alignment)
    // With --self flag and no filtering, expect 15k-20k alignments
    let line_count = count_lines(&output);
    assert!(
        line_count > 10000 && line_count < 25000,
        "Unexpected alignment count: {line_count}"
    );
}
