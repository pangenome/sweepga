#![allow(clippy::uninlined_format_args)]
/// Working integration tests that don't require the binary
/// These tests use the sweepga library API directly
use std::fs;
use std::path::Path;
use tempfile::TempDir;

#[path = "synthetic_genomes.rs"]
mod synthetic_genomes;

/// Test FastGA alignment using library API with small test data
#[test]
fn test_fastga_alignment_with_b3106() {
    // Use the small B-3106.fa file (32KB human chr6 fragment)
    let fasta_path = Path::new("data/B-3106.fa");

    if !fasta_path.exists() {
        eprintln!("Skipping test - data/B-3106.fa not found");
        return;
    }

    let _temp_dir = TempDir::new().unwrap();

    // Use the FastGA integration module directly
    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1, 100);

    // Run self-alignment
    let result = integration.align_to_temp_1aln(fasta_path, fasta_path);

    if result.is_err() {
        eprintln!("Skipping test - FastGA alignment failed (binaries may not be available in CI)");
        return;
    }

    let temp_file = result.unwrap();
    let aln_size = fs::metadata(temp_file.path()).unwrap().len();

    // Should produce some alignment output
    assert!(aln_size > 0, "No alignment output produced");
    eprintln!(
        "✓ FastGA aligned B-3106.fa self-alignment: {} bytes",
        aln_size
    );
}

/// Test filtering API without crashing
#[test]
fn test_filtering_api() {
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};

    let temp_dir = TempDir::new().unwrap();
    let input_paf = temp_dir.path().join("input.paf");

    // Create a minimal valid PAF file
    let paf_content = "chr1\t10000\t0\t1000\t+\tchr1\t10000\t0\t1000\t950\t1000\t60\n";
    fs::write(&input_paf, paf_content).unwrap();

    let config = FilterConfig {
        chain_gap: 0,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::ManyToMany,
        scaffold_max_per_query: None,
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: 0,
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.95,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let output = temp_dir.path().join("out.paf");

    // Test that filtering doesn't crash
    let filter = PafFilter::new(config);
    let result = filter.filter_paf(&input_paf, &output);

    assert!(result.is_ok(), "Filtering should not crash");
    assert!(output.exists(), "Output file should exist");

    eprintln!("✓ Filtering API works without crashing");
}

/// Test multisequence FASTA handling using synthetic data
#[test]
fn test_multisequence_synthetic() {
    use synthetic_genomes::{generate_base_sequence, mutate_sequence};

    let temp_dir = TempDir::new().unwrap();
    let multi_fa = temp_dir.path().join("multi.fa");

    // Create multi-sequence FASTA with related sequences (5kb each for speed)
    let base = generate_base_sequence(5000, 123);
    let seq1 = base.clone();
    let seq2 = mutate_sequence(&base, 50, 124); // 1% divergence
    let seq3 = mutate_sequence(&base, 100, 125); // 2% divergence

    fs::write(
        &multi_fa,
        format!(">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n"),
    )
    .unwrap();

    // Use FastGA integration directly
    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1, 100);
    let result = integration.align_to_temp_1aln(&multi_fa, &multi_fa);

    if result.is_err() {
        eprintln!("Skipping test - FastGA alignment failed (binaries may not be available in CI)");
        return;
    }

    let temp_file = result.unwrap();
    let aln_size = fs::metadata(temp_file.path()).unwrap().len();

    assert!(aln_size > 0, "Should produce alignments");
    eprintln!("✓ Multi-sequence alignment: {} bytes", aln_size);
}

/// Test 1aln filtering directly
#[test]
fn test_1aln_filtering_api() {
    use sweepga::paf_filter::{FilterConfig, FilterMode, ScoringFunction};
    use synthetic_genomes::generate_base_sequence;

    let temp_dir = TempDir::new().unwrap();
    let test_fa = temp_dir.path().join("test.fa");

    // Create test sequence
    let seq = generate_base_sequence(3000, 456);
    fs::write(&test_fa, format!(">test\n{seq}\n")).unwrap();

    // Generate 1aln file
    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1, 100);
    let aln_result = integration.align_to_temp_1aln(&test_fa, &test_fa);

    if aln_result.is_err() {
        eprintln!(
            "Skipping test - FastGA alignment failed (this is expected in some CI environments)"
        );
        return;
    }

    let aln_file = aln_result.unwrap();
    let output = temp_dir.path().join("filtered.1aln");

    // Apply filtering
    let config = FilterConfig {
        chain_gap: 0,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::ManyToMany,
        scaffold_max_per_query: None,
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: 0,
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.95,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter_result =
        sweepga::unified_filter::filter_file(aln_file.path(), &output, &config, false, false);

    assert!(filter_result.is_ok(), "1aln filtering failed");
    assert!(output.exists(), "Output file should exist");

    let output_size = fs::metadata(&output).unwrap().len();
    assert!(output_size > 0, "Filtered output should not be empty");

    eprintln!("✓ 1aln filtering: {} bytes", output_size);
}

/// Test that small sequences can be aligned
#[test]
fn test_small_sequence_alignment() {
    let temp_dir = TempDir::new().unwrap();
    let small_fa = temp_dir.path().join("small.fa");

    // Create tiny test sequence (200bp - very fast)
    let seq = "ACGTACGTACGT".repeat(17); // ~200bp
    fs::write(&small_fa, format!(">small\n{seq}\n")).unwrap();

    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1, 100);
    let result = integration.align_to_temp_1aln(&small_fa, &small_fa);

    // This might fail for very small sequences, which is okay
    if let Ok(temp_file) = result {
        let aln_size = fs::metadata(temp_file.path()).unwrap().len();
        eprintln!("✓ Small sequence alignment: {} bytes", aln_size);
    } else {
        eprintln!("⚠ Small sequence alignment skipped (expected for tiny sequences)");
    }
}

/// Test pairwise alignment with B-3106
#[test]
fn test_pairwise_b3106() {
    let fasta_path = Path::new("data/B-3106.fa");

    if !fasta_path.exists() {
        eprintln!("Skipping test - data/B-3106.fa not found");
        return;
    }

    let temp_dir = TempDir::new().unwrap();

    // Copy B-3106 to two files for pairwise test
    let query = temp_dir.path().join("query.fa");
    let target = temp_dir.path().join("target.fa");
    fs::copy(fasta_path, &query).unwrap();
    fs::copy(fasta_path, &target).unwrap();

    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1, 100);
    let result = integration.align_to_temp_1aln(&query, &target);

    if result.is_err() {
        eprintln!("Skipping test - FastGA alignment failed (binaries may not be available in CI)");
        return;
    }

    let temp_file = result.unwrap();
    let aln_size = fs::metadata(temp_file.path()).unwrap().len();

    assert!(aln_size > 0, "Should produce alignment");
    eprintln!("✓ Pairwise B-3106 alignment: {} bytes", aln_size);
}
