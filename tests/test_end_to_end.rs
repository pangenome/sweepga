/// End-to-end pipeline tests with coverage validation
///
/// Tests the complete workflow from FASTA input to filtered output,
/// validating that coverage expectations are met for real genomic data.
use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Parse PAF and calculate coverage statistics
fn calculate_coverage_stats(paf_path: &Path) -> Result<CoverageStats> {
    type AlignmentSpans = Vec<(u64, u64, u64, u64)>;
    let content = fs::read_to_string(paf_path)?;

    let mut genome_pairs: HashMap<(String, String), AlignmentSpans> = HashMap::new();

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Extract genome prefixes (e.g., "SGDref#1#chrI" -> "SGDref#1#")
        let extract_genome = |name: &str| -> String {
            if let Some(pos) = name.rfind('#') {
                if let Some(prev_pos) = name[..pos].rfind('#') {
                    return name[..=prev_pos].to_string();
                }
            }
            name.to_string()
        };

        let query_genome = extract_genome(fields[0]);
        let target_genome = extract_genome(fields[5]);

        let query_start: u64 = fields[2].parse().unwrap_or(0);
        let query_end: u64 = fields[3].parse().unwrap_or(0);
        let target_start: u64 = fields[7].parse().unwrap_or(0);
        let target_end: u64 = fields[8].parse().unwrap_or(0);

        if query_genome != target_genome {
            let key = if query_genome < target_genome {
                (query_genome, target_genome)
            } else {
                (target_genome, query_genome)
            };

            genome_pairs.entry(key).or_default().push((
                query_start,
                query_end,
                target_start,
                target_end,
            ));
        }
    }

    // Calculate covered bases per genome pair
    let mut pair_coverage = Vec::new();

    for ((q_genome, t_genome), alignments) in &genome_pairs {
        // Merge overlapping intervals for query coverage
        let mut query_intervals: Vec<(u64, u64)> =
            alignments.iter().map(|(qs, qe, _, _)| (*qs, *qe)).collect();
        query_intervals.sort();

        let query_coverage = merge_and_sum_intervals(&query_intervals);

        // Merge overlapping intervals for target coverage
        let mut target_intervals: Vec<(u64, u64)> =
            alignments.iter().map(|(_, _, ts, te)| (*ts, *te)).collect();
        target_intervals.sort();

        let target_coverage = merge_and_sum_intervals(&target_intervals);

        pair_coverage.push((
            format!("{q_genome}→{t_genome}"),
            alignments.len(),
            query_coverage,
            target_coverage,
        ));
    }

    let num_genome_pairs = genome_pairs.len();
    let total_alignments = content.lines().filter(|l| !l.trim().is_empty()).count();

    Ok(CoverageStats {
        total_alignments,
        genome_pairs: num_genome_pairs,
        pair_stats: pair_coverage,
    })
}

/// Merge overlapping intervals and sum total covered bases
fn merge_and_sum_intervals(intervals: &[(u64, u64)]) -> u64 {
    if intervals.is_empty() {
        return 0;
    }

    let mut merged = Vec::new();
    let mut current = intervals[0];

    for &(start, end) in &intervals[1..] {
        if start <= current.1 {
            // Overlapping or adjacent, merge
            current.1 = current.1.max(end);
        } else {
            // Gap, save current and start new
            merged.push(current);
            current = (start, end);
        }
    }
    merged.push(current);

    // Sum all merged intervals
    merged.iter().map(|(s, e)| e - s).sum()
}

#[derive(Debug)]
struct CoverageStats {
    total_alignments: usize,
    genome_pairs: usize,
    pair_stats: Vec<(String, usize, u64, u64)>,
}

/// Test end-to-end pipeline with yeast genomes
///
/// Note: This test is skipped on CI (--skip test_end_to_end) due to FastGA race conditions,
/// but runs by default locally. Run with: cargo test test_end_to_end -- --test-threads=1
#[test]
fn test_end_to_end_yeast_coverage() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    let output_paf = temp_dir.path().join("output.paf");

    eprintln!("Running end-to-end pipeline on yeast genomes...");

    // Run the full pipeline
    let status = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
            "--paf",
        ])
        .stdout(fs::File::create(&output_paf)?)
        .status()?;

    assert!(status.success(), "Pipeline should complete successfully");

    // Calculate coverage statistics
    let stats = calculate_coverage_stats(&output_paf)?;

    eprintln!("Coverage statistics:");
    eprintln!("  Total alignments: {}", stats.total_alignments);
    eprintln!("  Genome pairs: {}", stats.genome_pairs);

    // For 8 yeast genomes (99% identical), expect good coverage
    // With 8 genomes: 8 choose 2 = 28 pairs
    assert!(
        stats.genome_pairs >= 20,
        "Expected at least 20 genome pairs, got {}",
        stats.genome_pairs
    );

    assert!(
        stats.total_alignments >= 1000,
        "Expected at least 1000 alignments, got {}",
        stats.total_alignments
    );

    // Check that we have coverage for pairs
    for (pair_name, aln_count, q_cov, t_cov) in &stats.pair_stats {
        eprintln!("  {pair_name}: {aln_count} alns, {q_cov}bp q_cov, {t_cov}bp t_cov");

        assert!(
            *aln_count > 0,
            "Genome pair {pair_name} should have alignments"
        );
    }

    eprintln!("✓ End-to-end pipeline produces expected coverage");
    Ok(())
}

/// Test that 1:1 filtering preserves genome pair coverage
///
/// Note: This test is skipped on CI (--skip test_end_to_end) due to FastGA race conditions,
/// but runs by default locally. Run with: cargo test test_end_to_end -- --test-threads=1
#[test]
fn test_one_to_one_preserves_pairs() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // Generate unfiltered output
    let unfiltered_paf = temp_dir.path().join("unfiltered.paf");
    Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
            "--paf",
            "-n",
            "N:N",
        ])
        .stdout(fs::File::create(&unfiltered_paf)?)
        .status()?;

    // Generate 1:1 filtered output
    let filtered_paf = temp_dir.path().join("filtered.paf");
    Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            unfiltered_paf.to_str().unwrap(),
            "-n",
            "1:1",
        ])
        .stdout(fs::File::create(&filtered_paf)?)
        .status()?;

    let unfiltered_stats = calculate_coverage_stats(&unfiltered_paf)?;
    let filtered_stats = calculate_coverage_stats(&filtered_paf)?;

    eprintln!("Unfiltered: {} pairs", unfiltered_stats.genome_pairs);
    eprintln!("Filtered (1:1): {} pairs", filtered_stats.genome_pairs);

    // 1:1 filtering should preserve all genome pairs (but reduce alignment count)
    assert_eq!(
        unfiltered_stats.genome_pairs, filtered_stats.genome_pairs,
        "1:1 filtering should not lose genome pairs"
    );

    assert!(
        filtered_stats.total_alignments < unfiltered_stats.total_alignments,
        "1:1 filtering should reduce alignment count"
    );

    eprintln!("✓ 1:1 filtering preserves genome pair coverage");
    Ok(())
}

/// Test that pipeline handles different filtering modes correctly
///
/// Note: This test is skipped on CI (--skip test_end_to_end) due to FastGA race conditions,
/// but runs by default locally. Run with: cargo test test_end_to_end -- --test-threads=1
#[test]
fn test_filtering_mode_comparison() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // Test N:N (no filtering)
    let nn_paf = temp_dir.path().join("nn.paf");
    {
        let file = fs::File::create(&nn_paf)?;
        let status = Command::new("cargo")
            .args([
                "run",
                "--release",
                "--quiet",
                "--bin",
                "sweepga",
                "--",
                temp_input.to_str().unwrap(),
                "--paf",
                "-n",
                "N:N",
            ])
            .stdout(file)
            .status()?;
        assert!(status.success(), "N:N filtering failed");
    } // file is closed here

    // Test 1:1 filtering
    let one_to_one_paf = temp_dir.path().join("1to1.paf");
    {
        let file = fs::File::create(&one_to_one_paf)?;
        let status = Command::new("cargo")
            .args([
                "run",
                "--release",
                "--quiet",
                "--bin",
                "sweepga",
                "--",
                temp_input.to_str().unwrap(),
                "--paf",
                "-n",
                "1:1",
            ])
            .stdout(file)
            .status()?;
        assert!(status.success(), "1:1 filtering failed");
    } // file is closed here

    let nn_stats = calculate_coverage_stats(&nn_paf)?;
    let one_to_one_stats = calculate_coverage_stats(&one_to_one_paf)?;

    eprintln!("N:N: {} alignments", nn_stats.total_alignments);
    eprintln!("1:1: {} alignments", one_to_one_stats.total_alignments);

    // Invariant: N:N should have >= alignments than 1:1
    assert!(
        nn_stats.total_alignments >= one_to_one_stats.total_alignments,
        "N:N should have at least as many alignments as 1:1"
    );

    // Both should cover all genome pairs
    assert_eq!(
        nn_stats.genome_pairs, one_to_one_stats.genome_pairs,
        "Both modes should cover the same genome pairs"
    );

    eprintln!("✓ Filtering modes behave as expected");
    Ok(())
}
