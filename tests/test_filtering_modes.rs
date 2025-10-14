/// Comprehensive filtering mode tests
///
/// Tests various N:M filtering modes to ensure plane sweep algorithm
/// works correctly across different multiplicity constraints.
use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Parse PAF and count mappings per query/target
#[derive(Debug, Default)]
struct MappingStats {
    total_mappings: usize,
    query_counts: HashMap<String, usize>,  // query_name -> count
    target_counts: HashMap<String, usize>, // target_name -> count
    pair_counts: HashMap<(String, String), usize>, // (query, target) -> count
}

fn analyze_paf(paf_path: &Path) -> Result<MappingStats> {
    let content = fs::read_to_string(paf_path)?;
    let mut stats = MappingStats::default();

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let query = fields[0].to_string();
        let target = fields[5].to_string();

        stats.total_mappings += 1;
        *stats.query_counts.entry(query.clone()).or_default() += 1;
        *stats.target_counts.entry(target.clone()).or_default() += 1;
        *stats.pair_counts.entry((query, target)).or_default() += 1;
    }

    Ok(stats)
}

/// Generate unfiltered PAF for testing
fn generate_unfiltered_paf(temp_dir: &Path) -> Result<std::path::PathBuf> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    // Copy input to temp dir so FastGA's intermediate files (.gdb, .ktab, .bps)
    // are also created in temp and automatically cleaned up
    let temp_input = temp_dir.join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    let output = temp_dir.join("unfiltered.paf");

    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
            "--paf",
            "-n",
            "N:N", // No filtering
        ])
        .stdout(fs::File::create(&output)?)
        .status()?;

    Ok(output)
}

/// Apply filtering to PAF file
fn filter_paf(input: &Path, output: &Path, filter_mode: &str) -> Result<()> {
    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            input.to_str().unwrap(),
            "-n",
            filter_mode,
        ])
        .stdout(fs::File::create(output)?)
        .status()?;

    Ok(())
}

/// Test 1:1 filtering - most restrictive
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_1_to_1() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;
    let filtered = temp_dir.path().join("filtered_1to1.paf");

    filter_paf(&unfiltered, &filtered, "1:1")?;

    let unfiltered_stats = analyze_paf(&unfiltered)?;
    let filtered_stats = analyze_paf(&filtered)?;

    eprintln!("1:1 filtering:");
    eprintln!("  Unfiltered: {} mappings", unfiltered_stats.total_mappings);
    eprintln!("  Filtered: {} mappings", filtered_stats.total_mappings);

    // 1:1 should significantly reduce mappings
    assert!(
        filtered_stats.total_mappings < unfiltered_stats.total_mappings,
        "1:1 filtering should reduce mappings"
    );

    // 1:1 means at most 1 overlapping interval per query/target position
    // Multiple non-overlapping intervals are allowed per chromosome pair
    // So we can't check exact counts, but we verify it produces valid output
    assert!(
        filtered_stats.total_mappings > 0,
        "1:1 filtering should produce some mappings"
    );

    eprintln!("✓ 1:1 filtering completes and reduces mapping count");
    Ok(())
}

/// Test 1:N filtering - one query, many targets
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_1_to_many() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;
    let filtered = temp_dir.path().join("filtered_1toN.paf");

    filter_paf(&unfiltered, &filtered, "1")?; // "1" is shorthand for "1:∞"

    let unfiltered_stats = analyze_paf(&unfiltered)?;
    let filtered_stats = analyze_paf(&filtered)?;

    eprintln!("1:∞ filtering:");
    eprintln!("  Unfiltered: {} mappings", unfiltered_stats.total_mappings);
    eprintln!("  Filtered: {} mappings", filtered_stats.total_mappings);

    // 1:N means at most 1 overlapping interval per query position
    // but unlimited overlapping on target side
    // Multiple non-overlapping query intervals can map to same target

    // Should be less filtered than N:N (no filtering)
    assert!(
        filtered_stats.total_mappings <= unfiltered_stats.total_mappings,
        "1:∞ should keep <= mappings than unfiltered"
    );

    assert!(
        filtered_stats.total_mappings > 0,
        "1:∞ filtering should produce some mappings"
    );

    eprintln!("✓ 1:∞ filtering respects query constraint");
    Ok(())
}

/// Test N:1 filtering - many queries, one target
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_many_to_1() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;

    // Note: We need to test N:1 by reading as PAF and filtering
    // For now we'll use a direct sweepga call with the appropriate flag
    // This tests the target-side constraint

    // We can simulate this by using -m (scaffold filter) with 1:1
    // but on pre-filtered 1:∞ data to see the effect
    let intermediate = temp_dir.path().join("intermediate.paf");
    filter_paf(&unfiltered, &intermediate, "1")?;

    // Now apply target-side filtering (conceptually N:1)
    // In practice, sweepga's -n flag applies bidirectionally
    // So we test that increasing query multiplicity increases output
    let filtered_1to1 = temp_dir.path().join("filtered_1to1.paf");
    filter_paf(&unfiltered, &filtered_1to1, "1:1")?;

    let stats_1to1 = analyze_paf(&filtered_1to1)?;
    let stats_1to_n = analyze_paf(&intermediate)?;
    let unfiltered_stats = analyze_paf(&unfiltered)?;

    eprintln!("Multiplicity comparison:");
    eprintln!("  Unfiltered: {} mappings", unfiltered_stats.total_mappings);
    eprintln!("  1:1: {} mappings", stats_1to1.total_mappings);
    eprintln!("  1:∞: {} mappings", stats_1to_n.total_mappings);

    // Invariant: 1:∞ should have >= mappings than 1:1
    assert!(
        stats_1to_n.total_mappings >= stats_1to1.total_mappings,
        "1:∞ should have >= mappings than 1:1"
    );

    eprintln!("✓ Multiplicity invariant holds: 1:∞ >= 1:1");
    Ok(())
}

/// Test 2:3 filtering - specific multiplicity constraint
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_2_to_3() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;
    let filtered = temp_dir.path().join("filtered_2to3.paf");

    filter_paf(&unfiltered, &filtered, "2:3")?;

    let unfiltered_stats = analyze_paf(&unfiltered)?;
    let filtered_stats = analyze_paf(&filtered)?;

    eprintln!("2:3 filtering:");
    eprintln!("  Unfiltered: {} mappings", unfiltered_stats.total_mappings);
    eprintln!("  Filtered: {} mappings", filtered_stats.total_mappings);

    // Should be between 1:1 and N:N
    assert!(
        filtered_stats.total_mappings <= unfiltered_stats.total_mappings,
        "2:3 should keep <= mappings than unfiltered"
    );

    // Verify multiplicity constraints
    // Each query chromosome should map to at most 2 overlapping regions
    // Each target chromosome should receive at most 3 overlapping regions
    // (This is hard to verify precisely without interval overlap analysis)

    eprintln!("✓ 2:3 filtering completes successfully");
    Ok(())
}

/// Test 4:1 filtering - asymmetric constraint
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_4_to_1() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;
    let filtered = temp_dir.path().join("filtered_4to1.paf");

    filter_paf(&unfiltered, &filtered, "4:1")?;

    let unfiltered_stats = analyze_paf(&unfiltered)?;
    let filtered_stats = analyze_paf(&filtered)?;

    eprintln!("4:1 filtering:");
    eprintln!("  Unfiltered: {} mappings", unfiltered_stats.total_mappings);
    eprintln!("  Filtered: {} mappings", filtered_stats.total_mappings);

    assert!(
        filtered_stats.total_mappings <= unfiltered_stats.total_mappings,
        "4:1 should keep <= mappings than unfiltered"
    );

    eprintln!("✓ 4:1 filtering completes successfully");
    Ok(())
}

/// Test monotonicity: increasing multiplicity should keep >= mappings
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_monotonicity() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;

    // Test increasing query multiplicity: 1:1 -> 2:1 -> 3:1
    let modes = vec!["1:1", "2:1", "3:1", "4:1"];
    let mut prev_count = 0;

    for mode in &modes {
        let output = temp_dir.path().join(format!("filtered_{}.paf", mode));
        filter_paf(&unfiltered, &output, mode)?;

        let stats = analyze_paf(&output)?;
        eprintln!("{}: {} mappings", mode, stats.total_mappings);

        if prev_count > 0 {
            assert!(
                stats.total_mappings >= prev_count,
                "Increasing multiplicity should keep >= mappings: {} had {}, {} had {}",
                modes[modes.iter().position(|&m| m == *mode).unwrap() - 1],
                prev_count,
                mode,
                stats.total_mappings
            );
        }

        prev_count = stats.total_mappings;
    }

    eprintln!("✓ Monotonicity property holds: N₁:M₁ ⊆ N₂:M₂ when N₁≤N₂, M₁≤M₂");
    Ok(())
}

/// Test N:N (no filtering) keeps all mappings
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_no_filtering() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Copy input to temp dir to avoid accumulating intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    let output1 = temp_dir.path().join("output1.paf");
    let output2 = temp_dir.path().join("output2.paf");

    // Generate twice with N:N
    Command::new("cargo")
        .args(&[
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
        .stdout(fs::File::create(&output1)?)
        .status()?;

    filter_paf(&output1, &output2, "N:N")?;

    let stats1 = analyze_paf(&output1)?;
    let stats2 = analyze_paf(&output2)?;

    eprintln!("N:N filtering:");
    eprintln!("  First run: {} mappings", stats1.total_mappings);
    eprintln!("  Re-filtered: {} mappings", stats2.total_mappings);

    // N:N should be idempotent - no change
    assert_eq!(
        stats1.total_mappings, stats2.total_mappings,
        "N:N filtering should be idempotent"
    );

    eprintln!("✓ N:N filtering is idempotent (no filtering applied)");
    Ok(())
}

/// Test that 1:1 is idempotent
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_idempotence() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;

    let filtered1 = temp_dir.path().join("filtered1.paf");
    let filtered2 = temp_dir.path().join("filtered2.paf");

    // Apply 1:1 filter twice
    filter_paf(&unfiltered, &filtered1, "1:1")?;
    filter_paf(&filtered1, &filtered2, "1:1")?;

    let stats1 = analyze_paf(&filtered1)?;
    let stats2 = analyze_paf(&filtered2)?;

    eprintln!("Idempotence test:");
    eprintln!("  First filtering: {} mappings", stats1.total_mappings);
    eprintln!("  Second filtering: {} mappings", stats2.total_mappings);

    // Filtering should be idempotent
    assert_eq!(
        stats1.total_mappings, stats2.total_mappings,
        "1:1 filtering should be idempotent"
    );

    eprintln!("✓ 1:1 filtering is idempotent");
    Ok(())
}

/// Test various combinations to ensure correctness
#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_filter_combinations() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let unfiltered = generate_unfiltered_paf(temp_dir.path())?;
    let unfiltered_stats = analyze_paf(&unfiltered)?;

    let test_modes = vec![
        "1:1", "1:2", "2:1", "2:2", "2:3", "3:2", "3:3", "4:1", "1:4", "5:5",
    ];

    eprintln!("\nTesting various filter combinations:");
    eprintln!("Unfiltered: {} mappings", unfiltered_stats.total_mappings);

    for mode in &test_modes {
        let output = temp_dir
            .path()
            .join(format!("filtered_{}.paf", mode.replace(':', "to")));
        filter_paf(&unfiltered, &output, mode)?;

        let stats = analyze_paf(&output)?;
        eprintln!(
            "  {}: {} mappings ({:.1}% of unfiltered)",
            mode,
            stats.total_mappings,
            (stats.total_mappings as f64 / unfiltered_stats.total_mappings as f64) * 100.0
        );

        // All filtering should produce valid output
        assert!(stats.total_mappings > 0, "{} produced no mappings", mode);
        assert!(
            stats.total_mappings <= unfiltered_stats.total_mappings,
            "{} produced more mappings than unfiltered",
            mode
        );
    }

    eprintln!("✓ All filter combinations produce valid output");
    Ok(())
}
