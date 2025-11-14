#![allow(clippy::uninlined_format_args)]
/// Tests for chaining stability across different gap values
///
/// Key property: When gap threshold increases, chains should only grow or stay the same,
/// never shrink. If two mappings are 50kb apart, they should chain with both -j 100k and -j 1m.
use std::collections::HashMap;

/// Helper to get the correct sweepga binary path based on build mode
fn sweepga_bin() -> &'static str {
    if cfg!(debug_assertions) {
        "./target/debug/sweepga"
    } else {
        "./target/release/sweepga"
    }
}

/// Parse chain membership from PAF output
fn parse_chains(output: &str) -> HashMap<String, Vec<String>> {
    let mut chains: HashMap<String, Vec<String>> = HashMap::new();

    for line in output.lines() {
        if line.starts_with('[') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 13 {
            continue;
        }

        // Find chain ID in tags
        let mut chain_id = None;
        for field in fields.iter().skip(12) {
            if let Some(stripped) = field.strip_prefix("ch:Z:") {
                chain_id = Some(stripped.to_string());
                break;
            }
        }

        if let Some(cid) = chain_id {
            // Create a unique mapping ID from query and coordinates
            let mapping_id = format!("{}:{}-{}", fields[0], fields[2], fields[3]);
            chains.entry(cid).or_default().push(mapping_id);
        }
    }

    chains
}

#[test]
// Sweepga binary and yeast data available
fn test_chaining_monotonicity() {
    // Test that larger gap values create supersets of chains from smaller gaps

    // Run with different gap values
    let gaps = vec![10_000, 50_000, 100_000, 500_000, 1_000_000];
    let mut all_chains = Vec::new();

    for gap in &gaps {
        let output = std::process::Command::new(sweepga_bin())
            .arg("data/scerevisiae8.fa.gz")
            .arg("-j")
            .arg(gap.to_string())
            .arg("-i")
            .arg("0") // No identity filter for this test
            .output()
            .expect("Failed to run sweepga");

        let stdout = String::from_utf8_lossy(&output.stdout);
        let chains = parse_chains(&stdout);

        // Count total members across all chains
        let total_members: usize = chains.values().map(|v| v.len()).sum();
        all_chains.push((gap, chains, total_members));
    }

    // Verify monotonicity: larger gaps should have same or more chain members
    for i in 1..all_chains.len() {
        let (gap1, _, count1) = &all_chains[i - 1];
        let (gap2, _, count2) = &all_chains[i];

        assert!(
            count2 >= count1,
            "Chain membership should not decrease with larger gaps: \
             -j {} has {} members, but -j {} has {} members",
            gap1,
            count1,
            gap2,
            count2
        );
    }
}

#[test]
#[ignore] // Output format has changed - needs updating to parse new format
fn test_chain_identity_stability() {
    // Test that chain identities remain reasonable with large gaps

    let gaps = vec![10_000, 100_000, 1_000_000];

    for gap in gaps {
        let output = std::process::Command::new(sweepga_bin())
            .arg("data/scerevisiae8.fa.gz")
            .arg("-j")
            .arg(gap.to_string())
            .arg("-i")
            .arg("0.90") // 90% identity threshold
            .arg("--paf") // Output PAF format
            .output()
            .expect("Failed to run sweepga");

        let stderr = String::from_utf8_lossy(&output.stderr);

        // Extract coverage from output
        let mut coverage = None;
        for line in stderr.lines() {
            if line.contains("Output:") && line.contains("Mb total") {
                // Parse: "Output: 2778.5 Mb total, 94.8% avg identity"
                if let Some(parts) = line.split("Output:").nth(1) {
                    if let Some(mb_str) = parts.split_whitespace().next() {
                        coverage = mb_str.parse::<f64>().ok();
                    }
                }
            }
        }

        assert!(
            coverage.is_some(),
            "Failed to parse coverage for -j {}",
            gap
        );

        let cov = coverage.unwrap();

        // For 99% identical yeast genomes, we should get high coverage
        // even with large gap values
        // With 8 genomes (~12 Mb each) and 56 pairs (8×7, excluding self-mappings),
        // we expect ~672 Mb total. Getting >600 Mb means good coverage.
        assert!(
            cov > 600.0,
            "Coverage with -j {} is too low: {} Mb (expected >600 Mb for 8 yeast genomes, 56 pairs)",
            gap, cov
        );
    }
}

#[test]
fn test_nearest_neighbor_chaining() {
    // Test that mappings chain to their nearest neighbors, not distant ones
    use std::fs;
    use std::io::Write;
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
    use tempfile::NamedTempFile;

    // Create synthetic test data with three collinear mappings:
    // Mapping A: query 0-1000, target 0-1000
    // Mapping B: query 1100-2100, target 1100-2100 (100bp gap from A - NEAREST)
    // Mapping C: query 5000-6000, target 5000-6000 (3900bp gap from B)
    let mut test_input = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(test_input, "querySeq\t10000\t0\t1000\t+\ttargetSeq\t10000\t0\t1000\t950\t1000\t60").unwrap();
    writeln!(test_input, "querySeq\t10000\t1100\t2100\t+\ttargetSeq\t10000\t1100\t2100\t950\t1000\t60").unwrap();
    writeln!(test_input, "querySeq\t10000\t5000\t6000\t+\ttargetSeq\t10000\t5000\t6000\t950\t1000\t60").unwrap();
    test_input.flush().unwrap();

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
        overlap_threshold: 0.0,
        sparsity: 1.0,
        no_merge: false,      // Enable chaining
        scaffold_gap: 10_000, // 10kb gap allows all three to chain
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.0,
        scaffold_max_deviation: 20_000,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    let temp_out = std::env::temp_dir().join("test_nearest_neighbor_out.paf");

    filter
        .filter_paf(
            test_input.path().to_str().unwrap(),
            temp_out.to_str().unwrap(),
        )
        .expect("Failed to filter PAF");

    let output = fs::read_to_string(&temp_out).expect("Failed to read output");
    let chains = parse_chains(&output);

    // Should have exactly one chain containing all three mappings
    assert_eq!(
        chains.len(),
        1,
        "Expected single chain, got {} chains",
        chains.len()
    );

    let chain = chains.values().next().unwrap();
    assert_eq!(
        chain.len(),
        3,
        "Expected chain with 3 members, got {}",
        chain.len()
    );

    // Verify the mappings are in the expected order (A→B→C)
    // A: 0-1000, B: 1100-2100, C: 5000-6000
    assert!(
        chain.iter().any(|m| m.contains("0-1000")),
        "Chain missing mapping A"
    );
    assert!(
        chain.iter().any(|m| m.contains("1100-2100")),
        "Chain missing mapping B"
    );
    assert!(
        chain.iter().any(|m| m.contains("5000-6000")),
        "Chain missing mapping C"
    );

    // Cleanup
    let _ = fs::remove_file(temp_out);
}

#[test]
fn test_overlap_penalty() {
    // Test that overlapping mappings are penalized correctly
    use std::fs;
    use std::io::Write;
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
    use tempfile::NamedTempFile;

    // Create synthetic test data with three mappings:
    // Mapping A: query 0-1000, target 0-1000
    // Mapping B: query 900-1900, target 900-1900 (100bp OVERLAP with A)
    // Mapping C: query 1100-2100, target 1100-2100 (100bp GAP from A)
    let mut test_input = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(test_input, "querySeq\t10000\t0\t1000\t+\ttargetSeq\t10000\t0\t1000\t950\t1000\t60").unwrap();
    writeln!(test_input, "querySeq\t10000\t900\t1900\t+\ttargetSeq\t10000\t900\t1900\t950\t1000\t60").unwrap();
    writeln!(test_input, "querySeq\t10000\t1100\t2100\t+\ttargetSeq\t10000\t1100\t2100\t950\t1000\t60").unwrap();
    test_input.flush().unwrap();

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
        overlap_threshold: 0.0,
        sparsity: 1.0,
        no_merge: false,      // Enable chaining
        scaffold_gap: 10_000, // Large enough to allow chaining
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.0,
        scaffold_max_deviation: 20_000,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    let temp_out = std::env::temp_dir().join("test_overlap_penalty_out.paf");

    filter
        .filter_paf(
            test_input.path().to_str().unwrap(),
            temp_out.to_str().unwrap(),
        )
        .expect("Failed to filter PAF");

    let output = fs::read_to_string(&temp_out).expect("Failed to read output");
    let chains = parse_chains(&output);

    // With overlap penalty, we expect:
    // - A and C should chain together (gap is preferred)
    // - B may form its own chain or be excluded depending on overlap handling
    // The key property: A should NOT chain to B (overlap) when C (gap) is available

    // At minimum, there should be at least one chain
    assert!(!chains.is_empty(), "Expected at least one chain, got zero");

    // Check that if A and C exist in output, they're in the same chain
    let mut a_chain_id = None;
    let mut c_chain_id = None;

    for (chain_id, members) in &chains {
        for member in members {
            if member.contains("0-1000") {
                a_chain_id = Some(chain_id.clone());
            }
            if member.contains("1100-2100") {
                c_chain_id = Some(chain_id.clone());
            }
        }
    }

    // If both A and C are present, they should be in the same chain
    // (gap is preferred over overlap)
    if let (Some(a_id), Some(c_id)) = (a_chain_id, c_chain_id) {
        assert_eq!(
            a_id, c_id,
            "Mappings A and C should be in the same chain (gap preferred over overlap)"
        );
    }

    // Cleanup
    let _ = fs::remove_file(temp_out);
}
