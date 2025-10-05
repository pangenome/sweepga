/// Tests for chaining stability across different gap values
///
/// Key property: When gap threshold increases, chains should only grow or stay the same,
/// never shrink. If two mappings are 50kb apart, they should chain with both -j 100k and -j 1m.
use std::collections::HashMap;

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
            if field.starts_with("ch:Z:") {
                chain_id = Some(field[5..].to_string());
                break;
            }
        }

        if let Some(cid) = chain_id {
            // Create a unique mapping ID from query and coordinates
            let mapping_id = format!("{}:{}-{}", fields[0], fields[2], fields[3]);
            chains.entry(cid).or_insert_with(Vec::new).push(mapping_id);
        }
    }

    chains
}

#[test]
fn test_chaining_monotonicity() {
    // Test that larger gap values create supersets of chains from smaller gaps

    // Run with different gap values
    let gaps = vec![10_000, 50_000, 100_000, 500_000, 1_000_000];
    let mut all_chains = Vec::new();

    for gap in &gaps {
        let output = std::process::Command::new("./target/release/sweepga")
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
        let (gap1, _, count1) = &all_chains[i-1];
        let (gap2, _, count2) = &all_chains[i];

        assert!(
            count2 >= count1,
            "Chain membership should not decrease with larger gaps: \
             -j {} has {} members, but -j {} has {} members",
            gap1, count1, gap2, count2
        );
    }
}

#[test]
fn test_chain_identity_stability() {
    // Test that chain identities remain reasonable with large gaps

    let gaps = vec![10_000, 100_000, 1_000_000];

    for gap in gaps {
        let output = std::process::Command::new("./target/release/sweepga")
            .arg("data/scerevisiae8.fa.gz")
            .arg("-j")
            .arg(gap.to_string())
            .arg("-i")
            .arg("0.90") // 90% identity threshold
            .output()
            .expect("Failed to run sweepga");

        let stderr = String::from_utf8_lossy(&output.stderr);

        // Extract coverage from output
        let mut coverage = None;
        for line in stderr.lines() {
            if line.contains("Output:") && line.contains("Mb total") {
                // Parse: "Output: 2778.5 Mb total, 94.8% avg identity"
                if let Some(parts) = line.split("Output:").nth(1) {
                    if let Some(mb_str) = parts.trim().split_whitespace().next() {
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

    // Create a simple test case with three collinear mappings:
    // Mapping A: query 0-1000, target 0-1000
    // Mapping B: query 1100-2100, target 1100-2100 (100bp gap from A)
    // Mapping C: query 5000-6000, target 5000-6000 (3900bp gap from B)

    // With -j 10000, all three should be chainable, but:
    // - A should chain to B (nearest at 100bp)
    // - B should chain to C (next nearest at 3900bp)
    // Result: Single chain A-B-C

    // This is harder to test without creating custom PAF input
    // TODO: Create a test PAF file with known structure
}

#[test]
fn test_overlap_penalty() {
    // Test that overlapping mappings are penalized correctly

    // Create test case with:
    // Mapping A: query 0-1000, target 0-1000
    // Mapping B: query 900-1900, target 900-1900 (100bp overlap with A)
    // Mapping C: query 1100-2100, target 1100-2100 (100bp gap from A)

    // With overlap penalty, A should prefer to chain to C (gap) over B (overlap)

    // TODO: Create test PAF file with overlapping mappings
}