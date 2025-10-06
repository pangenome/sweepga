/// Test that plane sweep correctly filters overlapping scaffolds
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_overlapping_scaffolds_same_chromosome_pair() {
    // Two scaffolds on same chromosome pair, overlapping on query axis
    // Only the better one should survive
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: query chr1 10-20kb → target chr1 10-20kb, 95% identity
    writeln!(paf, "chr1\t100000\t10000\t15000\t+\ttarget_chr1\t100000\t10000\t15000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();
    writeln!(paf, "chr1\t100000\t15000\t20000\t+\ttarget_chr1\t100000\t15000\t20000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();

    // Scaffold 2: query chr1 12-22kb → target chr1 30-40kb, 98% identity (BETTER, overlaps scaffold 1)
    writeln!(paf, "chr1\t100000\t12000\t17000\t+\ttarget_chr1\t100000\t30000\t35000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();
    writeln!(paf, "chr1\t100000\t17000\t22000\t+\ttarget_chr1\t100000\t35000\t40000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("1000")
        .arg("-j")
        .arg("10000")
        .arg("-i")
        .arg("0")
        .arg("-m")
        .arg("1:1")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Check which scaffold survived
    let has_scaffold_1 = stdout.contains("10000\t15000") || stdout.contains("15000\t20000");
    let has_scaffold_2 = stdout.contains("12000\t17000") || stdout.contains("17000\t22000");

    eprintln!("Has scaffold 1 (95% identity): {}", has_scaffold_1);
    eprintln!("Has scaffold 2 (98% identity): {}", has_scaffold_2);

    // Only the better scaffold (98% identity) should survive
    assert!(has_scaffold_2, "Better scaffold (98%) should be kept");
    assert!(
        !has_scaffold_1,
        "Worse scaffold (95%) should be filtered out due to overlap"
    );
}

#[test]
fn test_overlapping_scaffolds_different_targets() {
    // Two scaffolds with SAME query region but DIFFERENT target chromosomes
    // With 1:1 filtering, only one should survive per query position
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: query chr1 10-20kb → target chr1, 95% identity
    writeln!(paf, "chr1\t100000\t10000\t15000\t+\ttarget_chr1\t100000\t10000\t15000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();
    writeln!(paf, "chr1\t100000\t15000\t20000\t+\ttarget_chr1\t100000\t15000\t20000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();

    // Scaffold 2: query chr1 10-20kb → target chr2 (DIFFERENT target), 98% identity (BETTER)
    writeln!(paf, "chr1\t100000\t10000\t15000\t+\ttarget_chr2\t100000\t10000\t15000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();
    writeln!(paf, "chr1\t100000\t15000\t20000\t+\ttarget_chr2\t100000\t15000\t20000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("1000")
        .arg("-j")
        .arg("10000")
        .arg("-i")
        .arg("0")
        .arg("-m")
        .arg("1:1")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Check which scaffolds survived
    let has_chr1 = stdout.contains("target_chr1");
    let has_chr2 = stdout.contains("target_chr2");

    eprintln!("Has mapping to target_chr1 (95%): {}", has_chr1);
    eprintln!("Has mapping to target_chr2 (98%): {}", has_chr2);

    // Count output lines
    let output_count = stdout
        .lines()
        .filter(|l| !l.starts_with('[') && !l.is_empty())
        .count();
    eprintln!("Total output alignments: {}", output_count);

    // With 1:1 filtering, overlapping scaffolds to different targets should both be filtered
    // Only the BEST one should survive
    assert!(has_chr2, "Better scaffold (98% to chr2) should be kept");

    // This is the key test: does 1:1 filter ACROSS target chromosomes?
    // If both survive, then we're doing "1:1 per chromosome pair" which is wrong
    // If only chr2 survives, then we're doing "1:1 globally" which is correct
    if has_chr1 && has_chr2 {
        panic!("BUG: Both scaffolds kept even though they overlap on query axis! 1:1 filtering should pick the best one.");
    }
}

#[test]
fn test_contained_scaffold_filtering() {
    // Scaffold 1 is completely CONTAINED within scaffold 2
    // Scaffold 2 should win
    let mut paf = NamedTempFile::new().unwrap();

    // Small scaffold: query chr1 15-18kb → target chr1, 98% identity
    writeln!(paf, "chr1\t100000\t15000\t18000\t+\ttarget_chr1\t100000\t15000\t18000\t2940\t3000\t60\tNM:i:60\tcg:Z:2940=60X").unwrap();

    // Large scaffold CONTAINING the small one: query chr1 10-25kb → target chr1, 95% identity
    writeln!(paf, "chr1\t100000\t10000\t17500\t+\ttarget_chr1\t100000\t10000\t17500\t7125\t7500\t60\tNM:i:375\tcg:Z:7125=375X").unwrap();
    writeln!(paf, "chr1\t100000\t17500\t25000\t+\ttarget_chr1\t100000\t17500\t25000\t7125\t7500\t60\tNM:i:375\tcg:Z:7125=375X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("1000")
        .arg("-j")
        .arg("10000")
        .arg("-i")
        .arg("0")
        .arg("-m")
        .arg("1:1")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Check which scaffold survived
    let has_small = stdout.contains("15000\t18000");
    let has_large = stdout.contains("10000\t17500") || stdout.contains("17500\t25000");

    eprintln!("Has small contained scaffold (98%): {}", has_small);
    eprintln!("Has large containing scaffold (95%): {}", has_large);

    // The large scaffold should win (more aligned bases) even though identity is slightly lower
    assert!(has_large, "Large containing scaffold should be kept");
    assert!(
        !has_small,
        "Small contained scaffold should be filtered out"
    );
}

#[test]
fn test_scaffolds_on_different_query_chromosomes() {
    // Two scaffolds on DIFFERENT query chromosomes but same target region
    // Both should be kept (no competition on target axis in current impl)
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: query_chr1 10-20kb → target chr1 10-20kb, 95% identity
    writeln!(paf, "query_chr1\t100000\t10000\t15000\t+\ttarget_chr1\t100000\t10000\t15000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();
    writeln!(paf, "query_chr1\t100000\t15000\t20000\t+\ttarget_chr1\t100000\t15000\t20000\t4750\t5000\t60\tNM:i:250\tcg:Z:4750=250X").unwrap();

    // Scaffold 2: query_chr2 10-20kb → target chr1 10-20kb (SAME target region!), 98% identity
    writeln!(paf, "query_chr2\t100000\t10000\t15000\t+\ttarget_chr1\t100000\t10000\t15000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();
    writeln!(paf, "query_chr2\t100000\t15000\t20000\t+\ttarget_chr1\t100000\t15000\t20000\t4900\t5000\t60\tNM:i:100\tcg:Z:4900=100X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("1000")
        .arg("-j")
        .arg("10000")
        .arg("-i")
        .arg("0")
        .arg("-m")
        .arg("1:1")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Check which scaffolds survived
    let has_chr1 = stdout.contains("query_chr1");
    let has_chr2 = stdout.contains("query_chr2");

    eprintln!("Has query_chr1: {}", has_chr1);
    eprintln!("Has query_chr2: {}", has_chr2);

    // For true 1:1 filtering, we filter on BOTH query and target axes
    // Two different queries mapping to the same target region should compete
    // The better one (98% identity) should win
    assert!(
        has_chr2,
        "Better scaffold (98% from query_chr2) should be kept"
    );
    assert!(
        !has_chr1,
        "Worse scaffold (95% from query_chr1) should be filtered - target axis competition"
    );
}
