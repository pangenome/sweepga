// Test to demonstrate the plane sweep grouping bug
use std::fs;
use std::process::Command;

#[test]
fn test_plane_sweep_grouping_bug() {
    // Create a PAF file with mappings from one query to multiple targets
    // The plane sweep should consider ALL mappings for each query sequence together
    let paf_content = "\
chrI_query\t10000\t1000\t2000\t+\tchrI_target1\t10000\t1000\t2000\t1000\t1000\t60\tcg:Z:1000M
chrI_query\t10000\t1000\t2000\t+\tchrII_target2\t10000\t2000\t3000\t1000\t1000\t60\tcg:Z:1000M
chrI_query\t10000\t1000\t2000\t+\tchrIII_target3\t10000\t3000\t4000\t1000\t1000\t60\tcg:Z:1000M
chrII_query\t15000\t2000\t3000\t+\tchrI_target1\t10000\t2000\t3000\t1000\t1000\t60\tcg:Z:1000M
chrII_query\t15000\t2000\t3000\t+\tchrII_target2\t10000\t4000\t5000\t1000\t1000\t60\tcg:Z:1000M
";

    // All mappings from chrI_query overlap on the query (1000-2000, 1500-2500, 1800-2800)
    // With correct plane sweep grouping by query only:
    // - chrI_query should keep only the best mapping (or all if same score)
    // With incorrect grouping by query-target pairs:
    // - All mappings would be kept (each is alone in its group)

    let temp_path = "/tmp/test_grouping_bug.paf";
    let output_path = "/tmp/test_grouping_bug_out.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run sweepga with plane sweep (n=0, keep only best)
    let output = Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--bin",
            "sweepga",
            "--",
            "-i",
            temp_path,
            "-o",
            output_path,
            "-n=1", // Keep only best per position
            "-s",
            "0",
        ]) // No scaffolding
        .output()
        .expect("Failed to run sweepga");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("sweepga failed");
    }

    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    // Count mappings for chrI_query
    let chr1_mappings = lines
        .iter()
        .filter(|line| line.starts_with("chrI_query"))
        .count();

    // Count mappings for chrII_query
    let chr2_mappings = lines
        .iter()
        .filter(|line| line.starts_with("chrII_query"))
        .count();

    eprintln!("chrI_query mappings: {}", chr1_mappings);
    eprintln!("chrII_query mappings: {}", chr2_mappings);
    eprintln!("Output:\n{}", result);

    // With CORRECT grouping by query only:
    // chrI_query has 3 overlapping mappings, should keep 1 (best)
    // chrII_query has 2 overlapping mappings, should keep 1 (best)

    // With INCORRECT grouping by query-target pairs:
    // Each mapping is alone in its group, so all 5 would be kept

    // This test will FAIL with the current bug (keeps all 5)
    // and PASS when fixed (keeps 2, one per query)
    assert_eq!(
        chr1_mappings, 1,
        "Should keep only best mapping for chrI_query"
    );
    assert_eq!(
        chr2_mappings, 1,
        "Should keep only best mapping for chrII_query"
    );
    assert_eq!(
        lines.len(),
        2,
        "Should keep exactly 2 mappings total (one per query)"
    );

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}

#[test]
fn test_multi_target_filtering() {
    // Another test: one query mapped to same position on multiple targets
    // These should compete in plane sweep
    let paf_content = "\
query1\t5000\t1000\t2000\t+\ttarget_A\t10000\t3000\t4000\t1000\t1000\t60\tcg:Z:1000M
query1\t5000\t1000\t2000\t+\ttarget_B\t10000\t5000\t6000\t1000\t1000\t60\tcg:Z:1000M
query1\t5000\t1000\t2000\t+\ttarget_C\t10000\t7000\t8000\t1000\t1000\t60\tcg:Z:1000M
query1\t5000\t1000\t2000\t+\ttarget_D\t10000\t1000\t2000\t1000\t1000\t60\tcg:Z:1000M
";

    // All 4 mappings have identical query coordinates (1000-2000)
    // They should compete in plane sweep, keeping only the best

    let temp_path = "/tmp/test_multi_target.paf";
    let output_path = "/tmp/test_multi_target_out.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    let output = Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--bin",
            "sweepga",
            "--",
            "-i",
            temp_path,
            "-o",
            output_path,
            "-n=1", // Keep only best
            "-s",
            "0",
        ]) // No scaffolding
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success(), "sweepga should succeed");

    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    eprintln!(
        "Multi-target test output ({} lines):\n{}",
        lines.len(),
        result
    );

    // With correct plane sweep: should keep only 1 mapping (all have same score, tie-breaking)
    // With incorrect grouping: would keep all 4 (each in different group)

    // Currently this will FAIL (keeps all 4) due to the bug
    assert_eq!(
        lines.len(),
        1,
        "Should keep only 1 mapping when all map to same query position"
    );

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}
