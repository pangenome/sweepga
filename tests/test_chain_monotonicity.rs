/// Test that chains are monotonic: larger gap thresholds should find supersets of smaller gaps
///
/// Create synthetic alignments with various gap patterns to try to reproduce
/// the coverage dropout observed with large -j values
use std::io::Write;
use tempfile::NamedTempFile;

/// Create a PAF with collinear high-identity alignments
fn create_simple_collinear_paf() -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();

    // 5 collinear high-identity (95%) alignments with varying gaps
    let alignments = vec![
        "query\t100000\t0\t1000\t+\ttarget\t100000\t0\t1000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
        "query\t100000\t2000\t3000\t+\ttarget\t100000\t2000\t3000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
        "query\t100000\t8000\t9000\t+\ttarget\t100000\t8000\t9000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
        "query\t100000\t20000\t21000\t+\ttarget\t100000\t20000\t21000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
        "query\t100000\t50000\t51000\t+\ttarget\t100000\t50000\t51000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
    ];

    for aln in alignments {
        writeln!(file, "{}", aln).unwrap();
    }
    file.flush().unwrap();
    file
}

/// Create a PAF with mixed identity alignments - high identity close together, lower identity far apart
/// This simulates what we might see in real data where distant regions have more divergence
fn create_mixed_identity_paf() -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();

    // 10 alignments: first 5 are high-identity and close, next 5 are lower-identity and far
    let alignments = vec![
        // Close together, 98% identity
        "query\t200000\t0\t1000\t+\ttarget\t200000\t0\t1000\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
        "query\t200000\t2000\t3000\t+\ttarget\t200000\t2000\t3000\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
        "query\t200000\t5000\t6000\t+\ttarget\t200000\t5000\t6000\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
        "query\t200000\t8000\t9000\t+\ttarget\t200000\t8000\t9000\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
        "query\t200000\t11000\t12000\t+\ttarget\t200000\t11000\t12000\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
        // Far away, 90% identity
        "query\t200000\t50000\t51000\t+\ttarget\t200000\t50000\t51000\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
        "query\t200000\t80000\t81000\t+\ttarget\t200000\t80000\t81000\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
        "query\t200000\t120000\t121000\t+\ttarget\t200000\t120000\t121000\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
        "query\t200000\t160000\t161000\t+\ttarget\t200000\t160000\t161000\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
        "query\t200000\t195000\t196000\t+\ttarget\t200000\t195000\t196000\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
    ];

    for aln in alignments {
        writeln!(file, "{}", aln).unwrap();
    }
    file.flush().unwrap();
    file
}

/// Create a PAF with many small fragments - tests if large gaps create bad chains
fn create_fragmented_paf() -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();

    // 20 small alignments with 95-97% identity, spaced 2kb apart
    for i in 0..20 {
        let qstart = i * 3000;
        let qend = qstart + 1000;
        let tstart = qstart;
        let tend = qend;
        let matches = 950 + (i % 3) * 10; // 95-97% identity
        let nm = 1000 - matches;

        writeln!(file, "query\t100000\t{}\t{}\t+\ttarget\t100000\t{}\t{}\t{}\t1000\t60\tNM:i:{}\tcg:Z:{}={}X",
                qstart, qend, tstart, tend, matches, nm, matches, nm).unwrap();
    }

    file.flush().unwrap();
    file
}

fn extract_chains(output: &str) -> Vec<Vec<String>> {
    use std::collections::HashMap;

    let mut chains: HashMap<String, Vec<String>> = HashMap::new();

    for line in output.lines() {
        if line.starts_with('[') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 13 {
            continue;
        }

        // Find chain ID
        let mut chain_id = None;
        for field in fields.iter().skip(12) {
            if field.starts_with("ch:Z:") {
                chain_id = Some(field[5..].to_string());
                break;
            }
        }

        if let Some(cid) = chain_id {
            // Create mapping ID from coordinates
            let mapping_id = format!("{}:{}-{}", fields[2], fields[3], fields[0]);
            chains.entry(cid).or_insert_with(Vec::new).push(mapping_id);
        }
    }

    chains.into_values().collect()
}

#[test]
fn test_simple_collinear_chaining() {
    let paf_file = create_simple_collinear_paf();
    let paf_path = paf_file.path();

    // Test with different gap values
    let gaps = vec![2_000, 10_000, 30_000, 100_000];

    for gap in &gaps {
        let output = std::process::Command::new("./target/release/sweepga")
            .arg("-i")
            .arg(paf_path)
            .arg("-j")
            .arg(gap.to_string())
            .arg("-Y")
            .arg("0.90") // 90% identity threshold
            .arg("-s")
            .arg("0") // No minimum scaffold length
            .output()
            .expect("Failed to run sweepga");

        let stdout = String::from_utf8_lossy(&output.stdout);
        let num_output = stdout.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();

        eprintln!("-j {} produced {} output mappings", gap, num_output);

        // All 5 input mappings should be in the output regardless of gap size
        // because they're all high-identity
        assert_eq!(num_output, 5,
                  "-j {} should keep all 5 mappings (95% identity)", gap);
    }
}

#[test]
fn test_mixed_identity_chaining() {
    let paf_file = create_mixed_identity_paf();
    let paf_path = paf_file.path();

    // Test with different gap and identity thresholds
    let test_cases = vec![
        (10_000, "0.95", 5),  // Gap 10k: chains 1-5 (98%) and 6-10 (90%) separate, only high-ID passes
        (100_000, "0.95", 0), // Gap 100k: all chain together (94% avg), fails 95% threshold
        (10_000, "0.85", 10), // Gap 10k: both chains pass 85% threshold
        (100_000, "0.85", 10), // Gap 100k: single chain (94%) passes 85% threshold
    ];

    for (gap, threshold, expected) in test_cases {
        let output = std::process::Command::new("./target/release/sweepga")
            .arg("-i")
            .arg(paf_path)
            .arg("-j")
            .arg(gap.to_string())
            .arg("-Y")
            .arg(threshold)
            .arg("-s")
            .arg("0")
            .output()
            .expect("Failed to run sweepga");

        let stdout = String::from_utf8_lossy(&output.stdout);
        let num_output = stdout.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();

        eprintln!("-j {}, -Y {} produced {} mappings (expected {})",
                 gap, threshold, num_output, expected);

        assert_eq!(num_output, expected,
                  "-j {}, -Y {} should keep {} mappings", gap, threshold, expected);
    }
}

#[test]
fn test_fragmented_chaining_coverage() {
    let paf_file = create_fragmented_paf();
    let paf_path = paf_file.path();

    // With 20 fragments of 1kb each (95-97% identity), we should get 20kb coverage
    // regardless of gap size
    let gaps = vec![5_000, 50_000, 500_000];

    for gap in &gaps {
        let output = std::process::Command::new("./target/release/sweepga")
            .arg("-i")
            .arg(paf_path)
            .arg("-j")
            .arg(gap.to_string())
            .arg("-Y")
            .arg("0.90")
            .arg("-s")
            .arg("0")
            .output()
            .expect("Failed to run sweepga");

        let stdout = String::from_utf8_lossy(&output.stdout);
        let num_output = stdout.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();

        eprintln!("-j {} produced {} mappings out of 20 input", gap, num_output);

        // Key test: coverage should NOT decrease with larger gap values
        // All 20 fragments have 95-97% identity, so they should all pass the 90% threshold
        assert_eq!(num_output, 20,
                  "-j {} should keep all 20 fragments", gap);
    }
}

/// Test that low-identity chains (like centromeric inversions) are correctly filtered
#[test]
fn test_centromere_inversion_filtering() {
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Create a synthetic centromere-like inversion with 76% identity
    let mut paf = NamedTempFile::new().unwrap();
    
    // Simulate multiple alignments that would chain together
    // Each has ~76% identity (like the real centromere data)
    let alignments = vec![
        "query\t200000000\t129000000\t130000000\t-\ttarget\t200000000\t132000000\t133000000\t760000\t1000000\t60\tNM:i:240000\tcg:Z:760000=240000X",
        "query\t200000000\t130000000\t131000000\t-\ttarget\t200000000\t133000000\t134000000\t760000\t1000000\t60\tNM:i:240000\tcg:Z:760000=240000X",
        "query\t200000000\t131000000\t132000000\t-\ttarget\t200000000\t134000000\t135000000\t760000\t1000000\t60\tNM:i:240000\tcg:Z:760000=240000X",
    ];
    
    for aln in alignments {
        writeln!(paf, "{}", aln).unwrap();
    }
    paf.flush().unwrap();
    
    // Test 1: With Y=0.80 (80%), chain should be filtered (76% < 80%)
    let output_80 = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-Y")
        .arg("0.80")
        .arg("-j")
        .arg("10000")
        .arg("-s")
        .arg("0")
        .output()
        .expect("Failed to run sweepga");
    
    let stdout_80 = String::from_utf8_lossy(&output_80.stdout);
    let stderr_80 = String::from_utf8_lossy(&output_80.stderr);
    eprintln!("Y=0.80 output:\n{}", stderr_80);
    
    // Should output 0 mappings (chain filtered due to 76% identity)
    let count_80 = stdout_80.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();
    assert_eq!(count_80, 0, "Chain with 76% identity should be filtered with Y=0.80");
    
    // Test 2: With Y=0.75 (75%), chain should pass (76% >= 75%)
    let output_75 = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-Y")
        .arg("0.75")
        .arg("-j")
        .arg("10000")
        .arg("-s")
        .arg("0")
        .output()
        .expect("Failed to run sweepga");
    
    let stdout_75 = String::from_utf8_lossy(&output_75.stdout);
    let stderr_75 = String::from_utf8_lossy(&output_75.stderr);
    eprintln!("Y=0.75 output:\n{}", stderr_75);
    
    // Should output 3 mappings (chain kept)
    let count_75 = stdout_75.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();
    assert!(count_75 > 0, "Chain with 76% identity should pass with Y=0.75");
    
    // Test 3: With Y=0 (no filter), chain should definitely pass
    let output_0 = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-Y")
        .arg("0")
        .arg("-j")
        .arg("10000")
        .arg("-s")
        .arg("0")
        .output()
        .expect("Failed to run sweepga");
    
    let stdout_0 = String::from_utf8_lossy(&output_0.stdout);
    assert!(stdout_0.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count() > 0,
            "Chain should pass with Y=0 (no filter)");
}
