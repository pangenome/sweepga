/// Test that scaffolds below the minimum length threshold are filtered out
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_scaffold_length_filtering() {
    // Create synthetic alignments that will form scaffolds of different lengths
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: Ten 1kb alignments (10kb aligned mass, should be kept with -s 10000)
    for i in 0..10 {
        let start = 10000 + i * 2000;
        let end = start + 1000;
        writeln!(paf, "query1\t100000\t{}\t{}\t+\ttarget\t100000\t{}\t{}\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X", start, end, start, end).unwrap();
    }

    // Scaffold 2: Five 1kb alignments (5kb aligned mass, should be dropped with -s 10000)
    for i in 0..5 {
        let start = 50000 + i * 2000;
        let end = start + 1000;
        writeln!(paf, "query2\t100000\t{}\t{}\t+\ttarget\t100000\t{}\t{}\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X", start, end, start, end).unwrap();
    }

    paf.flush().unwrap();

    // Test with -s 10000 (10kb minimum scaffold length)
    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("10000")
        .arg("-j")
        .arg("10000")  // Merge within 10kb
        .arg("-i")
        .arg("0")      // No identity filter
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output with -s 10000:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Count scaffolds by checking unique query sequences
    let query_seqs: std::collections::HashSet<&str> = stdout
        .lines()
        .filter(|l| !l.starts_with('[') && !l.is_empty())
        .map(|l| l.split('\t').next().unwrap())
        .collect();

    eprintln!("Query sequences in output: {:?}", query_seqs);

    // Should have query1 (10kb aligned mass) but not query2 (5kb aligned mass)
    assert!(query_seqs.contains("query1"), "Scaffold 1 (10kb aligned mass) should be kept");
    assert!(!query_seqs.contains("query2"), "Scaffold 2 (5kb aligned mass) should be filtered");

    // Count total output lines
    let output_count = stdout.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();

    // We should have 10 alignments in output (from scaffold 1)
    assert_eq!(output_count, 10, "Should keep only the 10 alignments from scaffold with 10kb aligned mass");
}

#[test]
fn test_scaffold_aligned_mass_filtering() {
    // Test that filtering uses aligned mass, not span
    let mut paf = NamedTempFile::new().unwrap();

    // Create a scaffold with:
    // - Span: 100kb (query 0-100000)
    // - Aligned mass: only 2kb total (2 × 1kb alignments)
    // - Large gap: 98kb between alignments
    // With -s 50000, should be FILTERED (2kb < 50kb)
    writeln!(paf, "query\t150000\t0\t1000\t+\ttarget\t150000\t0\t1000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X").unwrap();
    writeln!(paf, "query\t150000\t99000\t100000\t+\ttarget\t150000\t99000\t100000\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X").unwrap();

    paf.flush().unwrap();

    // Test with -s 50000 (50kb minimum aligned mass)
    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-s")
        .arg("50000")
        .arg("-j")
        .arg("100000")  // Merge within 100kb
        .arg("-i")
        .arg("0")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output with -s 50000:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    let output_count = stdout.lines().filter(|l| !l.starts_with('[') && !l.is_empty()).count();

    // Should filter out because aligned mass is only 2kb (< 50kb)
    // Even though span is 100kb
    assert_eq!(output_count, 0,
               "Scaffold with 100kb span but only 2kb aligned should be filtered with -s 50000");
}
