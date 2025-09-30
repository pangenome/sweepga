/// Test that plane sweep filters across target chromosomes
///
/// Bug scenario:
/// - genome1#chrA maps to genome2#chrA (full-length, high quality)
/// - genome1#chrA maps to genome2#chrB (overlapping region, lower quality)
///
/// Expected: The better mapping (chrA→chrA) should win, chrA→chrB should be filtered
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_same_query_different_targets_better_wins() {
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: genome1#chrA → genome2#chrA, 10-20kb, 98% identity, 10kb aligned
    for i in 0..10 {
        let start = 10000 + i * 1000;
        let end = start + 1000;
        writeln!(paf, "genome1#chrA\t100000\t{}\t{}\t+\tgenome2#chrA\t100000\t{}\t{}\t980\t1000\t60\tNM:i:20\tcg:Z:980=20X",
                 start, end, start, end).unwrap();
    }

    // Scaffold 2: genome1#chrA → genome2#chrB, 12-18kb (OVERLAPS!), 90% identity, 6kb aligned (WORSE)
    for i in 0..6 {
        let start = 12000 + i * 1000;
        let end = start + 1000;
        writeln!(paf, "genome1#chrA\t100000\t{}\t{}\t+\tgenome2#chrB\t100000\t{}\t{}\t900\t1000\t60\tNM:i:100\tcg:Z:900=100X",
                 start, end, start, end).unwrap();
    }

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-s")
        .arg("5000")  // 5kb minimum
        .arg("-j")
        .arg("2000")  // Merge within 2kb
        .arg("-Y")
        .arg("0")     // No identity filter
        .arg("-m")
        .arg("1:1")   // True 1:1 filtering
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    // Check which targets survived
    let has_chrA_target = stdout.lines().any(|l| l.contains("genome2#chrA"));
    let has_chrB_target = stdout.lines().any(|l| l.contains("genome2#chrB"));

    eprintln!("Has mapping to genome2#chrA (98%, 10kb): {}", has_chrA_target);
    eprintln!("Has mapping to genome2#chrB (90%, 6kb): {}", has_chrB_target);

    // The better scaffold (chrA→chrA) should be kept
    assert!(has_chrA_target, "Better scaffold (98% identity, 10kb to chrA) should be kept");

    // The worse scaffold (chrA→chrB) should be FILTERED because it overlaps on query axis
    assert!(!has_chrB_target,
            "BUG: Worse scaffold (90% identity, 6kb to chrB) should be filtered! \
             Overlapping scaffolds to different targets must compete in 1:1 mode.");
}

#[test]
fn test_non_overlapping_different_targets_both_kept() {
    let mut paf = NamedTempFile::new().unwrap();

    // Scaffold 1: genome1#chrA → genome2#chrA, 10-20kb
    for i in 0..10 {
        let start = 10000 + i * 1000;
        let end = start + 1000;
        writeln!(paf, "genome1#chrA\t100000\t{}\t{}\t+\tgenome2#chrA\t100000\t{}\t{}\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
                 start, end, start, end).unwrap();
    }

    // Scaffold 2: genome1#chrA → genome2#chrB, 50-60kb (NO OVERLAP)
    for i in 0..10 {
        let start = 50000 + i * 1000;
        let end = start + 1000;
        writeln!(paf, "genome1#chrA\t100000\t{}\t{}\t+\tgenome2#chrB\t100000\t{}\t{}\t950\t1000\t60\tNM:i:50\tcg:Z:950=50X",
                 start, end, start, end).unwrap();
    }

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-s")
        .arg("5000")
        .arg("-j")
        .arg("2000")
        .arg("-Y")
        .arg("0")
        .arg("-m")
        .arg("1:1")
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Output:\n{}", stderr);

    let has_chrA_target = stdout.lines().any(|l| l.contains("genome2#chrA"));
    let has_chrB_target = stdout.lines().any(|l| l.contains("genome2#chrB"));

    eprintln!("Has mapping to genome2#chrA: {}", has_chrA_target);
    eprintln!("Has mapping to genome2#chrB: {}", has_chrB_target);

    // Both should be kept since they don't overlap on query axis
    assert!(has_chrA_target, "Non-overlapping scaffold to chrA should be kept");
    assert!(has_chrB_target, "Non-overlapping scaffold to chrB should be kept");
}
