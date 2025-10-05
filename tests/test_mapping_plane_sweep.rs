/// Test that plane sweep for raw mappings (before scaffolding) works correctly
/// This tests the -n 1:1 flag, not -m 1:1
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_mapping_plane_sweep_across_targets() {
    // Create mappings where same query region maps to two different targets
    // With -n 1:1, only the better one should survive BEFORE scaffolding
    let mut paf = NamedTempFile::new().unwrap();

    // Mapping 1: genome1#chrA 10-20kb → genome2#chrA, 95% identity
    writeln!(paf, "genome1#chrA\t100000\t10000\t20000\t+\tgenome2#chrA\t100000\t10000\t20000\t9500\t10000\t60\tNM:i:500\tcg:Z:9500=500X").unwrap();

    // Mapping 2: genome1#chrA 12-18kb → genome2#chrB, 90% identity (overlaps mapping 1)
    writeln!(paf, "genome1#chrA\t100000\t12000\t18000\t+\tgenome2#chrB\t100000\t12000\t18000\t5400\t6000\t60\tNM:i:600\tcg:Z:5400=600X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-n")
        .arg("1:1")   // Pre-scaffold 1:1 filtering
        .arg("-j")
        .arg("0")     // No scaffolding
        .arg("-i")
        .arg("0")
        .arg("-o")
        .arg("0.5")   // Lower overlap threshold (default 0.95 is too high)
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Stderr:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    let has_chrA = stdout.lines().any(|l| l.contains("genome2#chrA"));
    let has_chrB = stdout.lines().any(|l| l.contains("genome2#chrB"));

    eprintln!("Has mapping to chrA (95%): {}", has_chrA);
    eprintln!("Has mapping to chrB (90%): {}", has_chrB);

    // The better mapping (chrA, 95%) should be kept
    assert!(has_chrA, "Better mapping (95% to chrA) should be kept");

    // The worse mapping (chrB, 90%, overlapping) should be filtered
    assert!(!has_chrB, "Worse mapping (90% to chrB) should be filtered by -n 1:1");
}

#[test]
fn test_mapping_plane_sweep_target_axis() {
    // Test target axis filtering: different queries to same target region
    let mut paf = NamedTempFile::new().unwrap();

    // Mapping 1: genome1#chrA → genome2#chrX 10-20kb, 95% identity
    writeln!(paf, "genome1#chrA\t100000\t10000\t20000\t+\tgenome2#chrX\t100000\t10000\t20000\t9500\t10000\t60\tNM:i:500\tcg:Z:9500=500X").unwrap();

    // Mapping 2: genome1#chrB → genome2#chrX 12-22kb, 98% identity, LONGER (overlaps on target)
    writeln!(paf, "genome1#chrB\t100000\t10000\t20000\t+\tgenome2#chrX\t100000\t12000\t22000\t9800\t10000\t60\tNM:i:200\tcg:Z:9800=200X").unwrap();

    paf.flush().unwrap();

    let output = std::process::Command::new("./target/release/sweepga")
        .arg(paf.path())
        .arg("-n")
        .arg("1:1")
        .arg("-j")
        .arg("0")
        .arg("-i")
        .arg("0")
        .arg("-o")
        .arg("0.5")   // Lower overlap threshold
        .output()
        .expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    eprintln!("Stderr:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);

    let has_chrA = stdout.lines().any(|l| l.contains("genome1#chrA"));
    let has_chrB = stdout.lines().any(|l| l.contains("genome1#chrB"));

    eprintln!("Has mapping from chrA (95%): {}", has_chrA);
    eprintln!("Has mapping from chrB (98%): {}", has_chrB);

    // The better mapping (chrB, 98%) should be kept
    assert!(has_chrB, "Better mapping (98% from chrB) should be kept");

    // The worse mapping (chrA, 95%) should be filtered on target axis
    assert!(!has_chrA, "Worse mapping (95% from chrA) should be filtered - target axis competition");
}
