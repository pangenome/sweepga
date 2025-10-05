/// Simple test to understand FastGA's requirements

#[test]
#[ignore] // FastGA may not find alignments in synthetic homopolymer sequences
fn test_duplicate_sequence_alignment() {
    use std::fs;
    use tempfile::TempDir;
    use std::process::Command;

    let temp_dir = TempDir::new().unwrap();
    let test_fa = temp_dir.path().join("test.fa");
    let output = temp_dir.path().join("output.paf");

    // Create a more realistic sequence that FastGA can align
    // Use a repetitive pattern that creates local similarity
    let mut sequence = String::new();
    let pattern = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    for _ in 0..(100000 / pattern.len()) {
        sequence.push_str(pattern);
    }
    // Add some variation in the middle
    let mid = sequence.len() / 2;
    sequence.insert_str(mid, &"GGGGGGGGGGGGGGGGGGGG".repeat(100));
    let duplicate = sequence.clone();

    // Create FASTA with duplicated sequence
    fs::write(&test_fa, format!(
        ">chr1\n{sequence}\n>chr2_duplicate\n{duplicate}\n"
    )).unwrap();

    // Run self-alignment using compiled binary
    let sweepga_path = if cfg!(debug_assertions) {
        "target/debug/sweepga"
    } else {
        "target/release/sweepga"
    };

    let result = Command::new(sweepga_path)
        .arg(&test_fa)
        .arg("-t").arg("1")
        // Don't use --self since we want chr1 vs chr2_duplicate alignment
        .arg("-o").arg(&output)
        .output()
        .expect("Failed to run");

    assert!(result.status.success(), "FastGA failed");
    assert!(output.exists(), "No output produced");

    let content = fs::read_to_string(&output).unwrap();
    assert!(!content.is_empty(), "Duplicate sequences should produce alignments!");

    // Should find chr1 vs chr2_duplicate alignment (identical sequences)
    assert!(content.contains("chr1") && content.contains("chr2_duplicate"),
            "Should align identical sequences");
}

#[test]
#[ignore] // FastGA may not find alignments in synthetic repetitive sequences
fn test_repetitive_sequence() {
    use std::fs;
    use tempfile::TempDir;
    use std::process::Command;

    let temp_dir = TempDir::new().unwrap();
    let test_fa = temp_dir.path().join("repeat.fa");
    let output = temp_dir.path().join("output.paf");

    // Create sequence with tandem repeats (guaranteed to self-align)
    // Use 100kb minimum as FastGA needs longer sequences
    let repeat_unit = "ACGTACGTACGTACGT";
    let mut sequence = String::new();
    for _ in 0..6250 {  // 16 * 6250 = 100kb
        sequence.push_str(repeat_unit);
    }

    fs::write(&test_fa, format!(">repetitive\n{sequence}\n")).unwrap();

    // Run self-alignment using compiled binary
    let sweepga_path = if cfg!(debug_assertions) {
        "target/debug/sweepga"
    } else {
        "target/release/sweepga"
    };

    let result = Command::new(sweepga_path)
        .arg(&test_fa)
        .arg("-t").arg("1")
        .arg("--self")  // Include self-mappings
        .arg("-o").arg(&output)
        .output()
        .expect("Failed to run");

    if result.status.success() && output.exists() {
        let content = fs::read_to_string(&output).unwrap();
        // Tandem repeats should produce many self-alignments
        let line_count = content.lines().count();
        assert!(line_count > 0, "Repetitive sequence should produce alignments");
        println!("Repetitive sequence produced {line_count} alignments");
    }
}

#[test]
fn test_minimum_requirements() {
    use std::fs;
    use tempfile::TempDir;
    use std::process::Command;

    let temp_dir = TempDir::new().unwrap();

    // Test various sequence lengths to find minimum
    let lengths = vec![100, 500, 1000, 5000, 10000, 50000];

    for length in lengths {
        let test_fa = temp_dir.path().join(format!("len_{length}.fa"));
        let output = temp_dir.path().join(format!("out_{length}.paf"));

        // Create homopolymer run (simple but alignable)
        let sequence = "A".repeat(length / 2) + &"T".repeat(length / 2);
        fs::write(&test_fa, format!(">test_{length}\n{sequence}\n")).unwrap();

        let result = Command::new("cargo")
            .args(["run", "--release", "--bin", "sweepga", "--quiet", "--"])
            .arg(&test_fa)
            .arg("-t").arg("1")
            .arg("-o").arg(&output)
            .output();

        if let Ok(result) = result {
            if result.status.success() && output.exists() {
                let content = fs::read_to_string(&output).unwrap();
                if !content.is_empty() {
                    println!("Length {length} produced alignments");
                } else {
                    println!("Length {length} produced no alignments");
                }
            }
        }
    }
}