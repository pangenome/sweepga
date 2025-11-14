#![allow(clippy::uninlined_format_args)]
/// Test that plane sweep respects genome pair grouping
///
/// When we have alignments like:
/// - GenomeA#1#chr1 -> GenomeB#1#chr1 (score 100)
/// - GenomeA#1#chr1 -> GenomeC#1#chr1 (score 90)
///
/// With 1:1 filtering, we should keep BOTH because they are different genome pairs
/// (GenomeA->GenomeB and GenomeA->GenomeC are separate pairs)
///
/// Currently BROKEN: we only keep the best one globally (score 100), losing the GenomeA->GenomeC pair

#[test]
fn test_plane_sweep_preserves_genome_pairs() {
    use std::fs;
    use tempfile::TempDir;

    let temp_dir = TempDir::new().unwrap();
    let input_paf = temp_dir.path().join("input.paf");
    let _output_paf = temp_dir.path().join("output.paf");

    // Create test PAF with 3 genomes: A, B, C
    // All alignments from A#1#chr1 to different targets
    let paf_content = "\
A#1#chr1\t1000\t0\t500\t+\tB#1#chr1\t1000\t0\t500\t450\t500\t60\tcg:Z:500M
A#1#chr1\t1000\t0\t500\t+\tC#1#chr1\t1000\t0\t500\t400\t500\t60\tcg:Z:500M
A#1#chr1\t1000\t0\t500\t+\tD#1#chr1\t1000\t0\t500\t350\t500\t60\tcg:Z:500M
";

    fs::write(&input_paf, paf_content).unwrap();

    // Run sweepga with 1:1 filtering (disable scaffolding for small test data)
    let result = std::process::Command::new("cargo")
        .args(["run", "--release", "--bin", "sweepga", "--quiet", "--"])
        .arg("-j")
        .arg("0") // Disable scaffolding to avoid filtering small test alignments
        .arg(&input_paf)
        .output()
        .expect("Failed to run sweepga");

    assert!(result.status.success(), "sweepga failed: {:?}", result);

    let output = String::from_utf8_lossy(&result.stdout);
    let lines: Vec<&str> = output.lines().collect();

    // With proper genome pair grouping, we should keep ALL 3 alignments
    // because they are from different genome pairs: A->B, A->C, A->D
    assert_eq!(
        lines.len(),
        3,
        "Expected 3 alignments (one per genome pair), got {}:\n{}",
        lines.len(),
        output
    );

    // Verify each genome pair is represented
    assert!(output.contains("B#1#chr1"), "Missing A->B genome pair");
    assert!(output.contains("C#1#chr1"), "Missing A->C genome pair");
    assert!(output.contains("D#1#chr1"), "Missing A->D genome pair");
}

#[test]
fn test_plane_sweep_within_genome_pair() {
    use std::fs;
    use tempfile::TempDir;

    let temp_dir = TempDir::new().unwrap();
    let input_paf = temp_dir.path().join("input.paf");

    // Create test PAF with SAME genome pair (A->B) but different chromosomes
    // With 1:1 filtering using intersection semantics:
    // - Query sweep: A#1#chr1 keeps best (B#1#chr1, score 450), A#1#chr2 keeps B#1#chr1 (350)
    // - Target sweep: B#1#chr1 keeps best (A#1#chr1, score 450), B#1#chr2 keeps A#1#chr1 (400)
    // - Intersection: Only A#1#chr1->B#1#chr1 passes both sweeps
    let paf_content = "\
A#1#chr1\t1000\t0\t500\t+\tB#1#chr1\t1000\t0\t500\t450\t500\t60\tcg:Z:500M
A#1#chr1\t1000\t0\t500\t+\tB#1#chr2\t1000\t0\t500\t400\t500\t60\tcg:Z:500M
A#1#chr2\t1000\t0\t500\t+\tB#1#chr1\t1000\t0\t500\t350\t500\t60\tcg:Z:500M
";

    fs::write(&input_paf, paf_content).unwrap();

    // Run sweepga with 1:1 filtering (disable scaffolding for small test data)
    let result = std::process::Command::new("cargo")
        .args(["run", "--release", "--bin", "sweepga", "--quiet", "--"])
        .arg("-j")
        .arg("0") // Disable scaffolding to avoid filtering small test alignments
        .arg(&input_paf)
        .output()
        .expect("Failed to run sweepga");

    assert!(result.status.success(), "sweepga failed");

    let output = String::from_utf8_lossy(&result.stdout);
    let lines: Vec<&str> = output.lines().collect();

    // Within the same genome pair (A->B), 1:1 filtering applies intersection logic:
    // Only keep alignments that pass BOTH query and target sweeps
    // Result: 1 alignment (A#1#chr1->B#1#chr1 is the only one that's best on both axes)
    assert_eq!(
        lines.len(),
        1,
        "Expected 1 alignment (best on both query and target axes), got {}:\n{}",
        lines.len(),
        output
    );

    // Should keep the alignment that's best on both axes
    assert!(
        output.contains("A#1#chr1\t1000\t0\t500\t+\tB#1#chr1"),
        "Missing A#1#chr1->B#1#chr1"
    );
}

#[test]
#[ignore] // Requires z.paf file which doesn't exist in CI
fn test_yeast_genome_pairs_preserved() {
    use std::collections::HashSet;
    use std::fs;

    // Read the unfiltered z.paf
    let input_content = fs::read_to_string("z.paf").expect("z.paf not found");

    // Extract genome pairs from input
    let mut input_pairs: HashSet<(String, String)> = HashSet::new();
    for line in input_content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 6 {
            let q_parts: Vec<&str> = fields[0].split('#').collect();
            let t_parts: Vec<&str> = fields[5].split('#').collect();
            if q_parts.len() >= 2 && t_parts.len() >= 2 {
                let q_genome = q_parts[0].to_string();
                let t_genome = t_parts[0].to_string();
                input_pairs.insert((q_genome, t_genome));
            }
        }
    }

    // Run sweepga with 1:1 filtering
    let result = std::process::Command::new("cargo")
        .args(["run", "--release", "--bin", "sweepga", "--quiet", "--"])
        .arg("-j")
        .arg("0") // Disable scaffolding for consistent comparison
        .arg("z.paf")
        .output()
        .expect("Failed to run sweepga");

    assert!(result.status.success(), "sweepga failed");

    let output = String::from_utf8_lossy(&result.stdout);

    // Extract genome pairs from output
    let mut output_pairs: HashSet<(String, String)> = HashSet::new();
    for line in output.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 6 {
            let q_parts: Vec<&str> = fields[0].split('#').collect();
            let t_parts: Vec<&str> = fields[5].split('#').collect();
            if q_parts.len() >= 2 && t_parts.len() >= 2 {
                let q_genome = q_parts[0].to_string();
                let t_genome = t_parts[0].to_string();
                output_pairs.insert((q_genome, t_genome));
            }
        }
    }

    // Check that all genome pairs are preserved
    let missing_pairs: Vec<_> = input_pairs.difference(&output_pairs).collect();

    assert!(
        missing_pairs.is_empty(),
        "1:1 filtering lost {} genome pairs: {:?}\nInput had {} pairs, output has {} pairs",
        missing_pairs.len(),
        missing_pairs,
        input_pairs.len(),
        output_pairs.len()
    );
}
