// Integration tests for sweepga with plane sweep
use std::process::Command;
use std::fs;
use std::io::Write;

#[test]
fn test_default_plane_sweep() {
    // Create a simple PAF file
    let paf_content = "\
query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1000\t800\t800\t60\tcg:Z:800M
query1\t1000\t150\t850\t+\ttarget1\t2000\t300\t1000\t700\t700\t60\tcg:Z:700M
query1\t1000\t200\t600\t+\ttarget1\t2000\t400\t800\t400\t400\t60\tcg:Z:400M
query2\t1500\t100\t1400\t+\ttarget2\t2500\t100\t1400\t1300\t1300\t60\tcg:Z:1300M
query2\t1500\t200\t1200\t+\ttarget2\t2500\t200\t1200\t1000\t1000\t60\tcg:Z:1000M
";

    // Write to temp file
    let temp_path = "/tmp/test_sweepga_input.paf";
    let output_path = "/tmp/test_sweepga_output.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run sweepga with default settings (plane sweep with n=0)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_path,
                "-s", "0"])  // No scaffolding, just plane sweep
        .output()
        .expect("Failed to run sweepga");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("sweepga failed");
    }

    // Read output
    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    // With n=0 (default), plane sweep should keep best mappings at each position
    // query1 has 3 mappings with overlapping query ranges:
    //   100-900 (800bp), 150-850 (700bp), 200-600 (400bp)
    // Only the longest (800bp) should be kept as it's best across all positions

    // query2 has 2 mappings with overlapping query ranges:
    //   100-1400 (1300bp), 200-1200 (1000bp)
    // Only the longest (1300bp) should be kept

    // Total expected: 2 mappings (one best per query sequence)
    assert_eq!(lines.len(), 2, "Should keep best mapping per query sequence");

    // Check that the longest mappings are kept
    assert!(result.contains("800\t800"), "Longest query1 mapping should be kept");
    assert!(result.contains("1300\t1300"), "Longest query2 mapping should be kept");

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}

#[test]
fn test_plane_sweep_with_secondaries() {
    // Create PAF with overlapping mappings
    let paf_content = "\
chr1\t10000\t1000\t2000\t+\tchr1_ref\t10000\t1000\t2000\t1000\t1000\t60\tcg:Z:1000M
chr1\t10000\t1000\t2000\t+\tchr1_ref\t10000\t3000\t4000\t1000\t1000\t60\tcg:Z:1000M
chr1\t10000\t1000\t2000\t+\tchr1_ref\t10000\t5000\t6000\t1000\t1000\t60\tcg:Z:1000M
chr1\t10000\t1000\t2000\t+\tchr1_ref\t10000\t7000\t8000\t1000\t1000\t60\tcg:Z:1000M
chr1\t10000\t1000\t2000\t+\tchr1_ref\t10000\t9000\t10000\t1000\t1000\t60\tcg:Z:1000M
";

    let temp_path = "/tmp/test_sweepga_secondaries.paf";
    let output_path = "/tmp/test_sweepga_secondaries_out.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run with n=2 (keep best + 2 secondaries = 3 total)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_path,
                "-n=2",  // Keep best + 2 secondaries
                "-s", "0"])  // No scaffolding
        .output()
        .expect("Failed to run sweepga");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("sweepga failed with n=2");
    }

    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    // All mappings have the same query range and same score
    // With identical scores, all should be kept
    assert_eq!(lines.len(), 5, "All mappings with identical scores should be kept");

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}

#[test]
fn test_plane_sweep_keep_all() {
    let paf_content = "\
read1\t5000\t500\t1500\t+\tref1\t10000\t1000\t2000\t1000\t1000\t60\tcg:Z:1000M
read1\t5000\t1000\t1800\t+\tref1\t10000\t2500\t3300\t800\t800\t60\tcg:Z:800M
read1\t5000\t2000\t2600\t+\tref1\t10000\t4000\t4600\t600\t600\t60\tcg:Z:600M
read1\t5000\t3000\t3400\t+\tref1\t10000\t5000\t5400\t400\t400\t60\tcg:Z:400M
";

    let temp_path = "/tmp/test_sweepga_all.paf";
    let output_path = "/tmp/test_sweepga_all_out.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run with n=-1 (keep all non-overlapping)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_path,
                "-n=-1",  // Keep all non-overlapping
                "-s", "0"])  // No scaffolding
        .output()
        .expect("Failed to run sweepga");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("sweepga failed with n=-1");
    }

    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    // All mappings should be kept as they're all best at different positions
    assert_eq!(lines.len(), 4, "All non-overlapping mappings should be kept with n=-1");

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}

#[test]
fn test_plane_sweep_with_overlap_filtering() {
    // Create mappings with high overlap
    let paf_content = "\
contig1\t8000\t1000\t3000\t+\tref1\t10000\t2000\t4000\t2000\t2000\t60\tcg:Z:2000M
contig1\t8000\t1100\t2900\t+\tref1\t10000\t5000\t6800\t1800\t1800\t60\tcg:Z:1800M
contig1\t8000\t1200\t2800\t+\tref1\t10000\t7000\t8600\t1600\t1600\t60\tcg:Z:1600M
contig1\t8000\t4000\t5000\t+\tref1\t10000\t4000\t5000\t1000\t1000\t60\tcg:Z:1000M
";

    let temp_path = "/tmp/test_sweepga_overlap.paf";
    let output_path = "/tmp/test_sweepga_overlap_out.paf";
    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run with overlap threshold
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_path,
                "-n=0",  // Keep only best
                "-p", "0.5",  // 50% overlap threshold
                "-s", "0"])  // No scaffolding
        .output()
        .expect("Failed to run sweepga");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("sweepga failed with overlap filtering");
    }

    let result = fs::read_to_string(output_path).expect("Failed to read output");
    let lines: Vec<&str> = result.lines().collect();

    // The first mapping (longest) should be kept
    // The second and third have high overlap and lower scores, might be filtered
    // The fourth has no overlap, should be kept
    assert!(lines.len() >= 2, "Should keep non-overlapping mappings");

    // Check that the non-overlapping mapping is kept
    let has_non_overlapping = result.contains("4000\t5000");
    assert!(has_non_overlapping, "Non-overlapping mapping should be kept");

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_path).ok();
}

#[test]
fn test_compare_with_no_filter() {
    // Create test PAF
    let paf_content = "\
seq1\t5000\t100\t1100\t+\tref1\t10000\t200\t1200\t1000\t1000\t60\tcg:Z:1000M
seq1\t5000\t500\t1300\t+\tref1\t10000\t1500\t2300\t800\t800\t60\tcg:Z:800M
seq1\t5000\t900\t1500\t+\tref1\t10000\t3000\t3600\t600\t600\t60\tcg:Z:600M
seq2\t4000\t200\t1200\t+\tref2\t8000\t300\t1300\t1000\t1000\t60\tcg:Z:1000M
";

    let temp_path = "/tmp/test_sweepga_compare.paf";
    let output_filtered = "/tmp/test_sweepga_filtered.paf";
    let output_unfiltered = "/tmp/test_sweepga_unfiltered.paf";

    fs::write(temp_path, paf_content).expect("Failed to write test PAF");

    // Run with plane sweep filtering (default)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_filtered,
                "-s", "0"])  // No scaffolding, just plane sweep
        .output()
        .expect("Failed to run sweepga with filtering");

    assert!(output.status.success(), "Filtered run failed");

    // Run without filtering
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "sweepga-old", "--",
                "-i", temp_path,
                "-o", output_unfiltered,
                "-f"])  // No filtering at all
        .output()
        .expect("Failed to run sweepga without filtering");

    assert!(output.status.success(), "Unfiltered run failed");

    let filtered = fs::read_to_string(output_filtered).expect("Failed to read filtered output");
    let unfiltered = fs::read_to_string(output_unfiltered).expect("Failed to read unfiltered output");

    let filtered_lines = filtered.lines().count();
    let unfiltered_lines = unfiltered.lines().count();

    // Plane sweep should keep all mappings (they're best at different positions)
    assert_eq!(filtered_lines, 4, "Plane sweep should keep all non-overlapping mappings");
    assert_eq!(unfiltered_lines, 4, "Unfiltered should have all mappings");

    // Clean up
    fs::remove_file(temp_path).ok();
    fs::remove_file(output_filtered).ok();
    fs::remove_file(output_unfiltered).ok();
}