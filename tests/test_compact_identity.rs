use std::fs::{self, File};
use std::io::{Write, BufReader, BufRead};
use std::process::Command;

#[test]
fn test_unfiltered_preserves_records() {
    // Create a small test PAF file
    let test_paf = "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1200\t800\t1000\t60\n\
                    query2\t500\t50\t450\t-\ttarget2\t600\t100\t500\t350\t400\t60\n\
                    query3\t800\t0\t800\t+\ttarget3\t800\t0\t800\t750\t800\t60\n";

    fs::write("test_input.paf", test_paf).expect("Failed to write test input");

    // Run sweepga in debug mode
    let output = Command::new("./target/release/sweepga")
        .args(&["-i", "test_input.paf", "-o", "test_output_debug.paf", "--debug"])
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success(), "sweepga failed: {:?}", String::from_utf8_lossy(&output.stderr));

    // Read output
    let output_content = fs::read_to_string("test_output_debug.paf")
        .expect("Failed to read output");

    // Check that all lines are preserved
    assert_eq!(test_paf, output_content, "Output doesn't match input in debug mode");

    // Cleanup
    fs::remove_file("test_input.paf").ok();
    fs::remove_file("test_output_debug.paf").ok();
}

#[test]
fn test_no_filter_mode_preserves_records() {
    // Create test PAF with varying block lengths
    let test_paf = "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1200\t800\t1000\t60\n\
                    query2\t500\t50\t450\t-\ttarget2\t600\t100\t500\t350\t400\t60\n";

    fs::write("test_input_n.paf", test_paf).expect("Failed to write test input");

    // Run with -m N (no filter)
    let output = Command::new("./target/release/sweepga")
        .args(&["-i", "test_input_n.paf", "-o", "test_output_n.paf", "-m", "N"])
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success());

    // Read output
    let output_content = fs::read_to_string("test_output_n.paf")
        .expect("Failed to read output");

    // Verify line count
    let input_lines: Vec<&str> = test_paf.lines().collect();
    let output_lines: Vec<&str> = output_content.lines().collect();

    assert_eq!(input_lines.len(), output_lines.len(),
               "Line count mismatch: input {} vs output {}",
               input_lines.len(), output_lines.len());

    // Cleanup
    fs::remove_file("test_input_n.paf").ok();
    fs::remove_file("test_output_n.paf").ok();
}

#[test]
fn test_min_block_length_filtering() {
    // Create test PAF with different block lengths
    let test_paf = "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1200\t800\t1000\t60\n\
                    query2\t500\t50\t450\t-\ttarget2\t600\t100\t500\t350\t50\t60\n\
                    query3\t800\t0\t800\t+\ttarget3\t800\t0\t800\t750\t800\t60\n";

    fs::write("test_block_filter.paf", test_paf).expect("Failed to write test input");

    // Run with min block length of 100
    let output = Command::new("./target/release/sweepga")
        .args(&["-i", "test_block_filter.paf", "-o", "test_block_output.paf",
                "--debug", "-b", "100"])
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success());

    // Read output and count lines
    let output_content = fs::read_to_string("test_block_output.paf")
        .expect("Failed to read output");

    let output_lines: Vec<&str> = output_content.lines().collect();

    // Should only have 2 lines (block lengths 1000 and 800, not 50)
    assert_eq!(output_lines.len(), 2,
               "Expected 2 records with block length >= 100, got {}",
               output_lines.len());

    // Verify the remaining records have correct block lengths
    for line in output_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 11 {
            let block_len: u32 = fields[10].parse().unwrap_or(0);
            assert!(block_len >= 100, "Found record with block length {} < 100", block_len);
        }
    }

    // Cleanup
    fs::remove_file("test_block_filter.paf").ok();
    fs::remove_file("test_block_output.paf").ok();
}

#[test]
fn test_file_offset_accuracy() {
    // Create a PAF file with known content
    let lines = vec![
        "q1\t100\t10\t90\t+\tt1\t200\t20\t180\t80\t160\t60\n",
        "q2\t200\t20\t180\t-\tt2\t300\t30\t270\t160\t240\t60\n",
        "q3\t300\t30\t270\t+\tt3\t400\t40\t360\t240\t320\t60\n",
    ];
    let test_paf = lines.join("");

    fs::write("test_offsets.paf", &test_paf).expect("Failed to write test input");

    // Run in debug mode
    let output = Command::new("./target/release/sweepga")
        .args(&["-i", "test_offsets.paf", "-o", "test_offsets_out.paf", "--debug"])
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success());

    // Read output
    let output_content = fs::read_to_string("test_offsets_out.paf")
        .expect("Failed to read output");

    // Verify each line is preserved correctly
    let output_lines: Vec<&str> = output_content.lines().collect();
    assert_eq!(lines.len(), output_lines.len(), "Line count mismatch");

    for (i, expected) in lines.iter().enumerate() {
        let expected_trimmed = expected.trim();
        assert_eq!(expected_trimmed, output_lines[i],
                  "Line {} mismatch", i + 1);
    }

    // Cleanup
    fs::remove_file("test_offsets.paf").ok();
    fs::remove_file("test_offsets_out.paf").ok();
}

#[test]
fn test_large_file_streaming() {
    // Generate a larger test file
    let mut large_paf = String::new();
    for i in 0..1000 {
        large_paf.push_str(&format!(
            "query{}\t1000\t100\t900\t+\ttarget{}\t2000\t200\t1200\t800\t1000\t60\n",
            i, i
        ));
    }

    fs::write("test_large.paf", &large_paf).expect("Failed to write test input");

    // Run in debug mode
    let output = Command::new("./target/release/sweepga")
        .args(&["-i", "test_large.paf", "-o", "test_large_out.paf", "--debug"])
        .output()
        .expect("Failed to run sweepga");

    assert!(output.status.success());

    // Count lines in input and output
    let input_lines = large_paf.lines().count();
    let output_content = fs::read_to_string("test_large_out.paf")
        .expect("Failed to read output");
    let output_lines = output_content.lines().count();

    assert_eq!(input_lines, output_lines,
               "Line count mismatch for large file: {} vs {}",
               input_lines, output_lines);

    // Cleanup
    fs::remove_file("test_large.paf").ok();
    fs::remove_file("test_large_out.paf").ok();
}