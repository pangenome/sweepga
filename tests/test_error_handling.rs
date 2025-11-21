/// Error handling tests for malformed and invalid inputs
///
/// Tests that the program fails gracefully with clear error messages
/// when given invalid or malformed input files.
use anyhow::Result;
use std::fs;
use std::process::Command;
use tempfile::TempDir;

/// Test handling of empty input file
#[test]
fn test_empty_file_error() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let empty_file = temp_dir.path().join("empty.paf");

    // Create empty file
    fs::write(&empty_file, "")?;

    // Try to process it
    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            empty_file.to_str().unwrap(),
        ])
        .output()?;

    // Should fail with error
    assert!(!output.status.success(), "Empty file should cause error");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Empty file") || stderr.contains("No alignments"),
        "Should mention empty file, got: {stderr}"
    );

    Ok(())
}

/// Test handling of malformed PAF lines
#[test]
fn test_malformed_paf_lines() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Test case 1: Too few fields (PAF requires 12+ fields)
    let malformed_paf = temp_dir.path().join("malformed.paf");
    fs::write(&malformed_paf, "seq1\t100\t200\n")?; // Only 3 fields

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            malformed_paf.to_str().unwrap(),
        ])
        .output()?;

    // Should either skip the line or fail with error
    let _stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Should not crash, and output should be empty or contain error
    assert!(
        !output.status.success() || stdout.is_empty(),
        "Malformed PAF should be rejected or produce no output"
    );

    Ok(())
}

/// Test handling of invalid numeric fields in PAF
#[test]
fn test_invalid_paf_numbers() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let bad_paf = temp_dir.path().join("bad_numbers.paf");

    // Write PAF line with invalid number (non-numeric in position field)
    fs::write(
        &bad_paf,
        "seq1\t1000\tNOT_A_NUMBER\t200\t+\tseq2\t2000\t100\t300\t150\t200\t60\n",
    )?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            bad_paf.to_str().unwrap(),
        ])
        .output()?;

    // Should fail or skip the line
    let stderr = String::from_utf8_lossy(&output.stderr);

    // Either it should fail, or produce no output
    assert!(
        !output.status.success() || output.stdout.is_empty(),
        "Invalid numbers should cause error or be skipped, stderr: {stderr}"
    );

    Ok(())
}

/// Test handling of missing input file
#[test]
fn test_missing_file_error() -> Result<()> {
    let nonexistent = "/tmp/this_file_definitely_does_not_exist_12345.paf";

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            nonexistent,
        ])
        .output()?;

    // Should fail with error
    assert!(!output.status.success(), "Missing file should cause error");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("No such file")
            || stderr.contains("not found")
            || stderr.contains("does not exist"),
        "Should mention file not found, got: {stderr}"
    );

    Ok(())
}

/// Test handling of invalid coordinate ranges
/// Note: Currently the program passes through invalid coordinates (start > end)
/// This documents the behavior but doesn't enforce validation
#[test]
fn test_invalid_coordinate_ranges() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let bad_coords = temp_dir.path().join("bad_coords.paf");

    // Write PAF with start > end (invalid, but currently accepted)
    fs::write(
        &bad_coords,
        "seq1\t1000\t500\t200\t+\tseq2\t2000\t100\t300\t150\t200\t60\n",
    )?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            bad_coords.to_str().unwrap(),
        ])
        .output()?;

    // Currently passes through without validation
    // This test documents that behavior - may want stricter validation in future
    assert!(
        output.status.success(),
        "Program should not crash on invalid coordinates (currently passes them through)"
    );

    Ok(())
}

/// Test handling of unsupported file format
#[test]
fn test_unsupported_format() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let binary_file = temp_dir.path().join("binary.bin");

    // Write some binary garbage
    fs::write(&binary_file, [0xFF, 0xFE, 0xFD, 0xFC, 0x00, 0x01])?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            binary_file.to_str().unwrap(),
        ])
        .output()?;

    // Should fail with format error
    assert!(!output.status.success(), "Binary file should cause error");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("format")
            || stderr.contains("parse")
            || stderr.contains("invalid")
            || stderr.contains("UTF")
            || stderr.contains("Empty file"),
        "Should mention format/parsing error, got: {stderr}"
    );

    Ok(())
}

/// Test handling of mixed valid/invalid lines
#[test]
fn test_partial_valid_paf() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let mixed_paf = temp_dir.path().join("mixed.paf");

    // Mix of valid and invalid lines
    let content = "\
seq1\t100\t0\t50\t+\tseq2\t200\t0\t50\t40\t50\t60
INVALID LINE WITH GARBAGE
seq3\t300\t0\t100\t+\tseq4\t400\t0\t100\t90\t100\t60
too\tfew\tfields
seq5\t500\t0\t150\t+\tseq6\t600\t0\t150\t140\t150\t60
";

    fs::write(&mixed_paf, content)?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            "-j",
            "0", // Disable scaffolding to avoid filtering small test alignments
            mixed_paf.to_str().unwrap(),
        ])
        .output()?;

    let stdout = String::from_utf8_lossy(&output.stdout);

    // Should process valid lines and skip/warn about invalid ones
    // Count output lines - should have at least some valid records
    let line_count = stdout.lines().filter(|l| !l.trim().is_empty()).count();

    // Should get at least 1 valid line (might filter some)
    assert!(
        line_count >= 1,
        "Should process at least some valid lines, got {line_count} lines"
    );

    Ok(())
}

/// Test that negative coordinates are handled
#[test]
fn test_negative_coordinates() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let neg_coords = temp_dir.path().join("negative.paf");

    // Negative start position (invalid in PAF)
    fs::write(
        &neg_coords,
        "seq1\t1000\t-100\t200\t+\tseq2\t2000\t100\t300\t150\t200\t60\n",
    )?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            neg_coords.to_str().unwrap(),
        ])
        .output()?;

    // Should fail to parse negative number or reject it
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(
        !output.status.success() || stdout.is_empty(),
        "Negative coordinates should be rejected"
    );

    Ok(())
}

/// Test handling of extremely large coordinates
#[test]
fn test_overflow_coordinates() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let overflow_paf = temp_dir.path().join("overflow.paf");

    // Coordinates that might overflow
    fs::write(
        &overflow_paf,
        "seq1\t999999999999999999999\t0\t100\t+\tseq2\t2000\t0\t100\t90\t100\t60\n",
    )?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            overflow_paf.to_str().unwrap(),
        ])
        .output()?;

    // Should fail to parse or reject
    assert!(
        !output.status.success() || output.stdout.is_empty(),
        "Overflow coordinates should be rejected"
    );

    Ok(())
}

/// Test handling of invalid strand field
#[test]
fn test_invalid_strand() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let bad_strand = temp_dir.path().join("bad_strand.paf");

    // Invalid strand (should be + or -)
    fs::write(
        &bad_strand,
        "seq1\t1000\t0\t100\tX\tseq2\t2000\t0\t100\t90\t100\t60\n",
    )?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            bad_strand.to_str().unwrap(),
        ])
        .output()?;

    // Program should not crash, regardless of whether it accepts or rejects
    // Currently passes through invalid strands
    assert!(
        output.status.success(),
        "Should handle invalid strand gracefully (currently passes through)"
    );

    Ok(())
}

/// Test permission denied error (if possible)
#[test]
#[cfg(unix)]
fn test_permission_denied() -> Result<()> {
    use std::os::unix::fs::PermissionsExt;

    let temp_dir = TempDir::new()?;
    let no_read = temp_dir.path().join("no_read.paf");

    // Create file and remove read permissions
    fs::write(
        &no_read,
        "seq1\t100\t0\t50\t+\tseq2\t200\t0\t50\t40\t50\t60\n",
    )?;
    fs::set_permissions(&no_read, fs::Permissions::from_mode(0o000))?;

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            no_read.to_str().unwrap(),
        ])
        .output()?;

    // Restore permissions for cleanup
    let _ = fs::set_permissions(&no_read, fs::Permissions::from_mode(0o644));

    // Should fail with permission error
    assert!(
        !output.status.success(),
        "Permission denied should cause error"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Permission denied")
            || stderr.contains("permission")
            || stderr.contains("access"),
        "Should mention permission error, got: {stderr}"
    );

    Ok(())
}
