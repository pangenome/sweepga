/// Golden file regression tests using SHA256 checksums
///
/// These tests lock down current behavior to prevent unintended changes.
/// ANY change to output format/coordinates/filtering will fail these tests.
use anyhow::Result;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Compute SHA256 checksum of a file
fn sha256sum(path: &Path) -> Result<String> {
    let output = Command::new("sha256sum").arg(path).output()?;

    let stdout = String::from_utf8(output.stdout)?;
    let checksum = stdout
        .split_whitespace()
        .next()
        .ok_or_else(|| anyhow::anyhow!("Failed to parse sha256sum output"))?;

    Ok(checksum.to_string())
}

/// Load expected checksums from golden_data/checksums.txt
fn load_golden_checksums() -> Result<std::collections::HashMap<String, String>> {
    let checksums_path = Path::new("tests/golden_data/checksums.txt");

    if !checksums_path.exists() {
        anyhow::bail!(
            "Golden checksums file not found: {:?}\n\
             Run: cd tests/golden_data && ./generate_golden.sh",
            checksums_path
        );
    }

    let content = fs::read_to_string(checksums_path)?;
    let mut checksums = std::collections::HashMap::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            checksums.insert(parts[1].to_string(), parts[0].to_string());
        }
    }

    Ok(checksums)
}

#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
#[ignore = ".1aln output is non-deterministic due to provenance containing temp file paths"]
fn test_golden_1aln_output() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_output.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_output.1aln checksum not found"))?;

    // Generate fresh output
    let temp_dir = TempDir::new()?;
    let output = temp_dir.path().join("test_output.1aln");

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
        ])
        .stdout(fs::File::create(&output)?)
        .status()?;

    let actual = sha256sum(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_output.1aln\n\
         Expected: {}\n\
         Got:      {}\n\
         \n\
         This means the .1aln output format has changed.\n\
         \n\
         If this change was INTENTIONAL:\n\
         1. Review the diff carefully\n\
         2. Run: cd tests/golden_data && ./generate_golden.sh\n\
         3. Commit with explanation in CHANGELOG.md\n\
         \n\
         If this change was UNINTENTIONAL:\n\
         1. You have introduced a regression\n\
         2. Revert your changes and fix the bug\n\
         ",
        expected, actual
    );

    eprintln!("✓ golden_output.1aln checksum matches");
    Ok(())
}

#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
fn test_golden_paf_output() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_output.paf")
        .ok_or_else(|| anyhow::anyhow!("golden_output.paf checksum not found"))?;

    let temp_dir = TempDir::new()?;
    let output = temp_dir.path().join("test_output.paf");

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
            "--paf",
        ])
        .stdout(fs::File::create(&output)?)
        .status()?;

    let actual = sha256sum(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_output.paf\n\
         Expected: {}\n\
         Got:      {}\n\
         \n\
         The PAF output format has changed. See above for next steps.\n\
         ",
        expected, actual
    );

    eprintln!("✓ golden_output.paf checksum matches");
    Ok(())
}

#[test]
#[cfg_attr(
    target_os = "macos",
    ignore = "FastGA ARM64 compatibility issue on macOS CI"
)]
#[ignore = ".1aln output is non-deterministic due to provenance containing temp file paths"]
fn test_golden_filtered_1to1() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_filtered_1to1.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_filtered_1to1.1aln checksum not found"))?;

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // First generate unfiltered
    let unfiltered = temp_dir.path().join("unfiltered.1aln");
    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
        ])
        .stdout(fs::File::create(&unfiltered)?)
        .status()?;

    // Then filter with 1:1
    let output = temp_dir.path().join("filtered.1aln");
    Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            unfiltered.to_str().unwrap(),
            "-n",
            "1:1",
        ])
        .stdout(fs::File::create(&output)?)
        .status()?;

    let actual = sha256sum(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_filtered_1to1.1aln\n\
         Expected: {}\n\
         Got:      {}\n\
         \n\
         The 1:1 filtering behavior has changed.\n\
         ",
        expected, actual
    );

    eprintln!("✓ golden_filtered_1to1.1aln checksum matches");
    Ok(())
}

/// Verify all required checksums are present
#[test]
fn test_golden_files_complete() -> Result<()> {
    let checksums = load_golden_checksums()?;

    let required_files = [
        "golden_output.1aln",
        "golden_output.paf",
        "golden_filtered_1to1.1aln",
    ];

    for filename in &required_files {
        assert!(
            checksums.contains_key(*filename),
            "Checksum missing for: {}\nRun: cd tests/golden_data && ./generate_golden.sh",
            filename
        );
    }

    eprintln!("✓ All required checksums present in checksums.txt");
    Ok(())
}
