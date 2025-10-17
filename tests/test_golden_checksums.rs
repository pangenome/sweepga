/// Golden file regression tests using SHA256 checksums
///
/// These tests lock down current behavior to prevent unintended changes.
/// ANY change to output format/coordinates/filtering will fail these tests.
///
/// IMPORTANT: These tests must run sequentially (not in parallel) because they
/// all use the same fixed temp directory path `/tmp/sweepga_golden_gen`.
/// The temp directory path affects .1aln binary format, so must match exactly.
///
/// Run with: `cargo test test_golden -- --test-threads=1`
use anyhow::{Context, Result};
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

/// Compute SHA256 checksum of a .1aln file with normalized output
/// Uses ONEview to convert to text, filters out provenance records ('!' lines),
/// and sorts the output to handle non-deterministic record ordering from FastGA multithreading
fn sha256sum_1aln_normalized(path: &Path) -> Result<String> {
    // Run: ONEview file.1aln | grep -v '^!' | sort | sha256sum
    // This pipeline:
    // 1. Converts binary .1aln to text
    // 2. Removes non-deterministic provenance lines
    // 3. Sorts to canonicalize record order (FastGA multithreading creates non-deterministic order)
    // 4. Computes checksum of normalized output

    let mut oneview = Command::new("ONEview")
        .arg(path)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn ONEview")?;

    let mut grep = Command::new("grep")
        .arg("-v")
        .arg("^!")
        .stdin(
            oneview
                .stdout
                .take()
                .ok_or_else(|| anyhow::anyhow!("Failed to pipe ONEview"))?,
        )
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn grep")?;

    let mut sort = Command::new("sort")
        .stdin(
            grep.stdout
                .take()
                .ok_or_else(|| anyhow::anyhow!("Failed to pipe grep"))?,
        )
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn sort")?;

    let sha256sum = Command::new("sha256sum")
        .stdin(
            sort.stdout
                .take()
                .ok_or_else(|| anyhow::anyhow!("Failed to pipe sort"))?,
        )
        .output()
        .context("Failed to run sha256sum")?;

    // Wait for all processes
    oneview.wait()?;
    grep.wait()?;
    sort.wait()?;

    let stdout = String::from_utf8(sha256sum.stdout)?;
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
fn test_golden_1aln_output() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_output.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_output.1aln checksum not found"))?;

    // Use EXACT same temp directory as generate_golden.sh (ensures deterministic output)
    // The temp directory PATH itself affects .1aln binary format, so must match exactly
    let temp_dir_path = Path::new("/tmp/sweepga_golden_gen");
    if temp_dir_path.exists() {
        fs::remove_dir_all(temp_dir_path)?;
    }
    fs::create_dir_all(temp_dir_path)?;
    let output = temp_dir_path.join("test_output.1aln");

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input_gz = temp_dir_path.join("test_input.fa.gz");
    fs::copy(input, &temp_input_gz)?;

    // Gunzip to match golden file generation (uncompressed FASTA)
    // Use gunzip (no -c flag) to gunzip in-place and remove .gz file
    // This is important because FastGA embeds the filename in .1aln format
    let temp_input = temp_dir_path.join("test_input.fa");
    Command::new("gunzip").arg(&temp_input_gz).status()?;

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

    // Use normalized checksum that filters out non-deterministic provenance records
    let actual = sha256sum_1aln_normalized(&output)?;

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

    // Use EXACT same temp directory as generate_golden.sh
    let temp_dir_path = Path::new("/tmp/sweepga_golden_gen");
    if temp_dir_path.exists() {
        fs::remove_dir_all(temp_dir_path)?;
    }
    fs::create_dir_all(temp_dir_path)?;
    let output = temp_dir_path.join("test_output.paf");

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input_gz = temp_dir_path.join("test_input.fa.gz");
    fs::copy(input, &temp_input_gz)?;

    // Gunzip to match golden file generation (uncompressed FASTA)
    // Use gunzip (no -c flag) to gunzip in-place and remove .gz file
    // This is important because FastGA embeds the filename in .1aln format
    let temp_input = temp_dir_path.join("test_input.fa");
    Command::new("gunzip").arg(&temp_input_gz).status()?;

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
fn test_golden_filtered_1to1() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_filtered_1to1.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_filtered_1to1.1aln checksum not found"))?;

    // Use EXACT same temp directory as generate_golden.sh (ensures deterministic output)
    let temp_dir_path = Path::new("/tmp/sweepga_golden_gen");
    if temp_dir_path.exists() {
        fs::remove_dir_all(temp_dir_path)?;
    }
    fs::create_dir_all(temp_dir_path)?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );
    let temp_input_gz = temp_dir_path.join("test_input.fa.gz");
    fs::copy(input, &temp_input_gz)?;

    // Gunzip to match golden file generation (uncompressed FASTA)
    // Use gunzip (no -c flag) to gunzip in-place and remove .gz file
    // This is important because FastGA embeds the filename in .1aln format
    let temp_input = temp_dir_path.join("test_input.fa");
    Command::new("gunzip").arg(&temp_input_gz).status()?;

    // First generate unfiltered
    let unfiltered = temp_dir_path.join("unfiltered.1aln");
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
    let output = temp_dir_path.join("filtered.1aln");
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

    // Use normalized checksum that filters out non-deterministic provenance records
    let actual = sha256sum_1aln_normalized(&output)?;

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
