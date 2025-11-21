/// Golden file regression tests using SHA256 checksums
///
/// These tests lock down current behavior to prevent unintended changes.
/// ANY change to output format/coordinates/filtering will fail these tests.
///
/// These tests now support parallel execution by filtering out path-specific
/// information from .1aln files (both `!` provenance lines and `<` header lines
/// containing directory paths).
///
/// Run with: `cargo test test_golden`
use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Find ONEview binary from fastga-rs build output
fn find_oneview() -> Result<PathBuf> {
    // ONEview is built by fastga-rs and placed in its OUT_DIR
    // The path pattern is: target/{profile}/build/fastga-rs-{hash}/out/ONEview
    let target_dir = Path::new("target");

    for profile in &["release", "debug"] {
        let build_dir = target_dir.join(profile).join("build");
        if !build_dir.exists() {
            continue;
        }

        for entry in fs::read_dir(&build_dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir()
                && entry
                    .file_name()
                    .to_string_lossy()
                    .starts_with("fastga-rs-")
            {
                let oneview = path.join("out").join("ONEview");
                if oneview.exists() {
                    return Ok(oneview);
                }
            }
        }
    }

    anyhow::bail!("ONEview not found in fastga-rs build output. Run: cargo build --release")
}

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
/// Uses ONEview to convert to text, filters out path-dependent lines ('!' and '<'),
/// and sorts the output to handle non-deterministic record ordering from FastGA multithreading
fn sha256sum_1aln_normalized(path: &Path) -> Result<String> {
    // Run: ONEview file.1aln | grep -v '^[!<]' | sort | sha256sum
    // This pipeline:
    // 1. Converts binary .1aln to text
    // 2. Removes non-deterministic provenance lines ('!') and path-containing headers ('<')
    // 3. Sorts to canonicalize record order (FastGA multithreading creates non-deterministic order)
    // 4. Computes checksum of normalized output

    let oneview_bin = find_oneview()?;

    let mut oneview = Command::new(&oneview_bin)
        .arg(path)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn ONEview")?;

    let mut grep = Command::new("grep")
        .arg("-v")
        .arg("^[!<]") // Filter out both '!' (provenance) and '<' (path-containing headers)
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
            "Golden checksums file not found: {checksums_path:?}\n\
             Run: cd tests/golden_data && ./generate_golden.sh"
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

/// Test that .1aln output matches golden checksums
///
/// Note: This test is skipped on CI (--skip test_golden) because ONEview is not available,
/// but runs by default locally. Run with: cd tests/golden_data && ./generate_golden.sh
#[test]
fn test_golden_1aln_output() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_output.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_output.1aln checksum not found"))?;

    // Use temporary directory (path-independent now that we filter out '<' lines)
    let temp_dir = TempDir::new()?;
    let temp_dir_path = temp_dir.path();
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

    // Use the release binary directly instead of cargo run to avoid subprocess issues
    let sweepga_bin = Path::new("target/release/sweepga");
    if !sweepga_bin.exists() {
        anyhow::bail!("Release binary not found. Run: cargo build --release");
    }

    let result = Command::new(sweepga_bin)
        .arg(temp_input.to_str().unwrap())
        .arg("--1aln")
        .output()?;

    if !result.status.success() {
        anyhow::bail!(
            "sweepga failed with status {:?}\nstderr: {}",
            result.status.code(),
            String::from_utf8_lossy(&result.stderr)
        );
    }

    fs::write(&output, result.stdout)?;

    // Use normalized checksum that filters out non-deterministic provenance records
    let actual = sha256sum_1aln_normalized(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_output.1aln\n\
         Expected: {expected}\n\
         Got:      {actual}\n\
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
         "
    );

    eprintln!("✓ golden_output.1aln checksum matches");
    Ok(())
}

/// Test that PAF output matches golden checksums
///
/// Note: This test is skipped on CI (--skip test_golden) because ONEview is not available,
/// but runs by default locally. Run with: cd tests/golden_data && ./generate_golden.sh
#[test]
fn test_golden_paf_output() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_output.paf")
        .ok_or_else(|| anyhow::anyhow!("golden_output.paf checksum not found"))?;

    // Use temporary directory
    let temp_dir = TempDir::new()?;
    let temp_dir_path = temp_dir.path();
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

    // Use the release binary directly instead of cargo run to avoid subprocess issues
    let sweepga_bin = Path::new("target/release/sweepga");
    if !sweepga_bin.exists() {
        anyhow::bail!("Release binary not found. Run: cargo build --release");
    }

    let result = Command::new(sweepga_bin)
        .args([temp_input.to_str().unwrap(), "--paf"])
        .output()?;

    if !result.status.success() {
        anyhow::bail!(
            "sweepga failed with status {:?}\nstderr: {}",
            result.status.code(),
            String::from_utf8_lossy(&result.stderr)
        );
    }

    fs::write(&output, result.stdout)?;

    let actual = sha256sum(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_output.paf\n\
         Expected: {expected}\n\
         Got:      {actual}\n\
         \n\
         The PAF output format has changed. See above for next steps.\n\
         "
    );

    eprintln!("✓ golden_output.paf checksum matches");
    Ok(())
}

/// Test that 1:1 filtered output matches golden checksums
///
/// Note: This test is skipped on CI (--skip test_golden) because ONEview is not available,
/// but runs by default locally. Run with: cd tests/golden_data && ./generate_golden.sh
#[test]
fn test_golden_filtered_1to1() -> Result<()> {
    let golden_checksums = load_golden_checksums()?;

    let expected = golden_checksums
        .get("golden_filtered_1to1.1aln")
        .ok_or_else(|| anyhow::anyhow!("golden_filtered_1to1.1aln checksum not found"))?;

    // Use temporary directory (path-independent now that we filter out '<' lines)
    let temp_dir = TempDir::new()?;
    let temp_dir_path = temp_dir.path();

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

    // Use the release binary directly instead of cargo run to avoid subprocess issues
    let sweepga_bin = Path::new("target/release/sweepga");
    if !sweepga_bin.exists() {
        anyhow::bail!("Release binary not found. Run: cargo build --release");
    }

    // First generate unfiltered
    let unfiltered = temp_dir_path.join("unfiltered.1aln");
    let result = Command::new(sweepga_bin)
        .arg(temp_input.to_str().unwrap())
        .arg("--1aln")
        .output()?;

    if !result.status.success() {
        anyhow::bail!(
            "sweepga failed with status {:?}\nstderr: {}",
            result.status.code(),
            String::from_utf8_lossy(&result.stderr)
        );
    }

    fs::write(&unfiltered, result.stdout)?;

    // Then filter with 1:1
    let output = temp_dir_path.join("filtered.1aln");
    let result = Command::new(sweepga_bin)
        .args([unfiltered.to_str().unwrap(), "-n", "1:1", "--1aln"])
        .output()?;

    if !result.status.success() {
        anyhow::bail!(
            "sweepga failed with status {:?}\nstderr: {}",
            result.status.code(),
            String::from_utf8_lossy(&result.stderr)
        );
    }

    fs::write(&output, result.stdout)?;

    // Use normalized checksum that filters out non-deterministic provenance records
    let actual = sha256sum_1aln_normalized(&output)?;

    assert_eq!(
        &actual, expected,
        "\n\n\
         ❌ GOLDEN FILE REGRESSION DETECTED!\n\
         \n\
         File: golden_filtered_1to1.1aln\n\
         Expected: {expected}\n\
         Got:      {actual}\n\
         \n\
         The 1:1 filtering behavior has changed.\n\
         "
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
            "Checksum missing for: {filename}\nRun: cd tests/golden_data && ./generate_golden.sh"
        );
    }

    eprintln!("✓ All required checksums present in checksums.txt");
    Ok(())
}
