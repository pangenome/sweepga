/// Golden file regression tests using SHA256 checksums
///
/// These tests lock down current behavior to prevent unintended changes.
/// ANY change to output format/coordinates/filtering will fail these tests.
///
/// Only the PAF output has a checksum. The `.1aln` output depends on the
/// machine (FastGA output changes with the thread count), so
/// `tests/test_1aln_output.rs` checks properties of it.
///
/// Run with: `cargo test test_golden`
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

/// Test that PAF output matches golden checksums
///
/// Regenerate with: cd tests/golden_data && ./generate_golden.sh
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

/// Verify all required checksums are present
#[test]
fn test_golden_files_complete() -> Result<()> {
    let checksums = load_golden_checksums()?;

    let required_files = [
        "golden_output.paf",
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
