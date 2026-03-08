/// CI integration tests using B-3106.fa (~32KB, 9 sequences of ~3.3kb)
///
/// These tests exercise both wfmash and FastGA aligners end-to-end.
/// wfmash segment/block lengths are adapted from the average sequence
/// length via the .fai index, so small sequences work correctly.
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

fn sweepga_bin() -> &'static str {
    if cfg!(debug_assertions) {
        "target/debug/sweepga"
    } else {
        "target/release/sweepga"
    }
}

fn count_paf_lines(path: &Path) -> usize {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count()
}

/// Copy B-3106.fa to a temp dir and create a .fai index via samtools.
fn setup_b3106(temp_dir: &Path) -> std::path::PathBuf {
    let src = Path::new("data/B-3106.fa");
    assert!(src.exists(), "data/B-3106.fa not found");
    let dst = temp_dir.join("B-3106.fa");
    fs::copy(src, &dst).unwrap();
    // Create .fai index
    let status = Command::new("samtools")
        .arg("faidx")
        .arg(&dst)
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools faidx failed");
    dst
}

/// Test wfmash aligner on B-3106.fa
#[test]
fn test_ci_wfmash_b3106() {
    let temp_dir = TempDir::new().unwrap();
    let input = setup_b3106(temp_dir.path());
    let output = temp_dir.path().join("wfmash.paf");

    let result = Command::new(sweepga_bin())
        .arg(&input)
        .arg("--aligner")
        .arg("wfmash")
        .arg("--paf")
        .arg("-t")
        .arg("1")
        .arg("-n")
        .arg("N:N")
        .arg("-j")
        .arg("0")
        .stdout(fs::File::create(&output).unwrap())
        .output()
        .expect("Failed to run sweepga");

    assert!(
        result.status.success(),
        "wfmash alignment failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let n = count_paf_lines(&output);
    assert!(n > 0, "wfmash produced no alignments");
    eprintln!("wfmash on B-3106.fa: {n} alignments");
}

/// Test FastGA aligner on B-3106.fa
#[test]
fn test_ci_fastga_b3106() {
    let temp_dir = TempDir::new().unwrap();
    let input = setup_b3106(temp_dir.path());
    let output = temp_dir.path().join("fastga.paf");

    let result = Command::new(sweepga_bin())
        .arg(&input)
        .arg("--aligner")
        .arg("fastga")
        .arg("--paf")
        .arg("-t")
        .arg("1")
        .arg("-n")
        .arg("N:N")
        .arg("-j")
        .arg("0")
        .stdout(fs::File::create(&output).unwrap())
        .output()
        .expect("Failed to run sweepga");

    assert!(
        result.status.success(),
        "FastGA alignment failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let n = count_paf_lines(&output);
    assert!(n > 0, "FastGA produced no alignments");
    eprintln!("FastGA on B-3106.fa: {n} alignments");
}
