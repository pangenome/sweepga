/// Determinism regression tests.
///
/// sweepga must produce byte-identical output across consecutive runs on the
/// same input. Two paths are exercised:
///
/// 1. `test_filter_determinism`: filter the same raw PAF twice. Fast (~2s).
///    Isolates non-determinism inside the filter (the known problem area).
///
/// 2. `test_full_pipeline_determinism`: run the full FASTA → PAF pipeline
///    twice. Slow (~60s). Ensures the aligner + filter combination is stable.
///    Marked `#[ignore]` so `cargo test` skips it by default; run with
///    `cargo test --release --test test_determinism -- --ignored`.
///
/// Both tests must FAIL on `main` until the HashMap-iteration-order bug is
/// fixed (see `union_find::get_sets` and HashMap iterations in `paf_filter`).
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn sweepga_bin() -> &'static str {
    if cfg!(debug_assertions) {
        "./target/debug/sweepga"
    } else {
        "./target/release/sweepga"
    }
}

fn require_binary() {
    let bin = Path::new(sweepga_bin());
    assert!(
        bin.exists(),
        "{} not found — run `cargo build --release` first",
        sweepga_bin()
    );
}

fn fixture() -> PathBuf {
    let p = PathBuf::from("data/scerevisiae8.fa.gz");
    assert!(p.exists(), "fixture missing: {}", p.display());
    p
}

/// Generate the raw PAF once. Returns its path inside `dir`.
fn generate_raw_paf(dir: &Path) -> PathBuf {
    let raw = dir.join("raw.paf");
    let out = Command::new(sweepga_bin())
        .arg("--no-filter")
        .arg(fixture())
        .output()
        .expect("failed to invoke sweepga --no-filter");
    assert!(
        out.status.success(),
        "sweepga --no-filter failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::write(&raw, out.stdout).expect("write raw.paf");
    raw
}

/// Filter `input` once, returning the bytes written to stdout.
fn filter_once(input: &Path) -> Vec<u8> {
    let out = Command::new(sweepga_bin())
        .arg(input)
        .output()
        .expect("failed to invoke sweepga (filter)");
    assert!(
        out.status.success(),
        "sweepga filter failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn first_diff(a: &[u8], b: &[u8]) -> Option<(usize, u8, u8)> {
    a.iter()
        .zip(b.iter())
        .enumerate()
        .find(|(_, (x, y))| x != y)
        .map(|(i, (x, y))| (i, *x, *y))
}

fn line_count(bytes: &[u8]) -> usize {
    bytes.iter().filter(|&&b| b == b'\n').count()
}

#[test]
fn test_filter_determinism() {
    require_binary();
    let dir = TempDir::new().expect("tempdir");
    let raw = generate_raw_paf(dir.path());

    let run1 = filter_once(&raw);
    let run2 = filter_once(&raw);

    if run1 != run2 {
        let n1 = line_count(&run1);
        let n2 = line_count(&run2);
        let diff = first_diff(&run1, &run2);
        panic!(
            "filter is non-deterministic: \
             run1 = {} bytes / {} lines, run2 = {} bytes / {} lines, \
             first byte diff at offset {:?}",
            run1.len(),
            n1,
            run2.len(),
            n2,
            diff
        );
    }
}

#[test]
#[ignore = "slow (~60s); run with --ignored"]
fn test_full_pipeline_determinism() {
    require_binary();
    let dir = TempDir::new().expect("tempdir");

    let run = |label: &str| -> Vec<u8> {
        let out = Command::new(sweepga_bin())
            .arg(fixture())
            .output()
            .unwrap_or_else(|e| panic!("{label}: spawn failed: {e}"));
        assert!(
            out.status.success(),
            "{label}: sweepga failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        out.stdout
    };

    let _scratch = dir.path(); // keep tempdir alive for any temp files sweepga writes

    let run1 = run("run1");
    let run2 = run("run2");

    if run1 != run2 {
        panic!(
            "full pipeline is non-deterministic: \
             run1 = {} bytes / {} lines, run2 = {} bytes / {} lines",
            run1.len(),
            line_count(&run1),
            run2.len(),
            line_count(&run2),
        );
    }
}
