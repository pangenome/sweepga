// Tests for ANI-based identity filtering

use std::fs;
use std::process::Command;
use tempfile::NamedTempFile;

/// Helper to run sweepga and capture output
fn run_sweepga_with_args(paf_content: &str, args: &[&str]) -> (String, String) {
    // Write PAF to temp file
    let mut temp_file = NamedTempFile::new().unwrap();
    fs::write(temp_file.path(), paf_content).unwrap();

    // Build command
    let mut cmd = Command::new("cargo");
    cmd.arg("run")
       .arg("--release")
       .arg("--quiet")
       .arg("--")
       .arg("-i")
       .arg(temp_file.path());

    for arg in args {
        cmd.arg(arg);
    }

    // Run and capture output
    let output = cmd.output().expect("Failed to run sweepga");

    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    (stdout, stderr)
}

#[test]
fn test_ani_calculation() {
    let paf = "\
genome1#chr1	10000	100	1100	+	genome2#chr1	20000	1000	2000	850	1000	60	dv:f:0.15
genome1#chr1	10000	2000	3000	+	genome2#chr1	20000	3000	4000	900	1000	60	dv:f:0.10
genome1#chr2	10000	100	1100	+	genome2#chr2	20000	1000	2000	800	1000	60	dv:f:0.20
genome2#chr1	20000	100	1100	+	genome3#chr1	15000	1000	2000	700	1000	60	dv:f:0.30
genome2#chr2	20000	100	1100	+	genome3#chr2	15000	1000	2000	750	1000	60	dv:f:0.25
genome1#chr1	10000	100	1100	+	genome3#chr1	15000	1000	2000	600	1000	60	dv:f:0.40
";

    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-y", "ani50", "-j", "0"]);

    // Check ANI statistics are calculated
    assert!(stderr.contains("ANI statistics from 3 genome pairs"));
    assert!(stderr.contains("Min: 60.0%"));
    assert!(stderr.contains("Median: 72.5%"));  // (60 + 72.5 + 85) / 3, median is middle value
    assert!(stderr.contains("Max: 85.0%"));
    assert!(stderr.contains("Using minimum identity threshold: 72.5%"));
}

#[test]
fn test_ani_offset_syntax() {
    let paf = "\
genome1#chr1	10000	100	1000	+	genome2#chr1	20000	1000	2000	900	1000	60	dv:f:0.10
genome2#chr1	20000	100	1000	+	genome3#chr1	15000	1000	2000	800	1000	60	dv:f:0.20
genome1#chr1	10000	100	1000	+	genome3#chr1	15000	1000	2000	700	1000	60	dv:f:0.30
";

    // Test ani50-5 (median minus 5%)
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-y", "ani50-5", "-j", "0"]);
    assert!(stderr.contains("Using minimum identity threshold: 75.0%"));  // 80% - 5%

    // Test ani50+10 (median plus 10%)
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-y", "ani50+10", "-j", "0"]);
    assert!(stderr.contains("Using minimum identity threshold: 90.0%"));  // 80% + 10%
}

#[test]
fn test_percentage_vs_fraction() {
    let paf = "seq1	10000	100	1000	+	chr1	20000	1000	2000	900	1000	60\n";

    // Test percentage (90 = 90%)
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-y", "90", "-j", "0"]);
    assert!(stderr.contains("Using minimum identity threshold: 90.0%"));

    // Test fraction (0.85 = 85%)
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-y", "0.85", "-j", "0"]);
    assert!(stderr.contains("Using minimum identity threshold: 85.0%"));
}

#[test]
fn test_scaffold_identity_filter() {
    let paf = "\
seq1	10000	100	1100	+	chr1	20000	1000	2000	900	1000	60	dv:f:0.20
seq1	10000	1200	2200	+	chr1	20000	2100	3100	950	1000	60	dv:f:0.15
seq2	10000	100	1100	+	chr2	20000	1000	2000	700	1000	60	dv:f:0.40
seq2	10000	1200	2200	+	chr2	20000	2100	3100	750	1000	60	dv:f:0.35
";

    // Test -Y filter (75% minimum for scaffolds)
    let (stdout, stderr) = run_sweepga_with_args(paf, &["-j", "500", "-s", "1000", "-Y", "75", "-m", "1:1", "-y", "0"]);

    // Should filter out seq2 chain (avg identity ~65%)
    assert!(stderr.contains("Length/identity filter: 1 chains kept, 1 removed"));

    // Output should only have seq1 mappings
    let output_lines: Vec<&str> = stdout.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(output_lines.len(), 2);
    assert!(output_lines.iter().all(|l| l.starts_with("seq1")));
}

#[test]
fn test_default_ani50() {
    let paf = "\
genome1#chr1	10000	100	1000	+	genome2#chr1	20000	1000	2000	850	1000	60	dv:f:0.15
genome2#chr1	20000	100	1000	+	genome3#chr1	15000	1000	2000	700	1000	60	dv:f:0.30
";

    // Test default (should be ani50)
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["-j", "0"]);

    // Should calculate ANI and use median
    assert!(stderr.contains("ANI statistics"));
    assert!(stderr.contains("Using minimum identity threshold"));
}

#[test]
fn test_matches_scoring() {
    let paf = "\
seq1	10000	100	1000	+	chr1	20000	1000	1900	850	900	60	cg:Z:850=50X
seq1	10000	100	1000	+	chr1	20000	1000	1900	700	900	60	cg:Z:700=200X
";

    // Test that matches scoring prefers more matches
    let (stdout, _stderr) = run_sweepga_with_args(paf, &["--scoring", "matches", "-n", "1", "-j", "0", "-y", "0"]);

    // Should keep only the first mapping with 850 matches
    let output_lines: Vec<&str> = stdout.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(output_lines.len(), 1);
    assert!(output_lines[0].contains("850=50X"));
}

#[test]
fn test_cigar_warning() {
    let paf = "seq1	10000	100	1000	+	chr1	20000	1000	2000	900	1000	60\n";

    // Test warning when using matches scoring without CIGAR
    let (_stdout, stderr) = run_sweepga_with_args(paf, &["--scoring", "matches", "-y", "0", "-j", "0"]);

    assert!(stderr.contains("WARNING: Using 'matches' scoring but input lacks CIGAR strings"));
    assert!(stderr.contains("Match counts will be estimated from PAF matches field"));
}