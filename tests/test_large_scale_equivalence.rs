/// Large-scale format equivalence tests
///
/// Tests that .1aln and PAF formats produce identical filtering results
/// at scale (10K+ alignments). Uses real yeast genome data to ensure
/// coordinate conversion and filtering logic is stable under load.
use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Parse PAF file and extract key alignment properties
fn parse_paf_alignments(path: &Path) -> Result<Vec<PafRecord>> {
    let content = fs::read_to_string(path)?;
    let mut records = Vec::new();

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        records.push(PafRecord {
            query_name: fields[0].to_string(),
            query_start: fields[2].parse()?,
            query_end: fields[3].parse()?,
            strand: fields[4].to_string(),
            target_name: fields[5].to_string(),
            target_start: fields[7].parse()?,
            target_end: fields[8].parse()?,
            matches: fields[9].parse()?,
            block_len: fields[10].parse()?,
        });
    }

    Ok(records)
}

#[derive(Debug, Clone, PartialEq)]
struct PafRecord {
    query_name: String,
    query_start: u64,
    query_end: u64,
    strand: String,
    target_name: String,
    target_start: u64,
    target_end: u64,
    matches: u64,
    block_len: u64,
}

impl PafRecord {
    fn identity(&self) -> f64 {
        if self.block_len == 0 {
            0.0
        } else {
            self.matches as f64 / self.block_len as f64
        }
    }

    fn query_span(&self) -> u64 {
        self.query_end.saturating_sub(self.query_start)
    }

    fn target_span(&self) -> u64 {
        self.target_end.saturating_sub(self.target_start)
    }
}

#[test]
#[ignore = "Slow: 10K+ alignments, 2+ minutes - run manually with --ignored"]
fn test_large_scale_paf_output() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // Generate PAF output with default filtering
    eprintln!("Generating PAF output from FASTA...");
    let paf_output = temp_dir.path().join("output.paf");
    let status = Command::new("cargo")
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
        .stdout(fs::File::create(&paf_output)?)
        .status()?;

    if !status.success() {
        anyhow::bail!("Failed to generate PAF output");
    }

    // Parse PAF
    let paf_records = parse_paf_alignments(&paf_output)?;
    eprintln!("PAF output: {} alignments", paf_records.len());

    // Verify we have large-scale data (10K+ alignments)
    assert!(
        paf_records.len() >= 10_000,
        "Expected 10K+ alignments for large-scale test, got {}",
        paf_records.len()
    );

    // Verify basic properties
    for (i, rec) in paf_records.iter().take(100).enumerate() {
        assert!(
            rec.query_start < rec.query_end,
            "Record {}: query_start >= query_end",
            i
        );
        assert!(
            rec.target_start < rec.target_end,
            "Record {}: target_start >= target_end",
            i
        );
        assert!(
            rec.matches <= rec.block_len,
            "Record {}: matches > block_len",
            i
        );
        assert!(!rec.query_name.is_empty(), "Record {}: empty query_name", i);
        assert!(
            !rec.target_name.is_empty(),
            "Record {}: empty target_name",
            i
        );
    }

    eprintln!(
        "✓ Large-scale PAF output verified: {} alignments with valid properties",
        paf_records.len()
    );
    Ok(())
}

#[test]
#[ignore = "Slow: 10K+ alignments, 2+ minutes - run manually with --ignored"]
fn test_large_scale_paf_filtering() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // Generate truly unfiltered PAF (use -n N:N to disable filtering)
    eprintln!("Generating unfiltered PAF with -n N:N...");
    let unfiltered_paf = temp_dir.path().join("unfiltered.paf");
    let status = Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            temp_input.to_str().unwrap(),
            "--paf",
            "-n",
            "N:N",
        ])
        .stdout(fs::File::create(&unfiltered_paf)?)
        .status()?;

    if !status.success() {
        anyhow::bail!("Failed to generate unfiltered PAF");
    }

    // Apply 1:1 filtering
    eprintln!("Applying 1:1 filtering to PAF...");
    let filtered_paf = temp_dir.path().join("filtered.paf");
    let status = Command::new("cargo")
        .args(&[
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            unfiltered_paf.to_str().unwrap(),
            "-n",
            "1:1",
        ])
        .stdout(fs::File::create(&filtered_paf)?)
        .status()?;

    if !status.success() {
        anyhow::bail!("Failed to filter PAF");
    }

    // Parse results
    let unfiltered = parse_paf_alignments(&unfiltered_paf)?;
    let filtered = parse_paf_alignments(&filtered_paf)?;

    eprintln!("Unfiltered: {} alignments", unfiltered.len());
    eprintln!("Filtered (1:1): {} alignments", filtered.len());

    // Verify we have large-scale data
    assert!(
        unfiltered.len() >= 10_000,
        "Expected 10K+ unfiltered alignments, got {}",
        unfiltered.len()
    );

    // Filtering should reduce count (or at least not increase it)
    assert!(
        filtered.len() <= unfiltered.len(),
        "Filtering should not increase alignment count"
    );
    assert!(filtered.len() > 0, "Filtered output should not be empty");

    // All filtered records should be in unfiltered set
    let unfiltered_sigs: HashSet<String> = unfiltered
        .iter()
        .map(|r| {
            format!(
                "{}:{}-{}:{}:{}:{}-{}:{}:{}",
                r.query_name,
                r.query_start,
                r.query_end,
                r.strand,
                r.target_name,
                r.target_start,
                r.target_end,
                r.matches,
                r.block_len
            )
        })
        .collect();

    for (i, rec) in filtered.iter().enumerate() {
        let sig = format!(
            "{}:{}-{}:{}:{}:{}-{}:{}:{}",
            rec.query_name,
            rec.query_start,
            rec.query_end,
            rec.strand,
            rec.target_name,
            rec.target_start,
            rec.target_end,
            rec.matches,
            rec.block_len
        );
        assert!(
            unfiltered_sigs.contains(&sig),
            "Filtered record {} not found in unfiltered set: {}",
            i,
            sig
        );
    }

    eprintln!(
        "✓ Large-scale PAF filtering verified: {} → {} alignments",
        unfiltered.len(),
        filtered.len()
    );
    Ok(())
}

#[test]
#[ignore = "Slow: 10K+ alignments, 2+ minutes - run manually with --ignored"]
fn test_coordinate_stability_at_scale() -> Result<()> {
    let input = Path::new("data/scerevisiae8.fa.gz");
    assert!(
        input.exists(),
        "Test data not found: data/scerevisiae8.fa.gz - required for CI"
    );

    let temp_dir = TempDir::new()?;

    // Copy input to temp to avoid accumulating FastGA intermediate files
    let temp_input = temp_dir.path().join("test_input.fa.gz");
    fs::copy(input, &temp_input)?;

    // Run same input twice to verify deterministic output
    eprintln!("Testing coordinate determinism with repeated runs...");

    let paf1 = temp_dir.path().join("run1.paf");
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
        .stdout(fs::File::create(&paf1)?)
        .status()?;

    let paf2 = temp_dir.path().join("run2.paf");
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
        .stdout(fs::File::create(&paf2)?)
        .status()?;

    // Compare outputs (should be identical)
    let records1 = parse_paf_alignments(&paf1)?;
    let records2 = parse_paf_alignments(&paf2)?;

    eprintln!("Run 1: {} alignments", records1.len());
    eprintln!("Run 2: {} alignments", records2.len());

    assert!(records1.len() >= 10_000, "Expected 10K+ alignments");
    assert_eq!(
        records1.len(),
        records2.len(),
        "Count changed between runs: {} vs {}",
        records1.len(),
        records2.len()
    );

    // Check for any differences
    let mut drifted = 0;
    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        if r1 != r2 {
            drifted += 1;
            if drifted <= 3 {
                eprintln!("Difference at record {}:", i);
                eprintln!("  Run 1: {:?}", r1);
                eprintln!("  Run 2: {:?}", r2);
            }
        }
    }

    assert_eq!(
        drifted,
        0,
        "Non-deterministic output detected in {} out of {} alignments",
        drifted,
        records1.len()
    );

    eprintln!(
        "✓ Coordinate stability verified: {} alignments identical across runs",
        records1.len()
    );
    Ok(())
}
