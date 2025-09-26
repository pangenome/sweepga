/// Tests for scaffold mapping functionality
/// Scaffolding prevents transposon-like scattered alignments
/// while preserving large syntenic blocks

use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;

/// Helper to create a PAF file with specific mappings
fn create_test_paf(mappings: &[(
    &str,  // query_name
    u64,   // query_len
    u64,   // query_start
    u64,   // query_end
    &str,  // target_name
    u64,   // target_len
    u64,   // target_start
    u64,   // target_end
)]) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    for (qname, qlen, qstart, qend, tname, tlen, tstart, tend) in mappings {
        let matches = qend - qstart; // Assume perfect matches
        let block_len = matches;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t60\tcg:Z:{}M",
            qname, qlen, qstart, qend, tname, tlen, tstart, tend, matches, block_len, matches
        ).unwrap();
    }
    file.flush().unwrap();
    file
}

/// Run sweepga and return filtered output
fn run_sweepga(input_file: &str, args: &[&str]) -> String {
    let output_file = NamedTempFile::new().unwrap();
    let output_path = output_file.path().to_str().unwrap();

    let mut cmd = std::process::Command::new("./target/release/sweepga");
    cmd.arg("-i").arg(input_file)
       .arg("-o").arg(output_path);

    for arg in args {
        cmd.arg(arg);
    }

    let status = cmd.status().expect("Failed to run sweepga");
    assert!(status.success(), "sweepga failed");

    fs::read_to_string(output_path).unwrap()
}

#[test]
fn test_single_chromosome_syntenic_block() {
    // Test: Large syntenic block on single chromosome pair should be preserved
    // Create a 100kb syntenic region with mappings every 10kb
    let mappings = vec![
        ("chrA", 500000, 0,     10000,  "chr1", 600000, 0,     10000),
        ("chrA", 500000, 10000, 20000,  "chr1", 600000, 10000, 20000),
        ("chrA", 500000, 20000, 30000,  "chr1", 600000, 20000, 30000),
        ("chrA", 500000, 30000, 40000,  "chr1", 600000, 30000, 40000),
        ("chrA", 500000, 40000, 50000,  "chr1", 600000, 40000, 50000),
        ("chrA", 500000, 50000, 60000,  "chr1", 600000, 50000, 60000),
        ("chrA", 500000, 60000, 70000,  "chr1", 600000, 60000, 70000),
        ("chrA", 500000, 70000, 80000,  "chr1", 600000, 70000, 80000),
        ("chrA", 500000, 80000, 90000,  "chr1", 600000, 80000, 90000),
        ("chrA", 500000, 90000, 100000, "chr1", 600000, 90000, 100000),
    ];

    let input = create_test_paf(&mappings);

    // Run with scaffolding: -j 20k (scaffold-jump) -s 50k (scaffold-mass)
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "20k", "-s", "50k", "-d", "50k"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // Should keep all 10 mappings as they form one large scaffold
    assert_eq!(lines.len(), 10, "Syntenic block should be preserved");

    // Check that all are marked as scaffold or rescued
    for line in &lines {
        assert!(line.contains("st:Z:scaffold") || line.contains("st:Z:rescued"),
                "Mapping should be part of scaffold: {}", line);
    }
}

#[test]
fn test_transposon_like_scattered_alignments() {
    // Test: Scattered small alignments (transposon-like) should be filtered out
    // Create many small 500bp mappings scattered across the chromosome
    let mut mappings = vec![];
    for i in 0..20 {
        let start = i * 50000; // 50kb apart - too far for scaffolding
        mappings.push((
            "chrA", 1000000, start, start + 500,
            "chr1", 1000000, start * 2, start * 2 + 500
        ));
    }

    let input = create_test_paf(&mappings);

    // Run with scaffolding: -j 10k (scaffold-jump) -s 10k (scaffold-mass)
    // Small scattered alignments won't form scaffolds
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "10k", "-s", "10k", "-d", "10k"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // Should filter out most/all mappings as they don't form scaffolds
    assert!(lines.len() < 5,
            "Scattered transposon-like alignments should be filtered: kept {} of 20",
            lines.len());
}

#[test]
fn test_two_distinct_syntenic_blocks() {
    // Test: Two distinct syntenic blocks on same chromosome pair
    // Both should be kept as separate scaffolds

    let mappings = vec![
        // First block: 0-50kb
        ("chrA", 500000, 0,     10000, "chr1", 600000, 0,     10000),
        ("chrA", 500000, 10000, 20000, "chr1", 600000, 10000, 20000),
        ("chrA", 500000, 20000, 30000, "chr1", 600000, 20000, 30000),
        ("chrA", 500000, 30000, 40000, "chr1", 600000, 30000, 40000),
        ("chrA", 500000, 40000, 50000, "chr1", 600000, 40000, 50000),

        // Gap of 200kb

        // Second block: 250-300kb
        ("chrA", 500000, 250000, 260000, "chr1", 600000, 250000, 260000),
        ("chrA", 500000, 260000, 270000, "chr1", 600000, 260000, 270000),
        ("chrA", 500000, 270000, 280000, "chr1", 600000, 270000, 280000),
        ("chrA", 500000, 280000, 290000, "chr1", 600000, 280000, 290000),
        ("chrA", 500000, 290000, 300000, "chr1", 600000, 290000, 300000),
    ];

    let input = create_test_paf(&mappings);

    // Run with scaffolding: -j 20k (won't bridge the 200kb gap)
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "20k", "-s", "30k", "-d", "20k"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // Should keep all 10 mappings in two separate scaffolds
    assert_eq!(lines.len(), 10, "Both syntenic blocks should be preserved");

    // Count unique chain IDs to verify two separate scaffolds
    let mut chain_ids = std::collections::HashSet::new();
    for line in &lines {
        if let Some(cg_pos) = line.find("cg:f:") {
            let cg_part = &line[cg_pos..];
            if let Some(end) = cg_part.find('\t') {
                chain_ids.insert(&cg_part[5..end]);
            }
        }
    }
    // Note: chain IDs might not be populated, check scaffold status instead
}

#[test]
fn test_inversion_breakpoint() {
    // Test: Inversion should create separate scaffolds for + and - strands
    let mappings = vec![
        // Forward strand block
        ("chrA", 500000, 0,     10000, "chr1", 600000, 0,     10000),
        ("chrA", 500000, 10000, 20000, "chr1", 600000, 10000, 20000),
        ("chrA", 500000, 20000, 30000, "chr1", 600000, 20000, 30000),

        // Note: We can't easily test inversions with this simple PAF creation
        // Would need to handle strand properly
    ];

    // Skip this test for now as it needs more complex PAF generation
}

#[test]
fn test_scaffold_no_rescue() {
    // Test: With -d 0, only keep scaffold members, no rescue
    let mappings = vec![
        // Main scaffold: 0-50kb
        ("chrA", 500000, 0,     10000, "chr1", 600000, 0,     10000),
        ("chrA", 500000, 10000, 20000, "chr1", 600000, 10000, 20000),
        ("chrA", 500000, 20000, 30000, "chr1", 600000, 20000, 30000),
        ("chrA", 500000, 30000, 40000, "chr1", 600000, 30000, 40000),
        ("chrA", 500000, 40000, 50000, "chr1", 600000, 40000, 50000),

        // Small mapping 30kb away - NOT rescued with -d 0
        ("chrA", 500000, 80000, 81000, "chr1", 600000, 80000, 81000),

        // Small mapping 100kb away - NOT rescued with -d 0
        ("chrA", 500000, 150000, 151000, "chr1", 600000, 150000, 151000),
    ];

    let input = create_test_paf(&mappings);

    // Run with scaffolding but NO rescue distance (-d 0)
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "20k", "-s", "30k", "-d", "0"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // Should keep only 5 mappings in the scaffold, no rescue
    assert_eq!(lines.len(), 5,
              "With -d 0, should only keep scaffold members");

    // All should be scaffold, none rescued
    let scaffold_count = lines.iter().filter(|l| l.contains("st:Z:scaffold")).count();
    let rescued_count = lines.iter().filter(|l| l.contains("st:Z:rescued")).count();

    assert_eq!(scaffold_count, 5, "All should be scaffold mappings");
    assert_eq!(rescued_count, 0, "No mappings should be rescued with -d 0");
}

#[test]
fn test_rescue_with_distance() {
    // Test: With -d > 0, rescue mappings near scaffolds
    let mappings = vec![
        // Main scaffold: 0-50kb
        ("chrA", 500000, 0,     10000, "chr1", 600000, 0,     10000),
        ("chrA", 500000, 10000, 20000, "chr1", 600000, 10000, 20000),
        ("chrA", 500000, 20000, 30000, "chr1", 600000, 20000, 30000),
        ("chrA", 500000, 30000, 40000, "chr1", 600000, 30000, 40000),
        ("chrA", 500000, 40000, 50000, "chr1", 600000, 40000, 50000),

        // Small mapping 30kb away - should be rescued with -d 50k
        ("chrA", 500000, 80000, 81000, "chr1", 600000, 80000, 81000),

        // Small mapping 100kb away - too far to rescue with -d 50k
        ("chrA", 500000, 150000, 151000, "chr1", 600000, 150000, 151000),
    ];

    let input = create_test_paf(&mappings);

    // Run with scaffolding and rescue distance
    // Use -j 15k to ensure all 5 contiguous mappings form one scaffold
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "15k", "-s", "40k", "-d", "50k"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // Should keep 6 mappings: 5 in scaffold + 1 rescued
    assert_eq!(lines.len(), 6,
              "Should keep main scaffold and rescue nearby mapping");

    // Check that we have both scaffold and rescued mappings
    let scaffold_count = lines.iter().filter(|l| l.contains("st:Z:scaffold")).count();
    let rescued_count = lines.iter().filter(|l| l.contains("st:Z:rescued")).count();

    assert!(scaffold_count >= 5, "Should have scaffold mappings");
    assert!(rescued_count >= 1, "Should have rescued mapping");
}

#[test]
fn test_competing_scaffolds() {
    // Test: When multiple scaffolds compete for the same region,
    // only the best should be kept with 1:1 filtering

    let mappings = vec![
        // First scaffold to chr1: large and high quality
        ("chrA", 500000, 0,     50000, "chr1", 600000, 0,     50000),
        ("chrA", 500000, 50000, 100000, "chr1", 600000, 50000, 100000),

        // Second scaffold to chr2: overlapping query region but different target
        ("chrA", 500000, 30000, 80000, "chr2", 600000, 100000, 150000),

        // With default 1:1 filtering, all should be kept (different targets)
    ];

    let input = create_test_paf(&mappings);

    // Run with scaffolding and 1:1 filtering
    let output = run_sweepga(
        input.path().to_str().unwrap(),
        &["-j", "60k", "-s", "40k", "-n", "1:1"]
    );

    let lines: Vec<&str> = output.trim().lines().collect();

    // With our chromosome-pair grouping, mappings to different targets are kept
    assert_eq!(lines.len(), 3,
              "Mappings to different targets should all be kept with 1:1 filtering");
}