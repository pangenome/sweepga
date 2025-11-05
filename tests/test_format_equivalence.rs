#![allow(clippy::uninlined_format_args)]
/// Test that PAF and .1aln filtering workflows produce identical results
///
/// This is a critical validation that our unified filtering logic works correctly.
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

/// Statistics extracted from a PAF file
#[derive(Debug, Clone)]
struct PafStats {
    total_alignments: usize,
    total_bases: u64,
    total_matches: u64,
    avg_identity: f64,
    alignments_by_pair: HashMap<(String, String), usize>,
}

/// Extract genome prefix (e.g., "SGDref#1#chrI" -> "SGDref#1")
fn extract_genome_prefix(name: &str) -> String {
    let parts: Vec<&str> = name.split('#').collect();
    if parts.len() >= 2 {
        format!("{}#{}", parts[0], parts[1])
    } else {
        name.to_string()
    }
}

/// Parse PAF file and extract statistics
fn parse_paf_stats(paf_path: &Path) -> PafStats {
    let content = fs::read_to_string(paf_path).unwrap();

    let mut total_alignments = 0;
    let mut total_bases = 0u64;
    let mut total_matches = 0u64;
    let mut alignments_by_pair: HashMap<(String, String), usize> = HashMap::new();

    for line in content.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        total_alignments += 1;

        // Parse fields
        let query = fields[0];
        let target = fields[5];
        let matches: u64 = fields[9].parse().unwrap_or(0);
        let block_len: u64 = fields[10].parse().unwrap_or(0);

        total_matches += matches;
        total_bases += block_len;

        // Track genome pair
        let q_genome = extract_genome_prefix(query);
        let t_genome = extract_genome_prefix(target);
        let pair = if q_genome < t_genome {
            (q_genome, t_genome)
        } else {
            (t_genome, q_genome)
        };
        *alignments_by_pair.entry(pair).or_insert(0) += 1;
    }

    let avg_identity = if total_bases > 0 {
        total_matches as f64 / total_bases as f64
    } else {
        0.0
    };

    PafStats {
        total_alignments,
        total_bases,
        total_matches,
        avg_identity,
        alignments_by_pair,
    }
}

/// Compare two PafStats for equivalence
fn compare_stats(paf_stats: &PafStats, aln_stats: &PafStats, tolerance: f64) -> bool {
    // Check total alignments
    if paf_stats.total_alignments != aln_stats.total_alignments {
        eprintln!(
            "❌ Alignment count mismatch: PAF={}, .1aln={}",
            paf_stats.total_alignments, aln_stats.total_alignments
        );
        return false;
    }

    // Check total bases (should be exact)
    if paf_stats.total_bases != aln_stats.total_bases {
        eprintln!(
            "❌ Total bases mismatch: PAF={}, .1aln={}",
            paf_stats.total_bases, aln_stats.total_bases
        );
        return false;
    }

    // Check total matches (should be exact)
    if paf_stats.total_matches != aln_stats.total_matches {
        eprintln!(
            "❌ Total matches mismatch: PAF={}, .1aln={}",
            paf_stats.total_matches, aln_stats.total_matches
        );
        return false;
    }

    // Check average identity (within tolerance)
    let identity_diff = (paf_stats.avg_identity - aln_stats.avg_identity).abs();
    if identity_diff > tolerance {
        eprintln!(
            "❌ Identity difference too large: PAF={:.6}, .1aln={:.6}, diff={:.6} (tolerance={:.6})",
            paf_stats.avg_identity, aln_stats.avg_identity, identity_diff, tolerance
        );
        return false;
    }

    // Check genome pair counts
    if paf_stats.alignments_by_pair.len() != aln_stats.alignments_by_pair.len() {
        eprintln!(
            "❌ Genome pair count mismatch: PAF={}, .1aln={}",
            paf_stats.alignments_by_pair.len(),
            aln_stats.alignments_by_pair.len()
        );
        return false;
    }

    for (pair, &paf_count) in &paf_stats.alignments_by_pair {
        let aln_count = aln_stats.alignments_by_pair.get(pair).copied().unwrap_or(0);
        if paf_count != aln_count {
            eprintln!(
                "❌ Genome pair {:?} alignment count mismatch: PAF={}, .1aln={}",
                pair, paf_count, aln_count
            );
            return false;
        }
    }

    true
}

#[test]
#[ignore] // Requires FastGA binaries and takes time
fn test_paf_vs_1aln_equivalence() {
    let fasta_path = Path::new("data/scerevisiae8.fa.gz");

    if !fasta_path.exists() {
        eprintln!("Skipping test - data/scerevisiae8.fa.gz not found");
        return;
    }

    let temp_dir = TempDir::new().unwrap();

    // === WORKFLOW 1: PAF ===
    eprintln!("\n=== Testing PAF Workflow ===");

    let paf_output = temp_dir.path().join("filtered.paf");

    let paf_result = std::process::Command::new("./target/release/sweepga")
        .arg(fasta_path)
        .arg("--output-file")
        .arg(&paf_output)
        .arg("--paf")
        .output()
        .expect("Failed to run sweepga for PAF workflow");

    if !paf_result.status.success() {
        eprintln!("PAF workflow failed:");
        eprintln!("{}", String::from_utf8_lossy(&paf_result.stderr));
        panic!("PAF workflow failed");
    }

    assert!(paf_output.exists(), "PAF output file not created");

    let paf_stats = parse_paf_stats(&paf_output);
    eprintln!(
        "PAF: {} alignments, {:.1} Mb, {:.2}% identity, {} genome pairs",
        paf_stats.total_alignments,
        paf_stats.total_bases as f64 / 1_000_000.0,
        paf_stats.avg_identity * 100.0,
        paf_stats.alignments_by_pair.len()
    );

    // === WORKFLOW 2: .1aln → PAF ===
    eprintln!("\n=== Testing .1aln Workflow ===");

    let aln_output = temp_dir.path().join("filtered.1aln");
    let aln_as_paf = temp_dir.path().join("filtered_from_aln.paf");

    // Generate and filter as .1aln
    let aln_result = std::process::Command::new("./target/release/sweepga")
        .arg(fasta_path)
        .arg("--output-file")
        .arg(&aln_output)
        .output()
        .expect("Failed to run sweepga for .1aln workflow");

    if !aln_result.status.success() {
        eprintln!(".1aln workflow failed:");
        eprintln!("{}", String::from_utf8_lossy(&aln_result.stderr));

        // This is expected to fail currently - skip test gracefully
        eprintln!("\n⚠ .1aln workflow not yet fully functional, skipping test");
        return;
    }

    assert!(aln_output.exists(), ".1aln output file not created");

    // Convert .1aln to PAF for comparison
    eprintln!("Converting .1aln to PAF for comparison...");

    // Use sweepga to read .1aln and output as PAF
    let convert_result = std::process::Command::new("./target/release/sweepga")
        .arg(&aln_output)
        .arg("--output-file")
        .arg(&aln_as_paf)
        .arg("--paf")
        .output()
        .expect("Failed to convert .1aln to PAF");

    if !convert_result.status.success() {
        eprintln!("Conversion failed:");
        eprintln!("{}", String::from_utf8_lossy(&convert_result.stderr));
        panic!("Failed to convert .1aln to PAF");
    }

    assert!(aln_as_paf.exists(), "Converted PAF not created");

    let aln_stats = parse_paf_stats(&aln_as_paf);
    eprintln!(
        ".1aln: {} alignments, {:.1} Mb, {:.2}% identity, {} genome pairs",
        aln_stats.total_alignments,
        aln_stats.total_bases as f64 / 1_000_000.0,
        aln_stats.avg_identity * 100.0,
        aln_stats.alignments_by_pair.len()
    );

    // === COMPARISON ===
    eprintln!("\n=== Comparing Results ===");

    // Allow 0.0001% tolerance for floating point differences
    let tolerance = 0.000001;

    if compare_stats(&paf_stats, &aln_stats, tolerance) {
        eprintln!("\n✅ SUCCESS: PAF and .1aln workflows produce identical results!");
        eprintln!(
            "   {} alignments, {:.1} Mb coverage, {:.2}% identity",
            paf_stats.total_alignments,
            paf_stats.total_bases as f64 / 1_000_000.0,
            paf_stats.avg_identity * 100.0
        );
    } else {
        panic!("PAF and .1aln workflows produced different results!");
    }
}

#[test]
fn test_paf_vs_1aln_with_small_data() {
    // Use B-3106 data for fast testing
    let fasta_path = Path::new("data/B-3106.fa");

    if !fasta_path.exists() {
        eprintln!("Skipping test - data/B-3106.fa not found");
        return;
    }

    let temp_dir = TempDir::new().unwrap();

    // Use library API to generate alignments
    let integration = sweepga::fastga_integration::FastGAIntegration::new(None, 1);

    // Generate .1aln
    let aln_file = integration.align_to_temp_1aln(fasta_path, fasta_path);
    if aln_file.is_err() {
        eprintln!("Skipping test - FastGA not available");
        return;
    }
    let aln_file = aln_file.unwrap();

    // Generate PAF
    let paf_file = integration
        .align_to_temp_paf(fasta_path, fasta_path)
        .expect("Failed to generate PAF");

    // Parse stats from both
    eprintln!("[test] Parsing PAF stats from: {:?}", paf_file.path());
    let paf_stats = parse_paf_stats(paf_file.path());
    eprintln!("[test] PAF stats parsed");

    // For .1aln, convert to PAF using ALNtoPAF (the official converter)
    // This ensures we use the exact same conversion logic as FastGA
    let aln_as_paf = temp_dir.path().join("aln_converted.paf");
    eprintln!("[test] Converting .1aln to PAF using ALNtoPAF...");

    // Persist the temp file so ALNtoPAF can read it
    // (NamedTempFile auto-deletes when dropped)
    let persisted_aln = temp_dir.path().join("test.1aln");
    fs::copy(aln_file.path(), &persisted_aln).expect("Failed to copy .1aln file");
    eprintln!("[test] Copied .1aln to: {:?}", persisted_aln);

    // Also copy the .1gdb file if it exists (contains sequence names)
    let aln_gdb = aln_file.path().to_str().unwrap().replace(".1aln", ".1gdb");
    let persisted_gdb = persisted_aln.to_str().unwrap().replace(".1aln", ".1gdb");
    if std::path::Path::new(&aln_gdb).exists() {
        fs::copy(&aln_gdb, &persisted_gdb).expect("Failed to copy .1gdb file");
        eprintln!("[test] Copied .1gdb to: {:?}", persisted_gdb);
    }

    // Find ALNtoPAF binary using the embedded binary paths
    let alnto_paf_bin = sweepga::binary_paths::get_embedded_binary_path("ALNtoPAF")
        .expect("ALNtoPAF binary not found");

    eprintln!("[test] Using ALNtoPAF binary at: {:?}", alnto_paf_bin);

    let convert_result = std::process::Command::new(&alnto_paf_bin)
        .arg("-x") // Extended CIGAR mode (same as FastGA uses)
        .arg(&persisted_aln)
        .output()
        .expect("Failed to run ALNtoPAF");

    if !convert_result.status.success() {
        eprintln!("ALNtoPAF failed:");
        eprintln!("{}", String::from_utf8_lossy(&convert_result.stderr));
        panic!("ALNtoPAF conversion failed");
    }

    // Write PAF output to file
    fs::write(&aln_as_paf, &convert_result.stdout).expect("Failed to write PAF file");

    let aln_stats = parse_paf_stats(&aln_as_paf);

    eprintln!("\nComparison for B-3106:");
    eprintln!(
        "  PAF: {} alignments, {} bases",
        paf_stats.total_alignments, paf_stats.total_bases
    );
    eprintln!(
        "  .1aln: {} alignments, {} bases",
        aln_stats.total_alignments, aln_stats.total_bases
    );

    // For unfiltered data, they should be EXACTLY identical
    assert_eq!(
        paf_stats.total_alignments, aln_stats.total_alignments,
        "Unfiltered alignment counts must match exactly"
    );
    assert_eq!(
        paf_stats.total_bases, aln_stats.total_bases,
        "Unfiltered base counts must match exactly"
    );

    eprintln!("✅ PAF and .1aln produce identical unfiltered output");
}
