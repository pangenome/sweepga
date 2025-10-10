#![allow(clippy::uninlined_format_args)]
/// Debug test to verify PAF and .1aln produce identical internal representations
///
/// This test:
/// 1. Generates alignments from B-3106.fa in both PAF and .1aln formats
/// 2. Reads them into RecordMeta structures
/// 3. Sorts and compares to verify they're identical
use std::path::Path;
use std::process::Command;

#[test]
fn test_paf_vs_1aln_internal_representation() {
    let fasta_path = Path::new("data/B-3106.fa");

    if !fasta_path.exists() {
        eprintln!("Skipping test - data/B-3106.fa not found");
        return;
    }

    println!("\n=== Testing PAF vs .1aln Internal Representation ===\n");

    // Find FastGA and ALNtoPAF binaries
    let fastga = find_binary("FastGA");
    let alnto_paf = find_binary("ALNtoPAF");

    if fastga.is_none() || alnto_paf.is_none() {
        eprintln!("Skipping test - FastGA or ALNtoPAF not found");
        return;
    }

    let fastga = fastga.unwrap();
    let alnto_paf = alnto_paf.unwrap();

    println!("Using FastGA: {}", fastga);
    println!("Using ALNtoPAF: {}\n", alnto_paf);

    // Generate .1aln
    println!("Generating .1aln from B-3106.fa...");
    let aln_path = "/tmp/test_format_reading.1aln";

    // Use absolute paths for data directory
    let data_dir = std::env::current_dir().unwrap().join("data");

    let output = Command::new(&fastga)
        .arg(format!("-1:{}", aln_path))
        .arg("-T1")
        .arg("B-3106.fa")
        .arg("B-3106.fa")
        .current_dir(&data_dir)
        .output();

    match output {
        Ok(output) if output.status.success() && Path::new(aln_path).exists() => {
            println!(
                "Generated: {} ({} bytes)",
                aln_path,
                std::fs::metadata(aln_path).unwrap().len()
            );
        }
        Ok(output) => {
            eprintln!("FastGA failed:");
            eprintln!("  stdout: {}", String::from_utf8_lossy(&output.stdout));
            eprintln!("  stderr: {}", String::from_utf8_lossy(&output.stderr));
            eprintln!("Skipping test - FastGA failed");
            return;
        }
        Err(e) => {
            eprintln!("Failed to execute FastGA: {}", e);
            eprintln!("Skipping test - FastGA execution error");
            return;
        }
    }

    // Generate PAF from .1aln
    println!("Converting .1aln to PAF...");
    let paf_path = "/tmp/test_format_reading.paf";
    let output = Command::new(&alnto_paf)
        .arg("-T1")
        .arg("/tmp/test_format_reading")
        .output();

    match output {
        Ok(output) if output.status.success() => {
            std::fs::write(paf_path, &output.stdout).expect("Failed to write PAF");
            println!(
                "Generated: {} ({} lines)\n",
                paf_path,
                String::from_utf8_lossy(&output.stdout).lines().count()
            );
        }
        Ok(output) => {
            eprintln!("ALNtoPAF failed:");
            eprintln!("  stdout: {}", String::from_utf8_lossy(&output.stdout));
            eprintln!("  stderr: {}", String::from_utf8_lossy(&output.stderr));
            eprintln!("Skipping test - ALNtoPAF failed");
            return;
        }
        Err(e) => {
            eprintln!("Failed to execute ALNtoPAF: {}", e);
            eprintln!("Skipping test - ALNtoPAF execution error");
            return;
        }
    }

    // Read PAF into RecordMeta
    println!("=== Reading PAF ===");
    let paf_records = read_paf_to_records(Path::new(paf_path));
    println!("Read {} records from PAF\n", paf_records.len());

    // Read .1aln into RecordMeta
    println!("=== Reading .1aln ===");
    let aln_records = read_1aln_to_records(Path::new(aln_path));
    println!("Read {} records from .1aln\n", aln_records.len());

    // Sort both by a consistent key (query, target, qstart, qend)
    let mut paf_sorted = paf_records;
    let mut aln_sorted = aln_records;

    paf_sorted.sort_by(|a, b| {
        (&a.query_name, &a.target_name, a.query_start, a.query_end).cmp(&(
            &b.query_name,
            &b.target_name,
            b.query_start,
            b.query_end,
        ))
    });

    aln_sorted.sort_by(|a, b| {
        (&a.query_name, &a.target_name, a.query_start, a.query_end).cmp(&(
            &b.query_name,
            &b.target_name,
            b.query_start,
            b.query_end,
        ))
    });

    // Compare counts
    println!("=== Comparison ===");
    println!("PAF records: {}", paf_sorted.len());
    println!(".1aln records: {}", aln_sorted.len());

    if paf_sorted.len() != aln_sorted.len() {
        eprintln!("\n❌ Record count mismatch!");
        eprintln!("  PAF: {}", paf_sorted.len());
        eprintln!("  .1aln: {}", aln_sorted.len());
        panic!("Record counts differ");
    }

    println!("✅ Record counts match: {}\n", paf_sorted.len());

    // Compare first 5 records in detail
    println!("=== First 5 Records (Detailed Comparison) ===\n");
    let compare_count = paf_sorted.len().min(5);

    for i in 0..compare_count {
        let paf = &paf_sorted[i];
        let aln = &aln_sorted[i];

        println!("--- Record {} ---", i + 1);
        println!(
            "PAF: {} → {} {} {}-{} → {}-{}",
            paf.query_name,
            paf.target_name,
            paf.strand,
            paf.query_start,
            paf.query_end,
            paf.target_start,
            paf.target_end
        );
        println!(
            "     matches={}, block_len={}, identity={:.4}",
            paf.matches, paf.block_length, paf.identity
        );

        println!(
            ".1aln: {} → {} {} {}-{} → {}-{}",
            aln.query_name,
            aln.target_name,
            aln.strand,
            aln.query_start,
            aln.query_end,
            aln.target_start,
            aln.target_end
        );
        println!(
            "       matches={}, block_len={}, identity={:.4}",
            aln.matches, aln.block_length, aln.identity
        );

        // Check if identical
        let matches = compare_records(paf, aln, i);
        if matches {
            println!("✅ IDENTICAL\n");
        } else {
            println!("❌ DIFFERENT\n");
        }
    }

    // Quick check all records
    println!("=== Checking All {} Records ===", paf_sorted.len());
    let mut mismatch_count = 0;
    for i in 0..paf_sorted.len() {
        if !compare_records(&paf_sorted[i], &aln_sorted[i], i) {
            mismatch_count += 1;
            if mismatch_count <= 3 {
                println!("Mismatch at record {}", i);
            }
        }
    }

    if mismatch_count == 0 {
        println!(
            "\n✅ SUCCESS: All {} records are IDENTICAL!",
            paf_sorted.len()
        );
        println!("   PAF and .1aln produce the same internal representation");
    } else {
        println!(
            "\n❌ FAILURE: {} / {} records differ",
            mismatch_count,
            paf_sorted.len()
        );
        println!(
            "   Percentage matching: {:.2}%",
            (paf_sorted.len() - mismatch_count) as f64 / paf_sorted.len() as f64 * 100.0
        );
        panic!("PAF and .1aln internal representations differ!");
    }
}

/// Compare two RecordMeta structures
fn compare_records(
    paf: &sweepga::paf_filter::RecordMeta,
    aln: &sweepga::paf_filter::RecordMeta,
    _idx: usize,
) -> bool {
    paf.query_name == aln.query_name
        && paf.target_name == aln.target_name
        && paf.query_start == aln.query_start
        && paf.query_end == aln.query_end
        && paf.target_start == aln.target_start
        && paf.target_end == aln.target_end
        && paf.strand == aln.strand
        && paf.matches == aln.matches
        && paf.block_length == aln.block_length
        && (paf.identity - aln.identity).abs() < 0.0001 // Allow tiny floating point diff
}

/// Read PAF file into RecordMeta structures (using existing PAF filter logic)
fn read_paf_to_records(paf_path: &Path) -> Vec<sweepga::paf_filter::RecordMeta> {
    // Use the public extract_metadata function
    let (metadata, _) =
        sweepga::paf_filter::extract_metadata(paf_path).expect("Failed to read PAF");
    metadata
}

/// Read .1aln file into RecordMeta structures (using unified filter logic)
fn read_1aln_to_records(aln_path: &Path) -> Vec<sweepga::paf_filter::RecordMeta> {
    let (metadata, _) =
        sweepga::unified_filter::extract_1aln_metadata(aln_path).expect("Failed to read .1aln");

    metadata
}

/// Find a binary in target/debug or target/release build directories
fn find_binary(name: &str) -> Option<String> {
    use std::env;
    use std::fs;

    // Get absolute path to project root
    let project_root = env::current_dir().ok()?;

    // Try both debug and release
    for build_type in &["debug", "release"] {
        let build_dir = project_root.join(format!("target/{}/build", build_type));
        if let Ok(entries) = fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                let out_dir = entry.path().join("out");
                let binary_path = out_dir.join(name);
                if binary_path.exists() {
                    return Some(binary_path.to_string_lossy().to_string());
                }
            }
        }
    }

    None
}
