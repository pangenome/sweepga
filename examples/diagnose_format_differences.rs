/// Diagnostic tool to understand why .1aln and PAF produce different filtering results
///
/// Strategy:
/// 1. Read .1aln directly → RecordMeta (our internal representation)
/// 2. Convert .1aln → PAF using ALNtoPAF → RecordMeta
/// 3. Compare field-by-field to find systematic differences
use anyhow::Result;
use std::path::Path;

fn main() -> Result<()> {
    // Use existing test file from tests
    let test_1aln = Path::new("/tmp/test_format_reading.1aln");

    if !test_1aln.exists() {
        eprintln!("Generating test .1aln file...");

        // Find FastGA binary
        let fastga = find_binary("FastGA").expect("FastGA not found");

        let output = std::process::Command::new(&fastga)
            .arg("-1:/tmp/test_format_reading.1aln")
            .arg("-T1")
            .arg("B-3106.fa")
            .arg("B-3106.fa")
            .current_dir("data")
            .output()?;

        if !output.status.success() {
            eprintln!("FastGA failed: {}", String::from_utf8_lossy(&output.stderr));
            return Ok(());
        }
    }

    println!("=== Diagnostic: .1aln vs PAF Internal Representation ===\n");

    // Step 1: Read .1aln directly into RecordMeta
    println!("Step 1: Reading .1aln directly...");
    let (aln_records, _) = sweepga::unified_filter::extract_1aln_metadata(test_1aln)?;
    println!("  Read {} records from .1aln\n", aln_records.len());

    // Step 2: Convert .1aln → PAF using ALNtoPAF
    println!("Step 2: Converting .1aln to PAF via ALNtoPAF...");
    let alnto_paf = find_binary("ALNtoPAF").expect("ALNtoPAF not found");

    let paf_output = std::process::Command::new(&alnto_paf)
        .arg("-T1")
        .arg("/tmp/test_format_reading") // ALNtoPAF auto-appends .1aln
        .output()?;

    if !paf_output.status.success() {
        eprintln!(
            "ALNtoPAF failed: {}",
            String::from_utf8_lossy(&paf_output.stderr)
        );
        return Ok(());
    }

    let paf_path = "/tmp/diagnose_test.paf";
    std::fs::write(paf_path, &paf_output.stdout)?;
    println!("  Converted to PAF ({} bytes)\n", paf_output.stdout.len());

    // Step 3: Read PAF into RecordMeta
    println!("Step 3: Reading PAF into RecordMeta...");
    let (paf_records, _) = sweepga::paf_filter::extract_metadata(paf_path)?;
    println!("  Read {} records from PAF\n", paf_records.len());

    // Step 4: Compare
    println!("=== Comparing Records ===\n");

    if aln_records.len() != paf_records.len() {
        println!(
            "❌ Record count mismatch: .1aln={}, PAF={}",
            aln_records.len(),
            paf_records.len()
        );
        return Ok(());
    }

    // Sort both by (query, target, qstart, qend)
    let mut aln_sorted = aln_records;
    let mut paf_sorted = paf_records;

    let sort_key = |r: &sweepga::paf_filter::RecordMeta| {
        (
            r.query_name.clone(),
            r.target_name.clone(),
            r.query_start,
            r.query_end,
        )
    };

    aln_sorted.sort_by_key(|r| sort_key(r));
    paf_sorted.sort_by_key(|r| sort_key(r));

    println!("First 10 records - detailed comparison:\n");

    let mut diffs = Vec::new();

    for i in 0..aln_sorted.len().min(10) {
        let aln = &aln_sorted[i];
        let paf = &paf_sorted[i];

        println!("--- Record {} ---", i + 1);
        println!(
            ".1aln: {}:{}-{} → {}:{}-{} {} bl={} id={:.4} m={}",
            aln.query_name,
            aln.query_start,
            aln.query_end,
            aln.target_name,
            aln.target_start,
            aln.target_end,
            aln.strand,
            aln.block_length,
            aln.identity,
            aln.matches
        );
        println!(
            "PAF:   {}:{}-{} → {}:{}-{} {} bl={} id={:.4} m={}",
            paf.query_name,
            paf.query_start,
            paf.query_end,
            paf.target_name,
            paf.target_start,
            paf.target_end,
            paf.strand,
            paf.block_length,
            paf.identity,
            paf.matches
        );

        // Check each field
        let mut field_diffs = Vec::new();

        if aln.query_name != paf.query_name {
            field_diffs.push(format!(
                "query_name: '{}' vs '{}'",
                aln.query_name, paf.query_name
            ));
        }
        if aln.target_name != paf.target_name {
            field_diffs.push(format!(
                "target_name: '{}' vs '{}'",
                aln.target_name, paf.target_name
            ));
        }
        if aln.query_start != paf.query_start {
            field_diffs.push(format!(
                "query_start: {} vs {}",
                aln.query_start, paf.query_start
            ));
        }
        if aln.query_end != paf.query_end {
            field_diffs.push(format!("query_end: {} vs {}", aln.query_end, paf.query_end));
        }
        if aln.target_start != paf.target_start {
            field_diffs.push(format!(
                "target_start: {} vs {}",
                aln.target_start, paf.target_start
            ));
        }
        if aln.target_end != paf.target_end {
            field_diffs.push(format!(
                "target_end: {} vs {}",
                aln.target_end, paf.target_end
            ));
        }
        if aln.strand != paf.strand {
            field_diffs.push(format!("strand: {} vs {}", aln.strand, paf.strand));
        }
        if aln.block_length != paf.block_length {
            field_diffs.push(format!(
                "block_length: {} vs {}",
                aln.block_length, paf.block_length
            ));
        }
        if aln.matches != paf.matches {
            field_diffs.push(format!("matches: {} vs {}", aln.matches, paf.matches));
        }
        if (aln.identity - paf.identity).abs() > 0.0001 {
            field_diffs.push(format!(
                "identity: {:.6} vs {:.6} (diff: {:.6})",
                aln.identity,
                paf.identity,
                (aln.identity - paf.identity).abs()
            ));
        }

        if field_diffs.is_empty() {
            println!("✅ IDENTICAL\n");
        } else {
            println!("❌ DIFFERENCES:");
            for diff in &field_diffs {
                println!("   {}", diff);
            }
            println!();
            diffs.push((i, field_diffs));
        }
    }

    // Summary statistics
    println!("\n=== Summary ===");

    let mut total_mismatches = 0;
    let mut coord_diffs = 0;
    let mut identity_diffs = 0;
    let mut block_len_diffs = 0;
    let mut matches_diffs = 0;

    for i in 0..aln_sorted.len() {
        let aln = &aln_sorted[i];
        let paf = &paf_sorted[i];

        let mut has_diff = false;

        if aln.query_start != paf.query_start
            || aln.query_end != paf.query_end
            || aln.target_start != paf.target_start
            || aln.target_end != paf.target_end
        {
            coord_diffs += 1;
            has_diff = true;
        }

        if (aln.identity - paf.identity).abs() > 0.0001 {
            identity_diffs += 1;
            has_diff = true;
        }

        if aln.block_length != paf.block_length {
            block_len_diffs += 1;
            has_diff = true;
        }

        if aln.matches != paf.matches {
            matches_diffs += 1;
            has_diff = true;
        }

        if has_diff {
            total_mismatches += 1;
        }
    }

    println!("Total records: {}", aln_sorted.len());
    println!(
        "Records with differences: {} ({:.1}%)",
        total_mismatches,
        100.0 * total_mismatches as f64 / aln_sorted.len() as f64
    );
    println!("\nDifference breakdown:");
    println!("  Coordinate differences: {}", coord_diffs);
    println!("  Identity differences: {}", identity_diffs);
    println!("  Block length differences: {}", block_len_diffs);
    println!("  Matches differences: {}", matches_diffs);

    Ok(())
}

fn find_binary(name: &str) -> Option<String> {
    use std::env;
    use std::fs;

    let project_root = env::current_dir().ok()?;

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
