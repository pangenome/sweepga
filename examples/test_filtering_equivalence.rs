/// Test that .1aln and PAF produce identical FILTERING results
/// (even if coordinates differ slightly)
use anyhow::Result;
use std::collections::HashSet;

fn main() -> Result<()> {
    println!("=== Testing Filtering Equivalence ===\n");

    // Use existing test file
    let test_1aln = "/tmp/test_format_reading.1aln";
    let test_paf = "/tmp/test_format_reading.paf";

    if !std::path::Path::new(test_1aln).exists() {
        eprintln!("Test file not found. Run the format_reading test first.");
        return Ok(());
    }

    // Convert .1aln to PAF for comparison
    let alnto_paf = find_binary("ALNtoPAF").expect("ALNtoPAF not found");
    let paf_output = std::process::Command::new(&alnto_paf)
        .arg("-T1")
        .arg("/tmp/test_format_reading")
        .output()?;

    std::fs::write(test_paf, &paf_output.stdout)?;

    println!("Input files:");
    println!("  .1aln: {test_1aln}");
    println!("  PAF:   {test_paf}\n");

    // Apply 1:1 filtering to both
    use sweepga::paf_filter::{FilterConfig, FilterMode, ScoringFunction};
    use sweepga::unified_filter::filter_file;

    let config = FilterConfig {
        chain_gap: 0,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::OneToOne,
        mapping_max_per_query: Some(1),
        mapping_max_per_target: Some(1),
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::ManyToMany,
        scaffold_max_per_query: None,
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: 0,
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.95,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    println!("Applying 1:1 filtering to both formats...\n");

    // Filter .1aln
    let filtered_1aln = tempfile::NamedTempFile::with_suffix(".1aln")?;
    filter_file(test_1aln, filtered_1aln.path(), &config, false, false)?;

    // Filter PAF
    let filtered_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    filter_file(test_paf, filtered_paf.path(), &config, false, false)?;

    // Read filtered results
    let (aln_filtered, _) = sweepga::unified_filter::extract_1aln_metadata(filtered_1aln.path())?;
    let (paf_filtered, _) = sweepga::paf_filter::extract_metadata(filtered_paf.path())?;

    println!("Filtered results:");
    println!("  .1aln: {} alignments", aln_filtered.len());
    println!("  PAF:   {} alignments\n", paf_filtered.len());

    // Compare by alignment signature (ignoring exact coordinates)
    // Signature: (query_name, target_name, strand, ~query_span, ~target_span)

    fn signature(r: &sweepga::paf_filter::RecordMeta) -> (String, String, char, u64, u64) {
        (
            r.query_name.clone(),
            r.target_name.clone(),
            r.strand,
            r.query_end - r.query_start,   // query span
            r.target_end - r.target_start, // target span
        )
    }

    let aln_sigs: HashSet<_> = aln_filtered.iter().map(signature).collect();
    let paf_sigs: HashSet<_> = paf_filtered.iter().map(signature).collect();

    let in_aln_not_paf: Vec<_> = aln_sigs.difference(&paf_sigs).collect();
    let in_paf_not_aln: Vec<_> = paf_sigs.difference(&aln_sigs).collect();

    println!("=== Results ===");
    if in_aln_not_paf.is_empty() && in_paf_not_aln.is_empty() {
        println!("✅ SUCCESS: Filtering produced IDENTICAL results!");
        println!(
            "   Both formats kept the same {} alignments",
            aln_filtered.len()
        );
        println!("   (Coordinate differences don't affect filtering decisions)");
    } else {
        println!("❌ FAILURE: Filtering produced different results");
        println!("   In .1aln but not PAF: {}", in_aln_not_paf.len());
        println!("   In PAF but not .1aln: {}", in_paf_not_aln.len());

        if !in_aln_not_paf.is_empty() {
            println!("\nFirst 3 in .1aln but not PAF:");
            for sig in in_aln_not_paf.iter().take(3) {
                println!("  {sig:?}");
            }
        }

        if !in_paf_not_aln.is_empty() {
            println!("\nFirst 3 in PAF but not .1aln:");
            for sig in in_paf_not_aln.iter().take(3) {
                println!("  {sig:?}");
            }
        }
    }

    Ok(())
}

fn find_binary(name: &str) -> Option<String> {
    use std::env;
    use std::fs;

    let project_root = env::current_dir().ok()?;

    for build_type in &["debug", "release"] {
        let build_dir = project_root.join(format!("target/{build_type}/build"));
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
