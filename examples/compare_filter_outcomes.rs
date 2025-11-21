/// Direct comparison of filtering outcomes from .1aln vs PAF inputs
/// Compares the RecordMeta BEFORE writing to files
use anyhow::Result;
use std::collections::HashSet;

fn main() -> Result<()> {
    println!("=== Comparing Filter Outcomes: .1aln vs PAF ===\n");

    let test_1aln = "/tmp/test_format_reading.1aln";

    // Convert .1aln to PAF
    let alnto_paf = find_binary("ALNtoPAF").expect("ALNtoPAF not found");
    let paf_output = std::process::Command::new(&alnto_paf)
        .arg("-T1")
        .arg("/tmp/test_format_reading")
        .output()?;

    let test_paf = "/tmp/test_format_reading.paf";
    std::fs::write(test_paf, &paf_output.stdout)?;

    println!("Input files:");
    println!("  .1aln: {test_1aln}");
    println!("  PAF:   {test_paf}\n");

    // Read both into RecordMeta
    let (aln_records, _) = sweepga::unified_filter::extract_1aln_metadata(test_1aln)?;
    let (paf_records, _) = sweepga::paf_filter::extract_metadata(test_paf)?;

    println!("Read {} records from .1aln", aln_records.len());
    println!("Read {} records from PAF\n", paf_records.len());

    // Apply same filtering config to both
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};

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

    let filter = PafFilter::new(config);

    println!("Applying 1:1 filtering...\n");
    let aln_filtered = filter.apply_filters(aln_records)?;
    let paf_filtered = filter.apply_filters(paf_records)?;

    println!("Filtered results:");
    println!("  .1aln: {} alignments kept", aln_filtered.len());
    println!("  PAF:   {} alignments kept\n", paf_filtered.len());

    // Compare by signature (ignore exact coordinates)
    fn signature(r: &sweepga::paf_filter::RecordMeta) -> (String, String, char, u64, u64, u64) {
        (
            r.query_name.clone(),
            r.target_name.clone(),
            r.strand,
            r.query_end - r.query_start,   // query span
            r.target_end - r.target_start, // target span
            r.matches,                     // matches count
        )
    }

    let aln_sigs: HashSet<_> = aln_filtered.values().map(signature).collect();
    let paf_sigs: HashSet<_> = paf_filtered.values().map(signature).collect();

    let in_aln_not_paf: Vec<_> = aln_sigs.difference(&paf_sigs).collect();
    let in_paf_not_aln: Vec<_> = paf_sigs.difference(&aln_sigs).collect();

    println!("=== Results ===");
    if in_aln_not_paf.is_empty() && in_paf_not_aln.is_empty() {
        println!("✅ SUCCESS: Both formats produced IDENTICAL filtering outcomes!");
        println!("   {} alignments kept in both cases", aln_filtered.len());
        println!("   Filtering logic is format-agnostic ✓");
    } else {
        println!("❌ FAILURE: Filtering outcomes differ");
        println!("   In .1aln but not PAF: {}", in_aln_not_paf.len());
        println!("   In PAF but not .1aln: {}", in_paf_not_aln.len());

        if !in_aln_not_paf.is_empty() {
            println!("\nFirst 5 in .1aln but not PAF:");
            for sig in in_aln_not_paf.iter().take(5) {
                println!("  {sig:?}");
            }
        }

        if !in_paf_not_aln.is_empty() {
            println!("\nFirst 5 in PAF but not .1aln:");
            for sig in in_paf_not_aln.iter().take(5) {
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
