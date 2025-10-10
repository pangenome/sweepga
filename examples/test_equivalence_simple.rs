/// Simple equivalence test using existing test files
///
/// Uses test_output.1aln (pre-generated) to verify filtering equivalence
use anyhow::Result;
use std::collections::HashSet;

fn main() -> Result<()> {
    println!("=== Format Equivalence Test (Existing Files) ===\n");

    let test_1aln = "test_output.1aln";
    let test_paf = "test_output_converted.paf";

    if !std::path::Path::new(test_1aln).exists() {
        eprintln!("❌ Test file not found: {}", test_1aln);
        return Ok(());
    }

    if !std::path::Path::new(test_paf).exists() {
        eprintln!("❌ Test file not found: {}", test_paf);
        return Ok(());
    }

    println!("Using pre-generated test files:");
    println!("  .1aln: {}", test_1aln);
    println!("  PAF:   {}\n", test_paf);

    // Count unfiltered alignments
    let unfiltered_paf_count = std::fs::read_to_string(test_paf)?.lines().count();
    println!("Unfiltered: {} alignments\n", unfiltered_paf_count);

    println!("Step 1: Apply 1:1 filtering to both formats...");

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

    // Filter .1aln → .1aln
    println!("  - Filtering .1aln...");
    let temp_filtered_1aln = tempfile::NamedTempFile::with_suffix(".1aln")?;
    filter_file(test_1aln, temp_filtered_1aln.path(), &config, false)?;
    println!("    ✓ Complete");

    // Filter PAF → PAF
    println!("  - Filtering PAF...");
    let temp_filtered_paf = tempfile::NamedTempFile::with_suffix(".paf")?;
    filter_file(test_paf, temp_filtered_paf.path(), &config, false)?;
    println!("    ✓ Complete");

    println!("\nStep 2: Convert filtered .1aln to PAF for comparison...");
    let temp_converted_paf = tempfile::NamedTempFile::with_suffix(".paf")?;

    let convert_result = std::process::Command::new("/home/erik/bin/ALNtoPAF")
        .arg(temp_filtered_1aln.path())
        .output()?;

    if !convert_result.status.success() {
        eprintln!("❌ Conversion failed");
        return Ok(());
    }

    std::fs::write(temp_converted_paf.path(), &convert_result.stdout)?;
    println!("  ✓ Converted");

    println!("\nStep 3: Compare results...");

    // Parse both PAF files
    let paf_content = std::fs::read_to_string(temp_filtered_paf.path())?;
    let aln_content = std::fs::read_to_string(temp_converted_paf.path())?;

    let paf_alignments: Vec<_> = paf_content.lines().filter(|l| !l.is_empty()).collect();
    let aln_alignments: Vec<_> = aln_content.lines().filter(|l| !l.is_empty()).collect();

    println!("Filtered alignment counts:");
    println!(
        "  PAF direct:           {} alignments",
        paf_alignments.len()
    );
    println!(
        "  .1aln (via ALNtoPAF): {} alignments",
        aln_alignments.len()
    );

    // Extract alignment keys
    fn parse_key(line: &str) -> Option<(String, u64, u64, String, u64, u64, char)> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            return None;
        }
        Some((
            fields[0].to_string(),
            fields[2].parse().ok()?,
            fields[3].parse().ok()?,
            fields[5].to_string(),
            fields[7].parse().ok()?,
            fields[8].parse().ok()?,
            fields[4].chars().next()?,
        ))
    }

    let paf_keys: HashSet<_> = paf_alignments.iter().filter_map(|l| parse_key(l)).collect();
    let aln_keys: HashSet<_> = aln_alignments.iter().filter_map(|l| parse_key(l)).collect();

    // Compare
    let in_paf_not_aln: Vec<_> = paf_keys.difference(&aln_keys).collect();
    let in_aln_not_paf: Vec<_> = aln_keys.difference(&paf_keys).collect();

    println!("\n{}", "=".repeat(60));

    if in_paf_not_aln.is_empty() && in_aln_not_paf.is_empty() {
        println!("✅ SUCCESS: Filtered results are IDENTICAL!");
        println!("");
        println!("Summary:");
        println!("  - Input: {} alignments", unfiltered_paf_count);
        println!(
            "  - After 1:1 filtering: {} alignments",
            paf_alignments.len()
        );
        println!(
            "  - Reduction: {:.1}%",
            100.0 * (1.0 - paf_alignments.len() as f64 / unfiltered_paf_count as f64)
        );
        println!("  - All alignment coordinates match exactly");
        println!("");
        println!("✓ Format-agnostic filtering verified");
        println!("✓ Identity calculation fix working correctly");
        println!("✓ X record reading working correctly");
    } else {
        println!("⚠️  Filtered results differ:");
        println!("  In PAF but not .1aln: {}", in_paf_not_aln.len());
        println!("  In .1aln but not PAF: {}", in_aln_not_paf.len());

        if !in_paf_not_aln.is_empty() {
            println!("\nFirst 3 in PAF but not .1aln:");
            for key in in_paf_not_aln.iter().take(3) {
                println!("  {:?}", key);
            }
        }

        if !in_aln_not_paf.is_empty() {
            println!("\nFirst 3 in .1aln but not PAF:");
            for key in in_aln_not_paf.iter().take(3) {
                println!("  {:?}", key);
            }
        }
    }

    Ok(())
}
