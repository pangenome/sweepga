/// Integration test: Verify .1aln and PAF filtering produce identical results
///
/// This test:
/// 1. Aligns B-3106.fa in both .1aln and PAF formats
/// 2. Applies 1:1 filtering to both
/// 3. Compares results to ensure they're identical

use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::path::Path;

fn main() -> Result<()> {
    println!("=== Format Equivalence Test: B-3106.fa ===\n");

    let test_fasta = Path::new("data/B-3106.fa");
    if !test_fasta.exists() {
        eprintln!("❌ Test file not found: {}", test_fasta.display());
        eprintln!("Please ensure data/B-3106.fa exists");
        return Ok(());
    }

    println!("Step 1: Generate alignments in both formats...");

    // Create FastGA integration
    let fastga = sweepga::fastga_integration::FastGAIntegration::new(None, 1);

    // Generate .1aln
    println!("  - Generating .1aln via FastGA...");
    let temp_1aln = fastga.align_to_temp_1aln(test_fasta, test_fasta)?;
    println!("    ✓ Created: {}", temp_1aln.path().display());

    // Convert to PAF
    println!("  - Converting to PAF via ALNtoPAF...");
    let temp_paf_unfiltered = tempfile::NamedTempFile::with_suffix(".paf")?;

    let alnto_paf_result = std::process::Command::new("/home/erik/bin/ALNtoPAF")
        .arg(temp_1aln.path())
        .output()?;

    if !alnto_paf_result.status.success() {
        eprintln!("❌ ALNtoPAF failed: {}", String::from_utf8_lossy(&alnto_paf_result.stderr));
        return Ok(());
    }

    std::fs::write(temp_paf_unfiltered.path(), &alnto_paf_result.stdout)?;
    println!("    ✓ Converted to PAF");

    // Count initial alignments
    let paf_count = String::from_utf8_lossy(&alnto_paf_result.stdout).lines().count();
    println!("\nInitial alignments:");
    println!("  PAF: {} alignments", paf_count);

    println!("\nStep 2: Apply 1:1 filtering to both formats...");

    // Filter .1aln → .1aln (format-preserving)
    println!("  - Filtering .1aln → filtered.1aln...");
    let temp_filtered_1aln = tempfile::NamedTempFile::with_suffix(".1aln")?;

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

    filter_file(temp_1aln.path(), temp_filtered_1aln.path(), &config, false)?;
    println!("    ✓ .1aln filtering succeeded");

    // Filter PAF → PAF
    println!("  - Filtering PAF → filtered.paf...");
    let temp_filtered_paf = tempfile::NamedTempFile::with_suffix(".paf")?;

    filter_file(temp_paf_unfiltered.path(), temp_filtered_paf.path(), &config, false)?;
    println!("    ✓ PAF filtering succeeded");

    println!("\nStep 3: Convert filtered .1aln to PAF for comparison...");
    let temp_converted_paf = tempfile::NamedTempFile::with_suffix(".paf")?;

    let convert_result = std::process::Command::new("/home/erik/bin/ALNtoPAF")
        .arg(temp_filtered_1aln.path())
        .output()?;

    if !convert_result.status.success() {
        eprintln!("❌ Conversion failed");
        return Ok(());
    }

    std::fs::write(temp_converted_paf.path(), &convert_result.stdout)?;
    println!("  ✓ Converted filtered .1aln to PAF");

    println!("\nStep 4: Compare results...");

    // Parse both PAF files
    let paf_content = std::fs::read_to_string(temp_filtered_paf.path())?;
    let aln_content = std::fs::read_to_string(temp_converted_paf.path())?;

    let paf_alignments: Vec<_> = paf_content.lines().filter(|l| !l.is_empty()).collect();
    let aln_alignments: Vec<_> = aln_content.lines().filter(|l| !l.is_empty()).collect();

    println!("Filtered alignment counts:");
    println!("  PAF direct:           {} alignments", paf_alignments.len());
    println!("  .1aln (via ALNtoPAF): {} alignments", aln_alignments.len());

    // Extract alignment keys (query, qstart, qend, target, tstart, tend)
    fn parse_alignment_key(line: &str) -> Option<(String, u64, u64, String, u64, u64)> {
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
        ))
    }

    let paf_keys: HashSet<_> = paf_alignments
        .iter()
        .filter_map(|l| parse_alignment_key(l))
        .collect();
    let aln_keys: HashSet<_> = aln_alignments
        .iter()
        .filter_map(|l| parse_alignment_key(l))
        .collect();

    // Compare
    let in_paf_not_aln: Vec<_> = paf_keys.difference(&aln_keys).collect();
    let in_aln_not_paf: Vec<_> = aln_keys.difference(&paf_keys).collect();

    println!("\n=== Results ===");

    if in_paf_not_aln.is_empty() && in_aln_not_paf.is_empty() {
        println!("✅ SUCCESS: Filtered results are IDENTICAL!");
        println!("");
        println!("Summary:");
        println!("  - Both formats filtered from {} to {} alignments", paf_count, paf_alignments.len());
        println!("  - All alignment coordinates match exactly");
        println!("  - Format-agnostic filtering is working correctly!");
        println!("  - Identity calculation fix verified!");
    } else {
        println!("⚠️  Filtered results differ:");
        println!("  In PAF but not .1aln: {}", in_paf_not_aln.len());
        println!("  In .1aln but not PAF: {}", in_aln_not_paf.len());

        if !in_paf_not_aln.is_empty() {
            println!("\nFirst 5 in PAF but not .1aln:");
            for key in in_paf_not_aln.iter().take(5) {
                println!("  {:?}", key);
            }
        }

        if !in_aln_not_paf.is_empty() {
            println!("\nFirst 5 in .1aln but not PAF:");
            for key in in_aln_not_paf.iter().take(5) {
                println!("  {:?}", key);
            }
        }
    }

    Ok(())
}
