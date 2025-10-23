/// Round-trip .1aln filtering test
///
/// Ensures that filtering .1aln → .1aln preserves data integrity:
/// 1. Read input.1aln → RecordMeta
/// 2. Filter → kept ranks
/// 3. Write filtered.1aln
/// 4. Read filtered.1aln → RecordMeta
/// 5. Verify filtered records match expected
///
/// This catches any data corruption in the write/read cycle
use anyhow::Result;
use sweepga::paf_filter::{FilterConfig, FilterMode, ScoringFunction};
use sweepga::unified_filter;
use tempfile::TempDir;

#[path = "synthetic_genomes.rs"]
mod synthetic_genomes;

#[test]
fn test_1aln_roundtrip_preserves_data() -> Result<()> {
    use synthetic_genomes::{generate_base_sequence, mutate_sequence};

    // Generate test data
    let temp_dir = TempDir::new()?;
    let input_fasta = temp_dir.path().join("test.fa");

    // Create 2 related sequences (3kb each for speed)
    let base = generate_base_sequence(3000, 777);
    let seq1 = base.clone();
    let seq2 = mutate_sequence(&base, 30, 778); // 1% divergence

    std::fs::write(&input_fasta, format!(">seq1\n{seq1}\n>seq2\n{seq2}\n"))?;

    // Align to create .1aln
    let fastga = sweepga::fastga_integration::FastGAIntegration::new(None, 1);
    let aln_result = fastga.align_to_temp_1aln(&input_fasta, &input_fasta);

    if aln_result.is_err() {
        eprintln!(
            "Skipping test - FastGA alignment failed (expected when running from temp directories)"
        );
        return Ok(());
    }

    let temp_1aln = aln_result.unwrap();
    let input_path = temp_1aln.path();

    // Read original metadata
    let (original_meta, _) = unified_filter::extract_1aln_metadata(input_path)?;
    let original_count = original_meta.len();

    assert!(original_count > 0, "Should have alignments");

    // Filter with permissive config (should keep most records)
    let config = FilterConfig {
        chain_gap: 0,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
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

    let output_path = temp_dir.path().join("filtered.1aln");
    unified_filter::filter_file(input_path, &output_path, &config, false, true)?; // keep_self=true

    // Read filtered metadata
    let (filtered_meta, _) = unified_filter::extract_1aln_metadata(&output_path)?;
    let filtered_count = filtered_meta.len();

    assert_eq!(
        filtered_count, original_count,
        "Permissive filter should keep all records"
    );

    // Verify each record matches
    for (i, (orig, filt)) in original_meta.iter().zip(filtered_meta.iter()).enumerate() {
        assert_eq!(
            orig.query_name, filt.query_name,
            "Record {} query name mismatch",
            i
        );
        assert_eq!(
            orig.target_name, filt.target_name,
            "Record {} target name mismatch",
            i
        );
        assert_eq!(
            orig.query_start, filt.query_start,
            "Record {} query_start mismatch",
            i
        );
        assert_eq!(
            orig.query_end, filt.query_end,
            "Record {} query_end mismatch",
            i
        );
        assert_eq!(
            orig.target_start, filt.target_start,
            "Record {} target_start mismatch",
            i
        );
        assert_eq!(
            orig.target_end, filt.target_end,
            "Record {} target_end mismatch",
            i
        );
        assert_eq!(orig.strand, filt.strand, "Record {} strand mismatch", i);
        assert_eq!(orig.matches, filt.matches, "Record {} matches mismatch", i);

        // Allow small floating point differences in identity
        let identity_diff = (orig.identity - filt.identity).abs();
        assert!(
            identity_diff < 0.001,
            "Record {} identity diff too large: {} vs {} (diff: {})",
            i,
            orig.identity,
            filt.identity,
            identity_diff
        );
    }

    Ok(())
}

#[test]
fn test_1aln_roundtrip_with_filtering() -> Result<()> {
    use synthetic_genomes::{generate_base_sequence, mutate_sequence};

    // Test that filtering actually removes records and remaining ones are intact
    let temp_dir = TempDir::new()?;
    let input_fasta = temp_dir.path().join("test.fa");

    // Create 3 related sequences for more alignments (3kb each)
    let base = generate_base_sequence(3000, 888);
    let seq1 = base.clone();
    let seq2 = mutate_sequence(&base, 30, 889);
    let seq3 = mutate_sequence(&base, 60, 890);

    std::fs::write(
        &input_fasta,
        format!(">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n"),
    )?;

    let fastga = sweepga::fastga_integration::FastGAIntegration::new(None, 1);
    let aln_result = fastga.align_to_temp_1aln(&input_fasta, &input_fasta);

    if aln_result.is_err() {
        eprintln!(
            "Skipping test - FastGA alignment failed (expected when running from temp directories)"
        );
        return Ok(());
    }

    let temp_1aln = aln_result.unwrap();
    let input_path = temp_1aln.path();

    let (original_meta, _) = unified_filter::extract_1aln_metadata(input_path)?;

    // Apply 1:1 filtering (should reduce records)
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

    let output_path = temp_dir.path().join("filtered.1aln");
    unified_filter::filter_file(input_path, &output_path, &config, false, true)?; // keep_self=true

    let (filtered_meta, _) = unified_filter::extract_1aln_metadata(&output_path)?;

    // Should have filtered some out
    assert!(
        filtered_meta.len() < original_meta.len(),
        "1:1 filter should reduce record count: {} >= {}",
        filtered_meta.len(),
        original_meta.len()
    );

    // Every filtered record should have a corresponding original with matching coordinates
    for filt in &filtered_meta {
        let matching_orig = original_meta.iter().find(|orig| {
            orig.query_name == filt.query_name
                && orig.target_name == filt.target_name
                && orig.query_start == filt.query_start
                && orig.query_end == filt.query_end
                && orig.target_start == filt.target_start
                && orig.target_end == filt.target_end
        });

        assert!(
            matching_orig.is_some(),
            "Filtered record not found in original: {} → {} [{}-{}, {}-{}]",
            filt.query_name,
            filt.target_name,
            filt.query_start,
            filt.query_end,
            filt.target_start,
            filt.target_end
        );
    }

    Ok(())
}
