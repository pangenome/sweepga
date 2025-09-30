use sweepga::plane_sweep_exact::{plane_sweep_query, plane_sweep_both, PlaneSweepMapping};
use sweepga::paf_filter::{PafFilter, FilterConfig, FilterMode, ScoringFunction};
use tempfile::NamedTempFile;
use std::io::Write;

#[test]
fn test_regular_overlap_threshold_50_percent() {
    // Test with 50% overlap threshold
    let mut mappings = vec![
        PlaneSweepMapping {
            idx: 0,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.95,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 1,
            query_start: 140,  // 60% overlap with first mapping
            query_end: 240,
            target_start: 500,
            target_end: 600,
            identity: 0.90,
            flags: 0,
        },
    ];

    // With 50% overlap threshold, the second mapping should be filtered out
    let kept = plane_sweep_query(&mut mappings, 2, 0.5, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 1, "With 50% overlap threshold, overlapping mapping should be filtered");
    assert_eq!(kept[0], 0, "First (better scoring) mapping should be kept");
}

#[test]
fn test_regular_overlap_threshold_95_percent() {
    // Test with 95% overlap threshold
    let mut mappings = vec![
        PlaneSweepMapping {
            idx: 2,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.95,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 3,
            query_start: 140,  // 60% overlap with first mapping
            query_end: 240,
            target_start: 500,
            target_end: 600,
            identity: 0.90,
            flags: 0,
        },
    ];

    // With 95% overlap threshold, both mappings should be kept (60% < 95%)
    let kept = plane_sweep_query(&mut mappings, 2, 0.95, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 2, "With 95% overlap threshold, both mappings should be kept");
}

#[test]
fn test_no_overlap_filtering() {
    // Test with overlap threshold = 1.0 (no overlap filtering)
    let mut mappings = vec![
        PlaneSweepMapping {
            idx: 4,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.95,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 5,
            query_start: 100,  // 100% overlap with first mapping
            query_end: 200,
            target_start: 500,
            target_end: 600,
            identity: 0.90,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 6,
            query_start: 100,  // 100% overlap with first mapping
            query_end: 200,
            target_start: 700,
            target_end: 800,
            identity: 0.85,
            flags: 0,
        },
    ];

    // With overlap threshold = 1.0, keep only the best scoring ones up to the limit
    let kept = plane_sweep_query(&mut mappings, 2, 1.0, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 2, "With no overlap filtering, should keep best 2");
    assert!(kept.contains(&4), "Best scoring mapping should be kept");
    assert!(kept.contains(&5), "Second best scoring mapping should be kept");
}

#[test]
fn test_overlap_calculation() {
    // Test that overlap is calculated correctly
    let mapping1 = PlaneSweepMapping {
        idx: 0,
        query_start: 100,
        query_end: 200,
        target_start: 300,
        target_end: 400,
        identity: 0.95,
        flags: 0,
    };

    let mapping2 = PlaneSweepMapping {
        idx: 1,
        query_start: 150,  // 50 bp overlap
        query_end: 250,
        target_start: 500,
        target_end: 600,
        identity: 0.90,
        flags: 0,
    };

    // Overlap should be 50bp / min(100bp, 100bp) = 0.5 (50%)
    let overlap = mapping1.query_overlap(&mapping2);
    assert!((overlap - 0.5).abs() < 0.001, "Overlap should be 50%");
}

#[test]
fn test_scaffold_overlap_filtering() {
    // Create a test PAF file with overlapping scaffold chains
    let mut temp_file = NamedTempFile::new().unwrap();
    writeln!(temp_file, "query1\t1000\t100\t600\t+\ttarget1\t1000\t100\t600\t450\t500\t60\tch:Z:chain_1").unwrap();
    writeln!(temp_file, "query1\t1000\t350\t850\t+\ttarget1\t1000\t350\t850\t450\t500\t60\tch:Z:chain_2").unwrap();
    writeln!(temp_file, "query1\t1000\t800\t1000\t+\ttarget1\t1000\t800\t1000\t180\t200\t60\tch:Z:chain_3").unwrap();

    let input_path = temp_file.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    // Test with 50% scaffold overlap threshold
    let config = FilterConfig {
        chain_gap: 0,  // No additional chaining
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 10000,
        min_scaffold_length: 100,
        scaffold_overlap_threshold: 0.5,  // 50% overlap threshold for scaffolds
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    let result = filter.filter_paf(input_path, output_path);
    assert!(result.is_ok(), "Scaffold filtering should succeed");
}

#[test]
fn test_both_axes_overlap() {
    // Test plane_sweep_both with overlap filtering on both query and target
    let mut mappings = vec![
        PlaneSweepMapping {
            idx: 7,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.95,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 8,
            query_start: 140,  // 60% query overlap
            query_end: 240,
            target_start: 340,  // 60% target overlap
            target_end: 440,
            identity: 0.90,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 9,
            query_start: 500,  // No overlap
            query_end: 600,
            target_start: 700,
            target_end: 800,
            identity: 0.85,
            flags: 0,
        },
    ];

    // With 50% overlap threshold on both axes, second mapping should be filtered
    let kept = plane_sweep_both(&mut mappings, 1, 1, 0.5, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 2, "First and third mappings should be kept (no overlap between them)");
    assert!(kept.contains(&7), "First mapping should be kept");
    assert!(kept.contains(&9), "Third mapping should be kept");
}

#[test]
fn test_edge_cases() {
    // Test with exact boundary overlap
    let mut mappings = vec![
        PlaneSweepMapping {
            idx: 10,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.95,
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 11,
            query_start: 200,  // Touching but not overlapping
            query_end: 300,
            target_start: 400,
            target_end: 500,
            identity: 0.90,
            flags: 0,
        },
    ];

    // Even with 0% overlap threshold, touching mappings shouldn't be filtered
    let kept = plane_sweep_query(&mut mappings, 2, 0.0, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 2, "Non-overlapping mappings should both be kept");
}

#[test]
fn test_different_scoring_functions_with_overlap() {
    // Test that overlap filtering works consistently with different scoring functions
    let mappings = vec![
        PlaneSweepMapping {
            idx: 12,
            query_start: 100,
            query_end: 200,
            target_start: 300,
            target_end: 400,
            identity: 0.80,  // Lower identity
            flags: 0,
        },
        PlaneSweepMapping {
            idx: 13,
            query_start: 140,  // 60% overlap
            query_end: 240,
            target_start: 500,
            target_end: 600,
            identity: 0.99,  // Higher identity
            flags: 0,
        },
    ];

    // With identity scoring and 50% overlap threshold
    let mut test_mappings = mappings.clone();
    let kept = plane_sweep_query(&mut test_mappings, 2, 0.5, ScoringFunction::Identity);
    assert_eq!(kept.len(), 1, "Only one mapping should be kept due to overlap");
    assert_eq!(kept[0], 1, "Higher identity mapping should be kept with identity scoring");

    // With length scoring and 50% overlap threshold
    let mut test_mappings = mappings.clone();
    let kept = plane_sweep_query(&mut test_mappings, 2, 0.5, ScoringFunction::Length);
    assert_eq!(kept.len(), 1, "Only one mapping should be kept due to overlap");
    // Both have same length, so first one wins due to position
    assert_eq!(kept[0], 12, "First mapping should be kept with length scoring when lengths are equal");
}