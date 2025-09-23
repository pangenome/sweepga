/// Test symmetry property: swapping query and target should give same results
use sweepga::plane_sweep_core::{plane_sweep, Interval};

/// Helper to create test intervals
fn make_interval(idx: usize, begin: u32, end: u32) -> Interval {
    Interval {
        idx,
        begin,
        end,
        score: (end - begin) as f64, // Use length as score
        flags: 0,
    }
}

#[test]
fn test_symmetry_simple() {
    // Create mappings with query and target coordinates
    let mappings = vec![
        (100, 200, 300, 400), // query: 100-200, target: 300-400
        (150, 250, 350, 450), // overlaps on both axes
        (300, 400, 100, 200), // non-overlapping
    ];

    // Test query axis filtering
    let mut query_intervals: Vec<Interval> = mappings
        .iter()
        .enumerate()
        .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
        .collect();

    // Test target axis filtering
    let mut target_intervals: Vec<Interval> = mappings
        .iter()
        .enumerate()
        .map(|(i, &(_, _, t_start, t_end))| make_interval(i, t_start, t_end))
        .collect();

    // Apply plane sweep with same parameters
    let query_kept = plane_sweep(&mut query_intervals, 1, 0.95);
    let target_kept = plane_sweep(&mut target_intervals, 1, 0.95);

    // Query axis: mappings 0 and 2 don't overlap (100-200 vs 300-400)
    // With n=1, we keep best at each position, so both are kept
    assert_eq!(query_kept.len(), 2, "Non-overlapping on query axis");

    // Target axis: mappings 0 and 2 don't overlap (300-400 vs 100-200)
    // With n=1, we keep best at each position, so both are kept
    assert_eq!(target_kept.len(), 2, "Non-overlapping on target axis");
}

#[test]
fn test_symmetry_transposed() {
    // Create mappings where we swap query/target
    let original = vec![
        (100, 500, 1000, 1400), // 400 length on both
        (200, 400, 1100, 1300), // 200 length on both
        (600, 900, 1500, 1800), // 300 length on both
    ];

    let transposed = vec![
        (1000, 1400, 100, 500), // Same mapping, swapped axes
        (1100, 1300, 200, 400),
        (1500, 1800, 600, 900),
    ];

    // Filter original on query axis
    let mut orig_query_intervals: Vec<Interval> = original
        .iter()
        .enumerate()
        .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
        .collect();

    // Filter transposed on "query" axis (which is actually original's target)
    let mut trans_query_intervals: Vec<Interval> = transposed
        .iter()
        .enumerate()
        .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
        .collect();

    let orig_kept = plane_sweep(&mut orig_query_intervals, 2, 0.95);
    let trans_kept = plane_sweep(&mut trans_query_intervals, 2, 0.95);

    // Should get same number of mappings kept
    assert_eq!(orig_kept.len(), trans_kept.len());
}

#[test]
fn test_symmetry_with_overlaps() {
    // Mappings with complex overlaps
    let mappings = vec![
        (0, 100, 0, 100),     // Square mapping
        (50, 150, 50, 150),   // Overlaps 50% with first
        (200, 300, 200, 300), // Non-overlapping square
        (250, 350, 250, 350), // Overlaps 50% with third
    ];

    // Test with different n values
    for n in [1, 2, 3, 4] {
        let mut query_intervals: Vec<Interval> = mappings
            .iter()
            .enumerate()
            .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
            .collect();

        let mut target_intervals: Vec<Interval> = mappings
            .iter()
            .enumerate()
            .map(|(i, &(_, _, t_start, t_end))| make_interval(i, t_start, t_end))
            .collect();

        let query_kept = plane_sweep(&mut query_intervals, n, 0.95);
        let target_kept = plane_sweep(&mut target_intervals, n, 0.95);

        // For symmetric mappings, should get same results
        assert_eq!(
            query_kept, target_kept,
            "Symmetry broken for n={}: query {:?} != target {:?}",
            n, query_kept, target_kept
        );
    }
}

#[test]
fn test_asymmetric_mappings() {
    // Mappings that are NOT symmetric - different behavior expected
    let mappings = vec![
        (0, 200, 0, 100),    // Longer on query axis
        (50, 150, 200, 400), // Longer on target axis
        (300, 500, 50, 150), // Longer on query axis
    ];

    let mut query_intervals: Vec<Interval> = mappings
        .iter()
        .enumerate()
        .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
        .collect();

    let mut target_intervals: Vec<Interval> = mappings
        .iter()
        .enumerate()
        .map(|(i, &(_, _, t_start, t_end))| make_interval(i, t_start, t_end))
        .collect();

    let query_kept = plane_sweep(&mut query_intervals, 1, 0.95);
    let target_kept = plane_sweep(&mut target_intervals, 1, 0.95);

    // Query axis should keep mapping 2 (longest query interval: 200)
    assert!(query_kept.contains(&2), "Query should keep mapping 2");

    // Target axis should keep mapping 1 (longest target interval: 200)
    assert!(target_kept.contains(&1), "Target should keep mapping 1");

    // They should be different due to asymmetry
    assert_ne!(
        query_kept, target_kept,
        "Asymmetric mappings should give different results"
    );
}

#[test]
fn test_perfect_symmetry_all_axes() {
    // Test that filtering on both axes with symmetric data produces symmetric results
    let symmetric_mappings = vec![
        (100, 300, 100, 300), // 200x200
        (400, 600, 400, 600), // 200x200
        (700, 900, 700, 900), // 200x200
    ];

    for n in [1, 2, 3, usize::MAX] {
        let mut query_intervals: Vec<Interval> = symmetric_mappings
            .iter()
            .enumerate()
            .map(|(i, &(q_start, q_end, _, _))| make_interval(i, q_start, q_end))
            .collect();

        let mut target_intervals: Vec<Interval> = symmetric_mappings
            .iter()
            .enumerate()
            .map(|(i, &(_, _, t_start, t_end))| make_interval(i, t_start, t_end))
            .collect();

        let query_kept = plane_sweep(&mut query_intervals, n, 0.95);
        let target_kept = plane_sweep(&mut target_intervals, n, 0.95);

        assert_eq!(
            query_kept, target_kept,
            "Perfect symmetry should give identical results for n={}",
            n
        );

        // All mappings are non-overlapping and equal score
        // So all should be kept regardless of n
        assert_eq!(
            query_kept.len(),
            3,
            "Should keep all 3 non-overlapping mappings"
        );
    }
}
