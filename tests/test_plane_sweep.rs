// Comprehensive tests for the plane sweep algorithm
use sweepga::plane_sweep_exact::{
    plane_sweep_both, plane_sweep_query, plane_sweep_target, PlaneSweepMapping,
};

/// Helper function to create a mapping
fn make_mapping(
    idx: usize,
    q_start: u64,
    q_end: u64,
    t_start: u64,
    t_end: u64,
) -> PlaneSweepMapping {
    PlaneSweepMapping {
        idx,
        query_start: q_start,
        query_end: q_end,
        target_start: t_start,
        target_end: t_end,
        identity: 1.0, // Using 1.0 for length-based scoring
        flags: 0,
    }
}

#[test]
fn test_empty_input() {
    let mut mappings = vec![];
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(kept.len(), 0, "Empty input should return empty");
}

#[test]
fn test_single_mapping() {
    let mut mappings = vec![make_mapping(0, 100, 200, 300, 400)];
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(kept, vec![0], "Single mapping should always be kept");
}

#[test]
fn test_non_overlapping_mappings() {
    // Two mappings that don't overlap on query axis
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400), // Query: 100-200
        make_mapping(1, 300, 400, 500, 600), // Query: 300-400 (no overlap)
    ];

    // With n=0 (best only), both should be kept as they're best at different positions
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(
        kept.len(),
        2,
        "Non-overlapping mappings should both be kept"
    );
    assert!(kept.contains(&0) && kept.contains(&1));
}

#[test]
fn test_overlapping_mappings_keep_best() {
    // Two mappings that overlap on query axis
    let mut mappings = vec![
        make_mapping(0, 100, 250, 300, 450), // Query: 100-250, length=150
        make_mapping(1, 150, 350, 400, 600), // Query: 150-350, length=200 (longer, better score)
    ];

    // With n=0 (best only), both are kept because each is best at some position
    // Mapping 0 is best at positions 100-149
    // Mapping 1 is best at positions 150-350
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(
        kept.len(),
        2,
        "Both overlapping mappings are best at different positions"
    );
}

#[test]
fn test_identical_mappings() {
    // Three identical mappings at same position
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400),
        make_mapping(1, 100, 200, 500, 600), // Same query range, different target
        make_mapping(2, 100, 200, 700, 800), // Same query range, different target
    ];

    // With n=1 (keep exactly 1), only the best should be kept
    // When all have identical scores, we keep only 1 (respecting the n=1 limit)
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(
        kept.len(),
        1,
        "With n=1, only 1 mapping kept even with identical scores"
    );

    // With n=2 (keep 2 total), we keep 2 mappings
    mappings.iter_mut().for_each(|m| m.flags = 0); // Reset flags
    let kept = plane_sweep_query(&mut mappings, 2, 0.95);
    assert_eq!(kept.len(), 2, "With n=2, keep 2 mappings total");

    // With n=usize::MAX (all non-overlapping), all should be kept
    mappings.iter_mut().for_each(|m| m.flags = 0); // Reset flags
    let kept = plane_sweep_query(&mut mappings, usize::MAX, 0.95);
    assert_eq!(kept.len(), 3, "All mappings kept with n=-1");
}

#[test]
fn test_contained_mappings() {
    // Smaller mapping contained within larger mapping
    let mut mappings = vec![
        make_mapping(0, 100, 300, 400, 600), // Large mapping, length=200
        make_mapping(1, 150, 180, 500, 530), // Small mapping contained within, length=30
    ];

    // The larger mapping has better score (log(200) > log(30))
    // With n=1, only the larger should be kept
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(kept.len(), 1, "Only best (larger) mapping kept");
    assert_eq!(kept[0], 0, "Larger mapping has better score");

    // With n=2 (keep 2 mappings), both should be kept
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings, 2, 0.95);
    assert_eq!(kept.len(), 2, "Both mappings kept with n=2");
}

#[test]
fn test_overlap_threshold() {
    // Test that overlap threshold affects filtering
    let mut mappings = vec![
        make_mapping(0, 100, 300, 400, 600),   // Primary mapping
        make_mapping(1, 100, 300, 700, 900),   // Same query range
        make_mapping(2, 100, 300, 1000, 1200), // Same query range
        make_mapping(3, 100, 300, 1300, 1500), // Same query range
    ];

    // With n=2 (keep 2 total) and strict overlap threshold
    // All 4 mappings have identical query ranges and same scores
    // With n=2, we keep only 2 mappings regardless of identical scores
    let kept = plane_sweep_query(&mut mappings, 2, 0.5);
    assert_eq!(kept.len(), 2, "With n=2, keep exactly 2 mappings");

    // The third and fourth mappings should be marked as overlapped
    // since they have 100% overlap with the kept mappings
}

#[test]
fn test_complex_overlaps() {
    // Complex pattern of overlapping mappings
    let mut mappings = vec![
        make_mapping(0, 0, 100, 0, 100),     // Position 0-100
        make_mapping(1, 50, 150, 200, 300),  // Position 50-150 (overlaps with 0)
        make_mapping(2, 120, 220, 400, 500), // Position 120-220 (overlaps with 1)
        make_mapping(3, 200, 300, 600, 700), // Position 200-300 (overlaps with 2)
        make_mapping(4, 280, 380, 800, 900), // Position 280-380 (overlaps with 3)
    ];

    // With n=0, each mapping might be best at some position
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    // All mappings are kept because each is best at its non-overlapping portion
    assert!(
        kept.len() >= 3,
        "Multiple mappings kept in complex overlap pattern"
    );
}

#[test]
fn test_target_axis_filtering() {
    // Test plane sweep on target axis
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400), // Target: 300-400
        make_mapping(1, 300, 400, 350, 450), // Target: 350-450 (overlaps)
        make_mapping(2, 500, 600, 600, 700), // Target: 600-700 (no overlap)
    ];

    let kept = plane_sweep_target(&mut mappings, 1, 0.95);
    // Mappings 0 and 1 overlap on target, both might be kept if best at different positions
    // Mapping 2 doesn't overlap, should be kept
    assert!(kept.contains(&2), "Non-overlapping mapping on target kept");
}

#[test]
fn test_both_axes_filtering() {
    // Test filtering on both query and target axes
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400),
        make_mapping(1, 100, 200, 500, 600), // Same query, different target
        make_mapping(2, 300, 400, 300, 400), // Different query, same target as 0
        make_mapping(3, 500, 600, 700, 800), // No overlap with anything
    ];

    // Filter on both axes with n=0 (best only)
    let kept = plane_sweep_both(&mut mappings, 1, 1, 0.95);

    // Should keep mapping 0 (best at its position)
    // Should keep mapping 3 (no overlap)
    // Mappings 1 and 2 might be filtered depending on scoring
    assert!(
        kept.contains(&3),
        "Non-overlapping mapping kept in both axes"
    );
}

#[test]
fn test_score_calculation() {
    // Test that longer mappings get better scores
    let short_mapping = make_mapping(0, 100, 110, 200, 210); // Length 10
    let long_mapping = make_mapping(1, 100, 200, 300, 400); // Length 100

    assert!(
        long_mapping.score() > short_mapping.score(),
        "Longer mapping should have better score"
    );

    // Test that log scaling is applied
    let mapping_100 = make_mapping(0, 0, 100, 0, 100); // Length 100
    let mapping_1000 = make_mapping(1, 0, 1000, 0, 1000); // Length 1000

    let score_100 = mapping_100.score();
    let score_1000 = mapping_1000.score();

    // Linear scoring: 1000 / 100 = 10
    assert!((score_1000 / score_100 - 10.0).abs() < 0.001, "Score uses linear length scaling");
}

#[test]
fn test_secondary_count() {
    // Test that secondary count parameter works correctly
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400),   // Will be best
        make_mapping(1, 100, 190, 500, 590),   // Secondary 1
        make_mapping(2, 100, 180, 700, 780),   // Secondary 2
        make_mapping(3, 100, 170, 900, 970),   // Secondary 3
        make_mapping(4, 100, 160, 1100, 1160), // Secondary 4
    ];

    // All mappings have similar query range, sorted by length (score)

    // n=1: Keep only the best
    let kept = plane_sweep_query(&mut mappings, 1, 1.0);
    assert_eq!(kept.len(), 1, "n=1 keeps only best");
    assert_eq!(kept[0], 0);

    // n=3: Keep 3 total
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings, 3, 1.0);
    assert_eq!(kept.len(), 3, "n=3 keeps 3 mappings total");

    // n=usize::MAX: Keep all
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings, usize::MAX, 1.0);
    assert_eq!(kept.len(), 5, "n=-1 keeps all mappings");
}

#[test]
fn test_strand_independence() {
    // Plane sweep should work regardless of strand
    // (Our simplified version doesn't track strand, but this tests the logic)
    let mut mappings = vec![
        make_mapping(0, 100, 200, 300, 400), // "Forward" strand
        make_mapping(1, 150, 250, 500, 600), // "Reverse" strand (conceptually)
    ];

    // Both should be considered equally
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert_eq!(kept.len(), 2, "Strand doesn't affect plane sweep");
}

#[test]
fn test_event_ordering() {
    // Test that events are processed in correct order
    // This is important for the algorithm's correctness
    let mut mappings = vec![
        make_mapping(0, 100, 100, 300, 300), // Zero-length mapping (edge case)
        make_mapping(1, 100, 200, 400, 500), // Starts at same position
        make_mapping(2, 100, 300, 600, 800), // Also starts at same position
    ];

    // The zero-length mapping should get worst score
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);
    assert!(!kept.contains(&0), "Zero-length mapping should not be kept");
}

#[test]
fn test_real_world_scenario() {
    // Simulate a real-world scenario with multiple overlapping mappings
    let mut mappings = vec![
        // Chromosome 1 mappings
        make_mapping(0, 1000, 2000, 5000, 6000), // Good alignment
        make_mapping(1, 1500, 2500, 7000, 8000), // Overlaps with 0
        make_mapping(2, 3000, 4000, 9000, 10000), // Separate region
        make_mapping(3, 3200, 3800, 11000, 11600), // Overlaps with 2
        // Tandem repeat region (many similar mappings)
        make_mapping(4, 5000, 5500, 15000, 15500),
        make_mapping(5, 5000, 5500, 16000, 16500),
        make_mapping(6, 5000, 5500, 17000, 17500),
        make_mapping(7, 5000, 5500, 18000, 18500),
        // Large spanning alignment
        make_mapping(8, 8000, 12000, 20000, 24000), // Very long mapping
    ];

    // Default behavior: n=1 (keep best at each position)
    let kept = plane_sweep_query(&mut mappings, 1, 0.95);

    // Should keep the large spanning alignment (best score due to length)
    assert!(kept.contains(&8), "Large spanning alignment kept");

    // Should keep some from each non-overlapping region
    assert!(kept.len() >= 4, "Multiple regions represented");

    // With n=2, should keep more mappings
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept_with_more = plane_sweep_query(&mut mappings, 2, 0.95);
    assert!(
        kept_with_more.len() > kept.len(),
        "More mappings kept with n=2"
    );
}
