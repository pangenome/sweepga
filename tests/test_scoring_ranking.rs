// Tests specifically for verifying correct scoring and ranking of mappings
use sweepga::paf_filter::ScoringFunction;
use sweepga::plane_sweep_exact::{plane_sweep_query, PlaneSweepMapping};

fn make_mapping_with_identity(
    idx: usize,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
    identity: f64,
) -> PlaneSweepMapping {
    PlaneSweepMapping {
        idx,
        query_start,
        query_end,
        target_start,
        target_end,
        identity,
        flags: 0,
    }
}

#[test]
fn test_identity_scoring_prefers_high_identity() {
    // Test that identity scoring correctly ranks by identity, ignoring length
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 500, 1000, 1400, 0.70), // Long but low identity
        make_mapping_with_identity(1, 100, 200, 2000, 2100, 0.99), // Short but high identity
        make_mapping_with_identity(2, 100, 300, 3000, 3200, 0.85), // Medium length, medium identity
    ];

    // With identity scoring and n=1, should keep the highest identity (idx 1)
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::Identity);
    assert_eq!(kept.len(), 1, "Should keep exactly 1 mapping");
    assert_eq!(
        kept[0], 1,
        "Should keep mapping with 99% identity, not the longest"
    );
}

#[test]
fn test_length_scoring_prefers_long_alignments() {
    // Test that length scoring correctly ranks by length, ignoring identity
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 200, 1000, 1100, 0.99), // Short but high identity
        make_mapping_with_identity(1, 100, 600, 2000, 2500, 0.50), // Long but low identity
        make_mapping_with_identity(2, 100, 350, 3000, 3250, 0.75), // Medium length, medium identity
    ];

    // With length scoring and n=1, should keep the longest (idx 1)
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::Length);
    assert_eq!(kept.len(), 1, "Should keep exactly 1 mapping");
    assert_eq!(
        kept[0], 1,
        "Should keep longest mapping (500bp), not highest identity"
    );
}

#[test]
fn test_length_identity_scoring_balances_both() {
    // Test that length*identity scoring balances both factors
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 200, 1000, 1100, 0.95), // 100 * 0.95 = 95
        make_mapping_with_identity(1, 100, 400, 2000, 2300, 0.60), // 300 * 0.60 = 180
        make_mapping_with_identity(2, 100, 300, 3000, 3200, 0.80), // 200 * 0.80 = 160
    ];

    // With length*identity scoring, should keep mapping 1 (highest product)
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::LengthIdentity);
    assert_eq!(kept.len(), 1, "Should keep exactly 1 mapping");
    assert_eq!(
        kept[0], 1,
        "Should keep mapping with highest length*identity product"
    );
}

#[test]
fn test_log_length_identity_scoring_dampens_length_effect() {
    // Test that log(length)*identity scoring dampens the effect of length
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 200, 1000, 1100, 0.95), // log(100) * 0.95 ≈ 4.37
        make_mapping_with_identity(1, 100, 1100, 2000, 3000, 0.60), // log(1000) * 0.60 ≈ 4.14
        make_mapping_with_identity(2, 100, 600, 3000, 3500, 0.75), // log(500) * 0.75 ≈ 4.65
    ];

    // With log(length)*identity scoring, should keep mapping 2
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::LogLengthIdentity);
    assert_eq!(kept.len(), 1, "Should keep exactly 1 mapping");
    assert_eq!(
        kept[0], 2,
        "Should keep mapping with highest log(length)*identity score"
    );
}

#[test]
fn test_ranking_with_identical_scores() {
    // Test behavior when multiple mappings have identical scores
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 300, 1000, 1200, 0.90), // 200 * 0.90 = 180
        make_mapping_with_identity(1, 100, 280, 2000, 2180, 1.00), // 180 * 1.00 = 180
        make_mapping_with_identity(2, 100, 460, 3000, 3360, 0.50), // 360 * 0.50 = 180
    ];

    // All have identical length*identity scores
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::LengthIdentity);
    assert_eq!(
        kept.len(),
        1,
        "Should keep exactly 1 mapping even with ties"
    );
    // Don't check which one is kept, as tie-breaking is implementation-dependent
}

#[test]
fn test_scoring_preserves_non_overlapping() {
    // Test that non-overlapping mappings are all kept regardless of scores
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 200, 1000, 1100, 0.50), // Low score
        make_mapping_with_identity(1, 300, 500, 2000, 2200, 0.99), // High score
        make_mapping_with_identity(2, 600, 700, 3000, 3100, 0.30), // Very low score
    ];

    // With n=1 (best per position), all should be kept as they don't overlap
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::Identity);
    assert_eq!(kept.len(), 3, "All non-overlapping mappings should be kept");
}

#[test]
fn test_overlapping_mappings_best_survives() {
    // Test that with overlapping mappings, the best scoring one survives
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 300, 1000, 1200, 0.85), // Overlaps all others
        make_mapping_with_identity(1, 150, 350, 2000, 2200, 0.90), // Better identity, overlaps
        make_mapping_with_identity(2, 200, 400, 3000, 3200, 0.95), // Best identity, overlaps
    ];

    // With identity scoring and n=1, mapping 2 should survive
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::Identity);

    // All three overlap significantly, but mapping 2 has the best identity
    // Due to the sweep algorithm, we expect the best scoring one to be kept
    assert!(
        kept.contains(&2),
        "Best identity mapping (95%) should be kept"
    );
}

#[test]
fn test_scoring_with_contained_mappings() {
    // Test scoring when one mapping is completely contained in another
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 500, 1000, 1400, 0.80), // Large, lower identity
        make_mapping_with_identity(1, 200, 300, 2000, 2100, 0.99), // Small, contained, high identity
    ];

    // Test with different scoring functions

    // Identity scoring should prefer the small high-identity mapping
    let kept = plane_sweep_query(&mut mappings.clone(), 1, 0.95, ScoringFunction::Identity);
    assert!(
        kept.contains(&1),
        "Identity scoring should keep high-identity mapping"
    );

    // Length scoring should prefer the large mapping
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings.clone(), 1, 0.95, ScoringFunction::Length);
    assert!(
        kept.contains(&0),
        "Length scoring should keep longer mapping"
    );

    // Log-length-identity might prefer either depending on exact values
    // log(400) * 0.80 ≈ 4.79 vs log(100) * 0.99 ≈ 4.55
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::LogLengthIdentity);
    assert!(
        kept.contains(&0),
        "Log-length-identity should keep mapping with better combined score"
    );
}

#[test]
fn test_ranking_order_with_multiple_mappings() {
    // Test that mappings are correctly ranked and filtered
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 200, 1000, 1100, 0.70), // Score: 100 * 0.70 = 70
        make_mapping_with_identity(1, 100, 250, 2000, 2150, 0.80), // Score: 150 * 0.80 = 120
        make_mapping_with_identity(2, 100, 300, 3000, 3200, 0.90), // Score: 200 * 0.90 = 180
        make_mapping_with_identity(3, 100, 180, 4000, 4080, 0.99), // Score: 80 * 0.99 = 79.2
        make_mapping_with_identity(4, 100, 220, 5000, 5120, 0.60), // Score: 120 * 0.60 = 72
    ];

    // With length*identity scoring and n=2, should keep top 2 (idx 2 and 1)
    let kept = plane_sweep_query(&mut mappings, 2, 0.95, ScoringFunction::LengthIdentity);
    assert_eq!(kept.len(), 2, "Should keep exactly 2 mappings");
    assert!(
        kept.contains(&2),
        "Should keep highest scoring mapping (180)"
    );
    assert!(
        kept.contains(&1),
        "Should keep second highest scoring mapping (120)"
    );
}

#[test]
fn test_extreme_values() {
    // Test with extreme identity and length values
    let mut mappings = vec![
        make_mapping_with_identity(0, 100, 101, 1000, 1001, 1.00), // Tiny but perfect
        make_mapping_with_identity(1, 100, 100100, 2000, 102000, 0.01), // Huge but terrible
        make_mapping_with_identity(2, 100, 1100, 3000, 4000, 0.50), // Medium both ways
    ];

    // Length scoring should strongly prefer the huge mapping
    let kept = plane_sweep_query(&mut mappings.clone(), 1, 0.95, ScoringFunction::Length);
    assert_eq!(
        kept[0], 1,
        "Length scoring should prefer 100kb mapping despite 1% identity"
    );

    // Identity scoring should strongly prefer the perfect match
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings.clone(), 1, 0.95, ScoringFunction::Identity);
    assert_eq!(
        kept[0], 0,
        "Identity scoring should prefer 100% identity despite 1bp length"
    );

    // Log-length-identity should be more balanced
    // log(1) * 1.0 = 0, log(100000) * 0.01 ≈ 0.115, log(1000) * 0.5 ≈ 3.45
    mappings.iter_mut().for_each(|m| m.flags = 0);
    let kept = plane_sweep_query(&mut mappings, 1, 0.95, ScoringFunction::LogLengthIdentity);
    assert_eq!(
        kept[0], 2,
        "Log-length-identity should prefer balanced mapping"
    );
}
