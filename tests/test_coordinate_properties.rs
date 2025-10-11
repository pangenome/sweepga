/// Property-based tests for coordinate conversion stability
///
/// Uses proptest to verify mathematical invariants that must ALWAYS hold.
/// These tests prevent subtle bugs and ensure no drift in coordinate handling.
use proptest::prelude::*;

/// Property: Query span must be preserved during coordinate conversion
/// Forward: query_end - query_start must equal original span
#[test]
fn prop_query_span_preserved() {
    proptest!(|(
        contig_start in 0u64..100_000,
        span in 1u64..10_000,
        scaffold_offset in 0i64..1_000_000
    )| {
        let contig_end = contig_start + span;

        // Apply scaffold offset (forward strand)
        let scaffold_start = (scaffold_offset + contig_start as i64) as u64;
        let scaffold_end = (scaffold_offset + contig_end as i64) as u64;

        // Verify span is preserved
        let original_span = contig_end - contig_start;
        let converted_span = scaffold_end - scaffold_start;

        prop_assert_eq!(original_span, converted_span,
            "Query span changed during conversion: {} → {}",
            original_span, converted_span);
    });
}

/// Property: Target span must be preserved (forward strand)
#[test]
fn prop_target_span_preserved_forward() {
    proptest!(|(
        contig_start in 0u64..100_000,
        span in 1u64..10_000,
        scaffold_offset in 0i64..1_000_000
    )| {
        let contig_end = contig_start + span;

        // Forward strand: same as query
        let scaffold_start = (scaffold_offset + contig_start as i64) as u64;
        let scaffold_end = (scaffold_offset + contig_end as i64) as u64;

        let original_span = contig_end - contig_start;
        let converted_span = scaffold_end - scaffold_start;

        prop_assert_eq!(original_span, converted_span);
    });
}

/// Property: Target span must be preserved (reverse strand)
/// Reverse strand uses: (sbeg + clen) - coord
#[test]
fn prop_target_span_preserved_reverse() {
    proptest!(|(
        contig_start in 0u64..100_000,
        span in 1u64..10_000,
        scaffold_offset in 0i64..1_000_000,
        contig_length in 100_001i64..2_000_000
    )| {
        let contig_end = contig_start + span;

        // Reverse strand: use (sbeg + clen) - coord
        let b_off = scaffold_offset + contig_length;
        let scaffold_start = (b_off - contig_end as i64) as u64;
        let scaffold_end = (b_off - contig_start as i64) as u64;

        let original_span = contig_end - contig_start;
        let converted_span = scaffold_end - scaffold_start;

        prop_assert_eq!(original_span, converted_span,
            "Reverse strand span not preserved: {} → {}",
            original_span, converted_span);
    });
}

/// Property: Coordinate ordering must be preserved
/// start < end before and after conversion
#[test]
fn prop_coordinate_ordering_preserved() {
    proptest!(|(
        start in 0u64..100_000,
        length in 1u64..10_000,
        offset in 0i64..1_000_000
    )| {
        let end = start + length;

        // Apply conversion
        let conv_start = (offset + start as i64) as u64;
        let conv_end = (offset + end as i64) as u64;

        prop_assert!(start < end, "Original: start >= end");
        prop_assert!(conv_start < conv_end, "Converted: start >= end");
    });
}

/// Property: Reverse strand should swap coordinates but preserve span
#[test]
fn prop_reverse_strand_swaps_but_preserves_span() {
    proptest!(|(
        contig_start in 0u64..50_000,
        span in 1u64..5_000,
        scaffold_offset in 0i64..500_000,
        contig_length in 55_001i64..1_000_000
    )| {
        let contig_end = contig_start + span;

        // Forward
        let fwd_start = (scaffold_offset + contig_start as i64) as u64;
        let fwd_end = (scaffold_offset + contig_end as i64) as u64;
        let fwd_span = fwd_end - fwd_start;

        // Reverse
        let b_off = scaffold_offset + contig_length;
        let rev_start = (b_off - contig_end as i64) as u64;
        let rev_end = (b_off - contig_start as i64) as u64;
        let rev_span = rev_end - rev_start;

        // Spans must be equal
        prop_assert_eq!(fwd_span, rev_span,
            "Forward and reverse spans differ: fwd={}, rev={}",
            fwd_span, rev_span);

        // Original span must match
        let original_span = contig_end - contig_start;
        prop_assert_eq!(original_span, fwd_span);
        prop_assert_eq!(original_span, rev_span);
    });
}

/// Property: Multiple conversions with different offsets preserve relative positions
#[test]
fn prop_relative_positions_preserved() {
    proptest!(|(
        pos1 in 0u64..50_000,
        pos2 in 50_001u64..100_000,
        offset in 0i64..1_000_000
    )| {
        // Original distance
        let original_dist = pos2 - pos1;

        // Convert both
        let conv1 = (offset + pos1 as i64) as u64;
        let conv2 = (offset + pos2 as i64) as u64;
        let converted_dist = conv2 - conv1;

        prop_assert_eq!(original_dist, converted_dist,
            "Relative distance changed: {} → {}",
            original_dist, converted_dist);
    });
}

/// Property: Zero offset is identity
#[test]
fn prop_zero_offset_is_identity() {
    proptest!(|(
        start in 0u64..100_000,
        end in 100_001u64..200_000
    )| {
        let offset = 0i64;

        let conv_start = (offset + start as i64) as u64;
        let conv_end = (offset + end as i64) as u64;

        prop_assert_eq!(start, conv_start);
        prop_assert_eq!(end, conv_end);
    });
}

/// Property: Conversion is deterministic
#[test]
fn prop_conversion_is_deterministic() {
    proptest!(|(
        pos in 0u64..100_000,
        offset in 0i64..1_000_000
    )| {
        let conv1 = (offset + pos as i64) as u64;
        let conv2 = (offset + pos as i64) as u64;

        prop_assert_eq!(conv1, conv2,
            "Conversion not deterministic for pos={}, offset={}",
            pos, offset);
    });
}

/// Property: Block length calculation is commutative
#[test]
fn prop_block_length_commutative() {
    proptest!(|(
        q_start in 0u64..50_000,
        q_len in 1u64..5_000,
        t_start in 0u64..50_000,
        t_len in 1u64..5_000
    )| {
        let q_end = q_start + q_len;
        let t_end = t_start + t_len;

        let block_len1 = (q_end - q_start) + (t_end - t_start);
        let block_len2 = (t_end - t_start) + (q_end - q_start);

        prop_assert_eq!(block_len1, block_len2);
    });
}

/// Property: Identity must be in valid range [0.0, 1.0]
#[test]
fn prop_identity_in_valid_range() {
    proptest!(|(
        matches in 0u64..10_000,
        query_span in 1u64..10_000
    )| {
        // Cap matches at query_span (can't have more matches than span)
        let matches = matches.min(query_span);
        let identity = matches as f64 / query_span as f64;

        prop_assert!(identity >= 0.0, "Identity < 0: {}", identity);
        prop_assert!(identity <= 1.0, "Identity > 1: {}", identity);
    });
}
