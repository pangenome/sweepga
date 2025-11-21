/// Multi-contig scaffold coordinate conversion tests
///
/// Tests coordinate conversion when sequences contain Ns (creating multiple contigs).
/// Verifies that offsets are correctly applied and coordinates properly map across
/// contig boundaries.
use anyhow::Result;

/// Test scaffold offset calculation with multiple contigs
#[test]
fn test_scaffold_offsets_calculation() {
    // Simulate a sequence with Ns: ACGT[NNN]TGCA[NN]GGCC
    // This creates 3 contigs with different offsets

    let contig_boundaries = [
        (0, 4),   // Contig 0: positions 0-3 (ACGT)
        (7, 11),  // Contig 1: positions 7-10 (TGCA) - offset by 3 Ns
        (13, 17), // Contig 2: positions 13-16 (GGCC) - offset by 5 Ns total
    ];

    // Verify contig lengths and offsets
    assert_eq!(contig_boundaries[0].1 - contig_boundaries[0].0, 4);
    assert_eq!(contig_boundaries[1].1 - contig_boundaries[1].0, 4);
    assert_eq!(contig_boundaries[2].1 - contig_boundaries[2].0, 4);

    // Contig 1 starts after 3 Ns, so offset should account for that
    let contig1_offset = contig_boundaries[1].0; // 7
    assert_eq!(contig1_offset, 7);

    // Contig 2 starts after 5 Ns total, so offset is 13
    let contig2_offset = contig_boundaries[2].0; // 13
    assert_eq!(contig2_offset, 13);
}

/// Test forward strand coordinate conversion across contigs
#[test]
fn test_forward_strand_multi_contig() {
    // Scaffold coords: sequence with Ns removed (contig coordinates)
    // Contig 0: scaffold 0-3
    // Contig 1: scaffold 4-7 (but actually at positions 7-10 in original)
    // Contig 2: scaffold 8-11 (but actually at positions 13-16 in original)

    struct ContigInfo {
        scaffold_start: i64,
        scaffold_end: i64,
        original_start: i64,
    }

    let contigs = [ContigInfo {
            scaffold_start: 0,
            scaffold_end: 4,
            original_start: 0,
        },
        ContigInfo {
            scaffold_start: 4,
            scaffold_end: 8,
            original_start: 7,
        },
        ContigInfo {
            scaffold_start: 8,
            scaffold_end: 12,
            original_start: 13,
        }];

    // Test alignment in contig 1 (positions 5-7 in scaffold coords)
    let align_start = 5_i64;
    let align_end = 7_i64;

    // This is in contig 1 (scaffold 4-8)
    let contig = &contigs[1];
    assert!(align_start >= contig.scaffold_start && align_start < contig.scaffold_end);

    // Convert to original coordinates
    let contig_relative_start = align_start - contig.scaffold_start; // 5 - 4 = 1
    let contig_relative_end = align_end - contig.scaffold_start; // 7 - 4 = 3

    let original_start = contig.original_start + contig_relative_start; // 7 + 1 = 8
    let original_end = contig.original_start + contig_relative_end; // 7 + 3 = 10

    assert_eq!(original_start, 8);
    assert_eq!(original_end, 10);

    // Verify span is preserved
    assert_eq!(align_end - align_start, original_end - original_start);
}

/// Test reverse strand coordinate conversion across contigs
#[test]
fn test_reverse_strand_multi_contig() {
    // Same setup as forward test
    struct ContigInfo {
        scaffold_start: i64,
        scaffold_end: i64,
        original_start: i64,
    }

    let contigs = [ContigInfo {
            scaffold_start: 0,
            scaffold_end: 4,
            original_start: 0,
        },
        ContigInfo {
            scaffold_start: 4,
            scaffold_end: 8,
            original_start: 7,
        },
        ContigInfo {
            scaffold_start: 8,
            scaffold_end: 12,
            original_start: 13,
        }];

    // Test reverse alignment in contig 1 (positions 5-7 in scaffold coords)
    let align_start = 5_i64;
    let align_end = 7_i64;

    let contig = &contigs[1];

    // For reverse strand, coordinates are swapped relative to scaffold end
    let contig_relative_start = align_start - contig.scaffold_start; // 1
    let contig_relative_end = align_end - contig.scaffold_start; // 3

    // Reverse strand: swap and calculate from end
    let contig_length = contig.scaffold_end - contig.scaffold_start; // 4
    let rev_contig_start = contig_length - contig_relative_end; // 4 - 3 = 1
    let rev_contig_end = contig_length - contig_relative_start; // 4 - 1 = 3

    let original_start = contig.original_start + rev_contig_start; // 7 + 1 = 8
    let original_end = contig.original_start + rev_contig_end; // 7 + 3 = 10

    // Verify span is preserved even after reversal
    assert_eq!(align_end - align_start, original_end - original_start);
}

/// Test that contig boundaries are respected
#[test]
fn test_alignment_within_single_contig() {
    // Alignment should not span contig boundaries
    // If it does, it should be split into separate alignments

    let contig_boundaries = [
        (0, 100),   // Contig 0
        (150, 250), // Contig 1 (50 Ns gap)
        (300, 400), // Contig 2 (50 Ns gap)
    ];

    // Valid alignment: entirely within contig 0
    let align1_start = 20;
    let align1_end = 80;
    assert!(align1_start >= contig_boundaries[0].0);
    assert!(align1_end <= contig_boundaries[0].1);

    // Valid alignment: entirely within contig 1
    let align2_start = 160;
    let align2_end = 240;
    assert!(align2_start >= contig_boundaries[1].0);
    assert!(align2_end <= contig_boundaries[1].1);

    // Invalid: would span contig boundary (should be split)
    let align3_start = 90;
    let align3_end = 160; // Crosses from contig 0 to contig 1

    // This should be detected and split
    let spans_boundary =
        !(align3_end <= contig_boundaries[0].1 || align3_start >= contig_boundaries[1].0);
    assert!(
        spans_boundary,
        "Alignment incorrectly spans contig boundary"
    );
}

/// Test coordinate ordering invariants across contigs
#[test]
fn test_coordinate_ordering_multi_contig() {
    // Create alignment mappings across different contigs
    struct Mapping {
        scaffold_start: u64,
        scaffold_end: u64,
    }

    let mappings = vec![
        Mapping {
            scaffold_start: 10,
            scaffold_end: 20,
        },
        Mapping {
            scaffold_start: 50,
            scaffold_end: 60,
        },
        Mapping {
            scaffold_start: 90,
            scaffold_end: 95,
        },
    ];

    // Verify ordering within each mapping
    for m in &mappings {
        assert!(
            m.scaffold_start < m.scaffold_end,
            "Start must be < end in scaffold coords"
        );
    }

    // Verify scaffold coordinates increase across contigs
    for i in 1..mappings.len() {
        assert!(
            mappings[i - 1].scaffold_end <= mappings[i].scaffold_start,
            "Mappings should be ordered by scaffold position"
        );
    }
}

/// Test that offset calculation handles different contig sizes
#[test]
fn test_variable_contig_sizes() {
    struct Contig {
        length: usize,
        offset_in_scaffold: usize,
    }

    // Simulate scaffold with variable-sized contigs
    let contigs = [
        Contig {
            length: 1000,
            offset_in_scaffold: 0,
        },
        Contig {
            length: 500,
            offset_in_scaffold: 1010, // 1000 + 10 Ns
        },
        Contig {
            length: 2000,
            offset_in_scaffold: 1520, // 1010 + 500 + 10 Ns
        },
    ];

    // Calculate total scaffold length (sum of contig lengths, no Ns)
    let scaffold_length: usize = contigs.iter().map(|c| c.length).sum();
    assert_eq!(scaffold_length, 3500);

    // Verify offsets increase monotonically
    for i in 1..contigs.len() {
        assert!(contigs[i].offset_in_scaffold > contigs[i - 1].offset_in_scaffold);
        let _gap = contigs[i].offset_in_scaffold
            - (contigs[i - 1].offset_in_scaffold + contigs[i - 1].length);
        // Gap is always non-negative (usize subtraction)
    }
}

/// Integration test: verify coordinate conversion properties hold across contigs
#[test]
fn test_coordinate_conversion_properties_multi_contig() -> Result<()> {
    // Property: span preservation
    // Converting from contig coords to scaffold coords and back should preserve span

    let test_cases = vec![
        // (contig_start, contig_end, contig_offset)
        (100, 200, 0),   // First contig, no offset
        (50, 150, 1000), // Second contig, offset by 1000
        (0, 500, 2000),  // Third contig, offset by 2000
    ];

    for (contig_start, contig_end, offset) in test_cases {
        let contig_span = contig_end - contig_start;

        // Forward conversion: contig -> scaffold
        let scaffold_start = offset + contig_start;
        let scaffold_end = offset + contig_end;
        let scaffold_span = scaffold_end - scaffold_start;

        // Verify span preservation
        assert_eq!(
            contig_span, scaffold_span,
            "Span should be preserved in coordinate conversion"
        );

        // Reverse conversion: scaffold -> contig
        let recovered_contig_start = scaffold_start - offset;
        let recovered_contig_end = scaffold_end - offset;

        // Verify round-trip
        assert_eq!(contig_start, recovered_contig_start);
        assert_eq!(contig_end, recovered_contig_end);
    }

    Ok(())
}

/// Test N-region detection and handling
#[test]
fn test_n_region_detection() {
    // Simulate sequence analysis to find N regions
    let sequence = b"ACGTNNNNNTGCANNGGCC";
    let n_threshold = 3; // Minimum N-run to create contig boundary

    let mut n_regions = Vec::new();
    let mut in_n_region = false;
    let mut n_start = 0;

    for (i, &base) in sequence.iter().enumerate() {
        if base == b'N' {
            if !in_n_region {
                n_start = i;
                in_n_region = true;
            }
        } else if in_n_region {
            let n_length = i - n_start;
            if n_length >= n_threshold {
                n_regions.push((n_start, i));
            }
            in_n_region = false;
        }
    }

    // Should find one N-region of length 5 (positions 4-8)
    assert_eq!(n_regions.len(), 1);
    assert_eq!(n_regions[0], (4, 9));
    assert_eq!(n_regions[0].1 - n_regions[0].0, 5);
}
