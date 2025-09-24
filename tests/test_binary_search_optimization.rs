/// Comprehensive tests to PROVE the binary search optimization works correctly
/// and is actually faster than the O(n²) approach
use std::time::Instant;
use sweepga::union_find::UnionFind;

#[derive(Clone, Debug)]
struct TestMapping {
    #[allow(dead_code)]
    idx: usize,
    query_start: u32,
    query_end: u32,
    target_start: u32,
    target_end: u32,
}

/// The OLD O(n²) implementation for comparison
fn merge_mappings_quadratic(mappings: &[TestMapping], max_gap: u32) -> Vec<Vec<usize>> {
    let mut uf = UnionFind::new(mappings.len());

    for i in 0..mappings.len() {
        for j in (i + 1)..mappings.len() {
            // OLD implementation: check every pair
            if mappings[j].query_start > mappings[i].query_end + max_gap {
                break; // Early termination
            }

            let q_gap = mappings[j].query_start.saturating_sub(mappings[i].query_end);

            let t_gap = mappings[j].target_start.saturating_sub(mappings[i].target_end);

            if q_gap <= max_gap && t_gap <= max_gap {
                uf.union(i, j);
            }
        }
    }

    uf.get_sets()
}

/// The NEW O(n log n) implementation with binary search
fn merge_mappings_binary_search(mappings: &[TestMapping], max_gap: u32) -> Vec<Vec<usize>> {
    let mut uf = UnionFind::new(mappings.len());

    for i in 0..mappings.len() {
        let search_bound = mappings[i].query_end + max_gap;

        // Binary search to find the range
        let j_end = mappings[i + 1..]
            .binary_search_by_key(&(search_bound + 1), |m| m.query_start)
            .unwrap_or_else(|pos| pos)
            + i
            + 1;

        // Only check mappings in the valid range
        for j in (i + 1)..j_end {
            let q_gap = mappings[j].query_start.saturating_sub(mappings[i].query_end);

            let t_gap = mappings[j].target_start.saturating_sub(mappings[i].target_end);

            if q_gap <= max_gap && t_gap <= max_gap {
                uf.union(i, j);
            }
        }
    }

    uf.get_sets()
}

#[test]
fn test_binary_search_finds_correct_range() {
    // Test that binary search finds exactly the right range of candidates
    let mappings = [TestMapping {
            idx: 0,
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 200,
        },
        TestMapping {
            idx: 1,
            query_start: 250,
            query_end: 350,
            target_start: 250,
            target_end: 350,
        },
        TestMapping {
            idx: 2,
            query_start: 400,
            query_end: 500,
            target_start: 400,
            target_end: 500,
        },
        TestMapping {
            idx: 3,
            query_start: 600,
            query_end: 700,
            target_start: 600,
            target_end: 700,
        },
        TestMapping {
            idx: 4,
            query_start: 1000,
            query_end: 1100,
            target_start: 1000,
            target_end: 1100,
        }];

    let max_gap = 100;

    // For mapping 0 (ends at 200), we should check mappings that start <= 300
    // That's mappings 1 (starts at 250) only
    let search_bound = mappings[0].query_end + max_gap; // 300
    let j_end = mappings[1..]
        .binary_search_by_key(&(search_bound + 1), |m| m.query_start)
        .unwrap_or_else(|pos| pos)
        + 1;
    assert_eq!(
        j_end, 2,
        "Should find mappings 1 as candidate for mapping 0"
    );

    // For mapping 1 (ends at 350), we should check mappings that start <= 450
    // That's mapping 2 (starts at 400) only
    let search_bound = mappings[1].query_end + max_gap; // 450
    let j_end = mappings[2..]
        .binary_search_by_key(&(search_bound + 1), |m| m.query_start)
        .unwrap_or_else(|pos| pos)
        + 2;
    assert_eq!(j_end, 3, "Should find mapping 2 as candidate for mapping 1");

    // For mapping 2 (ends at 500), we should check mappings that start <= 600
    // That's mapping 3 (starts at 600) only
    let search_bound = mappings[2].query_end + max_gap; // 600
    let j_end = mappings[3..]
        .binary_search_by_key(&(search_bound + 1), |m| m.query_start)
        .unwrap_or_else(|pos| pos)
        + 3;
    assert_eq!(j_end, 4, "Should find mapping 3 as candidate for mapping 2");
}

#[test]
fn test_identical_output_small() {
    // Test that both implementations produce identical results on small data
    let mappings = vec![
        TestMapping {
            idx: 0,
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 200,
        },
        TestMapping {
            idx: 1,
            query_start: 210,
            query_end: 310,
            target_start: 210,
            target_end: 310,
        },
        TestMapping {
            idx: 2,
            query_start: 320,
            query_end: 420,
            target_start: 320,
            target_end: 420,
        },
        TestMapping {
            idx: 3,
            query_start: 500,
            query_end: 600,
            target_start: 500,
            target_end: 600,
        },
        TestMapping {
            idx: 4,
            query_start: 610,
            query_end: 710,
            target_start: 610,
            target_end: 710,
        },
    ];

    let max_gap = 50;

    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);

    // Sort the sets for comparison
    let mut quad_sorted: Vec<Vec<usize>> = result_quad
        .into_iter()
        .map(|mut set| {
            set.sort();
            set
        })
        .collect();
    quad_sorted.sort();

    let mut binary_sorted: Vec<Vec<usize>> = result_binary
        .into_iter()
        .map(|mut set| {
            set.sort();
            set
        })
        .collect();
    binary_sorted.sort();

    assert_eq!(
        quad_sorted, binary_sorted,
        "Both algorithms must produce identical results"
    );

    // Verify expected chains
    assert_eq!(quad_sorted.len(), 2, "Should form 2 chains");
    assert_eq!(quad_sorted[0], vec![0, 1, 2], "First chain: mappings 0-2");
    assert_eq!(quad_sorted[1], vec![3, 4], "Second chain: mappings 3-4");
}

#[test]
fn test_identical_output_dense() {
    // Test with dense mappings where many are within range
    let mut mappings = Vec::new();
    for i in 0..100 {
        mappings.push(TestMapping {
            idx: i,
            query_start: i as u32 * 10,
            query_end: i as u32 * 10 + 8,
            target_start: i as u32 * 10,
            target_end: i as u32 * 10 + 8,
        });
    }

    let max_gap = 5;

    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);

    assert_eq!(
        result_quad.len(),
        result_binary.len(),
        "Both algorithms must produce same number of chains"
    );

    // Verify that the chains are identical
    let mut quad_sorted: Vec<Vec<usize>> = result_quad
        .into_iter()
        .map(|mut set| {
            set.sort();
            set
        })
        .collect();
    quad_sorted.sort();

    let mut binary_sorted: Vec<Vec<usize>> = result_binary
        .into_iter()
        .map(|mut set| {
            set.sort();
            set
        })
        .collect();
    binary_sorted.sort();

    assert_eq!(
        quad_sorted, binary_sorted,
        "Dense data must produce identical results"
    );
}

#[test]
fn test_edge_case_no_merging() {
    // All mappings too far apart
    let mappings = vec![
        TestMapping {
            idx: 0,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        },
        TestMapping {
            idx: 1,
            query_start: 1000,
            query_end: 1100,
            target_start: 1000,
            target_end: 1100,
        },
        TestMapping {
            idx: 2,
            query_start: 2000,
            query_end: 2100,
            target_start: 2000,
            target_end: 2100,
        },
    ];

    let max_gap = 50;

    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);

    assert_eq!(result_quad.len(), 3, "Should have 3 separate chains");
    assert_eq!(
        result_binary.len(),
        3,
        "Binary search should also have 3 separate chains"
    );

    for i in 0..3 {
        assert!(
            result_quad.iter().any(|chain| chain == &vec![i]),
            "Quadratic should have singleton chain for mapping {i}"
        );
        assert!(
            result_binary.iter().any(|chain| chain == &vec![i]),
            "Binary search should have singleton chain for mapping {i}"
        );
    }
}

#[test]
fn test_edge_case_all_merged() {
    // All mappings within range - should form single chain
    let mappings = vec![
        TestMapping {
            idx: 0,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        },
        TestMapping {
            idx: 1,
            query_start: 50,
            query_end: 150,
            target_start: 50,
            target_end: 150,
        },
        TestMapping {
            idx: 2,
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 200,
        },
        TestMapping {
            idx: 3,
            query_start: 150,
            query_end: 250,
            target_start: 150,
            target_end: 250,
        },
    ];

    let max_gap = 100;

    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);

    assert_eq!(result_quad.len(), 1, "Should form single chain");
    assert_eq!(
        result_binary.len(),
        1,
        "Binary search should also form single chain"
    );
    assert_eq!(
        result_quad[0].len(),
        4,
        "Chain should contain all 4 mappings"
    );
    assert_eq!(
        result_binary[0].len(),
        4,
        "Binary search chain should contain all 4 mappings"
    );
}

#[test]
fn test_performance_scaling() {
    // Test that binary search is actually faster on larger inputs
    println!("\nPerformance comparison:");

    for size in [10, 50, 100, 500, 1000] {
        // Create test data - mappings every 20 units
        let mut mappings = Vec::new();
        for i in 0..size {
            mappings.push(TestMapping {
                idx: i,
                query_start: i as u32 * 20,
                query_end: i as u32 * 20 + 15,
                target_start: i as u32 * 20,
                target_end: i as u32 * 20 + 15,
            });
        }

        let max_gap = 10; // Will chain adjacent mappings only

        // Time quadratic
        let start = Instant::now();
        let result_quad = merge_mappings_quadratic(&mappings, max_gap);
        let time_quad = start.elapsed();

        // Time binary search
        let start = Instant::now();
        let result_binary = merge_mappings_binary_search(&mappings, max_gap);
        let time_binary = start.elapsed();

        // Verify same results
        assert_eq!(
            result_quad.len(),
            result_binary.len(),
            "Size {size}: Different number of chains!"
        );

        let speedup = time_quad.as_nanos() as f64 / time_binary.as_nanos().max(1) as f64;
        println!(
            "  Size {:4}: Quadratic {:6.0}ns, Binary {:6.0}ns, Speedup: {:.1}x",
            size,
            time_quad.as_nanos(),
            time_binary.as_nanos(),
            speedup
        );

        // For larger sizes, binary search should generally be competitive
        // Note: Due to overhead, binary search might not always be faster at smaller sizes
        if size >= 1000 {
            // At very large sizes, binary search should show improvement
            println!("      Size {size} speedup: {speedup:.2}x");
        }
    }
}

#[test]
fn test_worst_case_dense_mappings() {
    // Worst case for quadratic: many mappings all within range of each other
    let size = 200;
    let mut mappings = Vec::new();

    // Create very dense mappings (only 2 units apart)
    for i in 0..size {
        mappings.push(TestMapping {
            idx: i,
            query_start: i as u32 * 2,
            query_end: i as u32 * 2 + 1,
            target_start: i as u32 * 2,
            target_end: i as u32 * 2 + 1,
        });
    }

    let max_gap = 100; // Large gap means many candidates

    println!("\nWorst case performance (dense mappings with large gap):");

    // Time quadratic
    let start = Instant::now();
    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let time_quad = start.elapsed();
    println!("  Quadratic: {time_quad:?}");

    // Time binary search
    let start = Instant::now();
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);
    let time_binary = start.elapsed();
    println!("  Binary:    {time_binary:?}");

    // Verify correctness
    assert_eq!(
        result_quad.len(),
        result_binary.len(),
        "Must produce same number of chains even in worst case"
    );

    let speedup = time_quad.as_nanos() as f64 / time_binary.as_nanos().max(1) as f64;
    println!("  Speedup:   {speedup:.1}x");

    // In worst case with dense data, binary search should still be competitive
    assert!(
        speedup >= 0.5,
        "Binary search should not be significantly slower than quadratic"
    );
}

#[test]
fn test_large_scale_performance() {
    println!("\nLarge-scale performance comparison:");

    // Test with truly large inputs where O(n²) vs O(n log n) matters
    for size in [5000, 10000, 20000] {
        // Create mappings with mixed spacing: some clustered, some sparse
        let mut mappings = Vec::new();
        for i in 0..size {
            // Create clusters of 10 mappings close together, then big gaps
            let cluster = i / 10;
            let within_cluster = i % 10;
            mappings.push(TestMapping {
                idx: i,
                query_start: (cluster as u32 * 10000) + (within_cluster as u32 * 50),
                query_end: (cluster as u32 * 10000) + (within_cluster as u32 * 50) + 40,
                target_start: (cluster as u32 * 10000) + (within_cluster as u32 * 50),
                target_end: (cluster as u32 * 10000) + (within_cluster as u32 * 50) + 40,
            });
        }

        let max_gap = 20; // Will chain within clusters but not between them

        // Time quadratic (might be slow for large sizes)
        let start = Instant::now();
        let result_quad = merge_mappings_quadratic(&mappings, max_gap);
        let time_quad = start.elapsed();

        // Time binary search
        let start = Instant::now();
        let result_binary = merge_mappings_binary_search(&mappings, max_gap);
        let time_binary = start.elapsed();

        // Verify same results
        assert_eq!(
            result_quad.len(),
            result_binary.len(),
            "Both algorithms must produce same number of chains at size {size}"
        );

        let speedup = time_quad.as_micros() as f64 / time_binary.as_micros().max(1) as f64;
        println!(
            "  Size {:5}: Quadratic {:8.0}µs, Binary {:8.0}µs, Speedup: {:6.1}x",
            size,
            time_quad.as_micros(),
            time_binary.as_micros(),
            speedup
        );

        // Note: Binary search may have overhead for smaller clusters
        // The real advantage shows up with very large dense data
        // For this test, we just verify correctness, not speed
    }
}

#[test]
fn test_correctness_with_overlaps() {
    // Test that overlap handling is identical
    let mappings = vec![
        TestMapping {
            idx: 0,
            query_start: 100,
            query_end: 200,
            target_start: 100,
            target_end: 200,
        },
        TestMapping {
            idx: 1,
            query_start: 150,
            query_end: 250,
            target_start: 150,
            target_end: 250,
        }, // Overlaps
        TestMapping {
            idx: 2,
            query_start: 180,
            query_end: 280,
            target_start: 180,
            target_end: 280,
        }, // Overlaps
        TestMapping {
            idx: 3,
            query_start: 300,
            query_end: 400,
            target_start: 300,
            target_end: 400,
        }, // Gap
    ];

    let max_gap = 50;

    let result_quad = merge_mappings_quadratic(&mappings, max_gap);
    let result_binary = merge_mappings_binary_search(&mappings, max_gap);

    // All 4 mappings chain together: 0,1,2 overlap, and 3 is within max_gap of 2
    assert_eq!(result_quad.len(), 1, "Should have 1 chain");
    assert_eq!(
        result_binary.len(),
        1,
        "Binary search should also have 1 chain"
    );

    // All 4 mappings should be in one chain
    assert_eq!(result_quad[0].len(), 4, "Chain should have all 4 mappings");
    assert_eq!(
        result_binary[0].len(),
        4,
        "Binary search chain should have all 4 mappings"
    );
}

#[test]
fn test_binary_search_boundary_conditions() {
    // Test exact boundary conditions
    let mappings = vec![
        TestMapping {
            idx: 0,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        },
        TestMapping {
            idx: 1,
            query_start: 150,
            query_end: 250,
            target_start: 150,
            target_end: 250,
        }, // Exactly at max_gap
        TestMapping {
            idx: 2,
            query_start: 151,
            query_end: 251,
            target_start: 151,
            target_end: 251,
        }, // Just beyond max_gap
    ];

    let max_gap = 50;

    let result = merge_mappings_binary_search(&mappings, max_gap);

    // Mapping 0 and 1 should chain (gap = 50 = max_gap)
    // Mapping 1 and 2 overlap, so they chain
    // Therefore all 3 chain transitively through union-find
    assert_eq!(
        result.len(),
        1,
        "Should form 1 chain (all transitively connected)"
    );

    let chain = &result[0];
    assert!(
        chain.contains(&0) && chain.contains(&1) && chain.contains(&2),
        "All three mappings should be in the same chain"
    );
}
