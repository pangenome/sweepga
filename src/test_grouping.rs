#[cfg(test)]
mod tests {
    use crate::compact_filter::{CompactMapping, CompactPafFilter};
    use crate::plane_sweep_exact::{plane_sweep_both, plane_sweep_query, PlaneSweepMapping};
    use std::collections::HashMap;

    /// Helper to create a compact mapping with names
    fn create_mapping(
        query_name: &str,
        target_name: &str,
        query_start: u32,
        query_end: u32,
        target_start: u32,
        target_end: u32,
        identity: f64,
    ) -> (CompactMapping, String, String) {
        let mapping = CompactMapping {
            query_id: 0,  // Will be set by test
            target_id: 0, // Will be set by test
            query_start,
            query_end,
            target_start,
            target_end,
            identity: (identity * 10000.0) as u16,
            strand: 0,
            flags: 0,
            file_offset: 0,
        };
        (mapping, query_name.to_string(), target_name.to_string())
    }

    /// Group mappings by prefix pairs
    fn group_by_prefix_pairs(
        mappings: &[(CompactMapping, String, String)],
        delimiter: char,
        use_prefix: bool,
    ) -> HashMap<(String, String), Vec<usize>> {
        let mut groups = HashMap::new();

        for (idx, (_mapping, query_name, target_name)) in mappings.iter().enumerate() {
            let query_prefix = if use_prefix {
                if let Some(pos) = query_name.rfind(delimiter) {
                    query_name[..=pos].to_string()
                } else {
                    query_name.clone()
                }
            } else {
                query_name.clone()
            };

            let target_prefix = if use_prefix {
                if let Some(pos) = target_name.rfind(delimiter) {
                    target_name[..=pos].to_string()
                } else {
                    target_name.clone()
                }
            } else {
                target_name.clone()
            };

            groups
                .entry((query_prefix, target_prefix))
                .or_insert_with(Vec::new)
                .push(idx);
        }

        groups
    }

    #[test]
    fn test_basic_grouping() {
        // Create mappings between different genome pairs
        let mappings = vec![
            create_mapping("genomeA#chr1", "genomeB#chr1", 100, 200, 100, 200, 0.95),
            create_mapping("genomeA#chr1", "genomeB#chr2", 300, 400, 300, 400, 0.90),
            create_mapping("genomeA#chr2", "genomeB#chr1", 500, 600, 500, 600, 0.92),
            create_mapping("genomeC#chr1", "genomeD#chr1", 100, 200, 100, 200, 0.98),
        ];

        // Group with prefix
        let groups = group_by_prefix_pairs(&mappings, '#', true);

        // Should have 2 groups: A->B and C->D
        // All A->B mappings go into one group regardless of chromosome
        assert_eq!(groups.len(), 2);

        // Check genomeA# -> genomeB# group
        let ab_group = groups.get(&("genomeA#".to_string(), "genomeB#".to_string()));
        assert!(ab_group.is_some());
        assert_eq!(ab_group.unwrap().len(), 3); // All 3 A->B mappings (chr1->chr1, chr1->chr2, chr2->chr1)

        // Check genomeC# -> genomeD# group
        let cd_group = groups.get(&("genomeC#".to_string(), "genomeD#".to_string()));
        assert!(cd_group.is_some());
        assert_eq!(cd_group.unwrap().len(), 1);
    }

    #[test]
    fn test_no_prefix_grouping() {
        // Test without prefix grouping - each chromosome pair is separate
        let mappings = vec![
            create_mapping("genomeA#chr1", "genomeB#chr1", 100, 200, 100, 200, 0.95),
            create_mapping("genomeA#chr1", "genomeB#chr2", 300, 400, 300, 400, 0.90),
            create_mapping("genomeA#chr2", "genomeB#chr1", 500, 600, 500, 600, 0.92),
        ];

        // Group without prefix (each chromosome is separate)
        let groups = group_by_prefix_pairs(&mappings, '#', false);

        // Should have 3 separate groups
        assert_eq!(groups.len(), 3);

        // Each group should have exactly 1 mapping
        for (_key, indices) in &groups {
            assert_eq!(indices.len(), 1);
        }
    }

    #[test]
    fn test_mixed_delimiter() {
        // Test with mixed delimiters and formats
        let mappings = vec![
            create_mapping(
                "human_haplotype1_chr1",
                "human_haplotype2_chr1",
                100,
                200,
                100,
                200,
                0.95,
            ),
            create_mapping(
                "human_haplotype1_chr2",
                "human_haplotype2_chr2",
                300,
                400,
                300,
                400,
                0.90,
            ),
            create_mapping(
                "mouse_strain1_chr1",
                "mouse_strain2_chr1",
                500,
                600,
                500,
                600,
                0.92,
            ),
        ];

        // Group with underscore delimiter (last occurrence)
        let groups = group_by_prefix_pairs(&mappings, '_', true);

        // Should have 2 groups: human_haplotype1_ -> human_haplotype2_ and mouse_strain1_ -> mouse_strain2_
        assert_eq!(groups.len(), 2);

        let human_group = groups.get(&(
            "human_haplotype1_".to_string(),
            "human_haplotype2_".to_string(),
        ));
        assert!(human_group.is_some());
        assert_eq!(human_group.unwrap().len(), 2);

        let mouse_group = groups.get(&("mouse_strain1_".to_string(), "mouse_strain2_".to_string()));
        assert!(mouse_group.is_some());
        assert_eq!(mouse_group.unwrap().len(), 1);
    }

    #[test]
    fn test_self_mappings() {
        // Test self-mappings within same genome
        let mappings = vec![
            create_mapping("genomeA#chr1", "genomeA#chr1", 100, 200, 300, 400, 0.99), // Self, different region
            create_mapping("genomeA#chr1", "genomeA#chr2", 100, 200, 100, 200, 0.95), // Same genome, diff chr
            create_mapping("genomeA#chr1", "genomeB#chr1", 100, 200, 100, 200, 0.90), // Different genome
        ];

        let groups = group_by_prefix_pairs(&mappings, '#', true);

        // Should have 2 groups: A->A and A->B
        assert_eq!(groups.len(), 2);

        let aa_group = groups.get(&("genomeA#".to_string(), "genomeA#".to_string()));
        assert!(aa_group.is_some());
        assert_eq!(aa_group.unwrap().len(), 2); // Both self-mappings

        let ab_group = groups.get(&("genomeA#".to_string(), "genomeB#".to_string()));
        assert!(ab_group.is_some());
        assert_eq!(ab_group.unwrap().len(), 1);
    }

    #[test]
    fn test_plane_sweep_within_groups() {
        // Test that plane sweep works correctly within groups
        let mappings = vec![
            // Group 1: genomeA# -> genomeB# (overlapping mappings)
            create_mapping("genomeA#chr1", "genomeB#chr1", 100, 500, 100, 500, 0.95),
            create_mapping("genomeA#chr1", "genomeB#chr1", 200, 600, 200, 600, 0.90), // Overlaps
            create_mapping("genomeA#chr1", "genomeB#chr1", 700, 800, 700, 800, 0.85), // Non-overlapping
            // Group 2: genomeC# -> genomeD# (should not affect Group 1)
            create_mapping("genomeC#chr1", "genomeD#chr1", 100, 500, 100, 500, 0.98),
            create_mapping("genomeC#chr1", "genomeD#chr1", 200, 600, 200, 600, 0.97), // Overlaps
        ];

        let groups = group_by_prefix_pairs(&mappings, '#', true);
        assert_eq!(groups.len(), 2);

        // Process each group with plane sweep (1:1 filtering)
        for ((_q_prefix, _t_prefix), indices) in groups {
            let mut group_mappings: Vec<PlaneSweepMapping> = indices
                .iter()
                .map(|&idx| {
                    let (m, _, _) = &mappings[idx];
                    PlaneSweepMapping {
                        idx,
                        query_start: m.query_start,
                        query_end: m.query_end,
                        target_start: m.target_start,
                        target_end: m.target_end,
                        identity: (m.identity as f64) / 10000.0,
                        flags: 0,
                    }
                })
                .collect();

            // Apply query-only plane sweep (0 secondaries = keep 1 best per position)
            let kept = plane_sweep_query(&mut group_mappings, 0, 0.95);

            // With overlap threshold 0.95, mappings are considered overlapping if >95% overlap
            // The mappings at 100-500 and 200-600 have (500-200)/(600-200) = 300/400 = 75% overlap
            // So they should NOT be filtered by overlap

            if indices.len() == 3 {
                // Group 1 (A->B): 3 mappings, two partially overlapping + one separate
                // With 0 secondaries and 95% overlap threshold:
                // - Best at each position is kept
                // Since 75% < 95%, both overlapping mappings may be kept
                // Plus the non-overlapping one
                // Actually, with proper plane sweep, at position 100-500, mapping 0 is best
                // At position 200-600, mapping 1 starts but 0 is still active and better
                // So mapping 1 should be discarded
                // At position 700-800, only mapping 2 is active
                assert!(kept.len() >= 2, "Group 1: should keep at least 2");
            } else if indices.len() == 2 {
                // Group 2 (C->D): 2 overlapping mappings
                // The overlap is 75%, which is < 95% threshold
                // But with 0 secondaries, only best at each position is kept
                // Since they overlap, the better one (0.98) should dominate
                assert!(kept.len() >= 1, "Group 2: should keep at least 1");
            }
        }
    }

    #[test]
    fn test_complex_pan_genome() {
        // Simulate a pan-genome with multiple haplotypes
        let mappings = vec![
            // SGDref to various strains
            create_mapping("SGDref#1#chrI", "DBVPG6044#1#chrI", 0, 1000, 0, 1000, 0.98),
            create_mapping(
                "SGDref#1#chrI",
                "DBVPG6044#1#chrI",
                2000,
                3000,
                2000,
                3000,
                0.96,
            ),
            create_mapping("SGDref#1#chrI", "SK1#1#chrI", 0, 1000, 0, 1000, 0.97),
            create_mapping("SGDref#1#chrI", "SK1#1#chrI", 500, 1500, 500, 1500, 0.95), // Overlaps
            // Different haplotype
            create_mapping("SGDref#2#chrI", "DBVPG6044#2#chrI", 0, 1000, 0, 1000, 0.99),
            // Different chromosome
            create_mapping(
                "SGDref#1#chrII",
                "DBVPG6044#1#chrII",
                0,
                1000,
                0,
                1000,
                0.94,
            ),
        ];

        // Group with double delimiter (PanSN format)
        let groups = group_by_prefix_pairs(&mappings, '#', true);

        // Count unique prefix pairs
        let mut unique_pairs = std::collections::HashSet::new();
        for (pair, _) in &groups {
            unique_pairs.insert(pair.clone());
        }

        // We should have groups for:
        // SGDref#1# -> DBVPG6044#1#
        // SGDref#1# -> SK1#1#
        // SGDref#2# -> DBVPG6044#2#
        assert_eq!(unique_pairs.len(), 3);

        // Test filtering within SGDref#1# -> DBVPG6044#1# group
        let sgd_dbvpg = groups.get(&("SGDref#1#".to_string(), "DBVPG6044#1#".to_string()));
        assert!(sgd_dbvpg.is_some());
        assert_eq!(sgd_dbvpg.unwrap().len(), 3); // chrI (2) + chrII (1)
    }

    #[test]
    fn test_no_delimiter_in_names() {
        // Test when names have no delimiter
        let mappings = vec![
            create_mapping("chr1", "chr2", 100, 200, 100, 200, 0.95),
            create_mapping("chr1", "chr3", 300, 400, 300, 400, 0.90),
            create_mapping("chr2", "chr3", 500, 600, 500, 600, 0.92),
        ];

        // Try to group with delimiter that doesn't exist
        let groups = group_by_prefix_pairs(&mappings, '#', true);

        // Each name becomes its own prefix (no delimiter found)
        assert_eq!(groups.len(), 3);

        // Each group should be a unique chromosome pair
        assert!(groups.contains_key(&("chr1".to_string(), "chr2".to_string())));
        assert!(groups.contains_key(&("chr1".to_string(), "chr3".to_string())));
        assert!(groups.contains_key(&("chr2".to_string(), "chr3".to_string())));
    }

    #[test]
    fn test_filtering_preserves_groups() {
        // Test that filtering respects group boundaries
        let mappings = vec![
            // Group A->B: High quality mappings
            create_mapping("genomeA#chr1", "genomeB#chr1", 100, 1000, 100, 1000, 0.99),
            create_mapping("genomeA#chr1", "genomeB#chr1", 2000, 3000, 2000, 3000, 0.98),
            // Group C->D: Lower quality mappings
            create_mapping("genomeC#chr1", "genomeD#chr1", 100, 1000, 100, 1000, 0.85),
            create_mapping("genomeC#chr1", "genomeD#chr1", 2000, 3000, 2000, 3000, 0.84),
        ];

        let groups = group_by_prefix_pairs(&mappings, '#', true);

        // Process each group with 1:1 filtering
        let mut total_kept = 0;
        for ((_q, _t), indices) in groups {
            let mut group_mappings: Vec<PlaneSweepMapping> = indices
                .iter()
                .map(|&idx| {
                    let (m, _, _) = &mappings[idx];
                    PlaneSweepMapping {
                        idx,
                        query_start: m.query_start,
                        query_end: m.query_end,
                        target_start: m.target_start,
                        target_end: m.target_end,
                        identity: (m.identity as f64) / 10000.0,
                        flags: 0,
                    }
                })
                .collect();

            let kept = plane_sweep_query(&mut group_mappings, 0, 0.95);
            total_kept += kept.len();
        }

        // Each group should keep its best mappings independently
        // Even though C->D has lower quality, they should still be kept within their group
        assert_eq!(total_kept, 4); // All mappings kept (non-overlapping within groups)
    }
}
