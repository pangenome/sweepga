use anyhow::Result;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Alignment {
    pub query_name: String,
    pub query_start: u32,
    pub query_end: u32,
    pub query_len: u32,
    pub target_name: String,
    pub target_start: u32,
    pub target_end: u32,
    pub target_len: u32,
    pub matches: u32,
    pub block_length: u32,
    pub strand: char,
    pub tags: Vec<String>,
}

/// Merge overlapping and adjacent alignments
pub fn merge_alignments(alignments: Vec<Alignment>) -> Vec<Alignment> {
    // Group by query-target-strand
    let mut groups: HashMap<(String, String, char), Vec<Alignment>> = HashMap::new();

    for aln in alignments {
        let key = (aln.query_name.clone(), aln.target_name.clone(), aln.strand);
        groups.entry(key).or_insert_with(Vec::new).push(aln);
    }

    let mut merged = Vec::new();

    for ((_q, _t, strand), mut group) in groups {
        // Sort by query start position
        group.sort_by_key(|a| a.query_start);

        let mut current = group[0].clone();

        for next in &group[1..] {
            // Check if alignments can be merged (overlapping or adjacent)
            let query_gap = if next.query_start > current.query_end {
                next.query_start - current.query_end
            } else {
                0
            };

            let target_gap = if strand == '+' {
                if next.target_start > current.target_end {
                    next.target_start - current.target_end
                } else {
                    0
                }
            } else {
                // For reverse strand, target coordinates go backwards
                if current.target_start > next.target_end {
                    current.target_start - next.target_end
                } else {
                    0
                }
            };

            // Merge if gaps are small enough (say < 100bp)
            if query_gap < 100 && target_gap < 100 {
                // Extend current alignment
                current.query_end = current.query_end.max(next.query_end);
                current.block_length = current.query_end - current.query_start;

                if strand == '+' {
                    current.target_end = current.target_end.max(next.target_end);
                } else {
                    current.target_start = current.target_start.max(next.target_start);
                }

                // Approximate match count
                current.matches += next.matches;
            } else {
                // Gap too large, output current and start new
                merged.push(current.clone());
                current = next.clone();
            }
        }

        merged.push(current);
    }

    merged
}

/// Chain alignments into scaffold chains
pub fn chain_alignments(alignments: Vec<Alignment>, max_gap: u32) -> Vec<Vec<Alignment>> {
    // Group by query-target pair
    let mut groups: HashMap<(String, String), Vec<Alignment>> = HashMap::new();

    for aln in alignments {
        let key = (aln.query_name.clone(), aln.target_name.clone());
        groups.entry(key).or_insert_with(Vec::new).push(aln);
    }

    let mut all_chains = Vec::new();

    for (_, mut group) in groups {
        // Sort by query position
        group.sort_by_key(|a| a.query_start);

        let mut chains: Vec<Vec<Alignment>> = Vec::new();
        let mut current_chain = vec![group[0].clone()];

        for next in &group[1..] {
            let last = current_chain.last().unwrap();
            let gap = if next.query_start > last.query_end {
                next.query_start - last.query_end
            } else {
                0
            };

            if gap <= max_gap {
                current_chain.push(next.clone());
            } else {
                if !current_chain.is_empty() {
                    chains.push(current_chain);
                }
                current_chain = vec![next.clone()];
            }
        }

        if !current_chain.is_empty() {
            chains.push(current_chain);
        }

        all_chains.extend(chains);
    }

    all_chains
}