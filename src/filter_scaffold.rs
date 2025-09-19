use anyhow::Result;
use crate::mapping::{Mapping, MappingAux, ChainStatus};
use crate::filter::FilterConfig;
use std::collections::HashMap;

/// Represents a merged scaffold chain with its boundaries
#[derive(Debug, Clone)]
struct ScaffoldChain {
    ref_seq_id: u32,
    query_start: u32,
    query_end: u32,
    ref_start: u32,
    ref_end: u32,
    total_length: u32,
    is_reverse: bool,
}

/// Scaffold filtering implementation matching wfmash's algorithm
pub fn filter_by_scaffolds(
    mappings: &[Mapping],
    aux_data: &[MappingAux],
    config: &FilterConfig,
) -> Result<(Vec<Mapping>, Vec<MappingAux>)> {
    if config.min_scaffold_length == 0 || mappings.is_empty() {
        return Ok((mappings.to_vec(), aux_data.to_vec()));
    }

    eprintln!("[SCAFFOLD_TRACE] Input: {} mappings", mappings.len());

    // Group mappings by query-ref pair
    let mut groups: HashMap<(i32, u32), Vec<usize>> = HashMap::new();
    for (i, mapping) in mappings.iter().enumerate() {
        let key = (aux_data[i].query_seq_id, mapping.ref_seq_id);
        groups.entry(key).or_insert_with(Vec::new).push(i);
    }

    let mut all_kept_indices = Vec::new();
    let mut all_kept_aux = Vec::new();

    // Process each query-ref group
    for ((query_id, ref_id), mut indices) in groups {
        eprintln!("[SCAFFOLD_TRACE] Processing query {} vs ref {}: {} mappings",
                  query_id, ref_id, indices.len());

        // Sort by position
        indices.sort_by_key(|&i| (mappings[i].ref_start_pos, mappings[i].query_start_pos));

        // Keep original mappings for this group
        let original_mappings: Vec<(usize, Mapping, MappingAux)> = indices.iter()
            .map(|&i| (i, mappings[i], aux_data[i].clone()))
            .collect();

        // Step 1: Build merged chains
        let merged_chains = merge_into_chains(&original_mappings, config.scaffold_gap);

        eprintln!("[SCAFFOLD_TRACE] After merging with gap={}: {} chains",
                  config.scaffold_gap, merged_chains.len());
        for (i, chain) in merged_chains.iter().enumerate() {
            eprintln!("[SCAFFOLD_TRACE]   Chain[{}]: ref={} q={}-{} r={}-{} len={}",
                      i, chain.ref_seq_id, chain.query_start, chain.query_end,
                      chain.ref_start, chain.ref_end, chain.total_length);
        }

        // Step 2: Filter by minimum scaffold length
        let scaffold_chains: Vec<_> = merged_chains.into_iter()
            .filter(|chain| chain.total_length >= config.min_scaffold_length)
            .collect();

        eprintln!("[SCAFFOLD_TRACE] After length filter (min={}): {} scaffolds",
                  config.min_scaffold_length, scaffold_chains.len());

        // Step 3: Apply plane sweep filter (simplified for now - TODO: implement properly)
        let filtered_scaffolds = scaffold_chains; // TODO: implement plane sweep

        eprintln!("[SCAFFOLD_TRACE] After plane sweep: {} chains", filtered_scaffolds.len());
        for (i, chain) in filtered_scaffolds.iter().enumerate() {
            eprintln!("[SCAFFOLD_TRACE]   Scaffold[{}]: ref={} q={}-{} r={}-{} len={}",
                      i, chain.ref_seq_id, chain.query_start, chain.query_end,
                      chain.ref_start, chain.ref_end, chain.total_length);
        }

        // Step 4: Identify anchors - original mappings within scaffold bounds
        let mut anchor_indices = Vec::new();
        for &(orig_idx, ref orig_mapping, _) in &original_mappings {
            for scaffold in &filtered_scaffolds {
                if orig_mapping.ref_seq_id == scaffold.ref_seq_id &&
                   orig_mapping.is_reverse() == scaffold.is_reverse &&
                   orig_mapping.query_start_pos >= scaffold.query_start &&
                   orig_mapping.query_end_pos() <= scaffold.query_end &&
                   orig_mapping.ref_start_pos >= scaffold.ref_start &&
                   orig_mapping.ref_end_pos() <= scaffold.ref_end {
                    anchor_indices.push(orig_idx);
                    break;
                }
            }
        }

        eprintln!("[SCAFFOLD_TRACE] Anchors identified: {} mappings", anchor_indices.len());

        // Step 5: Calculate distances and rescue nearby mappings
        let mut kept_indices = Vec::new();
        let mut kept_aux = Vec::new();

        // Add all anchors
        for &idx in &anchor_indices {
            kept_indices.push(idx);
            let mut aux = aux_data[idx].clone();
            aux.chain_status = ChainStatus::Scaffold;
            kept_aux.push(aux);
        }

        // Check non-anchors for rescue
        let mut rescued_count = 0;
        for &(orig_idx, ref orig_mapping, _) in &original_mappings {
            if anchor_indices.contains(&orig_idx) {
                continue; // Already an anchor
            }

            // Calculate Euclidean distance to nearest anchor
            let mut min_distance = f64::INFINITY;
            for &anchor_idx in &anchor_indices {
                let anchor = &mappings[anchor_idx];

                // Calculate midpoints
                let orig_q_mid = orig_mapping.query_start_pos as f64 +
                                (orig_mapping.block_length as f64 / 2.0);
                let orig_r_mid = orig_mapping.ref_start_pos as f64 +
                                (orig_mapping.block_length as f64 / 2.0);
                let anchor_q_mid = anchor.query_start_pos as f64 +
                                  (anchor.block_length as f64 / 2.0);
                let anchor_r_mid = anchor.ref_start_pos as f64 +
                                  (anchor.block_length as f64 / 2.0);

                // Euclidean distance
                let q_diff = orig_q_mid - anchor_q_mid;
                let r_diff = orig_r_mid - anchor_r_mid;
                let distance = (q_diff * q_diff + r_diff * r_diff).sqrt();

                min_distance = min_distance.min(distance);
            }

            // Rescue if within threshold
            if min_distance <= config.scaffold_max_deviation as f64 {
                kept_indices.push(orig_idx);
                let mut aux = aux_data[orig_idx].clone();
                aux.chain_status = ChainStatus::Rescued;
                kept_aux.push(aux);
                rescued_count += 1;

                eprintln!("[SCAFFOLD_TRACE] Mapping[{}] distance={:.0} rescued",
                          orig_idx, min_distance);
            }
        }

        eprintln!("[SCAFFOLD_TRACE] Final for group: {} -> {} (anchors={}, rescued={})",
                  original_mappings.len(), kept_indices.len(),
                  anchor_indices.len(), rescued_count);

        all_kept_indices.extend(kept_indices);
        all_kept_aux.extend(kept_aux);
    }

    // Build final result preserving original order
    let mut result_mappings = Vec::new();
    let mut result_aux = Vec::new();

    for (i, aux) in all_kept_aux.into_iter().enumerate() {
        let idx = all_kept_indices[i];
        result_mappings.push(mappings[idx]);
        result_aux.push(aux);
    }

    eprintln!("[SCAFFOLD_TRACE] Final: {} -> {} mappings",
              mappings.len(), result_mappings.len());

    Ok((result_mappings, result_aux))
}

/// Merge mappings into chains based on gap distance
fn merge_into_chains(
    mappings: &[(usize, Mapping, MappingAux)],
    max_gap: u32,
) -> Vec<ScaffoldChain> {
    if mappings.is_empty() {
        return Vec::new();
    }

    let mut chains = Vec::new();
    let mut current_chain: Option<ScaffoldChain> = None;

    for &(_, ref mapping, _) in mappings {
        match current_chain {
            None => {
                // Start new chain
                current_chain = Some(ScaffoldChain {
                    ref_seq_id: mapping.ref_seq_id,
                    query_start: mapping.query_start_pos,
                    query_end: mapping.query_end_pos(),
                    ref_start: mapping.ref_start_pos,
                    ref_end: mapping.ref_end_pos(),
                    total_length: mapping.block_length,
                    is_reverse: mapping.is_reverse(),
                });
            }
            Some(ref mut chain) => {
                // Check if can extend current chain
                let ref_gap = mapping.ref_start_pos.saturating_sub(chain.ref_end);
                let query_gap = mapping.query_start_pos.saturating_sub(chain.query_end);

                if mapping.ref_seq_id == chain.ref_seq_id &&
                   mapping.is_reverse() == chain.is_reverse &&
                   ref_gap <= max_gap &&
                   query_gap <= max_gap {
                    // Extend chain
                    chain.query_end = mapping.query_end_pos();
                    chain.ref_end = mapping.ref_end_pos();
                    chain.total_length = chain.ref_end - chain.ref_start;
                } else {
                    // Save current chain and start new one
                    chains.push(chain.clone());
                    *chain = ScaffoldChain {
                        ref_seq_id: mapping.ref_seq_id,
                        query_start: mapping.query_start_pos,
                        query_end: mapping.query_end_pos(),
                        ref_start: mapping.ref_start_pos,
                        ref_end: mapping.ref_end_pos(),
                        total_length: mapping.block_length,
                        is_reverse: mapping.is_reverse(),
                    };
                }
            }
        }
    }

    // Don't forget the last chain
    if let Some(chain) = current_chain {
        chains.push(chain);
    }

    chains
}