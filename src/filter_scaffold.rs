use anyhow::Result;
use crate::mapping::{Mapping, MappingAux, ChainStatus};
use crate::filter::{FilterConfig, FilterMode};
use crate::plane_sweep_scaffold::{ScaffoldLike, plane_sweep_scaffolds};
use crate::paf_filter::ScoringFunction;
use std::collections::HashMap;

/// Represents a merged scaffold chain with its boundaries
#[derive(Debug, Clone)]
struct ScaffoldChain {
    query_name: String,
    ref_name: String,
    ref_seq_id: u32,
    query_start: u32,
    query_end: u32,
    ref_start: u32,
    ref_end: u32,
    total_length: u32,
    identity: f64,
    is_reverse: bool,
}

impl ScaffoldLike for ScaffoldChain {
    fn query_name(&self) -> &str {
        &self.query_name
    }
    fn target_name(&self) -> &str {
        &self.ref_name
    }
    fn query_start(&self) -> u64 {
        self.query_start as u64
    }
    fn query_end(&self) -> u64 {
        self.query_end as u64
    }
    fn target_start(&self) -> u64 {
        self.ref_start as u64
    }
    fn target_end(&self) -> u64 {
        self.ref_end as u64
    }
    fn identity(&self) -> f64 {
        self.identity
    }
}

/// Scaffold filtering implementation matching wfmash's algorithm
/// This version identifies scaffolds per chromosome pair, not per prefix pair
pub fn filter_by_scaffolds(
    mappings: &[Mapping],
    aux_data: &[MappingAux],
    config: &FilterConfig,
) -> Result<(Vec<Mapping>, Vec<MappingAux>)> {
    // For backwards compatibility, use the same mappings for both anchor and rescue
    filter_by_scaffolds_with_rescue(mappings, aux_data, mappings, aux_data, config)
}

/// Scaffold filtering with separate anchor and rescue sets
/// Identifies scaffolds from filtered_mappings but rescues from raw_mappings
pub fn filter_by_scaffolds_with_rescue(
    filtered_mappings: &[Mapping],  // Use these to identify scaffolds
    filtered_aux: &[MappingAux],
    raw_mappings: &[Mapping],       // Rescue mappings from this set
    raw_aux: &[MappingAux],
    config: &FilterConfig,
) -> Result<(Vec<Mapping>, Vec<MappingAux>)> {
    if config.min_scaffold_length == 0 || filtered_mappings.is_empty() {
        return Ok((raw_mappings.to_vec(), raw_aux.to_vec()));
    }

    // eprintln!("[SCAFFOLD_TRACE] Identifying scaffolds from {} filtered mappings, rescuing from {} raw mappings",
    //           filtered_mappings.len(), raw_mappings.len());

    // Group FILTERED mappings by query-ref CHROMOSOME pair to identify scaffolds
    let mut groups: HashMap<(i32, u32), Vec<usize>> = HashMap::new();
    for (i, mapping) in filtered_mappings.iter().enumerate() {
        let key = (filtered_aux[i].query_seq_id, mapping.ref_seq_id);
        groups.entry(key).or_insert_with(Vec::new).push(i);
    }

    // eprintln!("[SCAFFOLD_TRACE] Grouped into {} chromosome pairs", groups.len());

    let mut all_kept_indices = Vec::new();
    let mut all_kept_aux = Vec::new();

    // Process each query-ref chromosome pair independently
    for ((query_id, ref_id), mut indices) in groups {
        // Get the actual chromosome names for debugging
        let query_name = if !indices.is_empty() {
            &filtered_aux[indices[0]].query_name
        } else {
            "unknown"
        };
        let ref_name = if !indices.is_empty() {
            &filtered_aux[indices[0]].ref_name
        } else {
            "unknown"
        };

        // eprintln!("[SCAFFOLD_TRACE] Processing {} (id={}) vs {} (id={}): {} filtered mappings",
        //           query_name, query_id, ref_name, ref_id, indices.len());

        // Sort by position
        indices.sort_by_key(|&i| (filtered_mappings[i].ref_start_pos, filtered_mappings[i].query_start_pos));

        // Collect FILTERED mappings for scaffold identification
        let filtered_for_scaffolds: Vec<(usize, Mapping, MappingAux)> = indices.iter()
            .map(|&i| (i, filtered_mappings[i], filtered_aux[i].clone()))
            .collect();

        // Step 1: Build merged chains FROM FILTERED MAPPINGS
        let merged_chains = merge_into_chains(&filtered_for_scaffolds, config.scaffold_gap);

        // eprintln!("[SCAFFOLD_TRACE] After merging with gap={}: {} chains",
        //           config.scaffold_gap, merged_chains.len());
        for (i, chain) in merged_chains.iter().enumerate() {
            // eprintln!("[SCAFFOLD_TRACE]   Chain[{}]: ref={} q={}-{} r={}-{} len={}",
            //           i, chain.ref_seq_id, chain.query_start, chain.query_end,
            //           chain.ref_start, chain.ref_end, chain.total_length);
        }

        // Step 2: Filter by minimum scaffold length
        let scaffold_chains: Vec<_> = merged_chains.into_iter()
            .filter(|chain| chain.total_length >= config.min_scaffold_length)
            .collect();

        // eprintln!("[SCAFFOLD_TRACE] After length filter (min={}): {} scaffolds",
        //           config.min_scaffold_length, scaffold_chains.len());

        // Step 3: Apply plane sweep filter to scaffolds
        let kept_indices = plane_sweep_scaffolds(
            &scaffold_chains,
            config.filter_mode,
            config.max_per_query,
            config.max_per_target,
            config.scaffold_overlap_threshold,
            ScoringFunction::LogLengthIdentity, // Use standard scoring
        )?;

        let filtered_scaffolds: Vec<ScaffoldChain> = kept_indices
            .iter()
            .map(|&idx| scaffold_chains[idx].clone())
            .collect();

        // eprintln!("[SCAFFOLD_TRACE] After plane sweep: {} -> {} chains",
        //           scaffold_chains.len(), filtered_scaffolds.len());
        for (i, chain) in filtered_scaffolds.iter().enumerate() {
            // eprintln!("[SCAFFOLD_TRACE]   Scaffold[{}]: ref={} q={}-{} r={}-{} len={}",
            //           i, chain.ref_seq_id, chain.query_start, chain.query_end,
            //           chain.ref_start, chain.ref_end, chain.total_length);
        }

        // Step 4: Identify anchors - FILTERED mappings within scaffold bounds
        let mut anchor_indices = Vec::new();
        for &(orig_idx, ref orig_mapping, _) in &filtered_for_scaffolds {
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

        // eprintln!("[SCAFFOLD_TRACE] Anchors identified: {} mappings from filtered set", anchor_indices.len());

        // Step 5: NOW RESCUE FROM RAW MAPPINGS
        // Collect all raw mappings for this chromosome pair
        let mut raw_indices_for_pair = Vec::new();
        for (i, mapping) in raw_mappings.iter().enumerate() {
            if raw_aux[i].query_seq_id == query_id && mapping.ref_seq_id == ref_id {
                raw_indices_for_pair.push(i);
            }
        }

        // eprintln!("[SCAFFOLD_TRACE] Found {} raw mappings for this chromosome pair to consider for rescue",
        //           raw_indices_for_pair.len());

        // Collect anchor mappings from FILTERED set for distance calculation
        let anchor_mappings: Vec<&Mapping> = anchor_indices.iter()
            .map(|&idx| &filtered_mappings[idx])
            .collect();

        // Calculate distances and rescue nearby mappings FROM RAW SET
        let mut kept_indices = Vec::new();
        let mut kept_aux = Vec::new();

        // Check each raw mapping
        let mut rescued_count = 0;
        let mut anchor_count = 0;

        for &raw_idx in &raw_indices_for_pair {
            let raw_mapping = &raw_mappings[raw_idx];

            // Check if this raw mapping is actually one of the anchors
            // (it might be in both filtered and raw sets)
            let is_anchor = anchor_mappings.iter().any(|anchor| {
                anchor.query_start_pos == raw_mapping.query_start_pos &&
                anchor.ref_start_pos == raw_mapping.ref_start_pos &&
                anchor.block_length == raw_mapping.block_length
            });

            if is_anchor {
                // This is an anchor
                kept_indices.push(raw_idx);
                let mut aux = raw_aux[raw_idx].clone();
                aux.chain_status = ChainStatus::Scaffold;
                kept_aux.push(aux);
                anchor_count += 1;
                continue;
            }

            // Calculate Euclidean distance to nearest anchor
            let mut min_distance = f64::INFINITY;
            for anchor in &anchor_mappings {

                // Calculate midpoints
                let raw_q_mid = raw_mapping.query_start_pos as f64 +
                               (raw_mapping.block_length as f64 / 2.0);
                let raw_r_mid = raw_mapping.ref_start_pos as f64 +
                               (raw_mapping.block_length as f64 / 2.0);
                let anchor_q_mid = anchor.query_start_pos as f64 +
                                  (anchor.block_length as f64 / 2.0);
                let anchor_r_mid = anchor.ref_start_pos as f64 +
                                  (anchor.block_length as f64 / 2.0);

                // Euclidean distance
                let q_diff = raw_q_mid - anchor_q_mid;
                let r_diff = raw_r_mid - anchor_r_mid;
                let distance = (q_diff * q_diff + r_diff * r_diff).sqrt();

                min_distance = min_distance.min(distance);
            }

            // Rescue if within threshold
            if min_distance <= config.scaffold_max_deviation as f64 {
                kept_indices.push(raw_idx);
                let mut aux = raw_aux[raw_idx].clone();
                aux.chain_status = ChainStatus::Rescued;
                kept_aux.push(aux);
                rescued_count += 1;

                // eprintln!("[SCAFFOLD_TRACE] Mapping[{}] distance={:.0} rescued",
                //           raw_idx, min_distance);
            }
        }

        // eprintln!("[SCAFFOLD_TRACE] Final for group: {} -> {} (anchors={}, rescued={})",
        //           raw_indices_for_pair.len(), kept_indices.len(),
        //           anchor_count, rescued_count);

        all_kept_indices.extend(kept_indices);
        all_kept_aux.extend(kept_aux);
    }

    // Build final result from RAW mappings
    let mut result_mappings = Vec::new();
    let mut result_aux = Vec::new();

    for (i, aux) in all_kept_aux.into_iter().enumerate() {
        let idx = all_kept_indices[i];
        result_mappings.push(raw_mappings[idx]);
        result_aux.push(aux);
    }

    // eprintln!("[SCAFFOLD_TRACE] Final: {} raw -> {} kept mappings",
    //           raw_mappings.len(), result_mappings.len());

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
    let mut current_matches = 0u64;
    let mut current_block_length = 0u64;

    for &(_, ref mapping, ref aux) in mappings {
        match current_chain {
            None => {
                // Start new chain
                current_matches = mapping.matches as u64;
                current_block_length = mapping.block_length as u64;
                let identity = if current_block_length > 0 {
                    current_matches as f64 / current_block_length as f64
                } else {
                    0.0
                };
                current_chain = Some(ScaffoldChain {
                    query_name: aux.query_name.clone(),
                    ref_name: aux.ref_name.clone(),
                    ref_seq_id: mapping.ref_seq_id,
                    query_start: mapping.query_start_pos,
                    query_end: mapping.query_end_pos(),
                    ref_start: mapping.ref_start_pos,
                    ref_end: mapping.ref_end_pos(),
                    total_length: mapping.block_length,
                    identity,
                    is_reverse: mapping.is_reverse(),
                });
            }
            Some(ref mut chain) => {
                // Check if can extend current chain
                // Calculate gaps (can be negative for overlaps)
                let ref_gap = mapping.ref_start_pos as i64 - chain.ref_end as i64;
                let query_gap = mapping.query_start_pos as i64 - chain.query_end as i64;

                // Allow small overlaps up to max_gap/5 (matching wfmash)
                let allowed_overlap = (max_gap / 5) as i64;

                if mapping.ref_seq_id == chain.ref_seq_id &&
                   mapping.is_reverse() == chain.is_reverse &&
                   ref_gap >= -allowed_overlap && ref_gap <= max_gap as i64 &&
                   query_gap >= -allowed_overlap && query_gap <= max_gap as i64 {
                    // Extend chain
                    current_matches += mapping.matches as u64;
                    current_block_length += mapping.block_length as u64;
                    chain.query_end = mapping.query_end_pos();
                    chain.ref_end = mapping.ref_end_pos();
                    chain.total_length = chain.ref_end - chain.ref_start;
                    chain.identity = if current_block_length > 0 {
                        current_matches as f64 / current_block_length as f64
                    } else {
                        0.0
                    };
                } else {
                    // Save current chain and start new one
                    chains.push(chain.clone());
                    current_matches = mapping.matches as u64;
                    current_block_length = mapping.block_length as u64;
                    let identity = if current_block_length > 0 {
                        current_matches as f64 / current_block_length as f64
                    } else {
                        0.0
                    };
                    *chain = ScaffoldChain {
                        query_name: aux.query_name.clone(),
                        ref_name: aux.ref_name.clone(),
                        ref_seq_id: mapping.ref_seq_id,
                        query_start: mapping.query_start_pos,
                        query_end: mapping.query_end_pos(),
                        ref_start: mapping.ref_start_pos,
                        ref_end: mapping.ref_end_pos(),
                        total_length: mapping.block_length,
                        identity,
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