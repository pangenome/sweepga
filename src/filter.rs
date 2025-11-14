use anyhow::Result;
use crate::mapping::{Mapping, MappingAux};
use std::collections::HashMap;
use crate::filter_scaffold;

/// Filtering mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterMode {
    OneToOne,   // 1:1 - best mapping per query AND per target
    OneToMany,  // 1:N - best mapping per query, N per target (controlled by -n)
    ManyToMany, // N:N - N mappings per query and per target (controlled by -n)
}

/// Filter configuration matching wfmash parameters
#[derive(Clone)]
pub struct FilterConfig {
    pub chain_gap: u32,           // -c/--chain-jump
    pub min_block_length: u32,    // -l/--block-length
    pub filter_mode: FilterMode,  // Filtering mode
    pub max_per_query: Option<usize>,  // Max mappings per query
    pub max_per_target: Option<usize>, // Max mappings per target
    pub overlap_threshold: f64,   // -O/--overlap
    pub sparsity: f64,            // -x/--sparsify
    pub no_merge: bool,           // -M/--no-merge
    pub scaffold_gap: u32,        // -j/--scaffold-jump
    pub min_scaffold_length: u32, // -S/--scaffold-mass
    pub scaffold_overlap_threshold: f64, // --scaffold-overlap
    pub scaffold_max_deviation: u32, // -D/--scaffold-dist
    pub prefix_delimiter: char,   // -Y/--group-prefix (default '#')
    pub skip_prefix: bool,         // Whether to group by prefix
}

impl Default for FilterConfig {
    fn default() -> Self {
        FilterConfig {
            chain_gap: 2000,
            min_block_length: 0,
            filter_mode: FilterMode::ManyToMany,
            max_per_query: None,
            max_per_target: None,
            overlap_threshold: 0.95,
            sparsity: 1.0,
            no_merge: false,
            scaffold_gap: 100000,
            min_scaffold_length: 10000,  // Default to scaffold filtering
            scaffold_overlap_threshold: 0.5,
            scaffold_max_deviation: 100000,  // wfmash default
            prefix_delimiter: '#',  // Default PanSN delimiter
            skip_prefix: true,      // Group by prefix by default
        }
    }
}

/// Extract prefix from sequence name using delimiter
/// For "genome#haplotype#contig", returns "genome#haplotype#"
fn extract_prefix(seq_name: &str, delimiter: char) -> String {
    // Find the last occurrence of the delimiter
    if let Some(last_pos) = seq_name.rfind(delimiter) {
        // Include the delimiter in the prefix
        seq_name[..=last_pos].to_string()
    } else {
        // No delimiter found, use the whole name as prefix
        seq_name.to_string()
    }
}

/// Main filtering engine
pub struct FilterEngine {
    pub mappings: Vec<Mapping>,
    pub aux_data: Vec<MappingAux>,
    pub config: FilterConfig,
}

impl FilterEngine {
    pub fn new(mappings: Vec<Mapping>, aux_data: Vec<MappingAux>, config: FilterConfig) -> Self {
        FilterEngine {
            mappings,
            aux_data,
            config,
        }
    }

    /// Group mappings by query-target prefix pairs
    fn group_by_prefix_pairs(&self) -> HashMap<(String, String), Vec<usize>> {
        let mut groups: HashMap<(String, String), Vec<usize>> = HashMap::new();

        for (i, (mapping, aux)) in self.mappings.iter().zip(self.aux_data.iter()).enumerate() {
            let query_prefix = if self.config.skip_prefix {
                extract_prefix(&aux.query_name, self.config.prefix_delimiter)
            } else {
                aux.query_name.clone()
            };

            let target_prefix = if self.config.skip_prefix {
                extract_prefix(&aux.ref_name, self.config.prefix_delimiter)
            } else {
                aux.ref_name.clone()
            };

            groups.entry((query_prefix, target_prefix))
                .or_insert_with(Vec::new)
                .push(i);
        }

        groups
    }

    /// Apply full filtering pipeline in order (matching wfmash)
    /// Now groups by query-target prefix pairs and filters each group independently
    pub fn apply_filters(&mut self) -> Result<()> {
        // Group mappings by query-target prefix pairs
        let groups = self.group_by_prefix_pairs();

        // Collect results from each group
        let mut all_mappings = Vec::new();
        let mut all_aux = Vec::new();

        for ((query_prefix, target_prefix), indices) in groups {
            // Extract mappings for this prefix pair
            let mut group_mappings: Vec<Mapping> = indices.iter()
                .map(|&i| self.mappings[i].clone())
                .collect();
            let mut group_aux: Vec<MappingAux> = indices.iter()
                .map(|&i| self.aux_data[i].clone())
                .collect();

            // Create a temporary filter engine for this group
            let mut group_engine = FilterEngine {
                mappings: group_mappings,
                aux_data: group_aux,
                config: self.config.clone(),
            };

            // Apply filters to this group
            group_engine.apply_filters_to_group()?;

            // Collect filtered results
            all_mappings.extend(group_engine.mappings);
            all_aux.extend(group_engine.aux_data);
        }

        // Replace with filtered results
        self.mappings = all_mappings;
        self.aux_data = all_aux;

        Ok(())
    }

    /// Apply filters to a single prefix pair group
    fn apply_filters_to_group(&mut self) -> Result<()> {
        // Keep a copy of the original raw mappings for scaffold rescue
        // This is crucial - wfmash identifies scaffolds from filtered mappings
        // but rescues from the original raw mappings
        let raw_mappings = self.mappings.clone();
        let raw_aux = self.aux_data.clone();

        // 1. Merge mappings (unless disabled)
        if !self.config.no_merge {
            self.merge_mappings()?;
        }

        // 2. Filter weak mappings
        self.filter_weak_mappings()?;

        // 3. Apply filtering mode BEFORE scaffold identification
        // This determines which mappings become scaffold anchors
        match self.config.filter_mode {
            FilterMode::OneToOne => {
                self.one_to_one_filter()?;
            }
            FilterMode::OneToMany => {
                // Keep best per query, limit per target
                self.filter_one_to_many()?;
            }
            FilterMode::ManyToMany => {
                // Limit both dimensions if specified
                if self.config.max_per_query.is_some() || self.config.max_per_target.is_some() {
                    self.filter_many_to_many()?;
                }
            }
        }

        // 4. Apply sparsification if requested
        if self.config.sparsity < 1.0 {
            self.sparsify_mappings()?;
        }

        // 5. Apply scaffold filtering using BOTH filtered (for anchors) and raw (for rescue)
        if self.config.min_scaffold_length > 0 {
            let (filtered_mappings, filtered_aux) = filter_scaffold::filter_by_scaffolds_with_rescue(
                &self.mappings,  // Use filtered mappings to identify scaffolds
                &self.aux_data,
                &raw_mappings,   // But rescue from raw mappings
                &raw_aux,
                &self.config,
            )?;
            self.mappings = filtered_mappings;
            self.aux_data = filtered_aux;
        }

        Ok(())
    }

    /// Merge mappings within chain gap distance
    fn merge_mappings(&mut self) -> Result<()> {
        if self.mappings.is_empty() {
            return Ok(());
        }

        // Sort mappings by query, ref, and position
        let mut indices: Vec<usize> = (0..self.mappings.len()).collect();
        indices.sort_by_key(|&i| {
            let m = &self.mappings[i];
            let a = &self.aux_data[i];
            (
                a.query_seq_id,
                m.ref_seq_id,
                m.ref_start_pos,
                m.query_start_pos
            )
        });

        let mut merged = Vec::new();
        let mut merged_aux = Vec::new();

        let mut i = 0;
        while i < indices.len() {
            let idx = indices[i];
            let mut current = self.mappings[idx];
            let mut current_aux = self.aux_data[idx].clone();
            let mut merged_count = 1u32;

            // Try to merge with subsequent mappings
            let mut j = i + 1;
            while j < indices.len() {
                let next_idx = indices[j];
                let next = &self.mappings[next_idx];
                let next_aux = &self.aux_data[next_idx];

                // Check if same query and ref
                if current_aux.query_seq_id != next_aux.query_seq_id ||
                   current.ref_seq_id != next.ref_seq_id ||
                   current.is_reverse() != next.is_reverse() {
                    break;
                }

                // Check gap distance
                let ref_gap = next.ref_start_pos.saturating_sub(current.ref_end_pos());
                let query_gap = next.query_start_pos.saturating_sub(current.query_end_pos());

                if ref_gap <= self.config.chain_gap && query_gap <= self.config.chain_gap {
                    // Merge the mapping
                    let new_ref_end = next.ref_end_pos();
                    let new_query_end = next.query_end_pos();

                    current.block_length = new_ref_end - current.ref_start_pos;
                    current.conserved_sketches += next.conserved_sketches;

                    // Update identity as weighted average
                    let curr_weight = current.n_merged as f64;
                    let next_weight = next.n_merged as f64;
                    let total_weight = curr_weight + next_weight;
                    let new_identity = (current.identity() * curr_weight +
                                       next.identity() * next_weight) / total_weight;
                    current.set_identity(new_identity);

                    current.n_merged += next.n_merged;
                    merged_count += 1;
                    j += 1;
                } else {
                    break;
                }
            }

            merged.push(current);
            merged_aux.push(current_aux);
            i = j;
        }

        self.mappings = merged;
        self.aux_data = merged_aux;
        Ok(())
    }

    /// Filter out weak mappings below minimum block length
    fn filter_weak_mappings(&mut self) -> Result<()> {
        if self.config.min_block_length == 0 {
            return Ok(());
        }

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            if mapping.block_length >= self.config.min_block_length {
                keep.push(*mapping);
                keep_aux.push(self.aux_data[i].clone());
            }
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Keep only top N mappings per query-ref group
    fn filter_by_group(&mut self, max_per_group: usize) -> Result<()> {
        // Group mappings by query-ref pair
        let mut groups: HashMap<(i32, u32), Vec<usize>> = HashMap::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let key = (self.aux_data[i].query_seq_id, mapping.ref_seq_id);
            groups.entry(key).or_insert_with(Vec::new).push(i);
        }

        let mut keep_indices = Vec::new();

        for (_key, mut indices) in groups {
            // Sort by block length (descending)
            indices.sort_by_key(|&i| std::cmp::Reverse(self.mappings[i].block_length));

            // Keep top N
            for &idx in indices.iter().take(max_per_group) {
                keep_indices.push(idx);
            }
        }

        keep_indices.sort_unstable();

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for idx in keep_indices {
            keep.push(self.mappings[idx]);
            keep_aux.push(self.aux_data[idx].clone());
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Sparsify mappings by keeping only a fraction
    fn sparsify_mappings(&mut self) -> Result<()> {
        if self.config.sparsity >= 1.0 {
            return Ok(());
        }

        let target_count = (self.mappings.len() as f64 * self.config.sparsity) as usize;
        if target_count >= self.mappings.len() {
            return Ok(());
        }

        // Sort by block length and keep top fraction
        let mut indices: Vec<usize> = (0..self.mappings.len()).collect();
        indices.sort_by_key(|&i| std::cmp::Reverse(self.mappings[i].block_length));

        indices.truncate(target_count);
        indices.sort_unstable();

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for idx in indices {
            keep.push(self.mappings[idx]);
            keep_aux.push(self.aux_data[idx].clone());
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Filter by scaffolds - group and chain mappings at scaffold level
    fn filter_by_scaffolds(&mut self) -> Result<()> {
        if self.config.min_scaffold_length == 0 {
            return Ok(());
        }

        // INSTRUMENTATION: Log input mappings
        // eprintln!("[SCAFFOLD_TRACE] Input: {} mappings", self.mappings.len());

        // Group mappings by query-ref pair at scaffold level
        let mut scaffold_groups: HashMap<(i32, u32), Vec<usize>> = HashMap::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let key = (self.aux_data[i].query_seq_id, mapping.ref_seq_id);
            scaffold_groups.entry(key).or_insert_with(Vec::new).push(i);
        }

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();
        let mut chain_counter = 0u32;

        for ((query_id, ref_id), mut indices) in scaffold_groups {
            // Sort by position
            indices.sort_by_key(|&i| (self.mappings[i].ref_start_pos,
                                      self.mappings[i].query_start_pos));

            // Build scaffold chains
            let mut chains: Vec<Vec<usize>> = Vec::new();

            for idx in indices {
                let mapping = &self.mappings[idx];
                let mut added = false;

                // Try to add to existing chain
                for chain in &mut chains {
                    if let Some(&last_idx) = chain.last() {
                        let last = &self.mappings[last_idx];

                        let ref_gap = mapping.ref_start_pos.saturating_sub(last.ref_end_pos());
                        let query_gap = mapping.query_start_pos.saturating_sub(last.query_end_pos());

                        if ref_gap <= self.config.scaffold_gap &&
                           query_gap <= self.config.scaffold_gap &&
                           mapping.is_reverse() == last.is_reverse() {
                            chain.push(idx);
                            added = true;
                            break;
                        }
                    }
                }

                if !added {
                    chains.push(vec![idx]);
                }
            }

            // INSTRUMENTATION: Log chains formed
            // eprintln!("[SCAFFOLD_TRACE] After merging with gap={}: {} chains",
            //           self.config.scaffold_gap, chains.len());
            for (i, chain) in chains.iter().enumerate() {
                if chain.is_empty() { continue; }
                let first = &self.mappings[chain[0]];
                let last = &self.mappings[*chain.last().unwrap()];
                let chain_len = if chain.len() > 1 {
                    last.ref_end_pos() - first.ref_start_pos
                } else {
                    first.block_length
                };
                // Copy packed fields to avoid alignment issues
                let ref_id = first.ref_seq_id;
                let q_start = first.query_start_pos;
                let q_end = last.query_end_pos();
                let r_start = first.ref_start_pos;
                let r_end = last.ref_end_pos();
                // eprintln!("[SCAFFOLD_TRACE]   Chain[{}]: ref={} q={}-{} r={}-{} len={} mappings={}",
                //           i, ref_id, q_start, q_end, r_start, r_end,
                //           chain_len, chain.len());
            }

            // Identify scaffold chains and rescued mappings
            let mut scaffold_chains = Vec::new();
            let mut rescued_mappings = Vec::new();

            let before_length_filter = chains.len();
            for (chain_idx, chain) in chains.iter().enumerate() {
                if chain.len() > 1 {
                    let first = &self.mappings[chain[0]];
                    let last = &self.mappings[*chain.last().unwrap()];
                    let scaffold_len = last.ref_end_pos() - first.ref_start_pos;

                    if scaffold_len >= self.config.min_scaffold_length {
                        scaffold_chains.push(chain_idx);
                    } else {
                        rescued_mappings.extend(chain.iter().copied());
                    }
                } else if self.mappings[chain[0]].block_length >= self.config.min_scaffold_length {
                    // Single mapping that meets threshold
                    scaffold_chains.push(chain_idx);
                } else {
                    rescued_mappings.push(chain[0]);
                }
            }

            // INSTRUMENTATION: Log scaffold identification
            // eprintln!("[SCAFFOLD_TRACE] After length filter (min={}): {} -> {} scaffolds",
            //           self.config.min_scaffold_length, before_length_filter, scaffold_chains.len());
            for &chain_idx in &scaffold_chains {
                let chain = &chains[chain_idx];
                let first = &self.mappings[chain[0]];
                let last = &self.mappings[*chain.last().unwrap()];
                let scaffold_len = if chain.len() > 1 {
                    last.ref_end_pos() - first.ref_start_pos
                } else {
                    first.block_length
                };
                // Copy packed fields to avoid alignment issues
                let ref_id = first.ref_seq_id;
                let q_start = first.query_start_pos;
                let q_end = last.query_end_pos();
                let r_start = first.ref_start_pos;
                let r_end = last.ref_end_pos();
                // eprintln!("[SCAFFOLD_TRACE]   Scaffold[{}]: ref={} q={}-{} r={}-{} len={}",
                //           chain_idx, ref_id, q_start, q_end, r_start, r_end,
                //           scaffold_len);
            }

            // Process scaffold chains
            let mut anchor_count = 0;
            for &chain_idx in &scaffold_chains {
                chain_counter += 1;
                let chain_id = format!("{}.{}.{}", query_id + 1, ref_id + 1, chain_counter);

                for &idx in &chains[chain_idx] {
                    let mut aux = self.aux_data[idx].clone();
                    aux.chain_id = chain_id.clone();
                    aux.chain_status = crate::mapping::ChainStatus::Scaffold;
                    keep.push(self.mappings[idx]);
                    keep_aux.push(aux);
                    anchor_count += 1;
                }
            }
            // eprintln!("[SCAFFOLD_TRACE] Anchors identified: {} mappings", anchor_count);

            // Process rescued mappings (calculate distance to nearest scaffold)
            let mut rescued_count = 0;
            for idx in rescued_mappings {
                chain_counter += 1;
                let chain_id = format!("{}.{}.{}", query_id + 1, ref_id + 1, chain_counter);

                let mut aux = self.aux_data[idx].clone();
                aux.chain_id = chain_id;

                // Calculate distance to nearest scaffold chain
                let mapping = &self.mappings[idx];
                let mut min_distance = u32::MAX;
                let mut nearest_scaffold = None;

                for &chain_idx in &scaffold_chains {
                    for &scaffold_idx in &chains[chain_idx] {
                        let scaffold_mapping = &self.mappings[scaffold_idx];
                        let dist = if mapping.ref_start_pos > scaffold_mapping.ref_end_pos() {
                            mapping.ref_start_pos - scaffold_mapping.ref_end_pos()
                        } else if scaffold_mapping.ref_start_pos > mapping.ref_end_pos() {
                            scaffold_mapping.ref_start_pos - mapping.ref_end_pos()
                        } else {
                            0  // Overlapping
                        };
                        if dist < min_distance {
                            min_distance = dist;
                            nearest_scaffold = Some(scaffold_idx);
                        }
                    }
                }

                // INSTRUMENTATION: Log rescue decisions
                if let Some(scaffold_idx) = nearest_scaffold {
                    let scaffold = &self.mappings[scaffold_idx];
                    // Copy packed fields to avoid alignment issues
                    let s_q_start = scaffold.query_start_pos;
                    let s_q_end = scaffold.query_end_pos();
                    let s_r_start = scaffold.ref_start_pos;
                    let s_r_end = scaffold.ref_end_pos();
                    let m_q_start = mapping.query_start_pos;
                    let m_q_end = mapping.query_end_pos();
                    let m_r_start = mapping.ref_start_pos;
                    let m_r_end = mapping.ref_end_pos();
                    // eprintln!("[SCAFFOLD_TRACE] Mapping[{}] distance={} to anchor q={}-{} r={}-{} (mapping: q={}-{} r={}-{})",
                    //           idx, min_distance,
                    //           s_q_start, s_q_end, s_r_start, s_r_end,
                    //           m_q_start, m_q_end, m_r_start, m_r_end);
                }

                aux.chain_status = crate::mapping::ChainStatus::Rescued;
                keep.push(self.mappings[idx]);
                keep_aux.push(aux);
                rescued_count += 1;
            }
            // eprintln!("[SCAFFOLD_TRACE] Rescued: {} mappings", rescued_count);
        }

        // INSTRUMENTATION: Final summary
        let anchors = keep_aux.iter().filter(|a| matches!(a.chain_status, crate::mapping::ChainStatus::Scaffold)).count();
        let rescued = keep_aux.iter().filter(|a| matches!(a.chain_status, crate::mapping::ChainStatus::Rescued)).count();
        // eprintln!("[SCAFFOLD_TRACE] Final: {} -> {} (anchors={}, rescued={})",
        //           self.mappings.len(), keep.len(), anchors, rescued);

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Apply one-to-one filtering (reciprocal best hits)
    fn one_to_one_filter(&mut self) -> Result<()> {
        // Find best mapping for each query
        let mut best_for_query: HashMap<i32, (usize, u32)> = HashMap::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let query_id = self.aux_data[i].query_seq_id;
            let score = mapping.block_length;

            match best_for_query.get(&query_id) {
                None => {
                    best_for_query.insert(query_id, (i, score));
                }
                Some(&(_, best_score)) if score > best_score => {
                    best_for_query.insert(query_id, (i, score));
                }
                _ => {}
            }
        }

        // Find best mapping for each reference
        let mut best_for_ref: HashMap<u32, (usize, u32)> = HashMap::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let ref_id = mapping.ref_seq_id;
            let score = mapping.block_length;

            match best_for_ref.get(&ref_id) {
                None => {
                    best_for_ref.insert(ref_id, (i, score));
                }
                Some(&(_, best_score)) if score > best_score => {
                    best_for_ref.insert(ref_id, (i, score));
                }
                _ => {}
            }
        }

        // Keep only reciprocal best hits
        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let query_id = self.aux_data[i].query_seq_id;
            let ref_id = mapping.ref_seq_id;

            if let (Some(&(best_q_idx, _)), Some(&(best_r_idx, _))) =
                (best_for_query.get(&query_id), best_for_ref.get(&ref_id)) {
                if best_q_idx == i && best_r_idx == i {
                    keep.push(*mapping);
                    keep_aux.push(self.aux_data[i].clone());
                }
            }
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Filter for 1:N mode - keep best per query, limit per target
    fn filter_one_to_many(&mut self) -> Result<()> {
        // First, keep only best mapping per query
        let mut best_per_query: HashMap<i32, (usize, u32)> = HashMap::new();

        for (i, mapping) in self.mappings.iter().enumerate() {
            let query_id = self.aux_data[i].query_seq_id;
            let score = mapping.block_length;

            match best_per_query.get(&query_id) {
                None => {
                    best_per_query.insert(query_id, (i, score));
                }
                Some(&(_, best_score)) if score > best_score => {
                    best_per_query.insert(query_id, (i, score));
                }
                _ => {}
            }
        }

        // Collect indices of best mappings per query
        let mut keep_indices: Vec<usize> = best_per_query.values().map(|&(i, _)| i).collect();

        // Now apply target limit if specified
        if let Some(max_per_target) = self.config.max_per_target {
            let mut target_counts: HashMap<u32, Vec<(usize, u32)>> = HashMap::new();

            for &idx in &keep_indices {
                let mapping = &self.mappings[idx];
                target_counts.entry(mapping.ref_seq_id)
                    .or_insert_with(Vec::new)
                    .push((idx, mapping.block_length));
            }

            keep_indices.clear();
            for (_, mut indices) in target_counts {
                indices.sort_by_key(|&(_, score)| std::cmp::Reverse(score));
                for &(idx, _) in indices.iter().take(max_per_target) {
                    keep_indices.push(idx);
                }
            }
        }

        keep_indices.sort_unstable();

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for idx in keep_indices {
            keep.push(self.mappings[idx]);
            keep_aux.push(self.aux_data[idx].clone());
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Filter for N:N mode - limit both query and target dimensions
    fn filter_many_to_many(&mut self) -> Result<()> {
        let mut keep_indices = Vec::new();

        // Apply per-query limit
        if let Some(max_per_query) = self.config.max_per_query {
            let mut query_groups: HashMap<i32, Vec<(usize, u32)>> = HashMap::new();

            for (i, mapping) in self.mappings.iter().enumerate() {
                let query_id = self.aux_data[i].query_seq_id;
                query_groups.entry(query_id)
                    .or_insert_with(Vec::new)
                    .push((i, mapping.block_length));
            }

            for (_, mut indices) in query_groups {
                indices.sort_by_key(|&(_, score)| std::cmp::Reverse(score));
                for &(idx, _) in indices.iter().take(max_per_query) {
                    keep_indices.push(idx);
                }
            }
        } else {
            keep_indices = (0..self.mappings.len()).collect();
        }

        // Apply per-target limit
        if let Some(max_per_target) = self.config.max_per_target {
            let mut target_groups: HashMap<u32, Vec<(usize, u32)>> = HashMap::new();

            for &idx in &keep_indices {
                let mapping = &self.mappings[idx];
                target_groups.entry(mapping.ref_seq_id)
                    .or_insert_with(Vec::new)
                    .push((idx, mapping.block_length));
            }

            keep_indices.clear();
            for (_, mut indices) in target_groups {
                indices.sort_by_key(|&(_, score)| std::cmp::Reverse(score));
                for &(idx, _) in indices.iter().take(max_per_target) {
                    keep_indices.push(idx);
                }
            }
        }

        keep_indices.sort_unstable();
        keep_indices.dedup();

        let mut keep = Vec::new();
        let mut keep_aux = Vec::new();

        for idx in keep_indices {
            keep.push(self.mappings[idx]);
            keep_aux.push(self.aux_data[idx].clone());
        }

        self.mappings = keep;
        self.aux_data = keep_aux;
        Ok(())
    }

    /// Get filtered mappings
    pub fn get_mappings(&self) -> &[Mapping] {
        &self.mappings
    }

    /// Get auxiliary data
    pub fn get_aux_data(&self) -> &[MappingAux] {
        &self.aux_data
    }
}