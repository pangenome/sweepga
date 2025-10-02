#![allow(dead_code)]

use sweepga::mapping::{Mapping, MappingAux};
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};

/// Event types for plane sweep
#[derive(Debug, Clone, Copy, PartialEq)]
enum EventType {
    Begin = 1,
    End = 2,
}

/// Event for plane sweep algorithm
#[derive(Debug, Clone)]
struct Event {
    position: u32,
    event_type: EventType,
    mapping_idx: usize,
}

impl Event {
    fn new(position: u32, event_type: EventType, mapping_idx: usize) -> Self {
        Event {
            position,
            event_type,
            mapping_idx,
        }
    }
}

/// Calculate score for a mapping: identity * log(block_length)
fn calculate_score(mapping: &Mapping) -> f64 {
    let identity = mapping.identity() / 100.0; // Convert from percent
    let block_len = mapping.block_length as f64;

    if block_len <= 0.0 || identity <= 0.0 {
        f64::NEG_INFINITY
    } else {
        identity * block_len.ln()
    }
}

/// Wrapper for mapping indices in BST
#[derive(Clone, Debug)]
struct MappingOrder {
    idx: usize,
    score: f64,
    start_pos: u32,
    ref_id: u32,
}

impl PartialEq for MappingOrder {
    fn eq(&self, other: &Self) -> bool {
        self.idx == other.idx
    }
}

impl Eq for MappingOrder {}

impl PartialOrd for MappingOrder {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MappingOrder {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order by score (descending), then start position, then ref_id
        other
            .score
            .partial_cmp(&self.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.start_pos.cmp(&other.start_pos))
            .then_with(|| self.ref_id.cmp(&other.ref_id))
    }
}

/// Plane sweep filter for query axis
pub fn filter_by_query(
    mappings: &mut Vec<Mapping>,
    aux_data: &mut Vec<MappingAux>,
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) {
    if mappings.is_empty() {
        return;
    }

    // Initially mark all mappings as discarded
    for mapping in mappings.iter_mut() {
        mapping.set_discard(true);
        mapping.set_overlapped(false);
    }

    // Create events for plane sweep
    let mut events = Vec::with_capacity(mappings.len() * 2);
    for (idx, mapping) in mappings.iter().enumerate() {
        events.push(Event::new(mapping.query_start_pos, EventType::Begin, idx));
        events.push(Event::new(mapping.query_end_pos(), EventType::End, idx));
    }

    // Sort events by position
    events.sort_by_key(|e| (e.position, e.event_type as u32));

    // Pre-calculate scores for all mappings
    let scores: Vec<f64> = mappings.iter().map(calculate_score).collect();

    // Plane sweep with BST for active mappings
    let mut active_mappings = BTreeSet::new();

    let mut i = 0;
    while i < events.len() {
        let current_pos = events[i].position;

        // Process all events at current position
        let mut j = i;
        while j < events.len() && events[j].position == current_pos {
            let event = &events[j];
            let mapping_order = MappingOrder {
                idx: event.mapping_idx,
                score: scores[event.mapping_idx],
                start_pos: mappings[event.mapping_idx].query_start_pos,
                ref_id: mappings[event.mapping_idx].ref_seq_id,
            };

            match event.event_type {
                EventType::Begin => {
                    active_mappings.insert(mapping_order);
                }
                EventType::End => {
                    active_mappings.remove(&mapping_order);
                }
            }
            j += 1;
        }

        // Mark best mappings as good
        mark_good_mappings(
            &mut active_mappings,
            mappings,
            secondaries_to_keep,
            overlap_threshold,
        );

        i = j;
    }

    // Remove discarded mappings
    let mut keep_indices = Vec::new();
    for (idx, mapping) in mappings.iter().enumerate() {
        if !mapping.is_discard() && !mapping.is_overlapped() {
            keep_indices.push(idx);
        }
    }

    // Filter both mappings and aux_data
    let new_mappings: Vec<_> = keep_indices.iter().map(|&idx| mappings[idx]).collect();
    let new_aux: Vec<_> = keep_indices
        .iter()
        .map(|&idx| aux_data[idx].clone())
        .collect();

    *mappings = new_mappings;
    *aux_data = new_aux;
}

/// Mark the best N mappings in the active set as good
fn mark_good_mappings(
    active: &mut BTreeSet<MappingOrder>,
    mappings: &mut [Mapping],
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) {
    let mut kept = 0;
    let to_keep = secondaries_to_keep + 1; // +1 for primary

    // First pass: mark best mappings
    let mut kept_indices = Vec::new();
    for mapping_order in active.iter() {
        if kept >= to_keep {
            break;
        }
        mappings[mapping_order.idx].set_discard(false);
        kept_indices.push(mapping_order.idx);
        kept += 1;
    }

    // Second pass: check for overlaps if threshold < 1.0
    if overlap_threshold < 1.0 {
        for mapping_order in active.iter().skip(kept) {
            let idx = mapping_order.idx;
            for &kept_idx in &kept_indices {
                let overlap = calculate_overlap(&mappings[idx], &mappings[kept_idx]);
                if overlap > overlap_threshold {
                    mappings[idx].set_overlapped(true);
                    mappings[idx].set_discard(true);
                    break;
                }
            }
        }
    }
}

/// Calculate overlap fraction between two mappings on query axis
fn calculate_overlap(a: &Mapping, b: &Mapping) -> f64 {
    let overlap_start = a.query_start_pos.max(b.query_start_pos);
    let overlap_end = a.query_end_pos().min(b.query_end_pos());
    let overlap_len = (overlap_end as i64 - overlap_start as i64).max(0) as u32;

    let a_len = a.query_end_pos() - a.query_start_pos;
    let b_len = b.query_end_pos() - b.query_start_pos;
    let min_len = a_len.min(b_len);

    if min_len > 0 {
        overlap_len as f64 / min_len as f64
    } else {
        0.0
    }
}

/// Plane sweep filter for reference/target axis
pub fn filter_by_reference(
    mappings: &mut Vec<Mapping>,
    aux_data: &mut Vec<MappingAux>,
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) {
    if mappings.is_empty() {
        return;
    }

    // Group mappings by reference sequence
    let mut ref_groups: Vec<Vec<usize>> = Vec::new();
    let mut ref_id_to_group: Vec<Option<usize>> = vec![None; mappings.len()];

    for (idx, mapping) in mappings.iter().enumerate() {
        let ref_id = mapping.ref_seq_id as usize;
        if ref_id >= ref_id_to_group.len() {
            ref_id_to_group.resize(ref_id + 1, None);
        }

        if let Some(group_idx) = ref_id_to_group[ref_id] {
            ref_groups[group_idx].push(idx);
        } else {
            let group_idx = ref_groups.len();
            ref_groups.push(vec![idx]);
            ref_id_to_group[ref_id] = Some(group_idx);
        }
    }

    // Mark all as discarded initially
    for mapping in mappings.iter_mut() {
        mapping.set_discard(true);
        mapping.set_overlapped(false);
    }

    // Process each reference sequence group
    for group_indices in ref_groups {
        if group_indices.is_empty() {
            continue;
        }

        // Create events for this reference
        let mut events = Vec::with_capacity(group_indices.len() * 2);
        for &idx in &group_indices {
            let mapping = &mappings[idx];
            events.push(Event::new(mapping.ref_start_pos, EventType::Begin, idx));
            events.push(Event::new(mapping.ref_end_pos(), EventType::End, idx));
        }

        // Sort events
        events.sort_by_key(|e| (e.position, e.event_type as u32));

        // Pre-calculate scores for this reference group
        let scores: Vec<f64> = group_indices
            .iter()
            .map(|&idx| calculate_score(&mappings[idx]))
            .collect();

        // Plane sweep
        let mut active_mappings = BTreeSet::new();

        let mut i = 0;
        while i < events.len() {
            let current_pos = events[i].position;

            // Process all events at current position
            let mut j = i;
            while j < events.len() && events[j].position == current_pos {
                let event = &events[j];
                // Find the score index for this mapping
                let score_idx = group_indices
                    .iter()
                    .position(|&idx| idx == event.mapping_idx)
                    .unwrap_or(0);
                let mapping_order = MappingOrder {
                    idx: event.mapping_idx,
                    score: scores[score_idx],
                    start_pos: mappings[event.mapping_idx].ref_start_pos,
                    ref_id: mappings[event.mapping_idx].ref_seq_id,
                };

                match event.event_type {
                    EventType::Begin => {
                        active_mappings.insert(mapping_order);
                    }
                    EventType::End => {
                        active_mappings.remove(&mapping_order);
                    }
                }
                j += 1;
            }

            // Mark best mappings on reference axis
            mark_good_mappings_ref(
                &mut active_mappings,
                mappings,
                secondaries_to_keep,
                overlap_threshold,
            );

            i = j;
        }
    }

    // Remove discarded mappings
    let mut keep_indices = Vec::new();
    for (idx, mapping) in mappings.iter().enumerate() {
        if !mapping.is_discard() && !mapping.is_overlapped() {
            keep_indices.push(idx);
        }
    }

    // Filter both mappings and aux_data
    let new_mappings: Vec<_> = keep_indices.iter().map(|&idx| mappings[idx]).collect();
    let new_aux: Vec<_> = keep_indices
        .iter()
        .map(|&idx| aux_data[idx].clone())
        .collect();

    *mappings = new_mappings;
    *aux_data = new_aux;
}

/// Mark the best N mappings in the active set as good (reference axis)
fn mark_good_mappings_ref(
    active: &mut BTreeSet<MappingOrder>,
    mappings: &mut [Mapping],
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) {
    let mut kept = 0;
    let to_keep = secondaries_to_keep + 1;

    // First pass: mark best mappings
    let mut kept_indices = Vec::new();
    for mapping_order in active.iter() {
        if kept >= to_keep {
            break;
        }
        mappings[mapping_order.idx].set_discard(false);
        kept_indices.push(mapping_order.idx);
        kept += 1;
    }

    // Second pass: check for overlaps on reference axis
    if overlap_threshold < 1.0 {
        for mapping_order in active.iter().skip(kept) {
            let idx = mapping_order.idx;
            for &kept_idx in &kept_indices {
                let overlap = calculate_overlap_ref(&mappings[idx], &mappings[kept_idx]);
                if overlap > overlap_threshold {
                    mappings[idx].set_overlapped(true);
                    mappings[idx].set_discard(true);
                    break;
                }
            }
        }
    }
}

/// Calculate overlap fraction between two mappings on reference axis
fn calculate_overlap_ref(a: &Mapping, b: &Mapping) -> f64 {
    let overlap_start = a.ref_start_pos.max(b.ref_start_pos);
    let overlap_end = a.ref_end_pos().min(b.ref_end_pos());
    let overlap_len = (overlap_end as i64 - overlap_start as i64).max(0) as u32;

    let a_len = a.ref_end_pos() - a.ref_start_pos;
    let b_len = b.ref_end_pos() - b.ref_start_pos;
    let min_len = a_len.min(b_len);

    if min_len > 0 {
        overlap_len as f64 / min_len as f64
    } else {
        0.0
    }
}

/// Extract prefix from sequence name using delimiter
/// For "genome#haplotype#contig", returns "genome#haplotype#"
fn extract_prefix(seq_name: &str, delimiter: char, skip_prefix: bool) -> String {
    if !skip_prefix {
        return seq_name.to_string();
    }

    // Find the last occurrence of the delimiter
    if let Some(last_pos) = seq_name.rfind(delimiter) {
        // Include the delimiter in the prefix
        seq_name[..=last_pos].to_string()
    } else {
        // No delimiter found, use the whole name as prefix
        seq_name.to_string()
    }
}

/// Apply plane sweep filtering based on filter specification
///
/// Filter specification format:
/// - "N" or "N:N" - No filtering
/// - "1:1" - Keep best mapping per query AND per target (within each genome pair)
/// - "1" or "1:∞" - Keep best mapping per query, all per target
/// - "∞:1" - Keep all per query, best per target
/// - "M:N" - Keep best M per query, best N per target
pub fn apply_plane_sweep(
    mappings: &mut Vec<Mapping>,
    aux_data: &mut Vec<MappingAux>,
    filter_spec: &str,
    overlap_threshold: f64,
    prefix_delimiter: char,
    skip_prefix: bool,
) {
    // Parse filter specification
    let (query_limit, target_limit) = parse_filter_spec(filter_spec);

    // If no filtering needed, return early
    if query_limit.is_none() && target_limit.is_none() {
        return;
    }

    // Group mappings by prefix pairs (genome pairs)
    let groups = group_by_prefix_pairs(mappings, aux_data, prefix_delimiter, skip_prefix);

    // Process each group independently
    let mut all_keep_indices = Vec::new();

    for group_indices in groups.values() {
        if group_indices.is_empty() {
            continue;
        }

        // Create temporary vectors for this group
        let mut group_mappings: Vec<Mapping> =
            group_indices.iter().map(|&idx| mappings[idx]).collect();
        let mut group_aux: Vec<MappingAux> = group_indices
            .iter()
            .map(|&idx| aux_data[idx].clone())
            .collect();

        // Apply query axis filtering if needed
        if let Some(limit) = query_limit {
            if limit > 0 {
                let secondaries = limit.saturating_sub(1);
                filter_by_query(
                    &mut group_mappings,
                    &mut group_aux,
                    secondaries,
                    overlap_threshold,
                );
            }
        }

        // Apply target axis filtering if needed
        if let Some(limit) = target_limit {
            if limit > 0 {
                let secondaries = limit.saturating_sub(1);
                filter_by_reference(
                    &mut group_mappings,
                    &mut group_aux,
                    secondaries,
                    overlap_threshold,
                );
            }
        }

        // Map back to original indices
        for (i, &orig_idx) in group_indices.iter().enumerate() {
            if i < group_mappings.len() {
                // This mapping survived filtering
                all_keep_indices.push(orig_idx);
            }
        }
    }

    // Sort indices to maintain order
    all_keep_indices.sort();

    // Filter both mappings and aux_data
    let new_mappings: Vec<_> = all_keep_indices.iter().map(|&idx| mappings[idx]).collect();
    let new_aux: Vec<_> = all_keep_indices
        .iter()
        .map(|&idx| aux_data[idx].clone())
        .collect();

    *mappings = new_mappings;
    *aux_data = new_aux;
}

/// Group mappings by query-target chromosome pairs
/// Groups by full sequence IDs (one per chromosome)
fn group_by_prefix_pairs(
    mappings: &[Mapping],
    aux_data: &[MappingAux],
    _prefix_delimiter: char,
    _skip_prefix: bool,
) -> HashMap<(String, String), Vec<usize>> {
    let mut groups: HashMap<(String, String), Vec<usize>> = HashMap::new();

    for (idx, (mapping, aux)) in mappings.iter().zip(aux_data.iter()).enumerate() {
        // Group by full chromosome IDs (query_seq_id, ref_seq_id)
        // Each chromosome gets its own group
        // This matches the behavior in paf_filter.rs where we group by
        // (query_name, target_name) which are full chromosome names
        let query_key = aux.query_seq_id.to_string();
        // Copy field first to avoid packed field reference
        let ref_id = mapping.ref_seq_id;
        let target_key = ref_id.to_string();

        groups
            .entry((query_key, target_key))
            .or_default()
            .push(idx);
    }

    groups
}

/// Parse filter specification string
fn parse_filter_spec(spec: &str) -> (Option<usize>, Option<usize>) {
    let spec = spec.trim();

    // Handle "N" or "N:N" - no filtering
    if spec == "N" || spec == "N:N" {
        return (None, None);
    }

    // Handle single number (query only)
    if !spec.contains(':') {
        if let Ok(n) = spec.parse::<usize>() {
            return (Some(n), None);
        } else if spec == "∞" {
            return (None, None);
        }
    }

    // Handle M:N format
    if let Some((query_part, target_part)) = spec.split_once(':') {
        let query_limit = if query_part == "N" || query_part == "∞" {
            None
        } else {
            query_part.parse::<usize>().ok()
        };

        let target_limit = if target_part == "N" || target_part == "∞" {
            None
        } else {
            target_part.parse::<usize>().ok()
        };

        return (query_limit, target_limit);
    }

    // Default: no filtering
    (None, None)
}
