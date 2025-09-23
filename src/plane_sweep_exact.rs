// Exact implementation of wfmash plane sweep algorithm
use std::collections::BTreeSet;
use std::cmp::Ordering;

/// Compact mapping for plane sweep (minimal fields needed)
#[derive(Debug, Clone, Copy)]
pub struct PlaneSweepMapping {
    pub idx: usize,           // Original index in input
    pub query_start: u32,
    pub query_end: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub identity: f64,        // 0.0 - 1.0
    pub flags: u8,            // bit 0: discard, bit 1: overlapped
}

impl PlaneSweepMapping {
    pub const FLAG_DISCARD: u8 = 0x01;
    pub const FLAG_OVERLAPPED: u8 = 0x02;

    pub fn score(&self) -> f64 {
        let length = (self.query_end - self.query_start) as f64;
        if length <= 0.0 || self.identity <= 0.0 {
            f64::NEG_INFINITY
        } else {
            self.identity * length.ln()
        }
    }

    pub fn is_discard(&self) -> bool {
        (self.flags & Self::FLAG_DISCARD) != 0
    }

    pub fn set_discard(&mut self, discard: bool) {
        if discard {
            self.flags |= Self::FLAG_DISCARD;
        } else {
            self.flags &= !Self::FLAG_DISCARD;
        }
    }

    pub fn is_overlapped(&self) -> bool {
        (self.flags & Self::FLAG_OVERLAPPED) != 0
    }

    pub fn set_overlapped(&mut self, overlapped: bool) {
        if overlapped {
            self.flags |= Self::FLAG_OVERLAPPED;
        } else {
            self.flags &= !Self::FLAG_OVERLAPPED;
        }
    }

    /// Calculate overlap fraction with another mapping (query axis)
    pub fn query_overlap(&self, other: &PlaneSweepMapping) -> f64 {
        let overlap_start = self.query_start.max(other.query_start);
        let overlap_end = self.query_end.min(other.query_end);
        let overlap_len = (overlap_end as i64 - overlap_start as i64).max(0) as f64;

        let self_len = (self.query_end - self.query_start) as f64;
        let other_len = (other.query_end - other.query_start) as f64;
        let min_len = self_len.min(other_len);

        if min_len > 0.0 {
            overlap_len / min_len
        } else {
            0.0
        }
    }

    /// Calculate overlap fraction with another mapping (target axis)
    pub fn target_overlap(&self, other: &PlaneSweepMapping) -> f64 {
        let overlap_start = self.target_start.max(other.target_start);
        let overlap_end = self.target_end.min(other.target_end);
        let overlap_len = (overlap_end as i64 - overlap_start as i64).max(0) as f64;

        let self_len = (self.target_end - self.target_start) as f64;
        let other_len = (other.target_end - other.target_start) as f64;
        let min_len = self_len.min(other_len);

        if min_len > 0.0 {
            overlap_len / min_len
        } else {
            0.0
        }
    }
}

/// Event for plane sweep
#[derive(Debug, Clone)]
struct Event {
    position: u32,
    event_type: EventType,
    mapping_idx: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum EventType {
    Begin = 0,
    End = 1,
}

/// BST order for mappings (matching wfmash Helper comparator)
#[derive(Clone, Debug)]
struct MappingOrder {
    idx: usize,
    score: f64,
    start_pos: u32,
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
        // Order by score (descending), then start position
        // This matches wfmash's Helper operator()
        other.score.partial_cmp(&self.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.start_pos.cmp(&other.start_pos))
            .then_with(|| self.idx.cmp(&other.idx))  // For determinism
    }
}

/// Mark good mappings (exact implementation of wfmash's markGood)
fn mark_good(
    bst: &BTreeSet<MappingOrder>,
    mappings: &mut [PlaneSweepMapping],
    secondaries_to_keep: usize,
    overlap_threshold: f64,
    axis: Axis,
) {
    if bst.is_empty() {
        return;
    }

    // Keep best secondaryToKeep + 1 mappings
    let mut kept = 0;
    let to_keep = secondaries_to_keep + 1;

    // First pass: mark best mappings as not discarded
    let mut kept_indices = Vec::new();
    let first_score = bst.iter().next().map(|m| m.score).unwrap_or(f64::NEG_INFINITY);

    for mapping_order in bst {
        // Check if we should stop keeping mappings
        // This matches the wfmash condition: (greater_score || !discard) && kept > secondaryToKeep
        if kept >= to_keep && mapping_order.score < first_score {
            break;
        }

        mappings[mapping_order.idx].set_discard(false);
        kept_indices.push(mapping_order.idx);
        kept += 1;

        if kept >= to_keep && !mappings[mapping_order.idx].is_discard() {
            break;
        }
    }

    // Second pass: check for overlaps if threshold < 1.0
    if overlap_threshold < 1.0 {
        let mut check_idx = 0;
        for mapping_order in bst {
            if check_idx < kept_indices.len() && mapping_order.idx == kept_indices[check_idx] {
                check_idx += 1;
                continue;  // Skip already kept mappings
            }

            let idx = mapping_order.idx;
            for &kept_idx in &kept_indices {
                let overlap = match axis {
                    Axis::Query => mappings[idx].query_overlap(&mappings[kept_idx]),
                    Axis::Target => mappings[idx].target_overlap(&mappings[kept_idx]),
                };

                if overlap > overlap_threshold {
                    mappings[idx].set_overlapped(true);
                    mappings[idx].set_discard(true);
                    break;
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    Query,
    Target,
}

/// Apply plane sweep on query axis (exact wfmash algorithm)
pub fn plane_sweep_query(
    mappings: &mut [PlaneSweepMapping],
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) -> Vec<usize> {
    if mappings.is_empty() || mappings.len() == 1 {
        return (0..mappings.len()).collect();
    }

    // Initially mark all as discarded
    for mapping in mappings.iter_mut() {
        mapping.set_discard(true);
        mapping.set_overlapped(false);
    }

    // Create event schedule
    let mut events = Vec::with_capacity(mappings.len() * 2);
    for (idx, mapping) in mappings.iter().enumerate() {
        events.push(Event {
            position: mapping.query_start,
            event_type: EventType::Begin,
            mapping_idx: idx,
        });
        events.push(Event {
            position: mapping.query_end,
            event_type: EventType::End,
            mapping_idx: idx,
        });
    }

    // Sort events by position, then type (BEGIN before END)
    events.sort_by_key(|e| (e.position, e.event_type));

    // Plane sweep with BST
    let mut bst = BTreeSet::new();
    let mut i = 0;

    while i < events.len() {
        let current_pos = events[i].position;

        // Find all events at current position
        let mut j = i;
        while j < events.len() && events[j].position == current_pos {
            j += 1;
        }

        // Update BST by processing all events at current position
        for k in i..j {
            let event = &events[k];
            let mapping_order = MappingOrder {
                idx: event.mapping_idx,
                score: mappings[event.mapping_idx].score(),
                start_pos: mappings[event.mapping_idx].query_start,
            };

            match event.event_type {
                EventType::Begin => {
                    bst.insert(mapping_order);
                }
                EventType::End => {
                    bst.remove(&mapping_order);
                }
            }
        }

        // Mark good mappings
        mark_good(&bst, mappings, secondaries_to_keep, overlap_threshold, Axis::Query);

        i = j;
    }

    // Return indices of kept mappings
    mappings.iter()
        .enumerate()
        .filter(|(_, m)| !m.is_discard() && !m.is_overlapped())
        .map(|(idx, _)| idx)
        .collect()
}

/// Apply plane sweep on target axis
pub fn plane_sweep_target(
    mappings: &mut [PlaneSweepMapping],
    secondaries_to_keep: usize,
    overlap_threshold: f64,
) -> Vec<usize> {
    if mappings.is_empty() || mappings.len() == 1 {
        return (0..mappings.len()).collect();
    }

    // Initially mark all as discarded
    for mapping in mappings.iter_mut() {
        mapping.set_discard(true);
        mapping.set_overlapped(false);
    }

    // Create event schedule
    let mut events = Vec::with_capacity(mappings.len() * 2);
    for (idx, mapping) in mappings.iter().enumerate() {
        events.push(Event {
            position: mapping.target_start,
            event_type: EventType::Begin,
            mapping_idx: idx,
        });
        events.push(Event {
            position: mapping.target_end,
            event_type: EventType::End,
            mapping_idx: idx,
        });
    }

    events.sort_by_key(|e| (e.position, e.event_type));

    let mut bst = BTreeSet::new();
    let mut i = 0;

    while i < events.len() {
        let current_pos = events[i].position;

        let mut j = i;
        while j < events.len() && events[j].position == current_pos {
            j += 1;
        }

        for k in i..j {
            let event = &events[k];
            let mapping_order = MappingOrder {
                idx: event.mapping_idx,
                score: mappings[event.mapping_idx].score(),
                start_pos: mappings[event.mapping_idx].target_start,
            };

            match event.event_type {
                EventType::Begin => {
                    bst.insert(mapping_order);
                }
                EventType::End => {
                    bst.remove(&mapping_order);
                }
            }
        }

        mark_good(&bst, mappings, secondaries_to_keep, overlap_threshold, Axis::Target);

        i = j;
    }

    mappings.iter()
        .enumerate()
        .filter(|(_, m)| !m.is_discard() && !m.is_overlapped())
        .map(|(idx, _)| idx)
        .collect()
}

/// Apply both query and target filtering
pub fn plane_sweep_both(
    mappings: &mut [PlaneSweepMapping],
    query_secondaries: usize,
    target_secondaries: usize,
    overlap_threshold: f64,
) -> Vec<usize> {
    // First apply query axis filtering
    let query_kept = plane_sweep_query(mappings, query_secondaries, overlap_threshold);

    // Create a filtered set for target sweep
    let mut filtered_mappings: Vec<PlaneSweepMapping> = query_kept.iter()
        .map(|&idx| mappings[idx])
        .collect();

    // Apply target axis filtering
    let target_kept = plane_sweep_target(&mut filtered_mappings, target_secondaries, overlap_threshold);

    // Map back to original indices
    target_kept.iter()
        .map(|&idx| query_kept[idx])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_input() {
        let mut mappings = vec![];
        let kept = plane_sweep_query(&mut mappings, 0, 0.95);
        assert_eq!(kept.len(), 0);
    }

    #[test]
    fn test_single_mapping() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 100,
                query_end: 200,
                target_start: 300,
                target_end: 400,
                identity: 0.95,
                flags: 0,
            }
        ];
        let kept = plane_sweep_query(&mut mappings, 0, 0.95);
        assert_eq!(kept, vec![0]);
    }

    #[test]
    fn test_non_overlapping_mappings() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 100,
                query_end: 200,
                target_start: 300,
                target_end: 400,
                identity: 0.95,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 1,
                query_start: 300,
                query_end: 400,
                target_start: 500,
                target_end: 600,
                identity: 0.90,
                flags: 0,
            },
        ];

        // Keep best 1 per position (0 secondaries)
        let kept = plane_sweep_query(&mut mappings, 0, 0.95);
        assert_eq!(kept.len(), 2);  // Both should be kept (non-overlapping)
    }

    #[test]
    fn test_overlapping_mappings() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 100,
                query_end: 200,
                target_start: 300,
                target_end: 400,
                identity: 0.95,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 1,
                query_start: 150,
                query_end: 250,
                target_start: 350,
                target_end: 450,
                identity: 0.90,
                flags: 0,
            },
        ];

        // Keep best 1 per position
        let kept = plane_sweep_query(&mut mappings, 0, 0.95);
        assert_eq!(kept, vec![0]);  // Should keep only the first (higher identity)
    }

    #[test]
    fn test_secondaries() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 100,
                query_end: 200,
                target_start: 300,
                target_end: 400,
                identity: 0.95,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 1,
                query_start: 100,
                query_end: 200,
                target_start: 500,
                target_end: 600,
                identity: 0.90,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 2,
                query_start: 100,
                query_end: 200,
                target_start: 700,
                target_end: 800,
                identity: 0.85,
                flags: 0,
            },
        ];

        // Keep best 2 (1 secondary)
        let kept = plane_sweep_query(&mut mappings, 1, 0.95);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&0));  // Best
        assert!(kept.contains(&1));  // Secondary
    }

    #[test]
    fn test_overlap_threshold() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 100,
                query_end: 200,
                target_start: 300,
                target_end: 400,
                identity: 0.95,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 1,
                query_start: 180,  // 20% overlap
                query_end: 280,
                target_start: 500,
                target_end: 600,
                identity: 0.90,
                flags: 0,
            },
        ];

        // With high overlap threshold (allow overlaps)
        let kept = plane_sweep_query(&mut mappings, 1, 1.0);
        assert_eq!(kept.len(), 2);  // Both kept

        // With low overlap threshold
        mappings[0].flags = 0;
        mappings[1].flags = 0;
        let kept = plane_sweep_query(&mut mappings, 1, 0.1);
        assert_eq!(kept.len(), 1);  // Second should be filtered
    }

    #[test]
    fn test_chromosome_boundaries() {
        let mut mappings = vec![
            PlaneSweepMapping {
                idx: 0,
                query_start: 0,
                query_end: 100,
                target_start: 0,
                target_end: 100,
                identity: 0.95,
                flags: 0,
            },
            PlaneSweepMapping {
                idx: 1,
                query_start: u32::MAX - 100,
                query_end: u32::MAX,
                target_start: 1000,
                target_end: 1100,
                identity: 0.90,
                flags: 0,
            },
        ];

        let kept = plane_sweep_query(&mut mappings, 0, 0.95);
        assert_eq!(kept.len(), 2);  // Both at extremes
    }
}