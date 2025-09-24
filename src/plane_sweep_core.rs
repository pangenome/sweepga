use std::cmp::Ordering;
/// Core plane sweep algorithm that works on intervals
/// This is the single implementation used for both query and target axes
use std::collections::BTreeSet;

/// An interval with scoring and metadata
pub struct Interval {
    pub idx: usize, // Original index in input array
    pub begin: u32, // Start position
    pub end: u32,   // End position
    pub score: f64, // Score for this interval
    pub flags: u32, // Status flags
}

impl Interval {
    pub fn length(&self) -> u32 {
        self.end - self.begin
    }

    pub fn overlaps(&self, other: &Interval, threshold: f64) -> bool {
        let overlap_start = self.begin.max(other.begin);
        let overlap_end = self.end.min(other.end);

        if overlap_start >= overlap_end {
            return false;
        }

        let overlap_len = overlap_end - overlap_start;
        let min_len = self.length().min(other.length());

        overlap_len as f64 / min_len as f64 > threshold
    }
}

/// Event type for sweep line
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EventType {
    Begin,
    End,
}

/// Event in the sweep line algorithm
struct Event {
    position: u32,
    event_type: EventType,
    interval_idx: usize,
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        // Sort by position, then by event type (Begin before End at same position)
        self.position
            .cmp(&other.position)
            .then_with(|| match (self.event_type, other.event_type) {
                (EventType::Begin, EventType::End) => Ordering::Less,
                (EventType::End, EventType::Begin) => Ordering::Greater,
                _ => Ordering::Equal,
            })
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Event {}
impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.event_type == other.event_type
    }
}

/// Core plane sweep algorithm
///
/// Takes intervals and returns indices of intervals to keep
/// This is the single source of truth for plane sweep filtering
pub fn plane_sweep(
    intervals: &mut [Interval],
    max_to_keep: usize,
    overlap_threshold: f64,
) -> Vec<usize> {
    if intervals.is_empty() {
        return vec![];
    }

    if intervals.len() == 1 {
        return vec![0];
    }

    // If keeping all, just return all indices
    if max_to_keep == usize::MAX {
        return (0..intervals.len()).collect();
    }

    // Create events
    let mut events = Vec::with_capacity(intervals.len() * 2);
    for (idx, interval) in intervals.iter().enumerate() {
        events.push(Event {
            position: interval.begin,
            event_type: EventType::Begin,
            interval_idx: idx,
        });
        events.push(Event {
            position: interval.end,
            event_type: EventType::End,
            interval_idx: idx,
        });
    }

    // Sort events
    events.sort_unstable();

    // Active set ordered by score (best first)
    let mut active: BTreeSet<(i64, usize)> = BTreeSet::new(); // (-score_bits, idx) for descending order
    let mut kept_indices = Vec::new();

    // Process events
    for event in events {
        match event.event_type {
            EventType::Begin => {
                let interval = &intervals[event.interval_idx];
                let score_bits = interval.score.to_bits() as i64;
                active.insert((-score_bits, event.interval_idx));

                // Mark the best intervals at this position
                mark_best(&active, intervals, max_to_keep, &mut kept_indices);
            }
            EventType::End => {
                let interval = &intervals[event.interval_idx];
                let score_bits = interval.score.to_bits() as i64;
                active.remove(&(-score_bits, event.interval_idx));
            }
        }
    }

    // Deduplicate kept indices
    kept_indices.sort_unstable();
    kept_indices.dedup();

    // Apply overlap filtering if threshold < 1.0
    if overlap_threshold < 1.0 {
        filter_by_overlap(&mut kept_indices, intervals, overlap_threshold);
    }

    kept_indices
}

/// Mark the best intervals at current position
fn mark_best(
    active: &BTreeSet<(i64, usize)>,
    _intervals: &[Interval],
    max_to_keep: usize,
    kept_indices: &mut Vec<usize>,
) {
    for (count, &(_, idx)) in active.iter().enumerate() {
        if count >= max_to_keep {
            break;
        }
        kept_indices.push(idx);
    }
}

/// Filter by overlap threshold
fn filter_by_overlap(
    kept_indices: &mut Vec<usize>,
    intervals: &[Interval],
    overlap_threshold: f64,
) {
    if kept_indices.len() <= 1 {
        return;
    }

    // Sort by score (best first) for overlap filtering
    kept_indices.sort_by(|&a, &b| {
        intervals[b]
            .score
            .partial_cmp(&intervals[a].score)
            .unwrap_or(Ordering::Equal)
    });

    // Keep only non-overlapping or below threshold
    let mut final_kept = vec![kept_indices[0]];

    for &idx in kept_indices.iter().skip(1) {
        let mut keep = true;
        for &kept_idx in &final_kept {
            if intervals[idx].overlaps(&intervals[kept_idx], overlap_threshold) {
                keep = false;
                break;
            }
        }
        if keep {
            final_kept.push(idx);
        }
    }

    *kept_indices = final_kept;
}
