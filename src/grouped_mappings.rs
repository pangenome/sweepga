#![allow(dead_code)]
/// Efficient grouped storage for mappings
use std::collections::HashMap;
use crate::compact_mapping::CompactRecordMeta;
use crate::sequence_index::SequenceIndex;

/// Key for grouping mappings by chromosome pair
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
pub struct ChromosomePair {
    pub query_id: u32,
    pub target_id: u32,
}

/// Storage for mappings pre-grouped by chromosome pair
#[derive(Debug)]
pub struct GroupedMappings {
    /// Groups of mappings indexed by chromosome pair
    groups: HashMap<ChromosomePair, Vec<CompactRecordMeta>>,
    /// Sequence name index
    seq_index: SequenceIndex,
}

impl GroupedMappings {
    /// Create new grouped mappings
    pub fn new() -> Self {
        Self {
            groups: HashMap::new(),
            seq_index: SequenceIndex::new(),
        }
    }

    /// Add a mapping to its group
    pub fn add_mapping(&mut self, meta: CompactRecordMeta) {
        let key = ChromosomePair {
            query_id: meta.query_id,
            target_id: meta.target_id,
        };
        self.groups.entry(key).or_default().push(meta);
    }

    /// Get the sequence index
    pub fn seq_index(&self) -> &SequenceIndex {
        &self.seq_index
    }

    /// Get mutable sequence index
    pub fn seq_index_mut(&mut self) -> &mut SequenceIndex {
        &mut self.seq_index
    }

    /// Get all groups
    pub fn groups(&self) -> &HashMap<ChromosomePair, Vec<CompactRecordMeta>> {
        &self.groups
    }

    /// Get mappings for a specific chromosome pair
    pub fn get_group(&self, query_id: u32, target_id: u32) -> Option<&Vec<CompactRecordMeta>> {
        let key = ChromosomePair { query_id, target_id };
        self.groups.get(&key)
    }

    /// Get mutable mappings for a specific chromosome pair
    pub fn get_group_mut(&mut self, query_id: u32, target_id: u32) -> Option<&mut Vec<CompactRecordMeta>> {
        let key = ChromosomePair { query_id, target_id };
        self.groups.get_mut(&key)
    }

    /// Number of groups
    pub fn num_groups(&self) -> usize {
        self.groups.len()
    }

    /// Total number of mappings across all groups
    pub fn num_mappings(&self) -> usize {
        self.groups.values().map(|v| v.len()).sum()
    }

    /// Clear all mappings but keep the sequence index
    pub fn clear_mappings(&mut self) {
        self.groups.clear();
    }

    /// Iterate over all groups
    pub fn iter(&self) -> impl Iterator<Item = (&ChromosomePair, &Vec<CompactRecordMeta>)> {
        self.groups.iter()
    }

    /// Iterate over all groups mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&ChromosomePair, &mut Vec<CompactRecordMeta>)> {
        self.groups.iter_mut()
    }

    /// Reserve capacity for efficiency
    pub fn reserve(&mut self, additional_groups: usize) {
        self.groups.reserve(additional_groups);
    }

    /// Build from a vector of mappings
    pub fn from_mappings(mappings: Vec<CompactRecordMeta>, seq_index: SequenceIndex) -> Self {
        let mut grouped = Self {
            groups: HashMap::new(),
            seq_index,
        };

        for meta in mappings {
            grouped.add_mapping(meta);
        }

        grouped
    }

    /// Apply a filter function to all groups
    pub fn filter_groups<F>(&mut self, mut filter_fn: F)
    where
        F: FnMut(&ChromosomePair, &mut Vec<CompactRecordMeta>),
    {
        for (key, mappings) in self.groups.iter_mut() {
            filter_fn(key, mappings);
        }
    }

    /// Get summary statistics
    pub fn stats(&self) -> GroupedMappingsStats {
        let group_sizes: Vec<usize> = self.groups.values().map(|v| v.len()).collect();
        let total = group_sizes.iter().sum();
        let min = group_sizes.iter().min().copied().unwrap_or(0);
        let max = group_sizes.iter().max().copied().unwrap_or(0);
        let mean = if !group_sizes.is_empty() {
            total as f64 / group_sizes.len() as f64
        } else {
            0.0
        };

        GroupedMappingsStats {
            num_groups: self.groups.len(),
            num_mappings: total,
            min_group_size: min,
            max_group_size: max,
            mean_group_size: mean,
        }
    }
}

/// Statistics about grouped mappings
#[derive(Debug)]
pub struct GroupedMappingsStats {
    pub num_groups: usize,
    pub num_mappings: usize,
    pub min_group_size: usize,
    pub max_group_size: usize,
    pub mean_group_size: f64,
}

impl Default for GroupedMappings {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapping::ChainStatus;

    fn make_test_mapping(query_id: u32, target_id: u32, start: u64) -> CompactRecordMeta {
        CompactRecordMeta {
            rank: 0,
            query_id,
            target_id,
            query_start: start,
            query_end: start + 1000,
            target_start: start,
            target_end: start + 1000,
            block_length: 1000,
            identity: 1.0,
            strand: '+',
            chain_id: None,
            chain_status: ChainStatus::Unassigned,
            discard: false,
            overlapped: false,
        }
    }

    #[test]
    fn test_grouped_mappings() {
        let mut grouped = GroupedMappings::new();

        // Add mappings to different groups
        grouped.add_mapping(make_test_mapping(0, 0, 1000));
        grouped.add_mapping(make_test_mapping(0, 0, 2000));
        grouped.add_mapping(make_test_mapping(0, 1, 3000));
        grouped.add_mapping(make_test_mapping(1, 0, 4000));

        assert_eq!(grouped.num_groups(), 3);
        assert_eq!(grouped.num_mappings(), 4);

        // Check specific group
        let group_0_0 = grouped.get_group(0, 0).unwrap();
        assert_eq!(group_0_0.len(), 2);

        let group_0_1 = grouped.get_group(0, 1).unwrap();
        assert_eq!(group_0_1.len(), 1);

        // Check stats
        let stats = grouped.stats();
        assert_eq!(stats.num_groups, 3);
        assert_eq!(stats.num_mappings, 4);
        assert_eq!(stats.min_group_size, 1);
        assert_eq!(stats.max_group_size, 2);
    }
}