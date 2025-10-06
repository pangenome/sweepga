#![allow(dead_code)]
/// Efficient sequence name indexing for compact storage
use std::collections::HashMap;

/// Maps sequence names to compact integer IDs
#[derive(Debug, Clone, Default)]
pub struct SequenceIndex {
    /// All unique sequence names (ID is the index in this vec)
    names: Vec<String>,
    /// Map from name to ID for fast lookup
    name_to_id: HashMap<String, u32>,
}

impl SequenceIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self::default()
    }

    /// Get or create an ID for a sequence name
    /// Returns the ID (creating a new one if the name is new)
    pub fn get_or_insert(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(name) {
            id
        } else {
            let id = self.names.len() as u32;
            self.names.push(name.to_string());
            self.name_to_id.insert(name.to_string(), id);
            id
        }
    }

    /// Get the ID for a name (returns None if not found)
    pub fn get_id(&self, name: &str) -> Option<u32> {
        self.name_to_id.get(name).copied()
    }

    /// Get the name for an ID
    pub fn get_name(&self, id: u32) -> Option<&str> {
        self.names.get(id as usize).map(|s| s.as_str())
    }

    /// Get the name for an ID (panics if invalid)
    pub fn name(&self, id: u32) -> &str {
        &self.names[id as usize]
    }

    /// Number of unique sequences
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Get all sequence names
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Clear the index
    pub fn clear(&mut self) {
        self.names.clear();
        self.name_to_id.clear();
    }

    /// Reserve capacity for efficiency
    pub fn reserve(&mut self, additional: usize) {
        self.names.reserve(additional);
        self.name_to_id.reserve(additional);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_index() {
        let mut index = SequenceIndex::new();

        // First insertion
        let id1 = index.get_or_insert("chr1");
        assert_eq!(id1, 0);

        // Same name returns same ID
        let id1_again = index.get_or_insert("chr1");
        assert_eq!(id1_again, 0);

        // Different name gets new ID
        let id2 = index.get_or_insert("chr2");
        assert_eq!(id2, 1);

        // Lookup by ID
        assert_eq!(index.get_name(0), Some("chr1"));
        assert_eq!(index.get_name(1), Some("chr2"));
        assert_eq!(index.get_name(999), None);

        // Lookup by name
        assert_eq!(index.get_id("chr1"), Some(0));
        assert_eq!(index.get_id("chr2"), Some(1));
        assert_eq!(index.get_id("chr3"), None);

        assert_eq!(index.len(), 2);
    }
}
