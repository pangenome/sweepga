/// Sequence name registry for mapping between integer IDs and string names
///
/// This is only used for PAF input/output. When working with .1aln files,
/// we work directly with integer IDs and skip string conversion entirely.

use std::collections::HashMap;

/// Registry that maps sequence IDs to names
#[derive(Debug, Default)]
pub struct SequenceRegistry {
    /// Map from reference ID to name
    ref_id_to_name: HashMap<u32, String>,

    /// Map from reference name to ID
    ref_name_to_id: HashMap<String, u32>,

    /// Map from query ID to name
    query_id_to_name: HashMap<i32, String>,

    /// Map from query name to ID
    query_name_to_id: HashMap<String, i32>,

    /// Next available reference ID
    next_ref_id: u32,

    /// Next available query ID
    next_query_id: i32,
}

impl SequenceRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get or assign a reference ID for a name
    pub fn get_or_assign_ref_id(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.ref_name_to_id.get(name) {
            return id;
        }

        let id = self.next_ref_id;
        self.next_ref_id += 1;
        self.ref_name_to_id.insert(name.to_string(), id);
        self.ref_id_to_name.insert(id, name.to_string());
        id
    }

    /// Get or assign a query ID for a name
    pub fn get_or_assign_query_id(&mut self, name: &str) -> i32 {
        if let Some(&id) = self.query_name_to_id.get(name) {
            return id;
        }

        let id = self.next_query_id;
        self.next_query_id += 1;
        self.query_name_to_id.insert(name.to_string(), id);
        self.query_id_to_name.insert(id, name.to_string());
        id
    }

    /// Get reference name by ID
    pub fn get_ref_name(&self, id: u32) -> Option<&str> {
        self.ref_id_to_name.get(&id).map(|s| s.as_str())
    }

    /// Get query name by ID
    pub fn get_query_name(&self, id: i32) -> Option<&str> {
        self.query_id_to_name.get(&id).map(|s| s.as_str())
    }

    /// Register a reference sequence (e.g., from .1aln GDB file)
    pub fn register_ref(&mut self, id: u32, name: String) {
        self.ref_id_to_name.insert(id, name.clone());
        self.ref_name_to_id.insert(name, id);
        if id >= self.next_ref_id {
            self.next_ref_id = id + 1;
        }
    }

    /// Register a query sequence (e.g., from .1aln GDB file)
    pub fn register_query(&mut self, id: i32, name: String) {
        self.query_id_to_name.insert(id, name.clone());
        self.query_name_to_id.insert(name, id);
        if id >= self.next_query_id {
            self.next_query_id = id + 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_registry() {
        let mut registry = SequenceRegistry::new();

        // Assign new IDs
        let id1 = registry.get_or_assign_ref_id("chr1");
        let id2 = registry.get_or_assign_ref_id("chr2");
        assert_eq!(id1, 0);
        assert_eq!(id2, 1);

        // Get existing ID
        let id1_again = registry.get_or_assign_ref_id("chr1");
        assert_eq!(id1_again, 0);

        // Lookup by ID
        assert_eq!(registry.get_ref_name(0), Some("chr1"));
        assert_eq!(registry.get_ref_name(1), Some("chr2"));
    }
}
