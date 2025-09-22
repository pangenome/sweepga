use std::collections::HashMap;
use anyhow::Result;

/// Compact representation of PAF records for efficient filtering
pub struct CompactPafIndex {
    /// String interning table: name -> ID
    sequence_names: Vec<String>,
    sequence_name_to_id: HashMap<String, u32>,

    /// Per-record data (all parallel arrays)
    byte_offsets: Vec<u64>,      // File position of record start
    byte_lengths: Vec<u32>,       // Length of record in bytes
    query_ids: Vec<u32>,          // Query sequence ID
    target_ids: Vec<u32>,         // Target sequence ID
    query_starts: Vec<u32>,       // Query start position
    query_ends: Vec<u32>,         // Query end position
    target_starts: Vec<u32>,      // Target start position
    target_ends: Vec<u32>,        // Target end position
    block_lengths: Vec<u32>,      // Alignment block length
    strands: Vec<bool>,           // true = forward, false = reverse
}

impl CompactPafIndex {
    pub fn new() -> Self {
        CompactPafIndex {
            sequence_names: Vec::new(),
            sequence_name_to_id: HashMap::new(),
            byte_offsets: Vec::new(),
            byte_lengths: Vec::new(),
            query_ids: Vec::new(),
            target_ids: Vec::new(),
            query_starts: Vec::new(),
            query_ends: Vec::new(),
            target_starts: Vec::new(),
            target_ends: Vec::new(),
            block_lengths: Vec::new(),
            strands: Vec::new(),
        }
    }

    /// Intern a sequence name, returning its ID
    pub fn intern_name(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.sequence_name_to_id.get(name) {
            return id;
        }

        let id = self.sequence_names.len() as u32;
        self.sequence_names.push(name.to_string());
        self.sequence_name_to_id.insert(name.to_string(), id);
        id
    }

    /// Get sequence name by ID
    pub fn get_name(&self, id: u32) -> Option<&str> {
        self.sequence_names.get(id as usize).map(|s| s.as_str())
    }

    /// Add a record to the index
    pub fn add_record(
        &mut self,
        byte_offset: u64,
        byte_length: u32,
        query_name: &str,
        target_name: &str,
        query_start: u32,
        query_end: u32,
        target_start: u32,
        target_end: u32,
        block_length: u32,
        strand: char,
    ) {
        let query_id = self.intern_name(query_name);
        let target_id = self.intern_name(target_name);

        self.byte_offsets.push(byte_offset);
        self.byte_lengths.push(byte_length);
        self.query_ids.push(query_id);
        self.target_ids.push(target_id);
        self.query_starts.push(query_start);
        self.query_ends.push(query_end);
        self.target_starts.push(target_start);
        self.target_ends.push(target_end);
        self.block_lengths.push(block_length);
        self.strands.push(strand == '+');
    }

    /// Get number of records
    pub fn len(&self) -> usize {
        self.byte_offsets.len()
    }

    /// Get record data by index (for filtering algorithms)
    pub fn get_record(&self, idx: usize) -> Option<CompactRecord> {
        if idx >= self.len() {
            return None;
        }

        Some(CompactRecord {
            idx,
            byte_offset: self.byte_offsets[idx],
            byte_length: self.byte_lengths[idx],
            query_id: self.query_ids[idx],
            target_id: self.target_ids[idx],
            query_start: self.query_starts[idx],
            query_end: self.query_ends[idx],
            target_start: self.target_starts[idx],
            target_end: self.target_ends[idx],
            block_length: self.block_lengths[idx],
            strand: self.strands[idx],
        })
    }
}

/// Temporary struct for accessing record data during filtering
#[derive(Debug, Clone, Copy)]
pub struct CompactRecord {
    pub idx: usize,
    pub byte_offset: u64,
    pub byte_length: u32,
    pub query_id: u32,
    pub target_id: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub block_length: u32,
    pub strand: bool,
}