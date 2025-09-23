use anyhow::{Result, bail};
use std::collections::{HashMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::cmp::Ordering;

/// Compact mapping record - only numerical data (32-40 bytes)
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CompactMapping {
    pub query_id: u32,        // 4 bytes - sequence ID from index
    pub target_id: u32,       // 4 bytes - sequence ID from index
    pub query_start: u32,     // 4 bytes
    pub query_end: u32,       // 4 bytes
    pub target_start: u32,    // 4 bytes
    pub target_end: u32,      // 4 bytes
    pub identity: u16,        // 2 bytes - scaled 0-10000 for 0-100%
    pub strand: u8,           // 1 byte - 0 for +, 1 for -
    pub flags: u8,            // 1 byte - for discard/overlapped flags
    pub file_offset: u64,     // 8 bytes - offset in PAF file
    // Total: 36 bytes
}

impl CompactMapping {
    pub fn score(&self) -> f64 {
        let identity = (self.identity as f64) / 10000.0;
        let length = ((self.query_end - self.query_start) as f64).max(1.0);
        identity * length.ln()
    }

    pub fn query_len(&self) -> u32 {
        self.query_end - self.query_start
    }

    pub fn target_len(&self) -> u32 {
        self.target_end - self.target_start
    }
}

/// Sequence name index for compact storage
pub struct SequenceIndex {
    name_to_id: HashMap<String, u32>,
    id_to_name: Vec<String>,
    next_id: u32,
}

impl SequenceIndex {
    pub fn new() -> Self {
        SequenceIndex {
            name_to_id: HashMap::new(),
            id_to_name: Vec::new(),
            next_id: 0,
        }
    }

    pub fn get_or_create_id(&mut self, name: &str) -> Result<u32> {
        if let Some(&id) = self.name_to_id.get(name) {
            return Ok(id);
        }

        // Check if we're about to overflow u32
        if self.next_id == u32::MAX {
            bail!("Too many unique sequence names (exceeded u32::MAX). This is a limitation of the current implementation.");
        }

        let id = self.next_id;
        self.name_to_id.insert(name.to_string(), id);
        self.id_to_name.push(name.to_string());
        self.next_id += 1;
        Ok(id)
    }

    pub fn get_name(&self, id: u32) -> Option<&str> {
        self.id_to_name.get(id as usize).map(|s| s.as_str())
    }

    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }
}

/// Memory-efficient PAF filter
pub struct CompactPafFilter {
    mappings: Vec<CompactMapping>,
    query_index: SequenceIndex,
    target_index: SequenceIndex,
    min_block_length: u32,
    keep_self: bool,
    prefix_delimiter: char,
    group_by_prefix: bool,
}

impl CompactPafFilter {
    pub fn new(min_block_length: u32, keep_self: bool, prefix_delimiter: char, group_by_prefix: bool) -> Self {
        CompactPafFilter {
            mappings: Vec::new(),
            query_index: SequenceIndex::new(),
            target_index: SequenceIndex::new(),
            min_block_length,
            keep_self,
            prefix_delimiter,
            group_by_prefix,
        }
    }

    /// Read PAF file and build compact representation
    pub fn read_paf(&mut self, input_path: &str) -> Result<()> {
        let file = File::open(input_path)?;
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let mut offset: u64 = 0;

        eprintln!("Reading PAF and building sequence indices...");

        while reader.read_line(&mut line)? > 0 {
            let line_start_offset = offset;
            offset += line.len() as u64;

            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 12 {
                line.clear();
                continue;
            }

            // Parse PAF fields
            let query_name = fields[0];
            let query_start = fields[2].parse::<u32>().unwrap_or(0);
            let query_end = fields[3].parse::<u32>().unwrap_or(0);
            let strand = if fields[4] == "+" { 0u8 } else { 1u8 };
            let target_name = fields[5];
            let target_start = fields[7].parse::<u32>().unwrap_or(0);
            let target_end = fields[8].parse::<u32>().unwrap_or(0);
            let matches = fields[9].parse::<u32>().unwrap_or(0);
            let block_len = fields[10].parse::<u32>().unwrap_or(1);

            // Skip if below minimum block length
            if block_len < self.min_block_length {
                line.clear();
                continue;
            }

            // Skip self-mappings unless requested
            if !self.keep_self && query_name == target_name {
                line.clear();
                continue;
            }

            // Calculate identity (check for divergence tag first)
            let mut identity = (matches as f64) / (block_len.max(1) as f64);

            // Look for divergence tag (dv:f:) or CIGAR string
            for field in &fields[11..] {
                if field.starts_with("dv:f:") {
                    if let Ok(div) = field[5..].parse::<f64>() {
                        identity = 1.0 - div;
                    }
                } else if field.starts_with("cg:Z:") {
                    // Could parse CIGAR here if needed
                    // For now, rely on divergence or matches/block_len
                }
            }

            // Get or create sequence IDs
            let query_id = self.query_index.get_or_create_id(query_name)?;
            let target_id = self.target_index.get_or_create_id(target_name)?;

            // Create compact mapping
            let mapping = CompactMapping {
                query_id,
                target_id,
                query_start,
                query_end,
                target_start,
                target_end,
                identity: (identity * 10000.0).round() as u16,
                strand,
                flags: 0,
                file_offset: line_start_offset,
            };

            self.mappings.push(mapping);
            line.clear();
        }

        eprintln!("Loaded {} mappings", self.mappings.len());
        eprintln!("Unique query sequences: {}", self.query_index.len());
        eprintln!("Unique target sequences: {}", self.target_index.len());

        // Check memory usage estimate
        let memory_mb = (self.mappings.len() * std::mem::size_of::<CompactMapping>()) as f64 / 1_048_576.0;
        eprintln!("Estimated memory usage for mappings: {:.1} MB", memory_mb);

        Ok(())
    }

    /// Extract prefix from name for grouping
    fn extract_prefix(&self, name: &str) -> String {
        if !self.group_by_prefix {
            return name.to_string();
        }

        if let Some(last_pos) = name.rfind(self.prefix_delimiter) {
            name[..=last_pos].to_string()
        } else {
            name.to_string()
        }
    }

    /// Group mappings by prefix pairs
    fn group_by_prefix_pairs(&self) -> HashMap<(String, String), Vec<usize>> {
        let mut groups = HashMap::new();

        for (idx, mapping) in self.mappings.iter().enumerate() {
            let query_name = self.query_index.get_name(mapping.query_id).unwrap_or("");
            let target_name = self.target_index.get_name(mapping.target_id).unwrap_or("");

            let query_prefix = self.extract_prefix(query_name);
            let target_prefix = self.extract_prefix(target_name);

            groups.entry((query_prefix, target_prefix))
                .or_insert_with(Vec::new)
                .push(idx);
        }

        groups
    }

    /// Apply plane sweep filtering
    pub fn apply_plane_sweep(&mut self, filter_spec: &str, overlap_threshold: f64) -> Result<Vec<u64>> {
        let (query_limit, target_limit) = parse_filter_spec(filter_spec);

        // If no filtering, return all offsets
        if query_limit.is_none() && target_limit.is_none() {
            return Ok(self.mappings.iter().map(|m| m.file_offset).collect());
        }

        // Group by prefix pairs
        let groups = self.group_by_prefix_pairs();
        eprintln!("Processing {} genome pair groups", groups.len());

        let mut kept_offsets = Vec::new();

        for ((query_prefix, target_prefix), group_indices) in groups {
            if group_indices.is_empty() {
                continue;
            }

            // Apply plane sweep within this group
            let group_mappings: Vec<&mut CompactMapping> = group_indices.iter()
                .map(|&idx| unsafe { &mut *(&mut self.mappings[idx] as *mut CompactMapping) })
                .collect();

            let kept = apply_plane_sweep_to_group(group_mappings, query_limit, target_limit, overlap_threshold);

            // Collect file offsets of kept mappings
            for &idx in &kept {
                kept_offsets.push(self.mappings[group_indices[idx]].file_offset);
            }
        }

        eprintln!("Kept {} mappings after plane sweep", kept_offsets.len());
        Ok(kept_offsets)
    }

    /// Get all file offsets (for debug mode)
    pub fn get_all_offsets(&self) -> Vec<u64> {
        self.mappings.iter().map(|m| m.file_offset).collect()
    }

    /// Write filtered output using file offsets
    pub fn write_filtered_output(&self, input_path: &str, output_path: &str, kept_offsets: Vec<u64>) -> Result<()> {
        let mut input = BufReader::new(File::open(input_path)?);
        let mut output = BufWriter::new(File::create(output_path)?);

        // Sort offsets for sequential reading
        let mut kept_offsets = kept_offsets;
        kept_offsets.sort_unstable();

        let mut line = String::new();
        let mut current_offset: u64 = 0;
        let mut offset_idx = 0;

        eprintln!("Writing {} filtered records...", kept_offsets.len());

        while input.read_line(&mut line)? > 0 {
            if offset_idx < kept_offsets.len() && current_offset == kept_offsets[offset_idx] {
                output.write_all(line.as_bytes())?;
                offset_idx += 1;
            }
            current_offset += line.len() as u64;
            line.clear();
        }

        output.flush()?;
        Ok(())
    }
}

/// Parse filter specification
fn parse_filter_spec(spec: &str) -> (Option<usize>, Option<usize>) {
    let spec = spec.trim();

    if spec == "N" || spec == "N:N" {
        return (None, None);
    }

    if !spec.contains(':') {
        if let Ok(n) = spec.parse::<usize>() {
            return (Some(n), None);
        }
    }

    if let Some((query_part, target_part)) = spec.split_once(':') {
        let query_limit = if query_part == "N" {
            None
        } else {
            query_part.parse::<usize>().ok()
        };

        let target_limit = if target_part == "N" {
            None
        } else {
            target_part.parse::<usize>().ok()
        };

        return (query_limit, target_limit);
    }

    (None, None)
}

/// Event for plane sweep
#[derive(Debug, Clone)]
struct Event {
    position: u32,
    event_type: u8,  // 0=begin, 1=end
    mapping_idx: usize,
}

/// Mapping order for BST
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
        other.score.partial_cmp(&self.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.start_pos.cmp(&other.start_pos))
    }
}

/// Apply plane sweep to a group of mappings
fn apply_plane_sweep_to_group(
    mut mappings: Vec<&mut CompactMapping>,
    query_limit: Option<usize>,
    target_limit: Option<usize>,
    overlap_threshold: f64,
) -> Vec<usize> {
    if mappings.is_empty() {
        return Vec::new();
    }

    // Mark all as discarded initially
    for mapping in &mut mappings {
        mapping.flags = 1;  // discard flag
    }

    // Apply query axis filtering
    if let Some(limit) = query_limit {
        if limit > 0 {
            plane_sweep_query(&mut mappings, limit, overlap_threshold);
        }
    }

    // Apply target axis filtering
    if let Some(limit) = target_limit {
        if limit > 0 {
            plane_sweep_target(&mut mappings, limit, overlap_threshold);
        }
    }

    // Collect indices of kept mappings
    mappings.iter()
        .enumerate()
        .filter(|(_, m)| m.flags == 0)  // not discarded
        .map(|(idx, _)| idx)
        .collect()
}

/// Plane sweep on query axis
fn plane_sweep_query(mappings: &mut [&mut CompactMapping], limit: usize, overlap_threshold: f64) {
    // Create events
    let mut events = Vec::with_capacity(mappings.len() * 2);
    for (idx, mapping) in mappings.iter().enumerate() {
        events.push(Event {
            position: mapping.query_start,
            event_type: 0,  // begin
            mapping_idx: idx,
        });
        events.push(Event {
            position: mapping.query_end,
            event_type: 1,  // end
            mapping_idx: idx,
        });
    }

    // Sort events
    events.sort_by_key(|e| (e.position, e.event_type));

    // Plane sweep with BST
    let mut active = BTreeSet::new();
    let mut i = 0;

    while i < events.len() {
        let current_pos = events[i].position;

        // Process all events at current position
        let mut j = i;
        while j < events.len() && events[j].position == current_pos {
            let event = &events[j];
            let mapping_order = MappingOrder {
                idx: event.mapping_idx,
                score: mappings[event.mapping_idx].score(),
                start_pos: mappings[event.mapping_idx].query_start,
            };

            if event.event_type == 0 {  // begin
                active.insert(mapping_order);
            } else {  // end
                active.remove(&mapping_order);
            }
            j += 1;
        }

        // Mark best mappings as good
        let mut kept = 0;
        for mapping_order in &active {
            if kept >= limit {
                break;
            }
            mappings[mapping_order.idx].flags = 0;  // not discarded
            kept += 1;
        }

        i = j;
    }
}

/// Plane sweep on target axis
fn plane_sweep_target(mappings: &mut [&mut CompactMapping], limit: usize, overlap_threshold: f64) {
    // Similar to query axis but using target coordinates
    let mut events = Vec::with_capacity(mappings.len() * 2);
    for (idx, mapping) in mappings.iter().enumerate() {
        events.push(Event {
            position: mapping.target_start,
            event_type: 0,
            mapping_idx: idx,
        });
        events.push(Event {
            position: mapping.target_end,
            event_type: 1,
            mapping_idx: idx,
        });
    }

    events.sort_by_key(|e| (e.position, e.event_type));

    let mut active = BTreeSet::new();
    let mut i = 0;

    while i < events.len() {
        let current_pos = events[i].position;

        let mut j = i;
        while j < events.len() && events[j].position == current_pos {
            let event = &events[j];
            let mapping_order = MappingOrder {
                idx: event.mapping_idx,
                score: mappings[event.mapping_idx].score(),
                start_pos: mappings[event.mapping_idx].target_start,
            };

            if event.event_type == 0 {
                active.insert(mapping_order);
            } else {
                active.remove(&mapping_order);
            }
            j += 1;
        }

        let mut kept = 0;
        for mapping_order in &active {
            if kept >= limit {
                break;
            }
            mappings[mapping_order.idx].flags = 0;
            kept += 1;
        }

        i = j;
    }
}