use anyhow::Result;
use std::fmt;

/// Compact 32-byte mapping structure matching wfmash's MappingResult
#[repr(C, packed)]
#[derive(Debug, Clone, Copy)]
pub struct Mapping {
    pub ref_seq_id: u32,         // 4 bytes - reference sequence ID
    pub ref_start_pos: u32,      // 4 bytes - reference start position
    pub query_start_pos: u32,    // 4 bytes - query start position
    pub block_length: u32,       // 4 bytes - mapping block length
    pub n_merged: u32,           // 4 bytes - number of merged segments
    pub conserved_sketches: u32, // 4 bytes - count of conserved sketches
    pub nuc_identity: u16,       // 2 bytes - scaled identity (0-10000 for 0.00-100.00%)
    pub flags: u8,               // 1 byte - bit-packed flags (strand, discard, overlapped)
    pub kmer_complexity: u8,     // 1 byte - scaled kmer complexity (0-100)
}

impl Mapping {
    pub const FLAG_STRAND: u8 = 0x01;
    pub const FLAG_DISCARD: u8 = 0x02;
    pub const FLAG_OVERLAPPED: u8 = 0x04;

    pub fn is_reverse(&self) -> bool {
        (self.flags & Self::FLAG_STRAND) != 0
    }

    pub fn set_reverse(&mut self, rev: bool) {
        if rev {
            self.flags |= Self::FLAG_STRAND;
        } else {
            self.flags &= !Self::FLAG_STRAND;
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

    pub fn identity(&self) -> f64 {
        // nuc_identity stores 0-10000 for 0.00-100.00%
        self.nuc_identity as f64 / 10000.0 * 100.0 // Return as percentage
    }

    pub fn set_identity(&mut self, identity: f64) {
        // identity is 0.0-1.0, store as 0-10000
        self.nuc_identity = (identity * 10000.0).round() as u16;
    }

    pub fn ref_end_pos(&self) -> u32 {
        self.ref_start_pos + self.block_length
    }

    pub fn query_end_pos(&self) -> u32 {
        self.query_start_pos + self.block_length
    }
}

/// Chain status for mappings
#[derive(Debug, Clone, PartialEq)]
pub enum ChainStatus {
    Scaffold,   // Part of a scaffold chain (is an anchor)
    Rescued,    // Rescued mapping (near an anchor)
    Unassigned, // Not yet assigned to a chain
}

/// Auxiliary data for mappings during filtering/merging
#[derive(Debug, Clone)]
pub struct MappingAux {
    pub split_mapping_id: u64,
    pub chain_pair_score: f64,
    pub chain_pair_id: i64,
    pub query_seq_id: i32,
    pub query_len: u32,
    pub ref_name: String,
    pub query_name: String,
    pub chain_id: String,          // Chain identifier (e.g., "1.1.1")
    pub chain_status: ChainStatus, // Scaffold or rescued
    pub ref_len: u32,              // Reference sequence length
}

impl Default for MappingAux {
    fn default() -> Self {
        MappingAux {
            split_mapping_id: 0,
            chain_pair_score: f64::MAX,
            chain_pair_id: i64::MIN,
            query_seq_id: -1,
            query_len: 0,
            ref_name: String::new(),
            query_name: String::new(),
            chain_id: String::new(),
            chain_status: ChainStatus::Unassigned,
            ref_len: 0,
        }
    }
}

/// PAF record structure
#[derive(Debug, Clone)]
pub struct PafRecord {
    pub query_name: String,
    pub query_len: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub strand: char,
    pub ref_name: String,
    pub ref_len: u32,
    pub ref_start: u32,
    pub ref_end: u32,
    pub matches: u32,
    pub block_len: u32,
    pub quality: u8,
    pub tags: Vec<(String, String)>,
    pub cigar: Option<String>,
}

impl PafRecord {
    /// Calculate identity from CIGAR string using mutually gapped metric
    /// Returns (matches, total_positions) where:
    /// - matches: number of '=' operations in CIGAR
    /// - total_positions: matches + mismatches + insertions + deletions
    pub fn calculate_cigar_identity(&self) -> Option<f64> {
        if let Some(ref cigar) = self.cigar {
            let mut matches = 0u32;
            let mut mismatches = 0u32;
            let mut insertions = 0u32;
            let mut deletions = 0u32;

            let mut num_str = String::new();

            for ch in cigar.chars() {
                if ch.is_ascii_digit() {
                    num_str.push(ch);
                } else {
                    let count = num_str.parse::<u32>().unwrap_or(1);
                    num_str.clear();

                    match ch {
                        '=' => matches += count,    // Match
                        'X' => mismatches += count, // Mismatch
                        'M' => matches += count, // Match or mismatch (treat as match if no = available)
                        'I' => insertions += count, // Insertion to reference
                        'D' => deletions += count, // Deletion from reference
                        'N' | 'S' | 'H' | 'P' => {} // Skip these operations
                        _ => {}
                    }
                }
            }

            // Calculate mutually gapped identity:
            // numerator = matches
            // denominator = matches + mismatches + insertions + deletions
            let total = matches + mismatches + insertions + deletions;
            if total > 0 {
                Some(matches as f64 / total as f64)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Convert PAF record to internal mapping structure
    pub fn to_mapping(&self, ref_id: u32, query_id: i32) -> Mapping {
        let mut mapping = Mapping {
            ref_seq_id: ref_id,
            ref_start_pos: self.ref_start,
            query_start_pos: self.query_start,
            block_length: self.ref_end - self.ref_start,
            n_merged: 1,
            conserved_sketches: 0,
            nuc_identity: 0,
            flags: 0,
            kmer_complexity: 100,
        };

        // Set strand
        if self.strand == '-' {
            mapping.set_reverse(true);
        }

        // Calculate identity - prefer CIGAR-based if available
        let identity = if let Some(cigar_identity) = self.calculate_cigar_identity() {
            cigar_identity
        } else if self.block_len > 0 {
            // Fall back to matches/block_len from PAF columns
            self.matches as f64 / self.block_len as f64
        } else {
            0.0
        };
        mapping.set_identity(identity);

        mapping
    }
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.strand,
            self.ref_name,
            self.ref_len,
            self.ref_start,
            self.ref_end,
            self.matches,
            self.block_len,
            self.quality
        )?;

        // Write CIGAR if present
        if let Some(ref cigar) = self.cigar {
            write!(f, "\tcg:Z:{}", cigar)?;
        }

        // Write other tags
        for (key, val) in &self.tags {
            if key != "cg" {
                write!(f, "\t{}:{}", key, val)?;
            }
        }

        Ok(())
    }
}
