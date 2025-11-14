//! .1aln reader for tree-based sparsification
//!
//! This module provides .1aln file reading functionality specifically for
//! tree-based sparsification (see tree_sparsify.rs).
//!
//! Key features:
//! - Reads .1aln files using fastga-rs AlnReader
//! - Extracts alignments with identity calculation for building genome-pair matrices
//! - Used exclusively by tree_sparsify module for identity-based filtering
//!
//! NOTE: For general .1aln filtering (applying plane sweep, scaffolding, etc.),
//! use the unified_filter module instead, which preserves .1aln format without conversion.

use anyhow::{Context, Result};
use std::path::Path;

/// Alignment record from .1aln with identity calculation
/// Used by tree_sparsify for building genome-pair identity matrices
#[derive(Debug, Clone)]
#[allow(dead_code)] // Some fields only used in tests or for future functionality
pub struct AlnAlignment {
    pub query_id: i64,
    pub query_name: String,
    pub query_start: i64,
    pub query_end: i64,
    pub query_len: i64,
    pub target_id: i64,
    pub target_name: String,
    pub target_start: i64,
    pub target_end: i64,
    pub target_len: i64,
    pub reverse: bool,
    pub matches: i64,
    pub diffs: i64,
    pub identity: f64, // Computed from X field
    pub mapping_quality: i64,
    // Store X values for potential re-writing
    pub x_values: Vec<i64>,
}

/// Reader for .1aln files that extracts X field for identity calculation
/// Uses fastga-rs native reader for compatibility
pub struct AlnFilterReader {
    reader: fastga_rs::AlnReader,
}

impl AlnFilterReader {
    /// Open a .1aln file for filtering using fastga-rs native reader
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_str().context("Invalid path")?;
        let reader = fastga_rs::AlnReader::open(path_str)
            .context(format!("Failed to open .1aln file: {path_str}"))?;

        Ok(AlnFilterReader { reader })
    }

    /// Read next alignment using fastga-rs native reader
    pub fn read_alignment(&mut self) -> Result<Option<AlnAlignment>> {
        // Use fastga-rs reader
        if let Some(rec) = self.reader.read_record()? {
            // Get sequence names
            let query_name = self.reader.get_seq_name(rec.query_id, 0)?;
            let target_name = self.reader.get_seq_name(rec.target_id, 1)?;

            // Calculate identity from diffs field
            let aln_len = (rec.query_end - rec.query_start) as usize;
            let matches = aln_len.saturating_sub(rec.diffs as usize);
            let identity = if aln_len > 0 {
                matches as f64 / aln_len as f64
            } else {
                0.0
            };

            let alignment = AlnAlignment {
                query_id: rec.query_id,
                query_name,
                query_start: rec.query_start,
                query_end: rec.query_end,
                query_len: rec.query_len,
                target_id: rec.target_id,
                target_name,
                target_start: rec.target_start,
                target_end: rec.target_end,
                target_len: rec.target_len,
                reverse: rec.reverse != 0,
                matches: matches as i64,
                diffs: rec.diffs as i64,
                identity,
                mapping_quality: 60,
                x_values: Vec::new(), // Not extracted anymore
            };

            Ok(Some(alignment))
        } else {
            Ok(None)
        }
    }

    /// Read all alignments from the file
    #[allow(dead_code)] // Kept for convenience, though not currently used
    pub fn read_all(&mut self) -> Result<Vec<AlnAlignment>> {
        let mut alignments = Vec::new();

        while let Some(alignment) = self.read_alignment()? {
            alignments.push(alignment);
        }

        Ok(alignments)
    }
}

impl AlnAlignment {
    /// Convert to PAF format line (with X-based identity in dv:f: tag)
    #[allow(dead_code)] // Kept for potential future use
    pub fn to_paf_line(&self) -> String {
        let divergence = 1.0 - self.identity;
        let strand = if self.reverse { '-' } else { '+' };
        let block_len = self.query_end - self.query_start;

        // PAF field 10 is matches - compute from identity and block_len
        // identity = matches / block_len, so matches = identity * block_len
        let matches = (self.identity * block_len as f64).round() as i64;

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tdv:f:{:.6}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            strand,
            self.target_name,
            self.target_len,
            self.target_start,
            self.target_end,
            matches,
            block_len,
            self.mapping_quality,
            divergence,
        )
    }
}

// NOTE: General .1aln filtering functions have been moved to unified_filter.rs
// This module now only contains the reader needed for tree_sparsify

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // Requires .1aln test file
    fn test_fastga_rs_reader() {
        // Test using fastga-rs's AlnReader directly
        use fastga_rs::onelib::AlnReader;

        let mut reader = AlnReader::open("test.1aln").unwrap();

        println!("Testing fastga-rs AlnReader...");
        if let Some(aln) = reader.read_alignment().unwrap() {
            println!(
                "First alignment: query='{}' target='{}'",
                aln.query_name, aln.target_name
            );
            println!("  matches={} block_len={}", aln.matches, aln.block_len);
        }
    }

    #[test]
    #[ignore] // Requires .1aln test file
    fn test_read_1aln_with_x_field() {
        let mut reader = AlnFilterReader::open("test.1aln").unwrap();

        let mut count = 0;
        while let Some(aln) = reader.read_alignment().unwrap() {
            println!(
                "Alignment {}: query='{}' target='{}' identity={:.4} (from {} X values)",
                count,
                aln.query_name,
                aln.target_name,
                aln.identity,
                aln.x_values.len()
            );

            // Verify identity is in valid range
            assert!(
                aln.identity >= 0.0 && aln.identity <= 1.0,
                "Invalid identity: {}",
                aln.identity
            );

            count += 1;
            if count >= 10 {
                break;
            }
        }

        assert!(count > 0, "Should have read at least one alignment");
    }

    // NOTE: Test for filtering pipeline removed - use unified_filter module for .1aln filtering

    #[test]
    #[ignore] // Requires .1aln test file
    fn test_identity_from_x_field() {
        let mut reader = AlnFilterReader::open("test.1aln").unwrap();

        // Read first few alignments and check X-based identity
        let mut count = 0;
        while let Some(aln) = reader.read_alignment().unwrap() {
            if !aln.x_values.is_empty() {
                let edit_distance: i64 = aln.x_values.iter().sum();
                let alignment_length = aln.query_end - aln.query_start;
                let expected_identity = if alignment_length > 0 {
                    (alignment_length - edit_distance) as f64 / alignment_length as f64
                } else {
                    0.0
                };

                // Identity should match what we calculated
                assert!(
                    (aln.identity - expected_identity).abs() < 0.001,
                    "Identity mismatch: got {}, expected {} (from X field)",
                    aln.identity,
                    expected_identity
                );

                println!("✓ Alignment {}: identity {:.4} correctly computed from X field (edit_distance={})",
                         count, aln.identity, edit_distance);
            }

            count += 1;
            if count >= 10 {
                break;
            }
        }

        assert!(
            count > 0,
            "Should have read at least one alignment with X values"
        );
    }
}
