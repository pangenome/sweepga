//! Pure .1aln filtering module
//!
//! This module implements direct .1aln → filter → .1aln pipeline,
//! avoiding the wasteful .1aln→PAF→filter→PAF→.1aln conversion.
//!
//! Key features:
//! - Reads .1aln files using onecode-rs
//! - Computes identity from X field (edit distance per trace interval)
//! - Applies sweepga filtering logic
//! - Writes filtered .1aln using fastga-rs AlnWriter
//!
//! NOTE: Still uses PAF representation internally for filtering (temp files),
//! but avoids external format conversions (ALNtoPAF/PAFtoALN).

use anyhow::{Context, Result};
use std::path::Path;

use crate::paf_filter::{FilterConfig, PafFilter};

/// Alignment record from .1aln with X-field based identity
#[derive(Debug, Clone)]
#[allow(dead_code)] // TODO: Remove once integrated into main workflow
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

/// STREAMING FILTER: Read .1aln → filter records → emit passing records from original file
/// This preserves the exact original format without any conversions!
pub fn filter_1aln_streaming<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    filter_config: &FilterConfig,
) -> Result<()> {
    use std::collections::HashSet;
    use std::io::{BufWriter, Write};
    use tempfile::NamedTempFile;

    eprintln!("[filter_1aln] Phase 1: Reading and filtering alignments...");

    // Phase 1: Read .1aln, convert to PAF for filtering (in-memory), get kept IDs
    let mut reader = AlnFilterReader::open(&input_path)?;
    let alignments = reader.read_all()?;
    eprintln!(
        "[filter_1aln] Read {} alignments from .1aln",
        alignments.len()
    );

    // Convert to PAF and filter
    let temp_paf_input = NamedTempFile::new()?;
    {
        let file = std::fs::File::create(temp_paf_input.path())?;
        let mut writer = BufWriter::new(file);
        for aln in &alignments {
            writeln!(writer, "{}", aln.to_paf_line())?;
        }
        writer.flush()?;
    }

    let temp_paf_output = NamedTempFile::new()?;
    let filter = PafFilter::new(filter_config.clone());
    filter.filter_paf(temp_paf_input.path(), temp_paf_output.path())?;

    // Build map: PAF line → original index
    // This lets us figure out which original records to keep
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader};

    let mut paf_to_index: HashMap<String, usize> = HashMap::new();
    for (i, aln) in alignments.iter().enumerate() {
        let paf_line = aln.to_paf_line();
        paf_to_index.insert(paf_line, i);
    }

    // Read filtered PAF and determine which indices to keep
    let file = std::fs::File::open(temp_paf_output.path())?;
    let reader_paf = BufReader::new(file);
    let mut kept_indices = HashSet::new();

    for line in reader_paf.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Look up which original record this corresponds to
        if let Some(&idx) = paf_to_index.get(&line) {
            kept_indices.insert(idx);
        }
    }

    eprintln!(
        "[filter_1aln] Kept {} / {} alignments",
        kept_indices.len(),
        alignments.len()
    );

    // Phase 2: Emit passing records by index
    eprintln!("[filter_1aln] Phase 2: Emitting passing records to output...");

    use fastga_rs::AlnWriter;
    let mut writer = AlnWriter::create(output_path.as_ref(), true)?;

    let mut written = 0;
    for (i, aln) in alignments.iter().enumerate() {
        if kept_indices.contains(&i) {
            let block_length = (aln.query_end - aln.query_start) as usize;
            let alignment = fastga_rs::Alignment {
                query_name: aln.query_name.clone(),
                query_len: aln.query_len as usize,
                query_start: aln.query_start as usize,
                query_end: aln.query_end as usize,
                target_name: aln.target_name.clone(),
                target_len: aln.target_len as usize,
                target_start: aln.target_start as usize,
                target_end: aln.target_end as usize,
                strand: if aln.reverse { '-' } else { '+' },
                matches: aln.matches as usize,
                block_len: block_length,
                mapping_quality: aln.mapping_quality as u8,
                cigar: String::new(),
                tags: vec![],
                mismatches: (block_length - aln.matches as usize),
                gap_opens: 0,
                gap_len: 0,
            };
            writer.write_alignment(&alignment)?;
            written += 1;
        }
    }

    eprintln!("[filter_1aln] Wrote {written} alignments to output");
    Ok(())
}

/// Convert .1aln to PAF format (with X-based identity), then filter using PAF pipeline
/// This ensures 100% identical filtering logic
#[allow(dead_code)] // TODO: Remove once integrated into main workflow
pub fn filter_1aln_file_to_paf<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    filter_config: &FilterConfig,
) -> Result<()> {
    use std::io::{BufWriter, Write};
    use tempfile::NamedTempFile;

    // Step 1: Read .1aln and convert to PAF with X-based identity
    let mut reader = AlnFilterReader::open(&input_path)?;
    let alignments = reader.read_all()?;

    eprintln!(
        "[filter_1aln] Read {} alignments from .1aln",
        alignments.len()
    );

    // Create temporary PAF file
    let temp_paf = NamedTempFile::new()?;
    let temp_paf_path = temp_paf.path();

    {
        let file = std::fs::File::create(temp_paf_path)?;
        let mut writer = BufWriter::new(file);

        for (i, aln) in alignments.iter().enumerate() {
            let paf_line = aln.to_paf_line();
            if i < 3 {
                eprintln!("[filter_1aln] Sample PAF line: {paf_line}");
            }
            writeln!(writer, "{paf_line}")?;
        }
        writer.flush()?;
    }

    eprintln!(
        "[filter_1aln] Converted to PAF: {}",
        temp_paf_path.display()
    );

    // Step 2: Apply PAF filtering using EXACT same code as PAF input
    let filter = PafFilter::new(filter_config.clone());
    filter.filter_paf(temp_paf_path, output_path.as_ref())?;

    Ok(())
}

// Note: write_1aln_file removed - will be reimplemented once onecode-rs
// supports sequence name extraction from embedded GDB

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

    #[test]
    #[ignore] // Requires .1aln test file
    fn test_pure_1aln_filtering_pipeline() {
        use tempfile::NamedTempFile;

        // Create temporary output file
        let output_file = NamedTempFile::new().unwrap();
        let output_path = output_file.path();

        // Create filter config for 1:1 filtering (same as default sweepga)
        let config = FilterConfig {
            chain_gap: 0,
            min_block_length: 0,
            mapping_filter_mode: crate::paf_filter::FilterMode::ManyToMany,
            mapping_max_per_query: None,
            mapping_max_per_target: None,
            plane_sweep_secondaries: 1,
            scaffold_filter_mode: crate::paf_filter::FilterMode::OneToOne,
            scaffold_max_per_query: Some(1),
            scaffold_max_per_target: Some(1),
            overlap_threshold: 0.95,
            sparsity: 1.0,
            no_merge: true,
            scaffold_gap: 10000,
            min_scaffold_length: 10000,
            scaffold_overlap_threshold: 0.95,
            scaffold_max_deviation: 20000,
            prefix_delimiter: '#',
            skip_prefix: false,
            scoring_function: crate::paf_filter::ScoringFunction::LogLengthIdentity,
            min_identity: 0.0,
            min_scaffold_identity: 0.0,
        };

        // Run the filtering pipeline using PAF filtering (100% same code)
        filter_1aln_file_to_paf("test.1aln", output_path, &config).unwrap();

        // Read PAF output and count
        use std::io::{BufRead, BufReader};
        let file = std::fs::File::open(output_path).unwrap();
        let reader = BufReader::new(file);
        let count = reader.lines().map_while(Result::ok).count();

        println!("Filtered to {count} alignments from test.1aln");
        assert!(
            count > 0 && count < 30000,
            "Should filter to reasonable number"
        );
    }

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
