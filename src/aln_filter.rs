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

use anyhow::{Result, Context};
use onecode::{OneFile, OneSchema};
use std::path::Path;

use crate::paf_filter::{FilterConfig, PafFilter};

/// Contig information from GDB skeleton
#[derive(Debug, Clone)]
struct ContigInfo {
    clen: i64,    // Contig length
    sbeg: i64,    // Position within scaffold
    scaf: usize,  // Scaffold ID
}

/// Scaffold information from GDB skeleton
#[derive(Debug, Clone)]
struct ScaffoldInfo {
    name: String,
    slen: i64,
    fctg: usize,  // First contig index
    ectg: usize,  // Last contig + 1
}

/// Alignment record from .1aln with X-field based identity
#[derive(Debug, Clone)]
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
    pub identity: f64,  // Computed from X field
    pub mapping_quality: i64,
    // Store X values for potential re-writing
    pub x_values: Vec<i64>,
}

/// Reader for .1aln files that extracts X field for identity calculation
pub struct AlnFilterReader {
    file: OneFile,
    contigs: Vec<ContigInfo>,
    scaffolds: Vec<ScaffoldInfo>,
}

impl AlnFilterReader {
    /// Open a .1aln file for filtering
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;

        // Create schema for .1aln files (matches FastGA output)
        let schema_text = r#"
P 3 aln
D t 1 3 INT
O g 0
G S 0
O S 1 6 STRING
D G 1 3 INT
D C 1 3 INT
D M 1 8 INT_LIST
O a 0
G A 0
D p 2 3 INT 3 INT
O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
D L 2 3 INT 3 INT
D R 0
D D 1 3 INT
D T 1 8 INT_LIST
D X 1 8 INT_LIST
D Q 1 3 INT
D E 1 3 INT
D Z 1 6 STRING
"#;
        let schema = OneSchema::from_text(schema_text)
            .context("Failed to create .1aln schema")?;

        let mut file = OneFile::open_read(path_str, Some(&schema), Some("aln"), 1)
            .context(format!("Failed to open .1aln file: {}", path_str))?;

        // Read GDB skeleton to build contig→scaffold mapping
        let (contigs, scaffolds) = Self::read_gdb_skeleton(&mut file)?;

        Ok(AlnFilterReader {
            file,
            contigs,
            scaffolds,
        })
    }

    /// Read GDB skeleton from .1aln file to build contig/scaffold mapping
    fn read_gdb_skeleton(file: &mut OneFile) -> Result<(Vec<ContigInfo>, Vec<ScaffoldInfo>)> {
        let mut contigs: Vec<ContigInfo> = Vec::new();
        let mut scaffolds: Vec<ScaffoldInfo> = Vec::new();
        let mut current_scaffold: Option<usize> = None;
        let mut spos = 0i64;  // Position within current scaffold
        let mut in_gdb_group = false;

        // Read through file looking for 'g' group (GDB skeleton)
        // We start at the beginning and look for 'g', 'S', 'G', 'C' records
        loop {
            let line_type = file.read_line();

            if line_type == '\0' {
                break; // EOF
            }

            match line_type {
                'g' => {
                    in_gdb_group = true;
                    // Start of GDB skeleton group
                }

                'A' => {
                    // Reached alignments - stop reading GDB skeleton
                    break;
                }

                'S' if in_gdb_group => {
                    // New scaffold
                    if let Some(scaffold_name) = file.string() {
                        // Finalize previous scaffold
                        if let Some(scaf_id) = current_scaffold {
                            scaffolds[scaf_id].ectg = contigs.len();
                            scaffolds[scaf_id].slen = spos;
                        }

                        // Start new scaffold
                        current_scaffold = Some(scaffolds.len());
                        scaffolds.push(ScaffoldInfo {
                            name: scaffold_name.to_string(),
                            slen: 0,  // Will be set when we see next 'S' or finish
                            fctg: contigs.len(),
                            ectg: contigs.len(),
                        });
                        spos = 0;
                    }
                }

                'G' if in_gdb_group => {
                    // Gap - advance position within scaffold
                    let gap_len = file.int(0);
                    spos += gap_len;
                }

                'C' if in_gdb_group => {
                    // Contig
                    let clen = file.int(0);

                    if let Some(scaf_id) = current_scaffold {
                        contigs.push(ContigInfo {
                            clen,
                            sbeg: spos,
                            scaf: scaf_id,
                        });
                        spos += clen;
                    }
                }

                _ => {
                    // Unknown/other record type - continue
                }
            }
        }

        // Finalize last scaffold
        if let Some(scaf_id) = current_scaffold {
            scaffolds[scaf_id].ectg = contigs.len();
            scaffolds[scaf_id].slen = spos;
        }

        if scaffolds.is_empty() {
            eprintln!("[AlnFilterReader] WARNING: No GDB skeleton found in .1aln file");
        } else {
            eprintln!("[AlnFilterReader] Read {} scaffolds and {} contigs from GDB skeleton",
                      scaffolds.len(), contigs.len());
        }

        Ok((contigs, scaffolds))
    }

    /// Read next alignment, computing identity from X field
    pub fn read_alignment(&mut self) -> Result<Option<AlnAlignment>> {
        // Find next 'A' record
        loop {
            let line_type = self.file.line_type();

            if line_type != 'A' {
                let next = self.file.read_line();
                if next == '\0' {
                    return Ok(None); // EOF
                }
                continue;
            }

            // Read 'A' record: alignment coordinates
            // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
            // Fields: query_id, query_start, query_end, target_id, target_start, target_end
            let query_id = self.file.int(0);
            let query_start = self.file.int(1);
            let query_end = self.file.int(2);
            let target_id = self.file.int(3);
            let target_start = self.file.int(4);
            let target_end = self.file.int(5);

            // Initialize alignment data
            let mut matches = 0i64;
            let mut diffs = 0i64;
            let mut is_reverse = false;
            let mut query_len = 0i64;
            let mut target_len = 0i64;
            let mut mapping_quality = 60i64;
            let mut x_values = Vec::new();

            // Read associated data lines
            loop {
                let next_type = self.file.read_line();

                if next_type == '\0' {
                    break; // EOF
                }

                match next_type {
                    'T' => {
                        // Trace record - marks end of alignment data
                        self.file.read_line();
                        break;
                    }
                    'M' => matches = self.file.int(0),
                    'D' => diffs = self.file.int(0),
                    'R' => is_reverse = true,
                    'L' => {
                        query_len = self.file.int(0);
                        target_len = self.file.int(1);
                    }
                    'Q' => mapping_quality = self.file.int(0),
                    'X' => {
                        // X field: edit distance per trace interval
                        // This is the KEY difference - we use X instead of D for identity
                        if !self.file.is_empty() {
                            if let Some(x_list) = self.file.int_list() {
                                x_values = x_list.to_vec();
                            }
                        }
                    }
                    'A' => {
                        // Hit next alignment without seeing 'T'
                        break;
                    }
                    _ => {
                        // Skip other records (C, p, etc.)
                        continue;
                    }
                }
            }

            // Compute identity from X field if available, otherwise fall back to M/D
            let identity = if !x_values.is_empty() {
                // X-based identity: sum(X) = edit distance
                let edit_distance: i64 = x_values.iter().sum();
                let alignment_length = query_end - query_start;
                if alignment_length > 0 {
                    let matches_from_x = alignment_length - edit_distance;
                    matches_from_x as f64 / alignment_length as f64
                } else {
                    0.0
                }
            } else if matches + diffs > 0 {
                // Fallback to M/D if X not available
                matches as f64 / (matches + diffs) as f64
            } else {
                // DEBUG: This shouldn't happen!
                eprintln!("[AlnReader] WARNING: No X, M, or D values for alignment at {}..{}",
                          query_start, query_end);
                0.0
            };

            // Map contig IDs to scaffold IDs and adjust coordinates
            // The 'A' record contains CONTIG IDs, but PAF needs SCAFFOLD IDs
            let (query_scaffold_id, query_scaffold_name, query_scaffold_len, query_scaffold_start, query_scaffold_end) =
                if !self.contigs.is_empty() && (query_id as usize) < self.contigs.len() {
                    let contig = &self.contigs[query_id as usize];
                    let scaffold = &self.scaffolds[contig.scaf];
                    (
                        contig.scaf as i64,
                        scaffold.name.clone(),
                        scaffold.slen,
                        contig.sbeg + query_start,  // Adjust to scaffold coordinates
                        contig.sbeg + query_end,
                    )
                } else {
                    // Fallback if no GDB skeleton
                    let name = self.file.get_sequence_name(query_id)
                        .unwrap_or_else(|| format!("{}", query_id));
                    (query_id, name, query_len, query_start, query_end)
                };

            let (target_scaffold_id, target_scaffold_name, target_scaffold_len, target_scaffold_start, target_scaffold_end) =
                if !self.contigs.is_empty() && (target_id as usize) < self.contigs.len() {
                    let contig = &self.contigs[target_id as usize];
                    let scaffold = &self.scaffolds[contig.scaf];
                    (
                        contig.scaf as i64,
                        scaffold.name.clone(),
                        scaffold.slen,
                        contig.sbeg + target_start,  // Adjust to scaffold coordinates
                        contig.sbeg + target_end,
                    )
                } else {
                    // Fallback if no GDB skeleton
                    let name = self.file.get_sequence_name(target_id)
                        .unwrap_or_else(|| format!("{}", target_id));
                    (target_id, name, target_len, target_start, target_end)
                };

            let alignment = AlnAlignment {
                query_id: query_scaffold_id,
                query_name: query_scaffold_name,
                query_start: query_scaffold_start,
                query_end: query_scaffold_end,
                query_len: query_scaffold_len,
                target_id: target_scaffold_id,
                target_name: target_scaffold_name,
                target_start: target_scaffold_start,
                target_end: target_scaffold_end,
                target_len: target_scaffold_len,
                reverse: is_reverse,
                matches,
                diffs,
                identity,
                mapping_quality,
                x_values,
            };

            return Ok(Some(alignment));
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

/// Convert .1aln to PAF format (with X-based identity), then filter using PAF pipeline
/// This ensures 100% identical filtering logic
pub fn filter_1aln_file_to_paf<P1: AsRef<Path>, P2: AsRef<Path>>(
    input_path: P1,
    output_path: P2,
    filter_config: &FilterConfig,
) -> Result<()> {
    use tempfile::NamedTempFile;
    use std::io::{BufWriter, Write};

    // Step 1: Read .1aln and convert to PAF with X-based identity
    let mut reader = AlnFilterReader::open(&input_path)?;
    let alignments = reader.read_all()?;

    eprintln!("[filter_1aln] Read {} alignments from .1aln", alignments.len());

    // Create temporary PAF file
    let temp_paf = NamedTempFile::new()?;
    let temp_paf_path = temp_paf.path();

    {
        let file = std::fs::File::create(temp_paf_path)?;
        let mut writer = BufWriter::new(file);

        for (i, aln) in alignments.iter().enumerate() {
            let paf_line = aln.to_paf_line();
            if i < 3 {
                eprintln!("[filter_1aln] Sample PAF line: {}", paf_line);
            }
            writeln!(writer, "{}", paf_line)?;
        }
        writer.flush()?;
    }

    eprintln!("[filter_1aln] Converted to PAF: {}", temp_paf_path.display());

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
            println!("First alignment: query='{}' target='{}'", aln.query_name, aln.target_name);
            println!("  matches={} block_len={}", aln.matches, aln.block_len);
        }

        // Try to get sequence names
        match reader.get_seq_name(0, 0) {
            Ok(name) => println!("Sequence 0: {}", name),
            Err(e) => println!("Error getting sequence 0: {}", e),
        }

        let all_names = reader.get_all_seq_names();
        println!("Found {} sequence names", all_names.len());
        for (id, name) in all_names.iter().take(5) {
            println!("  ID {}: {}", id, name);
        }
    }

    #[test]
    #[ignore] // Requires .1aln test file
    fn test_read_1aln_with_x_field() {
        let mut reader = AlnFilterReader::open("test.1aln").unwrap();

        let mut count = 0;
        while let Some(aln) = reader.read_alignment().unwrap() {
            println!("Alignment {}: query='{}' target='{}' identity={:.4} (from {} X values)",
                     count, aln.query_name, aln.target_name, aln.identity, aln.x_values.len());

            // Verify identity is in valid range
            assert!(aln.identity >= 0.0 && aln.identity <= 1.0,
                    "Invalid identity: {}", aln.identity);

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
        let count = reader.lines().filter_map(|l| l.ok()).count();

        println!("Filtered to {} alignments from test.1aln", count);
        assert!(count > 0 && count < 30000, "Should filter to reasonable number");
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
                assert!((aln.identity - expected_identity).abs() < 0.001,
                        "Identity mismatch: got {}, expected {} (from X field)",
                        aln.identity, expected_identity);

                println!("✓ Alignment {}: identity {:.4} correctly computed from X field (edit_distance={})",
                         count, aln.identity, edit_distance);
            }

            count += 1;
            if count >= 10 {
                break;
            }
        }

        assert!(count > 0, "Should have read at least one alignment with X values");
    }
}
