/// Merge multiple PAF files, taking the union of alignments
///
/// This module provides functionality to combine multiple filtered PAF files,
/// keeping an alignment if it appears in ANY of the input files (union operation).
use anyhow::Result;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Alignment signature for deduplication
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct AlignmentSignature {
    query_name: String,
    target_name: String,
    query_start: usize,
    target_start: usize,
    query_end: usize,
    target_end: usize,
}

impl AlignmentSignature {
    fn from_paf_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }

        Some(AlignmentSignature {
            query_name: fields[0].to_string(),
            target_name: fields[5].to_string(),
            query_start: fields[2].parse().ok()?,
            target_start: fields[7].parse().ok()?,
            query_end: fields[3].parse().ok()?,
            target_end: fields[8].parse().ok()?,
        })
    }
}

/// Merge multiple PAF files into one, taking the union of alignments
///
/// # Arguments
/// * `input_paths` - Paths to input PAF files
/// * `output_path` - Path to output merged PAF file
///
/// # Returns
/// Number of unique alignments in the merged output
pub fn merge_paf_files(input_paths: &[&str], output_path: &str) -> Result<usize> {
    // First pass: collect all unique alignment signatures and their first occurrence
    let mut seen_signatures = HashSet::new();
    let mut signature_to_line = std::collections::HashMap::new();

    for input_path in input_paths {
        let input_file = File::open(input_path)?;
        let reader = BufReader::new(input_file);

        for line in reader.lines() {
            let line = line?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            if let Some(sig) = AlignmentSignature::from_paf_line(&line) {
                if seen_signatures.insert(sig.clone()) {
                    // First time seeing this alignment - save it
                    signature_to_line.insert(sig, line);
                }
            }
        }
    }

    // Write all unique alignments to output
    let mut output_file = File::create(output_path)?;
    let count = signature_to_line.len();

    for (_sig, line) in signature_to_line.iter() {
        writeln!(output_file, "{}", line)?;
    }

    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_merge_identical_alignments() {
        // Create two PAF files with identical alignments
        let mut file1 = NamedTempFile::new().unwrap();
        writeln!(
            file1,
            "seq1\t100\t10\t90\t+\tseq2\t200\t20\t180\t80\t160\t60"
        )
        .unwrap();

        let mut file2 = NamedTempFile::new().unwrap();
        writeln!(
            file2,
            "seq1\t100\t10\t90\t+\tseq2\t200\t20\t180\t80\t160\t60"
        )
        .unwrap();

        let output = NamedTempFile::new().unwrap();

        let count = merge_paf_files(
            &[
                file1.path().to_str().unwrap(),
                file2.path().to_str().unwrap(),
            ],
            output.path().to_str().unwrap(),
        )
        .unwrap();

        assert_eq!(count, 1); // Should deduplicate to 1 alignment
    }

    #[test]
    fn test_merge_different_alignments() {
        let mut file1 = NamedTempFile::new().unwrap();
        writeln!(
            file1,
            "seq1\t100\t10\t90\t+\tseq2\t200\t20\t180\t80\t160\t60"
        )
        .unwrap();

        let mut file2 = NamedTempFile::new().unwrap();
        writeln!(
            file2,
            "seq3\t100\t10\t90\t+\tseq4\t200\t20\t180\t80\t160\t60"
        )
        .unwrap();

        let output = NamedTempFile::new().unwrap();

        let count = merge_paf_files(
            &[
                file1.path().to_str().unwrap(),
                file2.path().to_str().unwrap(),
            ],
            output.path().to_str().unwrap(),
        )
        .unwrap();

        assert_eq!(count, 2); // Should have 2 different alignments
    }
}
