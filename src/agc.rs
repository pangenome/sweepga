//! AGC (Assembled Genome Compressor) support for sweepga
//!
//! Provides selective extraction of genome sequences from AGC archives
//! for alignment with FastGA. Supports batch-based extraction to minimize
//! disk usage and memory footprint.

use anyhow::{Context, Result};
use ragc_core::{Decompressor, DecompressorConfig};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Wrapper around RAGC Decompressor for sweepga integration
pub struct AgcSource {
    #[allow(dead_code)]
    path: PathBuf,
    decompressor: Decompressor,
}

/// Information about a sample in an AGC archive
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SampleInfo {
    pub name: String,
    pub total_bp: u64,
    pub num_contigs: usize,
}

#[allow(dead_code)]
impl AgcSource {
    /// Open an AGC archive for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let path_str = path.to_string_lossy().to_string();

        let config = DecompressorConfig { verbosity: 0 };
        let decompressor = Decompressor::open(&path_str, config)
            .with_context(|| format!("Failed to open AGC archive: {}", path.display()))?;

        Ok(AgcSource { path, decompressor })
    }

    /// Get the path to the AGC archive
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// List all sample names in the archive
    pub fn list_samples(&self) -> Vec<String> {
        self.decompressor.list_samples()
    }

    /// List samples matching a prefix (e.g., "SGDref#" or "SGDref#1#chrI")
    pub fn list_samples_with_prefix(&self, prefix: &str) -> Vec<String> {
        self.decompressor.list_samples_with_prefix(prefix)
    }

    /// Get sizes (in bp) for all samples without extracting sequences
    ///
    /// This is efficient because it only reads segment metadata (raw_length fields),
    /// not the actual compressed sequence data.
    pub fn get_sample_sizes(&mut self) -> Result<HashMap<String, u64>> {
        let samples = self.list_samples();
        self.get_sample_sizes_for(&samples)
    }

    /// Get sizes for specific samples
    pub fn get_sample_sizes_for(&mut self, samples: &[String]) -> Result<HashMap<String, u64>> {
        let mut sizes = HashMap::new();

        for sample in samples {
            let size = self.get_sample_size(sample)?;
            sizes.insert(sample.clone(), size);
        }

        Ok(sizes)
    }

    /// Get the total size (in bp) for a single sample without extracting
    pub fn get_sample_size(&mut self, sample_name: &str) -> Result<u64> {
        // List contigs triggers loading of segment metadata
        let contigs = self.decompressor.list_contigs(sample_name)?;

        let mut total = 0u64;
        for contig_name in &contigs {
            let segments = self
                .decompressor
                .get_contig_segments_desc(sample_name, contig_name)?;

            for seg in &segments {
                total += seg.raw_length as u64;
            }
        }

        Ok(total)
    }

    /// Get detailed information about all samples
    pub fn get_sample_info(&mut self) -> Result<Vec<SampleInfo>> {
        let samples = self.list_samples();
        let mut info = Vec::with_capacity(samples.len());

        for name in samples {
            let contigs = self.decompressor.list_contigs(&name)?;
            let num_contigs = contigs.len();

            let mut total_bp = 0u64;
            for contig_name in &contigs {
                let segments = self
                    .decompressor
                    .get_contig_segments_desc(&name, contig_name)?;
                for seg in &segments {
                    total_bp += seg.raw_length as u64;
                }
            }

            info.push(SampleInfo {
                name,
                total_bp,
                num_contigs,
            });
        }

        Ok(info)
    }

    /// Extract specific samples to a FASTA file
    pub fn extract_samples_to_fasta<P: AsRef<Path>>(
        &mut self,
        samples: &[String],
        output_path: P,
    ) -> Result<()> {
        let output_path = output_path.as_ref();
        let mut writer = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

        for sample_name in samples {
            self.write_sample_to_writer(sample_name, &mut writer)?;
        }

        Ok(())
    }

    /// Extract a single sample to a FASTA file
    pub fn extract_sample_to_fasta<P: AsRef<Path>>(
        &mut self,
        sample_name: &str,
        output_path: P,
    ) -> Result<()> {
        let output_path = output_path.as_ref();
        let mut writer = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

        self.write_sample_to_writer(sample_name, &mut writer)
    }

    /// Extract samples to a temporary FASTA file, returning the temp file handle
    pub fn extract_samples_to_temp(
        &mut self,
        samples: &[String],
        temp_dir: Option<&Path>,
    ) -> Result<tempfile::NamedTempFile> {
        let temp_file = if let Some(dir) = temp_dir {
            tempfile::Builder::new()
                .prefix("sweepga_agc_")
                .suffix(".fa")
                .tempfile_in(dir)?
        } else {
            tempfile::Builder::new()
                .prefix("sweepga_agc_")
                .suffix(".fa")
                .tempfile()?
        };

        {
            let mut writer = std::io::BufWriter::new(temp_file.as_file());
            for sample_name in samples {
                self.write_sample_to_writer(sample_name, &mut writer)?;
            }
            writer.flush()?;
        }

        Ok(temp_file)
    }

    /// Write a single sample to a writer in FASTA format
    fn write_sample_to_writer<W: Write>(
        &mut self,
        sample_name: &str,
        writer: &mut W,
    ) -> Result<()> {
        // Get all contigs for this sample
        let sample_data = self
            .decompressor
            .get_sample(sample_name)
            .with_context(|| format!("Failed to extract sample: {}", sample_name))?;

        // Base conversion table (numeric to ASCII)
        const BASE_MAP: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

        for (contig_name, contig_data) in sample_data {
            // Write FASTA header
            writeln!(writer, ">{}", contig_name)?;

            // Convert numeric bases to ASCII and write in 80-char lines
            let mut line_buf = Vec::with_capacity(80);
            for &base in &contig_data {
                let ascii = if (base as usize) < BASE_MAP.len() {
                    BASE_MAP[base as usize]
                } else {
                    b'N'
                };
                line_buf.push(ascii);

                if line_buf.len() == 80 {
                    writer.write_all(&line_buf)?;
                    writeln!(writer)?;
                    line_buf.clear();
                }
            }

            // Write remaining bases
            if !line_buf.is_empty() {
                writer.write_all(&line_buf)?;
                writeln!(writer)?;
            }
        }

        Ok(())
    }

    /// Clone the decompressor for thread-safe parallel access
    pub fn clone_for_thread(&self) -> Result<Self> {
        let decompressor = self.decompressor.clone_for_thread()?;
        Ok(AgcSource {
            path: self.path.clone(),
            decompressor,
        })
    }

    /// Extract a sample's sequence directly to a Vec<u8> (ASCII bases)
    /// This avoids creating temporary files for operations like minhash sketching
    pub fn extract_sample_to_bytes(&mut self, sample_name: &str) -> Result<Vec<u8>> {
        let sample_data = self
            .decompressor
            .get_sample(sample_name)
            .with_context(|| format!("Failed to extract sample: {}", sample_name))?;

        // Base conversion table (numeric to ASCII)
        const BASE_MAP: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

        // Calculate total size for pre-allocation
        let total_len: usize = sample_data.iter().map(|(_, data)| data.len()).sum();
        let mut sequence = Vec::with_capacity(total_len);

        for (_contig_name, contig_data) in sample_data {
            for &base in &contig_data {
                let ascii = if (base as usize) < BASE_MAP.len() {
                    BASE_MAP[base as usize]
                } else {
                    b'N'
                };
                sequence.push(ascii);
            }
        }

        Ok(sequence)
    }
}

/// Parse a sample list specification
///
/// Supports:
/// - Comma-separated list: "sample1,sample2,sample3"
/// - File reference: "@samples.txt" (one sample per line)
pub fn parse_sample_list(spec: &str) -> Result<Vec<String>> {
    if let Some(file_path) = spec.strip_prefix('@') {
        // Read from file
        let content = std::fs::read_to_string(file_path)
            .with_context(|| format!("Failed to read sample list file: {}", file_path))?;

        Ok(content
            .lines()
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && !s.starts_with('#'))
            .map(|s| s.to_string())
            .collect())
    } else {
        // Comma-separated list
        Ok(spec
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect())
    }
}

/// Detect if a file is an AGC archive by extension
#[allow(dead_code)]
pub fn is_agc_file<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .map(|ext| ext.eq_ignore_ascii_case("agc"))
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sample_list_comma() {
        let result = parse_sample_list("sample1,sample2,sample3").unwrap();
        assert_eq!(result, vec!["sample1", "sample2", "sample3"]);
    }

    #[test]
    fn test_parse_sample_list_with_spaces() {
        let result = parse_sample_list("sample1, sample2 , sample3").unwrap();
        assert_eq!(result, vec!["sample1", "sample2", "sample3"]);
    }

    #[test]
    fn test_is_agc_file() {
        assert!(is_agc_file("genome.agc"));
        assert!(is_agc_file("path/to/genome.AGC"));
        assert!(!is_agc_file("genome.fa"));
        assert!(!is_agc_file("genome.paf"));
    }
}
