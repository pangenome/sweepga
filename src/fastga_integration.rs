use anyhow::{Result, Context};
use fastga_rs::{FastGA, Config};
use std::path::Path;
use tempfile::NamedTempFile;
use std::io::Write;

/// Simple FastGA integration using the fastga-rs FFI bindings
/// This version just runs FastGA and writes output to a temp PAF file
pub struct FastGAIntegration {
    config: Config,
}

impl FastGAIntegration {
    /// Create a new FastGA integration with minimal configuration
    pub fn new(num_threads: usize) -> Self {
        // Use FastGA defaults for most parameters, only set threads
        let config = Config::builder()
            .num_threads(num_threads)
            .build();

        FastGAIntegration { config }
    }

    /// Create with custom parameters (for future use)
    #[allow(dead_code)]
    pub fn new_with_params(min_identity: f64, min_alignment_length: u32, num_threads: usize) -> Self {
        let config = Config::builder()
            .min_identity(min_identity)
            .min_alignment_length(min_alignment_length as usize)
            .num_threads(num_threads)
            .build();

        FastGAIntegration { config }
    }

    /// Run FastGA alignment and write output to a temporary PAF file
    /// Returns the temporary file handle (which auto-deletes when dropped)
    pub fn align_to_temp_paf(
        &self,
        queries: &Path,
        targets: &Path,
    ) -> Result<NamedTempFile> {
        // Create FastGA aligner
        let aligner = FastGA::new(self.config.clone())
            .context("Failed to create FastGA aligner")?;

        // Run alignment
        let alignments = aligner.align_files(queries, targets)
            .context("Failed to run FastGA alignment")?;

        // Convert to PAF format with extended CIGAR
        let paf_output = alignments.to_paf()
            .context("Failed to convert alignments to PAF format")?;

        // Write to temporary file
        let mut temp_file = NamedTempFile::new()
            .context("Failed to create temporary PAF file")?;

        temp_file.write_all(paf_output.as_bytes())
            .context("Failed to write PAF output to temporary file")?;

        temp_file.flush()
            .context("Failed to flush temporary file")?;

        Ok(temp_file)
    }
}