use anyhow::{Result, Context};
use std::path::Path;
use std::process::Command;
use tempfile::NamedTempFile;

/// Simple FastGA integration that outputs to .1aln format
/// Uses direct command-line invocation for .1aln output since fastga-rs
/// doesn't yet support .1aln output format
pub struct FastGAIntegration {
    num_threads: usize,
}

impl FastGAIntegration {
    /// Create a new FastGA integration with minimal configuration
    pub fn new(num_threads: usize) -> Self {
        FastGAIntegration { num_threads }
    }

    /// Run FastGA alignment and write output to a temporary .1aln file
    /// Returns the temporary file handle (which auto-deletes when dropped)
    pub fn align_to_temp_1aln(
        &self,
        queries: &Path,
        targets: &Path,
    ) -> Result<NamedTempFile> {
        // Create temporary .1aln output file
        let mut temp_file = NamedTempFile::with_suffix(".1aln")
            .context("Failed to create temporary .1aln file")?;

        let temp_path = temp_file.path().to_str()
            .context("Invalid temp file path")?;

        // Find FastGA binary (from fastga-rs build)
        let fastga_bin = std::env::var("FASTGA_PATH")
            .unwrap_or_else(|_| {
                // Try to find in target/release/build
                let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").ok();
                if let Some(dir) = manifest_dir {
                    // Look for FastGA in build directory
                    let build_dir = format!("{}/target/release/build", dir);
                    if let Ok(entries) = std::fs::read_dir(&build_dir) {
                        for entry in entries.flatten() {
                            if entry.file_name().to_string_lossy().starts_with("fastga-rs-") {
                                let fastga_path = entry.path().join("out/FastGA");
                                if fastga_path.exists() {
                                    return fastga_path.to_string_lossy().to_string();
                                }
                            }
                        }
                    }
                }
                "FastGA".to_string() // Hope it's in PATH
            });

        eprintln!("[fastga] Running alignment");
        eprintln!("[fastga] Executing: FastGA -1:{} -T{} {:?} {:?}",
                  temp_path, self.num_threads, queries, targets);

        // Run FastGA with .1aln output
        let output = Command::new(&fastga_bin)
            .arg(format!("-1:{}", temp_path))
            .arg(format!("-T{}", self.num_threads))
            .arg(queries)
            .arg(targets)
            .output()
            .context("Failed to execute FastGA")?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!("FastGA failed: {}", stderr);
        }

        Ok(temp_file)
    }
}