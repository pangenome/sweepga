use anyhow::{Context, Result};
use fastga_rs::{Config, FastGA};
use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use tempfile::NamedTempFile;

/// FastGA parameter presets optimized for different ANI levels
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum FastGAPreset {
    Ani70,
    Ani80,
    Ani85,
    Ani90,
    Ani95,
    Ani99,
}

impl FastGAPreset {
    /// Convert preset to FastGA Config
    #[allow(clippy::wrong_self_convention, dead_code)]
    pub fn to_config(&self, num_threads: usize) -> Config {
        match self {
            FastGAPreset::Ani70 => Config::builder()
                .min_identity(0.70)
                .min_alignment_length(100)
                .chain_start_threshold(1000)
                .min_chain_coverage(0.85)
                .frequency(10)
                .chain_break(2000)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
            FastGAPreset::Ani80 => Config::builder()
                .min_identity(0.78)
                .min_alignment_length(200)
                .chain_start_threshold(1500)
                .min_chain_coverage(0.87)
                .frequency(10)
                .chain_break(2000)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
            FastGAPreset::Ani85 => Config::builder()
                .min_identity(0.83)
                .min_alignment_length(300)
                .chain_start_threshold(1800)
                .min_chain_coverage(0.88)
                .frequency(8)
                .chain_break(1800)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
            FastGAPreset::Ani90 => Config::builder()
                .min_identity(0.88)
                .min_alignment_length(400)
                .chain_start_threshold(2500)
                .min_chain_coverage(0.90)
                .frequency(8)
                .chain_break(1500)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
            FastGAPreset::Ani95 => Config::builder()
                .min_identity(0.93)
                .min_alignment_length(500)
                .chain_start_threshold(3000)
                .min_chain_coverage(0.92)
                .frequency(5)
                .chain_break(1200)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
            FastGAPreset::Ani99 => Config::builder()
                .min_identity(0.98)
                .min_alignment_length(1000)
                .chain_start_threshold(5000)
                .min_chain_coverage(0.95)
                .frequency(5)
                .chain_break(1000)
                .num_threads(num_threads)
                .verbose(true)
                .build(),
        }
    }

    /// List all available presets with descriptions
    #[allow(dead_code)]
    pub fn list_all() -> Vec<(&'static str, &'static str)> {
        vec![
            ("ani70", "Default - Distant relatives (20-30% divergence)"),
            ("ani80", "Moderately distant (15-20% divergence)"),
            ("ani85", "Related species (10-15% divergence)"),
            ("ani90", "Close relatives (5-10% divergence)"),
            (
                "ani95",
                "Closely related - Human-chimp level (2-5% divergence)",
            ),
            ("ani99", "Nearly identical - Human strains (<1% divergence)"),
        ]
    }
}

impl FromStr for FastGAPreset {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "ani70" | "70" => Ok(FastGAPreset::Ani70),
            "ani80" | "80" => Ok(FastGAPreset::Ani80),
            "ani85" | "85" => Ok(FastGAPreset::Ani85),
            "ani90" | "90" => Ok(FastGAPreset::Ani90),
            "ani95" | "95" => Ok(FastGAPreset::Ani95),
            "ani99" | "99" => Ok(FastGAPreset::Ani99),
            _ => anyhow::bail!(
                "Unknown preset '{}'. Valid presets: ani70, ani80, ani85, ani90, ani95, ani99",
                s
            ),
        }
    }
}

/// Simple FastGA integration using the fastga-rs FFI bindings
/// This version just runs FastGA and writes output to a temp PAF file
pub struct FastGAIntegration {
    config: Config,
}

impl FastGAIntegration {
    /// Create a new FastGA integration with optional frequency parameter
    pub fn new(frequency: Option<usize>, num_threads: usize) -> Self {
        let mut builder = Config::builder().num_threads(num_threads).verbose(true);

        if let Some(freq) = frequency {
            builder = builder.adaptive_seed_cutoff(freq);
        }

        let config = builder.build();
        FastGAIntegration { config }
    }

    /// Create with custom parameters (for future use)
    #[allow(dead_code)]
    pub fn new_with_params(
        min_identity: f64,
        min_alignment_length: u32,
        num_threads: usize,
    ) -> Self {
        let config = Config::builder()
            .min_identity(min_identity)
            .min_alignment_length(min_alignment_length as usize)
            .num_threads(num_threads)
            .verbose(true)
            .build();

        FastGAIntegration { config }
    }

    /// Run FastGA alignment and write output to a temporary PAF file
    /// Returns the temporary file handle (which auto-deletes when dropped)
    pub fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Create FastGA aligner
        let aligner =
            FastGA::new(self.config.clone()).context("Failed to create FastGA aligner")?;

        // Run alignment
        let alignments = aligner
            .align_files(queries, targets)
            .context("Failed to run FastGA alignment")?;

        eprintln!("[sweepga] Got {} alignments", alignments.len());

        // Convert to PAF format with extended CIGAR
        let paf_output = alignments
            .to_paf()
            .context("Failed to convert alignments to PAF format")?;

        eprintln!(
            "[sweepga] PAF output: {} bytes, {} lines",
            paf_output.len(),
            paf_output.lines().count()
        );

        // Write to temporary file
        let mut temp_file = NamedTempFile::new().context("Failed to create temporary PAF file")?;

        temp_file
            .write_all(paf_output.as_bytes())
            .context("Failed to write PAF output to temporary file")?;

        temp_file
            .flush()
            .context("Failed to flush temporary file")?;

        Ok(temp_file)
    }
}
