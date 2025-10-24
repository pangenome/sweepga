use anyhow::{Context, Result};
use fastga_rs::Config;
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

    /// Convert FASTA to GDB format (creates .gdb and .gix files)
    /// Returns the path to the .gdb file (without extension)
    /// Ensure fastga-rs can find FastGA binaries by adding them to PATH temporarily
    /// This is needed because fastga-rs uses which/PATH to find binaries
    fn ensure_fastga_in_path() -> Result<()> {
        // Try to find FastGA using our binary_paths module
        if let Ok(fastga_path) = crate::binary_paths::get_embedded_binary_path("FastGA") {
            eprintln!(
                "[FastGA] Found embedded binary at: {}",
                fastga_path.display()
            );
            if let Some(fastga_dir) = fastga_path.parent() {
                // Set ISOLATED PATH so FastGA can ONLY find its own utilities
                // This prevents FastGA from accidentally using system binaries
                eprintln!("[FastGA] Setting ISOLATED PATH to: {}", fastga_dir.display());
                std::env::set_var("PATH", fastga_dir.to_str().unwrap());
            }
        } else {
            eprintln!("[FastGA] WARNING: Could not find embedded FastGA binary");
        }
        Ok(())
    }

    pub fn prepare_gdb(&self, fasta_path: &Path) -> Result<String> {
        // Ensure fastga-rs can find binaries
        Self::ensure_fastga_in_path()?;

        // Create a temporary FastGA instance to access orchestrator methods
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.7),
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: std::env::var("TMPDIR").unwrap_or_else(|_| ".".to_string()),
        };

        // Prepare GDB
        let gdb_base = orchestrator
            .prepare_gdb(fasta_path)
            .map_err(|e| anyhow::anyhow!("Failed to prepare GDB: {}", e))?;

        // Create index
        orchestrator
            .create_index(
                &gdb_base,
                self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            )
            .map_err(|e| anyhow::anyhow!("Failed to create index: {}", e))?;

        Ok(gdb_base)
    }

    /// Run FastGA alignment and return .1aln file directly (NO PAF conversion)
    /// This is the native FastGA output format
    /// IMPORTANT: Also copies the .1gdb file alongside the .1aln to preserve sequence names
    pub fn align_to_temp_1aln(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Ensure fastga-rs can find binaries
        Self::ensure_fastga_in_path()?;

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0),
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: std::env::var("TMPDIR").unwrap_or_else(|_| ".".to_string()),
        };

        // Run alignment - FastGA creates BOTH .1aln and .1gdb files
        let aln_path = orchestrator
            .align_to_1aln(queries, targets)
            .map_err(|e| anyhow::anyhow!("Failed to run FastGA alignment: {}", e))?;

        // Derive .1gdb path from .1aln path
        let gdb_path = aln_path
            .strip_suffix(".1aln")
            .unwrap_or(&aln_path)
            .to_string()
            + ".1gdb";

        // Create temp file for .1aln (NamedTempFile ensures unique name)
        let temp_aln = NamedTempFile::new().context("Failed to create temp .1aln file")?;
        let temp_aln_path = temp_aln.path();

        // Derive temp .1gdb path (same base name as temp .1aln)
        let temp_gdb_path = temp_aln_path
            .to_str()
            .context("Invalid temp path")?
            .strip_suffix(".1aln")
            .or_else(|| temp_aln_path.to_str())
            .unwrap()
            .to_string()
            + ".1gdb";

        // Copy BOTH files to temp location
        std::fs::copy(&aln_path, temp_aln_path).context("Failed to copy .1aln to temp")?;

        if std::path::Path::new(&gdb_path).exists() {
            std::fs::copy(&gdb_path, &temp_gdb_path).context("Failed to copy .1gdb to temp")?;
            eprintln!(
                "[fastga] Preserved .1gdb with {} sequence names",
                std::fs::metadata(&temp_gdb_path)?.len()
            );
        } else {
            eprintln!("[fastga] WARNING: No .1gdb file found at {gdb_path}");
        }

        // Clean up originals
        let _ = std::fs::remove_file(&aln_path);
        let _ = std::fs::remove_file(&gdb_path);

        Ok(temp_aln)
    }

    /// Run FastGA alignment and write output to a temporary PAF file
    /// Assumes GDB/GIX indices already exist for the input files
    /// Returns the temporary file handle (which auto-deletes when dropped)
    pub fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Ensure fastga-rs can find binaries
        Self::ensure_fastga_in_path()?;

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0), // 0.0 means use FastGA default
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: std::env::var("TMPDIR").unwrap_or_else(|_| ".".to_string()),
        };

        // Run alignment with existing indices (returns PAF bytes directly)
        let paf_output = orchestrator
            .align_with_existing_indices(queries, targets)
            .map_err(|e| anyhow::anyhow!("Failed to run FastGA alignment: {}", e))?;

        eprintln!(
            "[sweepga] PAF output: {} bytes, {} lines",
            paf_output.len(),
            String::from_utf8_lossy(&paf_output).lines().count()
        );

        // Write to temporary file
        let mut temp_file = NamedTempFile::new().context("Failed to create temporary PAF file")?;

        temp_file
            .write_all(&paf_output)
            .context("Failed to write PAF output to temporary file")?;

        temp_file
            .flush()
            .context("Failed to flush temporary file")?;

        Ok(temp_file)
    }
}
