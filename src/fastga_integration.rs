//! FastGA integration module
//!
//! This module provides integration with FastGA for genome alignment.

use anyhow::{Context, Result};
use fastga_rs::Config;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use tempfile::NamedTempFile;

/// Get the preferred temp directory for FastGA operations.
/// Priority: TMPDIR env var > /dev/shm (if available) > current directory
fn get_temp_dir() -> String {
    // First check TMPDIR environment variable
    if let Ok(tmpdir) = std::env::var("TMPDIR") {
        if !tmpdir.is_empty() && std::path::Path::new(&tmpdir).is_dir() {
            return tmpdir;
        }
    }

    // Check if /dev/shm exists and is writable (preferred for performance)
    let dev_shm = std::path::Path::new("/dev/shm");
    if dev_shm.is_dir() {
        // Try to create a test file to verify write access
        let test_path = dev_shm.join(format!(".sweepga_test_{}", std::process::id()));
        if std::fs::write(&test_path, b"test").is_ok() {
            let _ = std::fs::remove_file(&test_path);
            return "/dev/shm".to_string();
        }
    }

    // Fallback to current directory
    ".".to_string()
}

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
                "Unknown preset '{s}'. Valid presets: ani70, ani80, ani85, ani90, ani95, ani99"
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
    /// If not provided, auto-detects from PanSN haplotype count during alignment
    pub fn new(frequency: Option<usize>, num_threads: usize, min_alignment_length: u64) -> Self {
        let mut builder = Config::builder()
            .num_threads(num_threads)
            .min_alignment_length(min_alignment_length as usize)
            .verbose(true);

        // Only set frequency if explicitly provided
        // Otherwise, alignment methods will auto-detect from PanSN haplotype count
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
            // eprintln!(
            //     "[FastGA] Found embedded binary at: {}",
            //     fastga_path.display()
            // );
            if let Some(fastga_dir) = fastga_path.parent() {
                // Set ISOLATED PATH so FastGA can ONLY find its own utilities
                // This prevents FastGA from accidentally using system binaries
                // eprintln!(
                //     "[FastGA] Setting ISOLATED PATH to: {}",
                //     fastga_dir.display()
                // );
                std::env::set_var("PATH", fastga_dir.to_str().unwrap());
            }
        } else {
            // eprintln!("[FastGA] WARNING: Could not find embedded FastGA binary");
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
            temp_dir: get_temp_dir(),
        };

        // Prepare GDB
        let gdb_base = orchestrator
            .prepare_gdb(fasta_path)
            .map_err(|e| anyhow::anyhow!("Failed to prepare GDB: {e}"))?;

        // Create index
        orchestrator
            .create_index(
                &gdb_base,
                self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            )
            .map_err(|e| anyhow::anyhow!("Failed to create index: {e}"))?;

        Ok(gdb_base)
    }

    /// Compress ktab index files using zstd seekable format
    /// This reduces disk usage by ~2x and can improve I/O performance
    /// Level: 1-19 (higher = smaller files but slower compression)
    pub fn compress_index(gdb_base: &str, level: u32) -> Result<()> {
        use std::process::Command;

        // Get the bin directory from FASTGA_BIN_DIR or PATH
        let gixpack_path = std::env::var("FASTGA_BIN_DIR")
            .map(|dir| format!("{}/GIXpack", dir))
            .unwrap_or_else(|_| "GIXpack".to_string());

        eprintln!("[FastGA] Compressing index with zstd (level {level})...");

        let output = Command::new(&gixpack_path)
            .arg(format!("-l{}", level))
            .arg(gdb_base)
            .output()
            .with_context(|| format!("Failed to run GIXpack: {}", gixpack_path))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!("GIXpack failed: {}", stderr);
        }

        // Remove uncompressed ktab files after successful compression
        let gix_path = format!("{}.gix", gdb_base);
        if let Ok(gix_file) = std::fs::File::open(&gix_path) {
            use std::io::Read;
            let mut header = [0u8; 8];
            let mut gix_reader = std::io::BufReader::new(gix_file);
            if gix_reader.read_exact(&mut header).is_ok() {
                let nthreads = i32::from_le_bytes([header[4], header[5], header[6], header[7]]);
                for p in 1..=nthreads {
                    let ktab_path = format!(
                        ".{}.ktab.{}",
                        std::path::Path::new(gdb_base)
                            .file_name()
                            .unwrap()
                            .to_str()
                            .unwrap(),
                        p
                    );
                    let full_path = std::path::Path::new(gdb_base)
                        .parent()
                        .map(|p| p.join(&ktab_path))
                        .unwrap_or_else(|| std::path::PathBuf::from(&ktab_path));
                    let _ = std::fs::remove_file(full_path);
                }
            }
        }

        eprintln!("[FastGA] Index compression complete");
        Ok(())
    }

    /// Count unique haplotypes in a FASTA file based on PanSN naming convention
    /// PanSN format: SAMPLE#HAPLOTYPE#CONTIG (e.g., HG01106#1#CM087962.1)
    /// Extracts SAMPLE#HAPLOTYPE prefix and counts unique haplotypes
    fn count_haplotypes(fasta_path: &Path) -> Result<usize> {
        use std::collections::HashSet;
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(fasta_path)?;

        // Handle both plain and gzipped FASTA files
        let reader: Box<dyn BufRead> =
            if fasta_path.extension().and_then(|s| s.to_str()) == Some("gz") {
                use flate2::read::MultiGzDecoder;
                Box::new(BufReader::new(MultiGzDecoder::new(file)))
            } else {
                Box::new(BufReader::new(file))
            };

        let mut haplotypes = HashSet::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Strip '>' and extract SAMPLE#HAPLOTYPE
                let header = line.trim_start_matches('>');

                // Split by '#' and take first two components
                let parts: Vec<&str> = header.split('#').collect();
                if parts.len() >= 2 {
                    // PanSN format: SAMPLE#HAPLOTYPE#...
                    let haplotype_id = format!("{}#{}", parts[0], parts[1]);
                    haplotypes.insert(haplotype_id);
                } else {
                    // Fallback: if not PanSN format, treat whole header as unique
                    haplotypes.insert(header.to_string());
                }
            }
        }

        Ok(haplotypes.len())
    }

    /// Run FastGA alignment and return .1aln file directly (NO PAF conversion)
    /// This is the native FastGA output format
    /// IMPORTANT: Also copies the .1gdb file alongside the .1aln to preserve sequence names
    pub fn align_to_temp_1aln(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Ensure fastga-rs can find binaries
        Self::ensure_fastga_in_path()?;

        // Determine kmer frequency threshold
        // Default to number of haplotypes in query file for pangenome workflows
        let kmer_freq = if let Some(freq) = self.config.adaptive_seed_cutoff {
            freq as i32
        } else {
            let num_haplotypes = Self::count_haplotypes(queries)?;
            // Use at least 10 (FastGA's default) to avoid over-filtering with small datasets
            std::cmp::max(num_haplotypes, 10) as i32
        };

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0),
            kmer_freq,
            temp_dir: get_temp_dir(),
        };

        // Run alignment - FastGA creates BOTH .1aln and .1gdb files
        let aln_path = orchestrator
            .align_to_1aln(queries, targets)
            .map_err(|e| anyhow::anyhow!("Failed to run FastGA alignment: {e}"))?;

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
            // eprintln!(
            //     "[fastga] Preserved .1gdb with {} sequence names",
            //     std::fs::metadata(&temp_gdb_path)?.len()
            // );
        } else {
            // eprintln!("[fastga] WARNING: No .1gdb file found at {gdb_path}");
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

        // Determine kmer frequency threshold
        // Default to number of haplotypes in query file for pangenome workflows
        let kmer_freq = if let Some(freq) = self.config.adaptive_seed_cutoff {
            freq as i32
        } else {
            let num_haplotypes = Self::count_haplotypes(queries)?;
            // Use at least 10 (FastGA's default) to avoid over-filtering with small datasets
            std::cmp::max(num_haplotypes, 10) as i32
        };

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0), // 0.0 means use FastGA default
            kmer_freq,
            temp_dir: get_temp_dir(),
        };

        // Run alignment with existing indices (returns PAF bytes directly)
        let paf_output = orchestrator
            .align_with_existing_indices(queries, targets)
            .map_err(|e| anyhow::anyhow!("Failed to run FastGA alignment: {e}"))?;

        // eprintln!(
        //     "[sweepga] PAF output: {} bytes, {} lines",
        //     paf_output.len(),
        //     String::from_utf8_lossy(&paf_output).lines().count()
        // );

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
