//! FastGA integration module
//!
//! This module provides integration with FastGA for genome alignment.

use anyhow::{Context, Result};
use fastga_rs::Config;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tempfile::NamedTempFile;

/// Error types specific to GIXmake index creation
#[derive(Debug)]
pub enum IndexCreationError {
    /// Index size likely exceeds GIXmake's internal limits (≥48MB sequence data)
    SizeLimitExceeded {
        batch_size_mb: u64,
        suggested_limit: u64,
        underlying_error: String
    },
    /// Other index creation failure
    Other(String),
}

impl std::fmt::Display for IndexCreationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IndexCreationError::SizeLimitExceeded { batch_size_mb, suggested_limit, underlying_error } => {
                write!(f,
                    "GIXmake index creation failed: batch size {}MB likely exceeds FastGA's index size limit. \
                     Try --batch-bytes {}M or smaller. Original error: {}",
                    batch_size_mb, suggested_limit, underlying_error
                )
            }
            IndexCreationError::Other(msg) => write!(f, "Index creation failed: {}", msg),
        }
    }
}

impl std::error::Error for IndexCreationError {}

/// Get the preferred temp directory for FastGA operations.
/// Priority: explicit override > TMPDIR env var > current directory.
///
/// The special value "ramdisk" maps to /dev/shm (Linux shared-memory tmpfs)
/// for users who want RAM-backed temp storage without knowing the path.
fn get_temp_dir(override_dir: Option<&str>) -> String {
    // First check explicit override (from --temp-dir CLI option)
    if let Some(dir) = override_dir {
        if !dir.is_empty() {
            // "ramdisk" is a convenience alias for /dev/shm
            if dir == "ramdisk" {
                let dev_shm = std::path::Path::new("/dev/shm");
                if dev_shm.is_dir() {
                    let test_path = dev_shm.join(format!(".sweepga_test_{}", std::process::id()));
                    if std::fs::write(&test_path, b"test").is_ok() {
                        let _ = std::fs::remove_file(&test_path);
                        return "/dev/shm".to_string();
                    }
                }
                log::info!("[sweepga] warning: --temp-dir ramdisk requested but /dev/shm is not available, falling back to current directory");
            } else if std::path::Path::new(dir).is_dir() {
                return dir.to_string();
            }
        }
    }

    // Then check TMPDIR environment variable
    if let Ok(tmpdir) = std::env::var("TMPDIR") {
        if !tmpdir.is_empty() && std::path::Path::new(&tmpdir).is_dir() {
            return tmpdir;
        }
    }

    // Default: current directory (avoids silently consuming RAM via /dev/shm)
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
    temp_dir: Option<String>,
}

impl FastGAIntegration {
    /// Create a new FastGA integration.
    pub fn new(
        frequency: usize,
        num_threads: usize,
        min_alignment_length: u64,
        temp_dir: Option<String>,
    ) -> Self {
        let config = Config::builder()
            .num_threads(num_threads)
            .min_alignment_length(min_alignment_length as usize)
            .verbose(true)
            .adaptive_seed_cutoff(frequency)
            .build();

        FastGAIntegration { config, temp_dir }
    }

    pub fn prepare_gdb(&self, fasta_path: &Path) -> Result<String> {
        // binary_paths::setup_binary_env() in main() already set up PATH

        // Set TMPDIR for FastGA's internal temp files (ONEcode uses this)
        if let Some(ref dir) = self.temp_dir {
            std::env::set_var("TMPDIR", dir);
        }

        // Create a temporary FastGA instance to access orchestrator methods
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.7),
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: get_temp_dir(self.temp_dir.as_deref()),
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

    /// Create GDB files only (no index) - used for batch alignment
    /// Returns the gdb_base path (without extension)
    pub fn create_gdb_only(&self, fasta_path: &Path) -> Result<String> {
        // binary_paths::setup_binary_env() in main() already set up PATH

        // Set TMPDIR for FastGA's internal temp files (ONEcode uses this)
        if let Some(ref dir) = self.temp_dir {
            std::env::set_var("TMPDIR", dir);
        }

        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.7),
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: get_temp_dir(self.temp_dir.as_deref()),
        };

        let gdb_base = orchestrator
            .prepare_gdb(fasta_path)
            .map_err(|e| anyhow::anyhow!("Failed to prepare GDB: {e}"))?;

        Ok(gdb_base)
    }

    /// Create index for an existing GDB - used for batch alignment
    pub fn create_index_only(&self, gdb_base: &str) -> Result<()> {
        // binary_paths::setup_binary_env() in main() already set up PATH

        // Set TMPDIR for FastGA's internal temp files (ONEcode uses this)
        if let Some(ref dir) = self.temp_dir {
            std::env::set_var("TMPDIR", dir);
        }

        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.7),
            kmer_freq: self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            temp_dir: get_temp_dir(self.temp_dir.as_deref()),
        };

        orchestrator
            .create_index(
                gdb_base,
                self.config.adaptive_seed_cutoff.unwrap_or(10) as i32,
            )
            .map_err(|e| {
                // Try to determine batch size from GDB file for better error messages
                let batch_size_mb = Self::estimate_batch_size_mb(gdb_base).unwrap_or(0);

                // Check if this looks like a GIXmake size limit failure
                if Self::is_likely_size_limit_error(&e.to_string(), batch_size_mb) {
                    anyhow::anyhow!(IndexCreationError::SizeLimitExceeded {
                        batch_size_mb,
                        suggested_limit: std::cmp::max(32, batch_size_mb * 3 / 4), // Suggest 25% reduction
                        underlying_error: e.to_string()
                    })
                } else {
                    anyhow::anyhow!(IndexCreationError::Other(e.to_string()))
                }
            })?;

        Ok(())
    }

    /// Estimate the batch size in MB from the GDB file size
    fn estimate_batch_size_mb(gdb_base: &str) -> Result<u64> {
        let gdb_path = format!("{}.1gdb", gdb_base);
        let metadata = std::fs::metadata(&gdb_path)?;
        // GDB file size is roughly proportional to sequence data size
        // Use a conservative multiplier to estimate original sequence size
        Ok((metadata.len() / (1024 * 1024)) * 3) // Rough approximation
    }

    /// Check if an error message indicates a likely GIXmake size limit failure
    fn is_likely_size_limit_error(error_msg: &str, batch_size_mb: u64) -> bool {
        // GIXmake size limit failures typically occur at ≥48MB sequence data
        if batch_size_mb >= 40 {
            // Look for typical patterns in GIXmake/FastGA failure messages
            error_msg.contains("GIXmake")
                || error_msg.contains("code None")
                || error_msg.contains("failed with code")
                || (error_msg.is_empty() && batch_size_mb >= 48) // Silent failures
        } else {
            false
        }
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

        log::info!("[FastGA] Compressing index with zstd (level {level})...");

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

        log::info!("[FastGA] Index compression complete");
        Ok(())
    }

    /// Run FastGA alignment with PAF output including CIGAR strings
    /// Uses .1aln intermediate format and ALNtoPAF -x for proper CIGAR generation
    pub fn align_direct_paf(&self, queries: &Path, targets: &Path) -> Result<Vec<u8>> {
        use std::process::Command;

        // Get binary paths
        let fastga_bin = std::env::var("FASTGA_BIN_DIR")
            .map(|dir| format!("{}/FastGA", dir))
            .unwrap_or_else(|_| "FastGA".to_string());
        let alnto_paf_bin = std::env::var("FASTGA_BIN_DIR")
            .map(|dir| format!("{}/ALNtoPAF", dir))
            .unwrap_or_else(|_| "ALNtoPAF".to_string());

        // Frequency is always set by `FastGAIntegration::new`
        let kmer_freq = self
            .config
            .adaptive_seed_cutoff
            .expect("FastGAIntegration::new always sets adaptive_seed_cutoff") as i32;

        // Use absolute paths to handle batch alignment where query and target
        // may be in different directories with the same filename
        let query_abs = queries
            .canonicalize()
            .unwrap_or_else(|_| queries.to_path_buf());
        let target_abs = targets
            .canonicalize()
            .unwrap_or_else(|_| targets.to_path_buf());

        // Set working directory to target's directory (where the index is)
        let working_dir = target_abs
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));

        // Create temporary .1aln file for FastGA output
        let temp_aln = working_dir.join(format!(
            "_tmp_{}_{}.1aln",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        let temp_aln_filename = temp_aln.file_name().unwrap();

        // Build FastGA command - output to .1aln format
        let mut cmd = Command::new(&fastga_bin);
        cmd.arg(format!("-1:{}", temp_aln_filename.to_string_lossy()))
            .arg(format!("-T{}", self.config.num_threads));

        if kmer_freq != 10 {
            cmd.arg(format!("-f{}", kmer_freq));
        }

        if self.config.min_alignment_length > 0 {
            cmd.arg(format!("-l{}", self.config.min_alignment_length));
        }

        if let Some(min_id) = self.config.min_identity {
            if min_id > 0.0 {
                cmd.arg(format!("-i{:.2}", min_id));
            }
        }

        // Use absolute paths so cross-directory batch alignments work correctly
        cmd.arg(&query_abs)
            .arg(&target_abs)
            .current_dir(working_dir);

        let output = cmd
            .output()
            .with_context(|| format!("Failed to run FastGA: {}", fastga_bin))?;

        if !output.status.success() {
            let _ = std::fs::remove_file(&temp_aln);
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!("FastGA failed: {}", stderr);
        }

        // Convert .1aln to PAF with CIGAR using ALNtoPAF -x
        let paf_output = Command::new(&alnto_paf_bin)
            .arg("-x") // Generate CIGAR with X/= operators
            .arg(temp_aln_filename)
            .current_dir(working_dir)
            .output();

        // Clean up temp .1aln file
        let _ = std::fs::remove_file(&temp_aln);

        match paf_output {
            Ok(output) if output.status.success() => Ok(output.stdout),
            Ok(output) => {
                // ALNtoPAF failed (possibly segfault) - fall back to direct PAF without CIGAR
                let stderr = String::from_utf8_lossy(&output.stderr);
                log::warn!(
                    "[FastGA] ALNtoPAF failed ({}); falling back to PAF without CIGAR",
                    stderr.trim()
                );

                // Re-run FastGA with direct PAF output (no CIGAR)
                let mut fallback_cmd = Command::new(&fastga_bin);
                fallback_cmd.arg(format!("-T{}", self.config.num_threads));
                if kmer_freq != 10 {
                    fallback_cmd.arg(format!("-f{}", kmer_freq));
                }
                if self.config.min_alignment_length > 0 {
                    fallback_cmd.arg(format!("-l{}", self.config.min_alignment_length));
                }
                if let Some(min_id) = self.config.min_identity {
                    if min_id > 0.0 {
                        fallback_cmd.arg(format!("-i{:.2}", min_id));
                    }
                }
                fallback_cmd
                    .arg(&query_abs)
                    .arg(&target_abs)
                    .current_dir(working_dir);

                let fallback_output = fallback_cmd
                    .output()
                    .with_context(|| "FastGA fallback failed")?;

                if !fallback_output.status.success() {
                    let stderr = String::from_utf8_lossy(&fallback_output.stderr);
                    anyhow::bail!("FastGA fallback failed: {}", stderr);
                }

                Ok(fallback_output.stdout)
            }
            Err(e) => {
                anyhow::bail!("Failed to run ALNtoPAF: {}", e);
            }
        }
    }

    /// Clean up GIX index files only (.gix, .ktab.*, .ktab.*.zst) for a given base path
    /// Keeps .1gdb and .bps files which are needed for sequence name lookups
    /// This frees most disk space while preserving ability to query
    pub fn cleanup_index(gdb_base: &str) -> Result<()> {
        let base_path = std::path::Path::new(gdb_base);
        let parent = base_path
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));
        let base_name = base_path.file_name().unwrap().to_str().unwrap();

        // Remove .gix file and associated .ktab.* files
        let gix_path = format!("{}.gix", gdb_base);
        if let Ok(gix_file) = std::fs::File::open(&gix_path) {
            use std::io::Read;
            let mut header = [0u8; 8];
            let mut gix_reader = std::io::BufReader::new(gix_file);
            if gix_reader.read_exact(&mut header).is_ok() {
                let nthreads = i32::from_le_bytes([header[4], header[5], header[6], header[7]]);
                // Remove .ktab.* and .ktab.*.zst files (the big ones)
                for p in 1..=nthreads {
                    let ktab_path = parent.join(format!(".{}.ktab.{}", base_name, p));
                    let ktab_zst_path = parent.join(format!(".{}.ktab.{}.zst", base_name, p));
                    let _ = std::fs::remove_file(&ktab_path);
                    let _ = std::fs::remove_file(&ktab_zst_path);
                }
            }
        }
        let _ = std::fs::remove_file(&gix_path);

        // NOTE: We keep .1gdb and .bps files - they're needed for sequence name resolution
        // when this batch is used as a QUERY (not target) later

        Ok(())
    }

    /// Completely remove all files (.gdb, .gix, .bps, .ktab.*) for a given base path
    pub fn cleanup_all(gdb_base: &str) -> Result<()> {
        let base_path = std::path::Path::new(gdb_base);
        let parent = base_path
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));
        let base_name = base_path.file_name().unwrap().to_str().unwrap();

        // First cleanup GIX
        Self::cleanup_index(gdb_base)?;

        // Then remove GDB files
        let _ = std::fs::remove_file(format!("{}.1gdb", gdb_base));
        let bps_path = parent.join(format!(".{}.bps", base_name));
        let _ = std::fs::remove_file(&bps_path);

        Ok(())
    }

    /// Run FastGA alignment and return .1aln file directly (NO PAF conversion)
    /// This is the native FastGA output format
    /// IMPORTANT: Also copies the .1gdb file alongside the .1aln to preserve sequence names
    pub fn align_to_temp_1aln(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Ensure fastga-rs can find binaries
        // binary_paths::setup_binary_env() in main() already set up PATH

        // Set TMPDIR for FastGA's internal temp files (ONEcode uses this)
        if let Some(ref dir) = self.temp_dir {
            std::env::set_var("TMPDIR", dir);
        }

        // Frequency is always set by `FastGAIntegration::new`
        let kmer_freq = self
            .config
            .adaptive_seed_cutoff
            .expect("FastGAIntegration::new always sets adaptive_seed_cutoff") as i32;

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0),
            kmer_freq,
            temp_dir: get_temp_dir(self.temp_dir.as_deref()),
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
            // log::info!(
            //     "[fastga] Preserved .1gdb with {} sequence names",
            //     std::fs::metadata(&temp_gdb_path)?.len()
            // );
        } else {
            // log::info!("[fastga] WARNING: No .1gdb file found at {gdb_path}");
        }

        // Clean up originals
        let _ = std::fs::remove_file(&aln_path);
        let _ = std::fs::remove_file(&gdb_path);

        Ok(temp_aln)
    }

    /// Run FastGA alignment and write output to a temporary PAF file
    /// Pre-builds GDB/GIX indices before invoking FastGA;
    /// Returns the temporary file handle (which auto-deletes when dropped)
    pub fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        // Ensure fastga-rs can find binaries
        // binary_paths::setup_binary_env() in main() already set up PATH

        // Pre-build GDB + GIX up front so FastGA doesn't have to build
        // them internally on-the-fly (which crashes with a buffer
        // overflow on small inputs in `batch_align.rs` for the working pattern).
        let q_gdb = self.create_gdb_only(queries)?;
        self.create_index_only(&q_gdb)?;
        let t_gdb = if targets == queries {
            None
        } else {
            let t = self.create_gdb_only(targets)?;
            self.create_index_only(&t)?;
            Some(t)
        };

        // Set TMPDIR for FastGA's internal temp files (ONEcode uses this)
        if let Some(ref dir) = self.temp_dir {
            std::env::set_var("TMPDIR", dir);
        }

        // Frequency is always set by `FastGAIntegration::new`
        let kmer_freq = self
            .config
            .adaptive_seed_cutoff
            .expect("FastGAIntegration::new always sets adaptive_seed_cutoff") as i32;

        // Create orchestrator to run FastGA binary directly
        let orchestrator = fastga_rs::orchestrator::FastGAOrchestrator {
            num_threads: self.config.num_threads as i32,
            min_length: self.config.min_alignment_length as i32,
            min_identity: self.config.min_identity.unwrap_or(0.0), // 0.0 means use FastGA default
            kmer_freq,
            temp_dir: get_temp_dir(self.temp_dir.as_deref()),
        };

        // Get the directory containing the input files (where indexes will be created)
        // Use canonicalize to get absolute path, falling back to current dir
        let index_dir = queries
            .parent()
            .and_then(|p| {
                if p.as_os_str().is_empty() {
                    None
                } else {
                    Some(p)
                }
            })
            .map(|p| p.to_path_buf())
            .or_else(|| std::env::current_dir().ok())
            .unwrap_or_else(|| Path::new(".").to_path_buf());

        // Start a background thread to monitor index file creation
        let stop_flag = Arc::new(AtomicBool::new(false));
        let stop_flag_clone = stop_flag.clone();
        let index_dir_clone = index_dir.clone();

        let monitor_handle = std::thread::spawn(move || {
            let mut last_size = 0u64;
            while !stop_flag_clone.load(Ordering::Relaxed) {
                // Scan for index files and track the delta
                if let Ok(current_size) =
                    crate::disk_usage::scan_fastga_index_files(&index_dir_clone)
                {
                    if current_size != last_size {
                        // Record the delta (positive = new files, negative = files deleted)
                        if current_size > last_size {
                            crate::disk_usage::add_bytes(current_size - last_size);
                        } else {
                            crate::disk_usage::remove_bytes(last_size - current_size);
                        }
                        last_size = current_size;
                    }
                }
                std::thread::sleep(std::time::Duration::from_millis(1000));
            }
            // Final cleanup - indexes are deleted after alignment
            if last_size > 0 {
                crate::disk_usage::remove_bytes(last_size);
            }
        });

        // Run alignment with existing indices (returns PAF bytes directly)
        let paf_output = orchestrator
            .align_with_existing_indices(queries, targets)
            .map_err(|e| anyhow::anyhow!("Failed to run FastGA alignment: {e}"))?;

        // Stop the monitor thread (it will cleanup the index bytes on exit)
        stop_flag.store(true, Ordering::Relaxed);
        let _ = monitor_handle.join();

        // Best-effort cleanup of the GDB/GIX artifacts we created upfront.
        let _ = Self::cleanup_all(&q_gdb);
        if let Some(t) = t_gdb {
            let _ = Self::cleanup_all(&t);
        }

        // Write to temporary file
        let mut temp_file = NamedTempFile::new().context("Failed to create temporary PAF file")?;

        temp_file
            .write_all(&paf_output)
            .context("Failed to write PAF output to temporary file")?;

        temp_file
            .flush()
            .context("Failed to flush temporary file")?;

        // Track disk usage of the temp PAF file
        crate::disk_usage::track_file_created(temp_file.path());

        Ok(temp_file)
    }
}

impl crate::aligner::Aligner for FastGAIntegration {
    fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        self.align_to_temp_paf(queries, targets)
    }

    fn align_direct_paf(&self, queries: &Path, targets: &Path) -> Result<Vec<u8>> {
        self.align_direct_paf(queries, targets)
    }

    fn align_to_temp_1aln(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        self.align_to_temp_1aln(queries, targets)
    }

    fn name(&self) -> &str {
        "fastga"
    }
}
