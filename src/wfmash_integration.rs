//! Wfmash integration module
//!
//! This module provides integration with wfmash for genome alignment.

use crate::aligner::Aligner;
use anyhow::Result;
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

/// Wfmash alignment backend
pub struct WfmashIntegration {
    wfmash: wfmash_rs::Wfmash,
}

impl WfmashIntegration {
    /// Create a new WfmashIntegration instance.
    ///
    /// `map_pct_identity` sets the minimum percent identity for mapping (-p).
    /// Accepts a float (e.g. "90") or an ANI preset (e.g. "ani50-2").
    /// When `None`, wfmash uses its default (~70%).
    pub fn new(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
    ) -> Result<Self> {
        Self::with_segment_length(num_threads, min_alignment_length, map_pct_identity, temp_dir, None)
    }

    /// Create a new WfmashIntegration with explicit segment (window) length.
    ///
    /// `segment_length` sets the wfmash `-w` parameter. When `None`, wfmash
    /// uses its default (1000). Use a smaller value (e.g. 500) when input
    /// sequences are shorter than the default window size.
    pub fn with_segment_length(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
        segment_length: Option<u64>,
    ) -> Result<Self> {
        let mut builder = wfmash_rs::Config::builder()
            .num_threads(num_threads)
            .no_filter(true);

        if let Some(len) = min_alignment_length {
            builder = builder.block_length(len);
        }

        if let Some(ref pct) = map_pct_identity {
            builder = builder.map_pct_identity(pct);
        }

        if let Some(ref dir) = temp_dir {
            builder = builder.temp_dir(PathBuf::from(dir));
        }

        if let Some(w) = segment_length {
            builder = builder.window_size(&w.to_string());
            // Also scale scaffold mass (default 10k) proportionally.
            // wfmash's -S flag filters out scaffolds shorter than this.
            // For short sequences, the default 10k is too large.
            // Use segment_length as scaffold mass (each mapping covers ~w bp).
            let scaffold_mass = w.clamp(100, 10000);
            builder = builder.extra_args(vec![
                format!("-S{}", scaffold_mass),
            ]);
        }

        let config = builder.build();
        let wfmash = wfmash_rs::Wfmash::new(config)
            .map_err(|e| anyhow::anyhow!("Failed to create wfmash instance: {e}"))?;

        Ok(WfmashIntegration { wfmash })
    }
}

impl Aligner for WfmashIntegration {
    fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile> {
        if queries == targets {
            self.wfmash
                .align_self_to_temp_paf(queries)
                .map_err(|e| anyhow::anyhow!("wfmash self-alignment failed: {e}"))
        } else {
            // wfmash-rs uses (target, query) argument order
            self.wfmash
                .align_to_temp_paf(targets, queries)
                .map_err(|e| anyhow::anyhow!("wfmash alignment failed: {e}"))
        }
    }

    fn align_direct_paf(&self, queries: &Path, targets: &Path) -> Result<Vec<u8>> {
        if queries == targets {
            self.wfmash
                .align_self(queries)
                .map_err(|e| anyhow::anyhow!("wfmash self-alignment failed: {e}"))
        } else {
            // wfmash-rs uses (target, query) argument order
            self.wfmash
                .align_files(targets, queries)
                .map_err(|e| anyhow::anyhow!("wfmash alignment failed: {e}"))
        }
    }

    fn name(&self) -> &str {
        "wfmash"
    }
}
