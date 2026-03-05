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
