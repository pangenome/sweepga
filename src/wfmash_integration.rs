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
        Self::with_segment_length(
            num_threads,
            min_alignment_length,
            map_pct_identity,
            temp_dir,
            None,
        )
    }

    /// Create a new WfmashIntegration with adaptive parameters for input sequence sizes.
    ///
    /// `segment_length` sets the wfmash segment length (`-s`). When `None`,
    /// wfmash uses its default.
    pub fn with_segment_length(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
        segment_length: Option<u64>,
    ) -> Result<Self> {
        Self::adaptive(
            num_threads,
            min_alignment_length,
            map_pct_identity,
            temp_dir,
            segment_length,
            None,
            None,
        )
    }

    /// Create a new WfmashIntegration with full control over parameters.
    ///
    /// `segment_length`: explicit wfmash segment length (`-s`). When `None`, adapts from `avg_seq_len`.
    /// `avg_seq_len`: used to compute adaptive segment length as `avg_seq_len / 2` (capped at 5000).
    /// `sparsify`: fraction of mappings to keep (wfmash `-x`). `None` or `Some(1.0)` = keep all.
    pub fn adaptive(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
        segment_length: Option<u64>,
        avg_seq_len: Option<u64>,
        sparsify: Option<f64>,
    ) -> Result<Self> {
        let mut builder = wfmash_rs::Config::builder()
            .num_threads(num_threads);

        if let Some(ref pct) = map_pct_identity {
            builder = builder.map_pct_identity(pct);
        }

        if let Some(ref dir) = temp_dir {
            builder = builder.temp_dir(PathBuf::from(dir));
        }

        // Set segment length (-s): use explicit value, or adapt to avg_seq_len/2, capped at 5000
        let effective_segment_length = segment_length.or_else(|| {
            avg_seq_len.map(|avg| (avg / 2).min(5000))
        });
        if let Some(s) = effective_segment_length {
            builder = builder.segment_length(s as usize);
        }

        // Set block length (-l): use explicit value, or adapt to min(s*3, avg/2)
        // For short sequences this gives l == s
        let effective_block_length = min_alignment_length.or_else(|| {
            match (effective_segment_length, avg_seq_len) {
                (Some(s), Some(avg)) => Some((s * 3).min(avg / 2)),
                (Some(s), None) => Some(s * 3),
                _ => None,
            }
        });
        if let Some(l) = effective_block_length {
            builder = builder.block_length(l);
        }

        let mut extra_args = Vec::new();

        // Sparsify mappings: keep only this fraction (-x flag).
        if let Some(frac) = sparsify {
            if frac < 1.0 {
                extra_args.push(format!("-x{}", frac));
            }
        }

        if !extra_args.is_empty() {
            builder = builder.extra_args(extra_args);
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
