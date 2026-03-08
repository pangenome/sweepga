//! Wfmash integration module
//!
//! This module provides integration with wfmash for genome alignment.

use crate::aligner::Aligner;
use anyhow::Result;
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

/// Round a value to a "nice" multiple based on its magnitude:
/// ≤500 → multiple of 50, ≤1000 → 100, ≤3000 → 200, >3000 → 500.
fn round_nice(v: u64) -> u64 {
    if v == 0 {
        return 0;
    }
    let step = if v <= 500 {
        50
    } else if v <= 1000 {
        100
    } else if v <= 3000 {
        200
    } else {
        500
    };
    ((v + step / 2) / step * step).max(step)
}

/// Wfmash alignment backend
pub struct WfmashIntegration {
    wfmash: wfmash_rs::Wfmash,
}

impl WfmashIntegration {
    /// Create a new WfmashIntegration with adaptive parameters.
    ///
    /// `segment_length`: explicit wfmash segment length (`-s`). When `None`, adapts from `avg_seq_len`.
    /// `avg_seq_len`: used to compute adaptive segment length as `avg_seq_len / 2` (capped at 5000).
    /// `sparsify`: fraction of mappings to keep (wfmash `-x`). `None` or `Some(1.0)` = keep all.
    /// `num_mappings`: number of mappings per segment pair (`-n`). `None` = wfmash default (1).
    pub fn adaptive(
        num_threads: usize,
        min_alignment_length: Option<u64>,
        map_pct_identity: Option<String>,
        temp_dir: Option<String>,
        segment_length: Option<u64>,
        avg_seq_len: Option<u64>,
        sparsify: Option<f64>,
        num_mappings: Option<usize>,
    ) -> Result<Self> {
        let mut builder = wfmash_rs::Config::builder()
            .num_threads(num_threads)
            .prefix_delimiter('#');

        if let Some(n) = num_mappings {
            builder = builder.num_mappings(n);
        }

        if let Some(ref pct) = map_pct_identity {
            builder = builder.map_pct_identity(pct);
        }

        if let Some(ref dir) = temp_dir {
            builder = builder.temp_dir(PathBuf::from(dir));
        }

        // Set segment length (-s): use explicit value, or adapt to avg_seq_len/2, capped at 5000
        let effective_segment_length = segment_length.or_else(|| {
            avg_seq_len.map(|avg| round_nice((avg / 2).min(5000)))
        });
        if let Some(s) = effective_segment_length {
            builder = builder.segment_length(s as usize);
        }

        // Set block length (-l): use explicit value, or adapt to min(s*3, avg/2)
        // For short sequences this gives l == s
        let effective_block_length = min_alignment_length.or_else(|| {
            match (effective_segment_length, avg_seq_len) {
                (Some(s), Some(avg)) => Some(round_nice((s * 3).min(avg / 2))),
                (Some(s), None) => Some(round_nice(s * 3)),
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
