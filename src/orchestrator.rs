//! High-level alignment orchestration combining pair selection, batching, and alignment.
//!
//! This module provides a unified entry point that resolves sparsification strategies,
//! selects pairs, partitions into batches when needed, runs alignment, and merges output.

use crate::knn_graph::SparsificationStrategy;
#[allow(unused_imports)]
use crate::paf_filter::FilterConfig;

/// Configuration for the alignment orchestrator.
pub struct OrchestrationConfig {
    /// Pair selection / mapping density strategy.
    pub strategy: SparsificationStrategy,
    /// Maximum sequence data per batch in bytes.
    pub batch_max_bytes: Option<u64>,
    /// Aligner backend name: "wfmash" or "fastga".
    pub aligner_name: String,
    /// K-mer frequency for fastga.
    pub kmer_frequency: usize,
    /// Number of threads.
    pub num_threads: usize,
    /// Minimum alignment length.
    pub min_alignment_length: u64,
    /// Temp directory override.
    pub temp_dir: Option<String>,
    /// Optional PAF filter config. When None, no filtering is applied.
    pub filter_config: Option<FilterConfig>,
    // wfmash-specific
    /// Mapping percent identity for wfmash (e.g. "90").
    pub map_pct_identity: Option<String>,
    /// Average sequence length (for adaptive wfmash parameters).
    pub avg_seq_len: Option<u64>,
    /// Wfmash segment length override.
    pub segment_length: Option<u64>,
}

/// Result of an orchestrated alignment run.
#[allow(dead_code)]
pub struct OrchestrationResult {
    /// Temp file containing merged PAF output.
    pub paf_file: tempfile::NamedTempFile,
    /// Number of sequence pairs that were aligned.
    pub num_pairs_aligned: usize,
    /// Total number of alignment records in the PAF.
    pub num_alignments: usize,
    /// Number of batches used (1 if no batching was needed).
    pub batches_used: usize,
}

/// Resolve the wfmash sparsify fraction from a `SparsificationStrategy`.
///
/// Returns `Some(fraction)` for `WfmashDensity` variants, `None` otherwise.
/// For `WfmashDensity(None)` (auto), computes `ln(n)/n*10` using `n_haps`.
pub fn resolve_wfmash_density(
    strategy: &SparsificationStrategy,
    n_haps: usize,
) -> Option<f64> {
    match strategy {
        SparsificationStrategy::WfmashDensity(Some(frac)) => Some(*frac),
        SparsificationStrategy::WfmashDensity(None) => {
            SparsificationStrategy::wfmash_auto_density(n_haps)
        }
        _ => None,
    }
}

/// Plan pair selection without executing alignment.
///
/// Returns `(pair_indices, pair_count)` where pair_indices are `(i, j)` indices
/// into the sequence list. Useful for joblist output.
pub fn plan_pairs(
    n_sequences: usize,
    strategy: &SparsificationStrategy,
) -> Vec<(usize, usize)> {
    crate::knn_graph::select_pairs(
        n_sequences,
        None,
        strategy,
        &crate::knn_graph::MashParams::default(),
    )
}

/// Plan pair selection using pre-computed mash sketches.
pub fn plan_pairs_from_sketches(
    sketches: &[crate::mash::KmerSketch],
    strategy: &SparsificationStrategy,
) -> Vec<(usize, usize)> {
    crate::knn_graph::select_pairs_from_sketches(sketches, strategy)
}

/// Validate that a strategy is compatible with the chosen aligner.
///
/// Returns an error message if the combination is invalid.
pub fn validate_strategy_aligner(
    strategy: &SparsificationStrategy,
    aligner_name: &str,
) -> Result<(), String> {
    if matches!(strategy, SparsificationStrategy::WfmashDensity(_)) && aligner_name != "wfmash" {
        return Err(format!(
            "Wfmash density sparsification ({}) requires --aligner wfmash, but '{}' was specified",
            strategy, aligner_name
        ));
    }
    Ok(())
}
