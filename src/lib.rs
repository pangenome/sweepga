// Library exports for sweepga
pub mod agc;
pub mod aligner;
pub mod batch_align;
pub mod binary_paths;
pub mod cli;
pub mod compact_mapping;
pub mod disk_usage;
pub mod fastga_integration;
pub mod filter_types;
pub mod grouped_mappings;
pub mod joblist;
pub mod knn_graph;
pub mod library_api;
pub mod mapping;
pub mod mash;
pub mod orchestrator;
pub mod paf;
pub mod pansn;
pub mod paf_filter;
pub mod plane_sweep_core;
pub mod plane_sweep_exact;
pub mod plane_sweep_scaffold;
pub mod seq_registry;
pub mod sequence_index;
pub mod unified_filter;
pub mod union_find;
pub mod wfmash_integration;

pub use cli::{parse_identity_value, parse_metric_number, AlnArgs};

use anyhow::Result;
use std::path::{Path, PathBuf};

/// Direct (non-batched) self-alignment: hands `fasta_path` to the
/// already-constructed `aligner` and returns its raw PAF temp file.
pub fn align_self_paf_direct(
    aligner: &dyn aligner::Aligner,
    fasta_path: &Path,
) -> Result<tempfile::NamedTempFile> {
    aligner.align_to_temp_paf(fasta_path, fasta_path)
}

/// Batched self-alignment: partitions genomes within `fasta_path` into
/// batches sized by `batch_bytes` (e.g. `"2G"`) and aligns each batch
/// pair, limiting per-batch resource usage.
///
/// `pairs_file`: optional path to a `query\ttarget` TSV (as produced by
/// pair-level sparsification upstream). When set AND the aligner is
/// `wfmash`, every per-batch wfmash invocation uses `--pairs-file` to
/// restrict to those pairs — so batching and pair sparsification
/// compose correctly. Ignored for FastGA (which doesn't support
/// pair-file filtering).
pub fn align_self_paf_batched(
    fasta_path: &Path,
    aligner_name: &str,
    kmer_frequency: usize,
    num_threads: usize,
    min_aln_length: u64,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
    wfmash_density: Option<f64>,
    pairs_file: Option<PathBuf>,
    batch_bytes: &str,
    quiet: bool,
) -> Result<tempfile::NamedTempFile> {
    let max_bytes = crate::cli::parse_metric_number(batch_bytes)
        .map_err(|e| anyhow::anyhow!("Invalid --batch-bytes: {e}"))?;

    let batch_aligner: Box<dyn batch_align::BatchAligner> = match aligner_name {
        "wfmash" => Box::new(batch_align::WfmashBatchAligner::new(
            num_threads,
            if min_aln_length > 0 { Some(min_aln_length) } else { None },
            map_pct_identity,
            temp_dir.clone(),
            wfmash_density,
            pairs_file,
        )),
        _ => Box::new(batch_align::FastGABatchAligner::new(
            kmer_frequency,
            num_threads,
            min_aln_length,
            temp_dir.clone(),
            false, // zstd_compress
            3,     // zstd_level
        )),
    };

    let batch_config = batch_align::BatchAlignConfig {
        keep_self: true,
        quiet,
    };

    let fasta_path_str = fasta_path.to_string_lossy().to_string();
    batch_align::run_batch_alignment_generic(
        &[fasta_path_str],
        max_bytes,
        batch_aligner.as_ref(),
        &batch_config,
        temp_dir.as_deref(),
    )
}
