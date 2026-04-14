// Library exports for sweepga
pub mod agc;
pub mod aligner;
pub mod batch_align;
pub mod binary_paths;
pub mod compact_mapping;
pub mod det_map;
pub mod disk_usage;
pub mod fastga_integration;
pub mod filter_types;
pub mod grouped_mappings;
pub mod knn_graph;
pub mod mapping;
pub mod mash;
pub mod orchestrator;
pub mod paf;
pub mod paf_filter;
pub mod plane_sweep_core;
pub mod plane_sweep_exact;
pub mod plane_sweep_scaffold;
pub mod seq_registry;
pub mod sequence_index;
pub mod unified_filter;
pub mod union_find;
pub mod wfmash_integration;

use anyhow::Result;
use std::path::Path;

/// Self-alignment with optional batching.
///
/// If `batch_bytes` is `None`, uses `aligner.align_to_temp_paf()` directly.
/// If `batch_bytes` is `Some("2G")` etc., partitions genomes into batches and
/// aligns each batch pair, limiting per-batch resource usage.
///
/// The `aligner` parameter is the already-constructed aligner (used when not batching).
/// When batching, the function creates the appropriate `BatchAligner` from the
/// provided parameters.
pub fn align_self_paf(
    fasta_path: &Path,
    aligner: &dyn aligner::Aligner,
    aligner_name: &str,
    kmer_frequency: usize,
    num_threads: usize,
    min_aln_length: u64,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
    batch_bytes: Option<&str>,
    quiet: bool,
) -> Result<tempfile::NamedTempFile> {
    match batch_bytes {
        None => {
            // Direct alignment, no batching
            aligner.align_to_temp_paf(fasta_path, fasta_path)
        }
        Some(batch_str) => {
            let max_bytes = batch_align::parse_size_string(batch_str)?;

            let batch_aligner: Box<dyn batch_align::BatchAligner> = match aligner_name {
                "wfmash" => Box::new(batch_align::WfmashBatchAligner::new(
                    num_threads,
                    if min_aln_length > 0 { Some(min_aln_length) } else { None },
                    map_pct_identity,
                    temp_dir.clone(),
                )),
                _ => Box::new(batch_align::FastGABatchAligner::new(
                    Some(kmer_frequency),
                    num_threads,
                    min_aln_length,
                    temp_dir.clone(),
                    false, // zstd_compress
                    3,     // zstd_level
                )),
            };

            let batch_config = batch_align::BatchAlignConfig {
                frequency: Some(kmer_frequency),
                threads: num_threads,
                min_alignment_length: min_aln_length,
                zstd_compress: false,
                zstd_level: 3,
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
    }
}
