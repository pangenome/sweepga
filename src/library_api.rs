//! High-level alignment orchestration API consumed by external crates
//! (originally by impg). Runs sweepga against in-memory named sequences,
//! optionally restricted to a pair list produced by sparsification
//! (`SparsificationStrategy`), then applies sweepga's PAF filter.
//!
//! Exists because downstream consumers often build their sequence set
//! in memory (e.g. pangenome-window excerpts) and can't hand a FASTA
//! path to `sweepga::align_self_paf`. This module accepts
//! `&[(name, &[u8])]` slices and does the FASTA/pairs-file marshaling.
//!
//! Public surface:
//! - [`SweepgaAlignConfig`] — flags forwarded to the alignment call.
//! - [`sweepga_align`] — dispatches between all-vs-all and pairwise modes.
//! - [`generate_pairs_for_sequences`] — PanSN-aware pair selection helper.
//! - [`apply_paf_filter`] — runs sweepga's `PafFilter` on a PAF file.
//! - [`create_aligner_adaptive`] — factory for the boxed `Aligner` trait
//!   object; picks FastGA or wfmash based on the `aligner_name` string.

use anyhow::Result;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::aligner::Aligner;
use crate::fastga_integration::FastGAIntegration;
use crate::knn_graph::SparsificationStrategy;
use crate::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
use crate::wfmash_integration::WfmashIntegration;

/// Parse filter-mode string like `"1:1"`, `"many:many"`, `"5:3"` into
/// `(FilterMode, max_per_query, max_per_target)`.
pub fn parse_filter_mode(s: &str) -> (FilterMode, Option<usize>, Option<usize>) {
    let s_lower = s.to_lowercase();
    if s_lower == "many:many" || s_lower == "n:n" {
        return (FilterMode::ManyToMany, None, None);
    }

    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return (FilterMode::OneToOne, Some(1), Some(1));
    }

    let query_max = if parts[0] == "many" || parts[0] == "n" {
        None
    } else {
        parts[0].parse().ok()
    };

    let target_max = if parts[1] == "many" || parts[1] == "n" {
        None
    } else {
        parts[1].parse().ok()
    };

    match (query_max, target_max) {
        (Some(1), Some(1)) => (FilterMode::OneToOne, Some(1), Some(1)),
        (Some(1), _) => (FilterMode::OneToMany, Some(1), target_max),
        (_, Some(1)) => (FilterMode::OneToMany, query_max, Some(1)),
        (Some(q), Some(t)) => (FilterMode::ManyToMany, Some(q), Some(t)),
        (Some(q), None) => (FilterMode::ManyToMany, Some(q), None),
        (None, Some(t)) => (FilterMode::ManyToMany, None, Some(t)),
        (None, None) => (FilterMode::ManyToMany, None, None),
    }
}

/// Construct a boxed aligner (FastGA or wfmash) with adaptive-by-avg-seq-len
/// parameters suitable for short-sequence workflows (pangenome excerpts).
#[allow(clippy::too_many_arguments)]
pub fn create_aligner_adaptive(
    aligner_name: &str,
    kmer_frequency: usize,
    num_threads: usize,
    min_aln_length: u64,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
    segment_length: Option<u64>,
    avg_seq_len: Option<u64>,
    sparsify: Option<f64>,
    num_mappings: Option<usize>,
    pairs_file: Option<PathBuf>,
) -> Result<Box<dyn Aligner>> {
    match aligner_name {
        "wfmash" => {
            let block_len = if min_aln_length > 0 {
                Some(min_aln_length)
            } else {
                None
            };
            let wfmash = WfmashIntegration::adaptive(
                num_threads,
                block_len,
                map_pct_identity,
                temp_dir,
                segment_length,
                avg_seq_len,
                sparsify,
                num_mappings,
                pairs_file,
            )
            .map_err(|e| anyhow::anyhow!("Failed to create wfmash aligner: {e}"))?;
            Ok(Box::new(wfmash))
        }
        "fastga" => Ok(Box::new(FastGAIntegration::new(
            kmer_frequency,
            num_threads,
            min_aln_length,
            temp_dir,
        ))),
        other => anyhow::bail!(
            "Unknown aligner: {other}. Valid options: wfmash, fastga"
        ),
    }
}

/// Configuration for in-memory sweepga alignment.
pub struct SweepgaAlignConfig {
    /// Number of threads for alignment.
    pub num_threads: usize,
    /// K-mer frequency for FastGA (pre-resolved via
    /// `pansn::resolve_fastga_frequency`).
    pub kmer_frequency: usize,
    /// Minimum alignment length.
    pub min_aln_length: u64,
    /// Skip filtering entirely.
    pub no_filter: bool,
    /// Filter: n:m-best mappings.
    pub num_mappings: String,
    /// Filter: scaffold jump distance.
    pub scaffold_jump: u64,
    /// Filter: minimum scaffold chain length.
    pub scaffold_mass: u64,
    /// Filter: scaffold filter mode.
    pub scaffold_filter: String,
    /// Filter: max overlap ratio.
    pub overlap: f64,
    /// Filter: minimum identity fraction (0..=1).
    pub min_identity: f64,
    /// Filter: max scaffold deviation.
    pub scaffold_dist: u64,
    /// Filter: minimum mapping length.
    pub min_map_length: u64,
    /// Optional temp directory.
    pub temp_dir: Option<String>,
    /// Unified sparsification strategy (pair selection + mapping density).
    pub sparsify: SparsificationStrategy,
    /// Mash distance parameters for sparsification sketching.
    pub mash_params: crate::knn_graph::MashParams,
    /// Aligner backend: `"wfmash"` or `"fastga"`.
    pub aligner: String,
    /// Minimum mapping identity for wfmash (e.g. `"70"`). `None` = wfmash auto-estimates.
    pub map_pct_identity: Option<String>,
    /// Batch alignment: max resource usage per batch (e.g. `"2G"`).
    pub batch_bytes: Option<String>,
}

impl Default for SweepgaAlignConfig {
    fn default() -> Self {
        SweepgaAlignConfig {
            num_threads: 4,
            kmer_frequency: 10,
            min_aln_length: 0,
            no_filter: false,
            num_mappings: "many:many".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "many:many".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_map_length: 0,
            temp_dir: None,
            sparsify: SparsificationStrategy::None,
            mash_params: crate::knn_graph::MashParams::default(),
            aligner: "wfmash".to_string(),
            map_pct_identity: None,
            batch_bytes: None,
        }
    }
}

/// Generate alignment pairs from named in-memory sequences. For simple
/// strategies (None/Random/WfmashDensity) works from the count alone;
/// other strategies (tree/giant/connectivity) compute mash sketches.
pub fn generate_pairs_for_sequences(
    sequences: &[(String, &[u8])],
    strategy: &SparsificationStrategy,
    mash_params: &crate::knn_graph::MashParams,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n <= 1 {
        return vec![];
    }

    match strategy {
        SparsificationStrategy::None
        | SparsificationStrategy::Random(_)
        | SparsificationStrategy::WfmashDensity(_) => {
            crate::knn_graph::select_pairs(n, None, strategy, mash_params)
        }
        _ => {
            let raw_seqs: Vec<Vec<u8>> = sequences.iter().map(|(_, s)| s.to_vec()).collect();
            let sketches = crate::mash::compute_sketches_parallel(
                &raw_seqs,
                mash_params.kmer_size,
                mash_params.sketch_size,
            );
            crate::knn_graph::select_pairs_from_sketches(&sketches, strategy)
        }
    }
}

/// Build a `FilterConfig` with adaptive scaffold parameters from a
/// `SweepgaAlignConfig`.
///
/// Defaults (`scaffold_mass=10kb`, `scaffold_jump=50kb`) are tuned for
/// whole-genome alignments. For short sequences (e.g. kilobase-scale
/// pangenome excerpts), these thresholds would filter out every
/// alignment, so we clamp them based on `avg_seq_len` via
/// `pansn::clamp_scaffold_params`.
///
/// Exposed so external callers (impg) can derive a `FilterConfig` from
/// the alignment config without re-implementing the field mapping or the
/// short-sequence clamping logic.
pub fn filter_config_from_align_cfg(cfg: &SweepgaAlignConfig, avg_seq_len: u64) -> FilterConfig {
    let (mapping_mode, mapping_per_query, mapping_per_target) =
        parse_filter_mode(&cfg.num_mappings);
    let (scaffold_mode, scaffold_per_query, scaffold_per_target) =
        parse_filter_mode(&cfg.scaffold_filter);

    let (scaffold_jump, scaffold_mass) = crate::pansn::clamp_scaffold_params(
        cfg.scaffold_jump,
        cfg.scaffold_mass,
        if avg_seq_len > 0 { Some(avg_seq_len) } else { None },
        true, // adaptive default: impg's historical behavior
    );

    FilterConfig {
        chain_gap: 0,
        min_block_length: cfg.min_map_length,
        mapping_filter_mode: mapping_mode,
        mapping_max_per_query: mapping_per_query,
        mapping_max_per_target: mapping_per_target,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: scaffold_mode,
        scaffold_max_per_query: scaffold_per_query,
        scaffold_max_per_target: scaffold_per_target,
        overlap_threshold: cfg.overlap,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: scaffold_jump,
        min_scaffold_length: scaffold_mass,
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: cfg.scaffold_dist,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: cfg.min_identity,
        min_scaffold_identity: cfg.min_identity,
    }
}

/// Apply sweepga filtering to a PAF temp file, returning a new filtered
/// temp file. Exposed so external orchestrators can filter PAFs produced
/// by other means.
///
/// Callers that have a `SweepgaAlignConfig` can build the `FilterConfig`
/// via [`filter_config_from_align_cfg`].
pub fn apply_paf_filter(
    paf_temp: tempfile::NamedTempFile,
    filter_config: FilterConfig,
) -> Result<tempfile::NamedTempFile> {
    let filtered_paf_file = tempfile::Builder::new()
        .suffix(".filtered.paf")
        .tempfile()?;

    let filter = PafFilter::new(filter_config).with_keep_self(false);
    filter
        .filter_paf(paf_temp.path(), filtered_paf_file.path())
        .map_err(|e| anyhow::anyhow!("PAF filtering failed: {}", e))?;

    Ok(filtered_paf_file)
}

/// Run sweepga alignment on in-memory named sequences with optional
/// sparsification.
///
/// When `sparsify` is `SparsificationStrategy::None`, runs all-vs-all
/// alignment. Otherwise, generates selected pairs via the strategy and
/// runs pairwise alignments only for those pairs.
///
/// Returns a `NamedTempFile` with the (optionally filtered) PAF.
pub fn sweepga_align(
    sequences: &[(String, &[u8])],
    config: &SweepgaAlignConfig,
) -> Result<tempfile::NamedTempFile> {
    if sequences.len() < 2 {
        // Nothing to align — return empty PAF
        let temp = tempfile::Builder::new().suffix(".paf").tempfile()?;
        return Ok(temp);
    }
    let pairs = generate_pairs_for_sequences(sequences, &config.sparsify, &config.mash_params);
    let total_possible = sequences.len() * (sequences.len() - 1) / 2;

    if pairs.is_empty() {
        let temp = tempfile::Builder::new().suffix(".paf").tempfile()?;
        return Ok(temp);
    }

    // If all pairs are selected (or None strategy), use fast all-vs-all
    let paf_temp = if pairs.len() == total_possible {
        sweepga_align_all_vs_all(sequences, config)?
    } else {
        log::info!(
            "sweepga: sparsification selected {} of {} pairs ({:.1}%)",
            pairs.len(),
            total_possible,
            (pairs.len() as f64 / total_possible as f64) * 100.0
        );
        sweepga_align_pairwise(sequences, &pairs, config)?
    };

    if config.no_filter {
        Ok(paf_temp)
    } else {
        let avg_seq_len = if !sequences.is_empty() {
            sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64
        } else {
            0
        };
        let filter_config = filter_config_from_align_cfg(config, avg_seq_len);
        apply_paf_filter(paf_temp, filter_config)
    }
}

/// All-vs-all alignment: write all sequences to one FASTA, align against
/// itself via `sweepga::align_self_paf`.
fn sweepga_align_all_vs_all(
    sequences: &[(String, &[u8])],
    config: &SweepgaAlignConfig,
) -> Result<tempfile::NamedTempFile> {
    let mut combined_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
    {
        let mut writer = BufWriter::new(&mut combined_fasta);
        for (name, seq) in sequences {
            writeln!(writer, ">{}", name)?;
            writer.write_all(seq)?;
            writeln!(writer)?;
        }
        writer.flush()?;
    }

    // Create FASTA index (.fai) — required by wfmash
    if config.aligner == "wfmash" {
        rust_htslib::faidx::Reader::from_path(combined_fasta.path())
            .map_err(|e| anyhow::anyhow!("Failed to create FASTA index: {e}"))?;
    }

    let avg_len = if !sequences.is_empty() {
        Some(sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64)
    } else {
        None
    };

    let wfmash_density =
        crate::orchestrator::resolve_wfmash_density(&config.sparsify, sequences.len());

    // Build the per-call aligner only on the non-batched branch; the
    // batched path constructs its own `BatchAligner` internally and would
    // drop this one untouched.
    let result = match config.batch_bytes.as_deref() {
        None => {
            let aligner: Box<dyn Aligner> = create_aligner_adaptive(
                &config.aligner,
                config.kmer_frequency,
                config.num_threads,
                config.min_aln_length,
                config.map_pct_identity.clone(),
                config.temp_dir.clone(),
                None, // segment_length: adapt from avg_len
                avg_len,
                wfmash_density,
                None, // num_mappings: wfmash default (-n 1)
                None, // pairs_file
            )?;
            crate::align_self_paf_direct(aligner.as_ref(), combined_fasta.path())
        }
        Some(batch_bytes) => crate::align_self_paf_batched(
            combined_fasta.path(),
            &config.aligner,
            config.kmer_frequency,
            config.num_threads,
            config.min_aln_length,
            config.map_pct_identity.clone(),
            config.temp_dir.clone(),
            wfmash_density,
            batch_bytes,
            true, // quiet
        ),
    };
    result.map_err(|e| anyhow::anyhow!("{} alignment failed: {}", config.aligner, e))
}

/// Pairwise alignment dispatcher.
fn sweepga_align_pairwise(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
) -> Result<tempfile::NamedTempFile> {
    let avg_len = if !sequences.is_empty() {
        Some(sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64)
    } else {
        None
    };
    let wfmash_density =
        crate::orchestrator::resolve_wfmash_density(&config.sparsify, sequences.len());

    if config.aligner == "wfmash" {
        sweepga_align_pairwise_wfmash(sequences, pairs, config, avg_len, wfmash_density)
    } else {
        sweepga_align_pairwise_generic(sequences, pairs, config, avg_len, wfmash_density)
    }
}

/// Optimized pairwise alignment using wfmash `--pairs-file`: single
/// wfmash invocation over a combined FASTA plus a pairs TSV.
fn sweepga_align_pairwise_wfmash(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
    avg_len: Option<u64>,
    wfmash_density: Option<f64>,
) -> Result<tempfile::NamedTempFile> {
    let mut combined_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
    {
        let mut writer = BufWriter::new(&mut combined_fasta);
        for (name, seq) in sequences {
            writeln!(writer, ">{}", name)?;
            writer.write_all(seq)?;
            writeln!(writer)?;
        }
        writer.flush()?;
    }

    rust_htslib::faidx::Reader::from_path(combined_fasta.path())
        .map_err(|e| anyhow::anyhow!("Failed to create FASTA index: {e}"))?;

    let mut pairs_tsv = tempfile::Builder::new().suffix(".pairs.tsv").tempfile()?;
    {
        let mut writer = BufWriter::new(&mut pairs_tsv);
        writeln!(writer, "# query_name\ttarget_name")?;
        for &(i, j) in pairs {
            // Write both directions so wfmash aligns A→B and B→A.
            writeln!(writer, "{}\t{}", sequences[i].0, sequences[j].0)?;
            writeln!(writer, "{}\t{}", sequences[j].0, sequences[i].0)?;
        }
        writer.flush()?;
    }

    log::info!(
        "sweepga: wfmash batch alignment: {} sequences, {} pairs, pairs file: {}",
        sequences.len(),
        pairs.len(),
        pairs_tsv.path().display()
    );

    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        config.kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        config.map_pct_identity.clone(),
        config.temp_dir.clone(),
        None,
        avg_len,
        wfmash_density,
        None,
        Some(pairs_tsv.path().to_path_buf()),
    )?;

    aligner
        .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
        .map_err(|e| anyhow::anyhow!("wfmash batch alignment failed: {}", e))
}

/// Fallback per-pair alignment for aligners without `--pairs-file`
/// support (currently FastGA). Pre-writes one FASTA per unique index
/// to avoid O(pairs) tempfile churn.
fn sweepga_align_pairwise_generic(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
    avg_len: Option<u64>,
    wfmash_density: Option<f64>,
) -> Result<tempfile::NamedTempFile> {
    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        config.kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        config.map_pct_identity.clone(),
        config.temp_dir.clone(),
        None,
        avg_len,
        wfmash_density,
        None,
        None,
    )?;

    let mut combined_paf = tempfile::Builder::new().suffix(".paf").tempfile()?;

    let unique_indices: std::collections::BTreeSet<usize> =
        pairs.iter().flat_map(|&(i, j)| [i, j]).collect();
    let mut fasta_files: std::collections::HashMap<usize, tempfile::NamedTempFile> =
        std::collections::HashMap::new();
    for idx in unique_indices {
        let mut fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
        {
            let mut writer = BufWriter::new(&mut fasta);
            writeln!(writer, ">{}", sequences[idx].0)?;
            writer.write_all(sequences[idx].1)?;
            writeln!(writer)?;
            writer.flush()?;
        }
        if config.aligner == "wfmash" {
            rust_htslib::faidx::Reader::from_path(fasta.path())
                .map_err(|e| anyhow::anyhow!("Failed to create FASTA index: {e}"))?;
        }
        fasta_files.insert(idx, fasta);
    }

    for &(i, j) in pairs {
        match aligner.align_to_temp_paf(fasta_files[&i].path(), fasta_files[&j].path()) {
            Ok(pair_paf) => {
                let contents = std::fs::read(pair_paf.path())?;
                if !contents.is_empty() {
                    combined_paf.write_all(&contents)?;
                }
            }
            Err(e) => {
                log::warn!(
                    "sweepga: pairwise alignment failed for {} vs {}: {}",
                    sequences[i].0,
                    sequences[j].0,
                    e
                );
            }
        }
    }

    combined_paf.flush()?;
    Ok(combined_paf)
}
