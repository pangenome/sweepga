//! Shared alignment-options clap struct.
//!
//! `AlnArgs` is the single source of truth for every flag that controls an
//! alignment run. The sweepga binary flattens it into its own `Args`, and
//! downstream consumers (e.g. impg) flatten it into their own clap structs
//! with `#[clap(flatten)] aln: sweepga::AlnArgs`. **There is no other
//! definition of these flags anywhere.**
//!
//! Design notes:
//! - No short flags live on `AlnArgs`. Short flags (`-i`, `-b`, `-o`, etc.)
//!   collide with impg's per-command shorts; dropping them is the price of
//!   a single flattenable struct.
//! - `--threads`, `--quiet`, the positional input files, and the output
//!   path stay on the binary's own `Args`, NOT on `AlnArgs`. They are
//!   cross-cutting / output-layer concerns, not alignment parameters.

use anyhow::Result;
use clap::Args;

use crate::knn_graph::SparsificationStrategy;

/// Parse a number that may have a metric suffix (k/K=1e3, m/M=1e6, g/G=1e9).
///
/// Accepted: integer or float bases, optional single-char suffix.
/// Returns a descriptive error string on invalid input.
pub fn parse_metric_number(s: &str) -> Result<u64, String> {
    if s.is_empty() {
        return Err("Empty string".to_string());
    }

    let (num_part, suffix) = if s.ends_with(|c: char| c.is_ascii_alphabetic()) {
        let last_char = s.chars().last().unwrap();
        (&s[..s.len() - last_char.len_utf8()], Some(last_char))
    } else {
        (s, None)
    };

    let base: f64 = num_part
        .parse()
        .map_err(|e| format!("Invalid number: {e}"))?;

    let multiplier = match suffix {
        Some('k') | Some('K') => 1000.0,
        Some('m') | Some('M') => 1_000_000.0,
        Some('g') | Some('G') => 1_000_000_000.0,
        Some(c) => {
            return Err(format!(
                "Unknown suffix '{c}'. Use k/K (1000), m/M (1e6), or g/G (1e9)"
            ))
        }
        None => 1.0,
    };

    let result = base * multiplier;

    if result > u64::MAX as f64 {
        return Err(format!("Value {result} too large for u64"));
    }

    Ok(result as u64)
}

/// Parse an identity threshold string into a 0..=1 fraction.
///
/// Accepts three forms:
/// - A fraction in [0, 1] (e.g. `"0.9"`).
/// - A percentage > 1 (e.g. `"90"` → 0.9).
/// - An `aniN[+/-offset]` preset (e.g. `"ani50"`, `"ani50-2"`). Requires
///   `ani_percentile` to be `Some(median_identity)`; typically populated
///   by a first-pass ANI survey of the input alignments. Only the
///   median (`ani50`) percentile is currently honored; other values are
///   silently mapped to the median.
///
/// Returns an error for unparseable inputs or ANI-based values when
/// `ani_percentile` is None.
pub fn parse_identity_value(value: &str, ani_percentile: Option<f64>) -> Result<f64> {
    let lower = value.to_lowercase();

    if let Some(remainder) = lower.strip_prefix("ani") {
        // Parse aniN, aniN+X, or aniN-X
        if let Some(ani_value) = ani_percentile {
            if remainder.is_empty() {
                // Default "ani" means ani50
                return Ok(ani_value);
            }

            // Parse percentile number and optional offset
            let (percentile_str, offset_part) = if let Some(plus_pos) = remainder.find('+') {
                (
                    &remainder[..plus_pos],
                    Some(('+', &remainder[plus_pos + 1..])),
                )
            } else if let Some(minus_pos) = remainder.find('-') {
                (
                    &remainder[..minus_pos],
                    Some(('-', &remainder[minus_pos + 1..])),
                )
            } else {
                (remainder, None)
            };

            // Only ani50 (median) honored today; other percentiles silently
            // fall back to it. TODO: extend for ani25/ani75 presets.
            let _ = percentile_str;

            if let Some((sign, offset_str)) = offset_part {
                let offset: f64 = offset_str
                    .parse()
                    .map_err(|_| anyhow::anyhow!("Invalid ANI offset: {offset_str}"))?;
                match sign {
                    '+' => Ok((ani_value + offset / 100.0).min(1.0)),
                    '-' => Ok((ani_value - offset / 100.0).max(0.0)),
                    _ => Ok(ani_value),
                }
            } else {
                Ok(ani_value)
            }
        } else {
            anyhow::bail!("Cannot use ANI-based threshold without input alignments");
        }
    } else if let Ok(val) = value.parse::<f64>() {
        if val > 1.0 {
            Ok(val / 100.0) // percentage → fraction
        } else {
            Ok(val) // already a fraction
        }
    } else {
        anyhow::bail!("Invalid identity value: {value}");
    }
}

/// All alignment-related CLI flags, as a single flattenable `clap::Args`.
///
/// See the module docstring for invariants.
#[derive(Args, Debug, Clone)]
pub struct AlnArgs {
    // ========================================================================
    // Alignment (FASTA input only)
    // ========================================================================
    /// Aligner to use for FASTA input
    #[clap(long = "aligner", default_value = "fastga", value_parser = ["fastga", "wfmash"],
           help_heading = "Alignment options")]
    pub aligner: String,

    /// Use FastGA aligner (shorthand for --aligner fastga)
    #[clap(long = "fastga", help_heading = "Alignment options")]
    pub use_fastga: bool,

    /// Use wfmash aligner (shorthand for --aligner wfmash)
    #[clap(long = "wfmash", help_heading = "Alignment options")]
    pub use_wfmash: bool,

    /// FastGA k-mer frequency threshold (use k-mers occurring ≤ N times).
    /// Overrides --fastga-frequency-multiplier when set.
    #[clap(long = "fastga-frequency", help_heading = "Alignment options")]
    pub frequency: Option<usize>,

    /// FastGA k-mer frequency multiplier: when --fastga-frequency is unset,
    /// the effective frequency = `num_pansn_haplotypes * multiplier`
    /// (clamped to at least 10). Default 1 reproduces FastGA's historical
    /// behavior of one k-mer hit per haplotype.
    #[clap(long = "fastga-frequency-multiplier", default_value = "1",
           help_heading = "Alignment options")]
    pub fastga_frequency_multiplier: usize,

    /// Minimum percent identity for wfmash mapping (e.g. "90" or ANI preset "ani50-2")
    #[clap(long = "map-pct-identity", help_heading = "Alignment options")]
    pub map_pct_identity: Option<String>,

    /// Align all genome pairs separately (slower, uses more memory, but handles many genomes)
    #[clap(long = "all-pairs", help_heading = "Alignment options")]
    pub all_pairs: bool,

    /// Maximum sequence data per batch (e.g., "50M", "2G"). Partitions
    /// genomes into batches. Stored as a raw string so library consumers
    /// (like impg's `sweepga::align_self_paf`) can pass it verbatim;
    /// parse via `crate::cli::parse_metric_number` at use sites when you
    /// need a numeric byte count.
    #[clap(long = "batch-bytes", help_heading = "Alignment options")]
    pub batch_bytes: Option<String>,

    /// Explicit number of genomes per batch (manual override). Overrides --batch-bytes and --max-disk.
    #[clap(long = "batch-size", help_heading = "Alignment options")]
    pub batch_size: Option<usize>,

    /// Maximum disk space for temporary files during alignment (e.g., "100G", "500M").
    /// Computes batch size automatically to stay within budget. Overrides --batch-bytes.
    #[clap(long = "max-disk", value_parser = parse_metric_number,
           help_heading = "Alignment options")]
    pub max_disk: Option<u64>,

    /// Compress k-mer index with zstd for ~2x disk savings and faster I/O
    #[clap(long = "zstd", help_heading = "Alignment options")]
    pub zstd_compress: bool,

    /// Zstd compression level (1-19, higher = smaller files but slower). Default: 3
    #[clap(long = "zstd-level", default_value = "3", help_heading = "Alignment options")]
    pub zstd_level: u32,

    // ========================================================================
    // Basic Filtering
    // ========================================================================
    /// Minimum alignment block length (omit to use aligner default)
    #[clap(long = "min-aln-length", value_parser = parse_metric_number,
           help_heading = "Basic filtering")]
    pub block_length: Option<u64>,

    /// n:m-best mappings kept in query:target dimensions. 1:1 (orthogonal),
    /// use ∞/many for unbounded. Default matches impg's historical
    /// pangenome setting: no mapping-axis filter before scaffolding.
    #[clap(long = "num-mappings", default_value = "many:many", help_heading = "Basic filtering")]
    pub num_mappings: String,

    /// Maximum overlap ratio for plane sweep filtering
    #[clap(long = "overlap", default_value = "0.95", help_heading = "Basic filtering")]
    pub overlap: f64,

    /// Scoring function for plane sweep
    #[clap(long = "scoring", default_value = "log-length-ani",
           value_parser = ["ani", "length", "length-ani", "log-length-ani", "matches"],
           help_heading = "Basic filtering")]
    pub scoring: String,

    /// Minimum per-mapping identity filter applied to PAF records after
    /// alignment (0-1 fraction, 1-100%, or "aniN" for Nth percentile).
    /// This is a post-alignment PAF filter, not an alignment-stage knob —
    /// for wfmash's alignment-stage identity cutoff, use `--map-pct-identity`.
    #[clap(long = "min-aln-identity", default_value = "0", help_heading = "Basic filtering")]
    pub min_identity: String,

    /// Keep self-mappings (excluded by default)
    #[clap(long = "self", help_heading = "Basic filtering")]
    pub keep_self: bool,

    /// Disable all filtering
    #[clap(long = "no-filter", help_heading = "Basic filtering")]
    pub no_filter: bool,

    // ========================================================================
    // Scaffolding and Chaining
    // ========================================================================
    /// Scaffold jump (gap) distance. 0 = disable scaffolding, >0 = enable
    /// (accepts k/m/g suffix). Default matches impg's historical 50kb —
    /// pass `--scaffold-jump 0` to turn scaffolding off.
    #[clap(long = "scaffold-jump", default_value = "50k", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    pub scaffold_jump: u64,

    /// Minimum scaffold chain length (accepts k/m/g suffix)
    #[clap(long = "scaffold-mass", default_value = "10k", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    pub scaffold_mass: u64,

    /// Scaffold filter mode: "1:1" (best), "M:N" (M per query, N per target), "many" (unbounded)
    #[clap(long = "scaffold-filter", default_value = "many:many",
           help_heading = "Scaffolding and chaining")]
    pub scaffold_filter: String,

    /// Scaffold chain overlap threshold for plane sweep filtering
    #[clap(long = "scaffold-overlap", default_value = "0.5",
           help_heading = "Scaffolding and chaining")]
    pub scaffold_overlap: f64,

    /// Maximum Euclidean distance from scaffold anchor for rescue (0 = no rescue, accepts k/m/g suffix)
    #[clap(long = "scaffold-dist", default_value = "0", value_parser = parse_metric_number,
           help_heading = "Scaffolding and chaining")]
    pub scaffold_dist: u64,

    /// Minimum scaffold identity threshold (0-1 fraction, 1-100%, "aniN", or defaults to --min-aln-identity)
    #[clap(long = "min-scaffold-identity", default_value = "0",
           help_heading = "Scaffolding and chaining")]
    pub min_scaffold_identity: String,

    /// Output scaffold chains only (for debugging)
    #[clap(long = "scaffolds-only", help_heading = "Scaffolding and chaining")]
    pub scaffolds_only: bool,

    /// Disable scaffold-parameter adaptation to input sequence length.
    ///
    /// By default, when the FASTA average sequence length is known,
    /// `--scaffold-mass` is clamped to `avg_seq_len * 3/5` (nice-rounded)
    /// and `--scaffold-jump` is clamped to `avg_seq_len * 10`. This is a
    /// no-op for whole-genome inputs (user values are already smaller
    /// than the clamps) but prevents default thresholds from filtering
    /// out everything when inputs are short (e.g. pangenome-window
    /// excerpts). Pass `--no-adaptive-scaffolds` to turn the clamping off.
    #[clap(long = "no-adaptive-scaffolds", help_heading = "Scaffolding and chaining")]
    pub no_adaptive_scaffolds: bool,

    // ========================================================================
    // Advanced Filtering
    // ========================================================================
    /// Sparsification strategy: `none` (default = all pairs), `auto`,
    /// `<frac>` / `random:<frac>`, `giant:<prob>` / `connectivity:<prob>`,
    /// `tree:<near>[:<far>[:<random>]]` / `knn:…`,
    /// `wfmash:auto` / `wfmash:<frac>`. Parsed once by clap at
    /// CLI-parse time via `SparsificationStrategy::from_str`.
    #[clap(
        long = "sparsify",
        default_value = "none",
        help_heading = "Advanced filtering",
        value_parser = |s: &str| s.parse::<SparsificationStrategy>()
    )]
    pub sparsify: SparsificationStrategy,

    /// K-mer size used by the mash sketches driving tree/giant sparsification
    #[clap(long = "mash-kmer-size", default_value_t = crate::mash::DEFAULT_KMER_SIZE,
           help_heading = "Advanced filtering")]
    pub mash_kmer_size: usize,

    /// Sketch size (number of minimizers) for the mash sketches driving tree/giant sparsification
    #[clap(long = "mash-sketch-size", default_value_t = crate::mash::DEFAULT_SKETCH_SIZE,
           help_heading = "Advanced filtering")]
    pub mash_sketch_size: usize,

    /// Instead of running alignments, emit one shell command per genome
    /// pair to stdout (or `--output-file`) and exit. Useful for cluster /
    /// grid dispatch: each emitted line is a standalone `sweepga` invocation
    /// that aligns one pair. Requires multiple FASTA inputs.
    #[clap(long = "joblist", help_heading = "Advanced filtering")]
    pub joblist: bool,

    /// Output directory for per-pair PAFs emitted by `--joblist`. Defaults
    /// to the current working directory. Ignored when `--joblist` is unset.
    #[clap(long = "joblist-output-dir", help_heading = "Advanced filtering")]
    pub joblist_output_dir: Option<String>,

    /// Method for calculating ANI: all, orthogonal, nX[-sort] (e.g. n50, n90-identity, n100-score)
    #[clap(long = "ani-method", default_value = "n100", help_heading = "Advanced filtering")]
    pub ani_method: String,

    // ========================================================================
    // Temp directory (cross-cutting but alignment-scoped)
    // ========================================================================
    /// Temporary directory for intermediate files (defaults to TMPDIR env var, then current directory; use "ramdisk" for /dev/shm)
    #[clap(long = "temp-dir", help_heading = "Alignment options")]
    pub tempdir: Option<String>,

    // ========================================================================
    // AGC Archive Options
    // ========================================================================
    /// Extract only samples matching this prefix (AGC input)
    #[clap(long = "agc-prefix", help_heading = "AGC archive options")]
    pub agc_prefix: Option<String>,

    /// Extract only these samples (comma-separated or @file with one per line)
    #[clap(long = "agc-samples", help_heading = "AGC archive options")]
    pub agc_samples: Option<String>,

    /// Query samples for asymmetric alignment (prefix or comma-separated/@file list)
    #[clap(long = "agc-queries", help_heading = "AGC archive options")]
    pub agc_queries: Option<String>,

    /// Target samples for asymmetric alignment (prefix or comma-separated/@file list)
    #[clap(long = "agc-targets", help_heading = "AGC archive options")]
    pub agc_targets: Option<String>,

    /// Temp directory for AGC extraction (defaults to --temp-dir)
    #[clap(long = "agc-temp-dir", help_heading = "AGC archive options")]
    pub agc_tempdir: Option<String>,

    // ========================================================================
    // Pair Selection (for incremental/targeted alignment)
    // ========================================================================
    /// File of sample pairs to align (TSV: query<tab>target, one pair per line)
    #[clap(long = "pairs", help_heading = "Pair selection")]
    pub pairs_file: Option<String>,

    /// File to append completed pairs (for resume capability)
    #[clap(long = "pairs-done", help_heading = "Pair selection")]
    pub pairs_done: Option<String>,

    /// File to write remaining pairs after run (for resume capability)
    #[clap(long = "pairs-remaining", help_heading = "Pair selection")]
    pub pairs_remaining: Option<String>,

    /// List all sample pairs and exit (for generating pairs files)
    #[clap(long = "list-pairs", help_heading = "Pair selection")]
    pub list_pairs: bool,

    /// Shuffle pair order (for load balancing across runs)
    #[clap(long = "shuffle-pairs", help_heading = "Pair selection")]
    pub shuffle_pairs: bool,

    /// Seed for pair shuffling (for reproducible ordering)
    #[clap(long = "shuffle-seed", help_heading = "Pair selection")]
    pub shuffle_seed: Option<u64>,

    /// Maximum number of pairs to process in this run (0 = unlimited)
    #[clap(long = "max-pairs", default_value = "0", help_heading = "Pair selection")]
    pub max_pairs: usize,

    /// Start index for pair range (0-based, use with shuffle-seed for stable ordering)
    #[clap(long = "pair-start", default_value = "0", help_heading = "Pair selection")]
    pub pair_start: usize,

    /// Sparsification strategy for pair selection (none, auto, random:<frac>, giant:<prob>, tree:<near>:<far>:<random>[:<kmer>])
    #[clap(long = "sparsify-pairs", default_value = "none", help_heading = "Pair selection")]
    pub sparsify_pairs: String,
}
