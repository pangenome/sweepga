/// Shared filtering types used across multiple modules
///
/// This module contains common types to avoid circular dependencies between
/// paf_filter, plane_sweep_scaffold, and other filtering modules.
///
/// Scoring function for plane sweep
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScoringFunction {
    Identity,          // Identity only
    Length,            // Length only
    LengthIdentity,    // Length * Identity
    LogLengthIdentity, // log(Length) * Identity (default)
    Matches,           // Total matches only (gap-neutral)
}

/// Filtering mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterMode {
    OneToOne,   // 1:1 - best mapping per query AND per target
    OneToMany,  // 1:N - best mapping per query, N per target
    ManyToMany, // N:N - N mappings per query and per target
}

/// How sequences are grouped into genomes for filtering.
///
/// Grouping decides which sequences' mappings compete with each other during
/// the plane sweep, which is what removes cross-chromosome noise. It is never
/// inferred from contig names: without a real signal the safe default is to
/// not group (`None`), which reproduces the historical per-sequence behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenomeGrouping {
    /// PanSN prefix if present, else one genome per input file, else `None`.
    Auto,
    /// PanSN `#` prefix (genome#haplotype#).
    PanSn,
    /// One genome per input file.
    File,
    /// All query sequences are one genome; all target sequences are another.
    Pairwise,
    /// One genome per sequence (no cross-sequence competition).
    None,
}

impl Default for GenomeGrouping {
    fn default() -> Self {
        GenomeGrouping::Auto
    }
}
