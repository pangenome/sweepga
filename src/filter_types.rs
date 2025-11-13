/// Shared filtering types used across multiple modules
///
/// This module contains common types to avoid circular dependencies between
/// paf_filter, plane_sweep_scaffold, and other filtering modules.

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
