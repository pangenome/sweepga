//! Alignment orchestration helpers shared between the sweepga binary and
//! external library consumers (impg).

use crate::knn_graph::SparsificationStrategy;

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

