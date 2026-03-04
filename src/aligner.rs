//! Aligner trait for pluggable alignment backends.

use anyhow::Result;
use std::path::Path;
use tempfile::NamedTempFile;

/// Trait for alignment backends (FastGA, wfmash, etc.)
#[allow(dead_code)]
pub trait Aligner: Send + Sync {
    /// Run alignment and write output to a temporary PAF file.
    fn align_to_temp_paf(&self, queries: &Path, targets: &Path) -> Result<NamedTempFile>;

    /// Run alignment and return PAF bytes directly (no temp file).
    fn align_direct_paf(&self, queries: &Path, targets: &Path) -> Result<Vec<u8>>;

    /// Run alignment and return .1aln output (FastGA only).
    fn align_to_temp_1aln(&self, _queries: &Path, _targets: &Path) -> Result<NamedTempFile> {
        anyhow::bail!("This aligner does not support .1aln output")
    }

    /// Name of this aligner backend.
    fn name(&self) -> &str;
}
