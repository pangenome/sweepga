//! PanSN (Pangenome Sequence Naming) helpers.
//!
//! PanSN format is `SAMPLE#HAPLOTYPE#CONTIG`. For pangenome alignment, the
//! relevant "genome count" is the number of distinct `SAMPLE#HAPLOTYPE`
//! prefixes across all input FASTAs — not the number of sequences, which
//! can be orders of magnitude larger (one per contig).
//!
//! FastGA uses this count as its default k-mer frequency threshold: a k-mer
//! appearing in all genomes is likely non-informative. With
//! `--fastga-frequency-multiplier M`, the threshold becomes `genomes * M`,
//! letting users tune for diverged or repetitive pangenomes.

use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Open a FASTA file, transparently handling `.gz`.
///
/// Returns a boxed `BufRead` over the decompressed stream. Uses
/// `MultiGzDecoder` so concatenated gzip members are all consumed (bgzipped
/// FASTAs are stored this way).
fn open_fasta(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open FASTA '{}'", path.display()))?;

    let is_gz = path.extension().and_then(|s| s.to_str()) == Some("gz");
    if is_gz {
        use flate2::read::MultiGzDecoder;
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Grouping level for PanSN keys. Directly usable as a `clap::ValueEnum`
/// for CLI parsing; the `Sequence` variant is skipped from CLI choices
/// because "no PanSN aggregation" is expressed by omitting the flag
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum PanSnLevel {
    /// Whole name. Use for per-contig/per-sequence
    /// counting where PanSN grouping is not desired.
    #[value(skip)]
    Sequence,
    /// First `#`-separated segment (`SAMPLE`). Falls back to the whole
    /// name for non-PanSN inputs.
    Sample,
    /// First two `#`-separated segments (`SAMPLE#HAPLOTYPE`). Falls back
    /// to just `SAMPLE` if the second segment is missing/empty, which in
    /// turn falls back to the whole name for non-PanSN inputs.
    Haplotype,
}

/// Extract a grouping key from a sequence name at the requested level.
///
/// Returns `None` when the normalized name is empty.
pub fn extract_pansn_key(name: &str, level: PanSnLevel) -> Option<String> {
    let name = name.trim_start_matches('>').trim();
    let name = name.split_whitespace().next().unwrap_or("");
    let base = name.split(':').next().unwrap_or(name);
    if base.is_empty() {
        return None;
    }
    match level {
        PanSnLevel::Sequence => Some(base.to_string()),
        PanSnLevel::Sample => {
            let sample = base.split('#').next().unwrap_or(base);
            if sample.is_empty() {
                None
            } else {
                Some(sample.to_string())
            }
        }
        PanSnLevel::Haplotype => {
            let mut parts = base.split('#');
            let sample = parts.next().unwrap_or("");
            if sample.is_empty() {
                return None;
            }
            match parts.next() {
                Some(h) if !h.is_empty() => Some(format!("{sample}#{h}")),
                _ => Some(sample.to_string()),
            }
        }
    }
}

/// Count unique PanSN keys at the requested level across an iterator of names.
pub fn count_pansn_keys<'a, I>(names: I, level: PanSnLevel) -> usize
where
    I: IntoIterator<Item = &'a str>,
{
    let set: HashSet<String> = names
        .into_iter()
        .filter_map(|n| extract_pansn_key(n, level))
        .collect();
    set.len().max(1)
}

/// Group sequence indices by their PanSN key at the requested level.
///
/// Returns a `Vec` of groups; each group is the list of input indices that
/// share a key. Non-PanSN names (or names with no extractable key) fall back
/// to the whole name, so every such contig becomes its own group — callers
/// can therefore use this uniformly and treat "one group per name" as
/// "input is not PanSN".
///
/// Ordering is stable: groups are sorted by their PanSN key, which makes the
/// output reproducible across runs. Indices within each group keep their
/// input order.
pub fn group_indices_by_pansn<'a, I>(names: I, level: PanSnLevel) -> Vec<Vec<usize>>
where
    I: IntoIterator<Item = &'a str>,
{
    use std::collections::BTreeMap;
    let mut map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, name) in names.into_iter().enumerate() {
        let key = extract_pansn_key(name, level).unwrap_or_else(|| name.to_string());
        map.entry(key).or_default().push(i);
    }
    map.into_values().collect()
}

/// Count unique PanSN haplotypes across one or more FASTA files.
pub fn count_haplotypes<P: AsRef<Path>>(fasta_paths: &[P]) -> Result<usize> {
    let mut haplotypes: HashSet<String> = HashSet::new();

    for path in fasta_paths {
        let path = path.as_ref();
        let reader = open_fasta(path)?;

        for line in reader.lines() {
            let line = line.with_context(|| {
                format!("Failed to read a line from '{}'", path.display())
            })?;
            if line.starts_with('>') {
                if let Some(key) = extract_pansn_key(&line, PanSnLevel::Haplotype) {
                    haplotypes.insert(key);
                }
            }
        }
    }

    Ok(haplotypes.len().max(1))
}

/// Resolve the effective FastGA k-mer frequency threshold.
///
/// - If `explicit` is `Some(n)`, returns `n` (user override wins).
/// - Otherwise counts PanSN haplotypes across `fasta_paths` and returns
///   `num_haplotypes * multiplier` (no lower bound). This matches impg's
///   historical behavior of using the raw haplotype count for small
///   pangenomes — FastGA's internal floor of 10 is not reimposed here
///   because users explicitly requesting `multiplier * num_haplotypes`
///   may want a tight threshold on small inputs.
///
/// Callers that already have a haplotype count should pass `explicit`
/// rather than re-counting.
pub fn resolve_fastga_frequency<P: AsRef<Path>>(
    explicit: Option<usize>,
    multiplier: usize,
    fasta_paths: &[P],
) -> Result<usize> {
    if let Some(n) = explicit {
        return Ok(n);
    }
    let num_haplotypes = count_haplotypes(fasta_paths)?;
    Ok(num_haplotypes.saturating_mul(multiplier.max(1)))
}

/// Round a value to a "nice" multiple based on its magnitude:
/// ≤500 → multiple of 50, ≤1000 → 100, ≤3000 → 200, >3000 → 500.
/// Zero stays zero. Non-zero values are floored at the chosen step so
/// rounding never collapses a small positive number to 0.
pub fn round_nice(v: u64) -> u64 {
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

/// Clamp scaffold parameters based on average sequence length.
///
/// Default scaffold thresholds (scaffold_mass=10kb, scaffold_jump=50kb) are
/// tuned for whole-genome alignment. For short sequences (e.g. 1-2 kb query
/// excerpts from `impg query -o fasta`), these thresholds filter out every
/// alignment. When `adaptive` is true and `avg_seq_len > 0`:
///
/// - `scaffold_mass` clamps to `avg_seq_len * 3/5` (rounded via `round_nice`)
/// - `scaffold_jump` clamps to `avg_seq_len * 10`
///
/// When `adaptive` is false or `avg_seq_len` is `None`/`Some(0)`, the user
/// values pass through unchanged — callers can wire this helper at every
/// `FilterConfig` construction site without worrying about over-eagerly
/// shrinking thresholds.
#[allow(dead_code)] // Part of the public library API consumed by impg.
pub fn clamp_scaffold_params(
    user_jump: u64,
    user_mass: u64,
    avg_seq_len: Option<u64>,
    adaptive: bool,
) -> (u64, u64) {
    if !adaptive {
        return (user_jump, user_mass);
    }
    let Some(avg) = avg_seq_len else {
        return (user_jump, user_mass);
    };
    if avg == 0 {
        return (user_jump, user_mass);
    }
    let jump = user_jump.min(avg.saturating_mul(10));
    let mass = round_nice(user_mass.min(avg.saturating_mul(3) / 5));
    (jump, mass)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_pansn_key_haplotype() {
        assert_eq!(
            extract_pansn_key(">HG01106#1#CM087962.1", PanSnLevel::Haplotype),
            Some("HG01106#1".to_string())
        );
        assert_eq!(
            extract_pansn_key(">HG01106#1#chr1 extra annotation", PanSnLevel::Haplotype),
            Some("HG01106#1".to_string())
        );
    }

    #[test]
    fn extract_pansn_key_sample() {
        assert_eq!(
            extract_pansn_key("HG01106#1#CM087962.1", PanSnLevel::Sample),
            Some("HG01106".to_string())
        );
    }

    #[test]
    fn extract_pansn_key_sequence_strips_interval() {
        // impg query-output names carry `:start-end` — interval suffix is stripped.
        assert_eq!(
            extract_pansn_key("HG01106#1#chr1:29000000-29003000", PanSnLevel::Sequence),
            Some("HG01106#1#chr1".to_string())
        );
        assert_eq!(
            extract_pansn_key("HG01106#1#chr1:29000000-29003000", PanSnLevel::Haplotype),
            Some("HG01106#1".to_string())
        );
    }

    #[test]
    fn extract_pansn_key_non_pansn_returns_whole_name() {
        // Fewer than 2 `#` segments: fall back to whole name for haplotype mode.
        assert_eq!(
            extract_pansn_key("chr1", PanSnLevel::Haplotype),
            Some("chr1".to_string())
        );
        assert_eq!(
            extract_pansn_key(">chr1 extra", PanSnLevel::Haplotype),
            Some("chr1".to_string())
        );
    }

    #[test]
    fn count_pansn_keys_basic() {
        let names = [
            "HG01#1#chr1",
            "HG01#1#chr2",
            "HG01#2#chr1",
            "HG02#1#chr1",
            "chm13#0#chr6",
        ];
        assert_eq!(count_pansn_keys(names.iter().copied(), PanSnLevel::Sample), 3);
        assert_eq!(count_pansn_keys(names.iter().copied(), PanSnLevel::Haplotype), 4);
        assert_eq!(count_pansn_keys(names.iter().copied(), PanSnLevel::Sequence), 5);
    }

    #[test]
    fn count_pansn_keys_empty_returns_one() {
        let names: [&str; 0] = [];
        // `.max(1)` invariant: never return 0 for `num * multiplier` math.
        assert_eq!(count_pansn_keys(names.iter().copied(), PanSnLevel::Haplotype), 1);
    }

    #[test]
    fn resolve_frequency_uses_explicit_override() {
        let empty: [&Path; 0] = [];
        let f = resolve_fastga_frequency(Some(42), 99, &empty).unwrap();
        assert_eq!(f, 42);
    }

    #[test]
    fn round_nice_steps() {
        assert_eq!(round_nice(0), 0);
        assert_eq!(round_nice(120), 100);
        assert_eq!(round_nice(480), 500);
        assert_eq!(round_nice(950), 1000);
        assert_eq!(round_nice(2900), 3000);
        assert_eq!(round_nice(7200), 7000);
    }

    #[test]
    fn clamp_scaffolds_noop_when_disabled() {
        let (j, m) = clamp_scaffold_params(50_000, 10_000, Some(1000), false);
        assert_eq!((j, m), (50_000, 10_000));
    }

    #[test]
    fn clamp_scaffolds_noop_without_avg() {
        let (j, m) = clamp_scaffold_params(50_000, 10_000, None, true);
        assert_eq!((j, m), (50_000, 10_000));
    }

    #[test]
    fn clamp_scaffolds_shrinks_for_short_inputs() {
        // 1 kb avg seq: jump clamped to 10 kb, mass clamped to 600 bp → round_nice → 600.
        let (j, m) = clamp_scaffold_params(50_000, 10_000, Some(1000), true);
        assert_eq!(j, 10_000);
        assert_eq!(m, 600);
    }

    #[test]
    fn clamp_scaffolds_preserves_small_user_values() {
        // If the user already passed tighter values, clamping never loosens them.
        let (j, m) = clamp_scaffold_params(5_000, 3_000, Some(1_000_000), true);
        assert_eq!(j, 5_000);
        // 3_000 clamped to min(3_000, 600_000) == 3_000 → round_nice(3000) == 3000.
        assert_eq!(m, 3_000);
    }
}
