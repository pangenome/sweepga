//! `--joblist` output: emit one shell command per genome pair instead of
//! running alignments in-process.
//!
//! Intended for cluster / grid submission: the user dispatches the list to
//! their scheduler, which runs each pair independently. Downstream
//! filtering / scaffolding happens either per-pair (the user passes
//! filter flags on the emitted command line) or in a post-processing
//! merge step.
//!
//! Kept intentionally small: one free function, no trait zoo. The
//! template is a `sweepga` self-invocation, so the emitted commands
//! reproduce the exact flag set the user is running with.

use anyhow::{Context, Result};
use std::io::Write;
use std::path::Path;

/// Configuration for a joblist emission. Only the fields actually
/// reproduced into the shell command belong here — keep this lean.
pub struct JoblistEmitConfig<'a> {
    /// Path to the sweepga binary to invoke. Normally the current
    /// executable's path (`std::env::current_exe`).
    pub sweepga_bin: &'a Path,
    /// Output directory for per-pair PAFs. Each emitted command writes
    /// `<output_dir>/<query_stem>_vs_<target_stem>.paf`.
    pub output_dir: &'a Path,
    /// Number of threads to pass via `--threads`.
    pub threads: usize,
    /// Extra CLI flags to append to every emitted command. Typically the
    /// filter/scaffolding args the user invoked sweepga with, so each
    /// dispatched job reproduces the same filtering behavior.
    pub extra_flags: &'a [String],
}

/// Emit one `sweepga QUERY TARGET --output-file …` line per pair.
///
/// `pairs` are `(query_path, target_path)` tuples. The output filename is
/// derived from the file stems, so pairs with the same stems will collide
/// — callers should pre-deduplicate and canonicalize paths.
pub fn write_pair_commands<W: Write>(
    pairs: &[(String, String)],
    cfg: &JoblistEmitConfig,
    writer: &mut W,
) -> Result<()> {
    for (query, target) in pairs {
        let query_stem = Path::new(query)
            .file_stem()
            .and_then(|s| s.to_str())
            .context("Query path has no file stem")?;
        let target_stem = Path::new(target)
            .file_stem()
            .and_then(|s| s.to_str())
            .context("Target path has no file stem")?;
        let output = cfg
            .output_dir
            .join(format!("{query_stem}_vs_{target_stem}.paf"));

        writeln!(
            writer,
            "{bin} {query} {target} --output-file {out} --paf --threads {t}{extra}",
            bin = cfg.sweepga_bin.display(),
            query = query,
            target = target,
            out = output.display(),
            t = cfg.threads,
            extra = if cfg.extra_flags.is_empty() {
                String::new()
            } else {
                format!(" {}", cfg.extra_flags.join(" "))
            },
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_one_line_per_pair() {
        let pairs = vec![
            ("a.fa".to_string(), "b.fa".to_string()),
            ("c.fa".to_string(), "d.fa".to_string()),
        ];
        let cfg = JoblistEmitConfig {
            sweepga_bin: Path::new("/usr/bin/sweepga"),
            output_dir: Path::new("/out"),
            threads: 4,
            extra_flags: &["--scaffold-jump".to_string(), "50k".to_string()],
        };
        let mut buf = Vec::new();
        write_pair_commands(&pairs, &cfg, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            "/usr/bin/sweepga a.fa b.fa --output-file /out/a_vs_b.paf --paf --threads 4 --scaffold-jump 50k"
        );
        assert_eq!(
            lines[1],
            "/usr/bin/sweepga c.fa d.fa --output-file /out/c_vs_d.paf --paf --threads 4 --scaffold-jump 50k"
        );
    }

    #[test]
    fn empty_pairs_produces_empty_output() {
        let cfg = JoblistEmitConfig {
            sweepga_bin: Path::new("sweepga"),
            output_dir: Path::new("."),
            threads: 1,
            extra_flags: &[],
        };
        let mut buf = Vec::new();
        write_pair_commands(&[], &cfg, &mut buf).unwrap();
        assert!(buf.is_empty());
    }
}
