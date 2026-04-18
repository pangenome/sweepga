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

/// Configuration for emitting a wfmash-per-haplotype-pair joblist.
///
/// Targets the "one PanSN-named FASTA in, one wfmash invocation per
/// selected `SAMPLE#HAPLOTYPE` pair out" flow: wfmash filters queries and
/// targets to those haplotype prefixes via `-Q` / `-T`, aligning only the
/// contigs that match, without the user having to pre-split the input
/// into per-haplotype FASTAs.
pub struct WfmashPansnEmitConfig<'a> {
    /// Path to the PanSN-named FASTA input (used for both query and
    /// target sides of each wfmash invocation).
    pub fasta: &'a Path,
    /// Output directory for per-pair PAFs.
    pub output_dir: &'a Path,
    /// Number of threads to pass via `-t`.
    pub threads: usize,
    /// Optional block-length filter (`-l`). 0 means "omit the flag".
    pub block_length: u64,
}

/// Replace filesystem-hostile characters (notably PanSN's `#`) so a
/// `SAMPLE#HAPLOTYPE` key is safe to embed in a filename.
fn sanitize_for_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '/' | '\\' | '#' | ':' | ' ' | '\t' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            other => other,
        })
        .collect()
}

/// Emit one `wfmash -T <target_hap> -Q <query_hap>` command per PanSN
/// haplotype pair.
///
/// `hap_pairs` are `(target_haplotype_key, query_haplotype_key)` string
/// tuples — e.g. `("HG01106#1", "HG00733#2")`. The output filename for
/// each pair is `<output_dir>/<target>_vs_<query>.paf`, with `#` replaced
/// by `_`. When `block_length` is non-zero the emitted command includes
/// `-l <block_length>`; this is the only filtering flag that wfmash's
/// `-T`/`-Q` workflow honors.
pub fn write_wfmash_pansn_commands<W: Write>(
    hap_pairs: &[(String, String)],
    cfg: &WfmashPansnEmitConfig,
    writer: &mut W,
) -> Result<()> {
    for (target_hap, query_hap) in hap_pairs {
        let output = cfg.output_dir.join(format!(
            "{}_vs_{}.paf",
            sanitize_for_filename(target_hap),
            sanitize_for_filename(query_hap),
        ));
        let mut cmd = format!("wfmash -t {}", cfg.threads);
        if cfg.block_length > 0 {
            cmd.push_str(&format!(" -l {}", cfg.block_length));
        }
        cmd.push_str(&format!(" -T {} -Q {}", target_hap, query_hap));
        cmd.push_str(&format!(" {}", cfg.fasta.display()));
        cmd.push_str(&format!(" > {}", output.display()));
        writeln!(writer, "{}", cmd)?;
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

    #[test]
    fn wfmash_pansn_emits_T_Q_per_hap_pair() {
        let pairs = vec![
            ("HG01106#1".to_string(), "HG00733#2".to_string()),
            ("SGDref#0".to_string(), "HG01106#1".to_string()),
        ];
        let cfg = WfmashPansnEmitConfig {
            fasta: Path::new("/data/pangenome.fa"),
            output_dir: Path::new("/out"),
            threads: 8,
            block_length: 0,
        };
        let mut buf = Vec::new();
        write_wfmash_pansn_commands(&pairs, &cfg, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            "wfmash -t 8 -T HG01106#1 -Q HG00733#2 /data/pangenome.fa > /out/HG01106_1_vs_HG00733_2.paf"
        );
        assert_eq!(
            lines[1],
            "wfmash -t 8 -T SGDref#0 -Q HG01106#1 /data/pangenome.fa > /out/SGDref_0_vs_HG01106_1.paf"
        );
    }

    #[test]
    fn wfmash_pansn_includes_block_length_when_set() {
        let pairs = vec![("A#0".to_string(), "B#0".to_string())];
        let cfg = WfmashPansnEmitConfig {
            fasta: Path::new("pg.fa"),
            output_dir: Path::new("."),
            threads: 4,
            block_length: 20_000,
        };
        let mut buf = Vec::new();
        write_wfmash_pansn_commands(&pairs, &cfg, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("-l 20000"));
        assert!(s.contains("-T A#0 -Q B#0"));
    }
}
