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

/// One PanSN haplotype-pair wfmash job: the `(target_hap, query_hap)`
/// keys and the target / query FASTA paths the command should run on.
///
/// Single-file callers set `target_fasta == query_fasta`; multi-file
/// callers point each side at whichever FASTA contains the relevant
/// haplotype's contigs. When the two paths are equal the emitted command
/// omits the query file (wfmash self-map form).
pub struct WfmashPansnJob {
    pub target_hap: String,
    pub query_hap: String,
    pub target_fasta: std::path::PathBuf,
    pub query_fasta: std::path::PathBuf,
}

/// Shared emission config for [`write_wfmash_pansn_commands`].
pub struct WfmashPansnEmitConfig<'a> {
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

/// Emit one `wfmash -T <target_hap> -Q <query_hap> target.fa [query.fa]`
/// command per PanSN haplotype-pair job.
///
/// Works for both single-FASTA PanSN inputs (each job's target and query
/// FASTA are the same physical file; wfmash self-maps with `-T`/`-Q`
/// filtering) and multi-FASTA PanSN inputs (jobs point at whichever
/// FASTA holds the haplotype's contigs; the second file is passed
/// through when it differs from the first).
///
/// Output filename is `<output_dir>/<target>_vs_<query>.paf`, with
/// PanSN `#` replaced by `_` for filesystem safety.
pub fn write_wfmash_pansn_commands<W: Write>(
    jobs: &[WfmashPansnJob],
    cfg: &WfmashPansnEmitConfig,
    writer: &mut W,
) -> Result<()> {
    for job in jobs {
        let output = cfg.output_dir.join(format!(
            "{}_vs_{}.paf",
            sanitize_for_filename(&job.target_hap),
            sanitize_for_filename(&job.query_hap),
        ));
        let mut cmd = format!("wfmash -t {}", cfg.threads);
        if cfg.block_length > 0 {
            cmd.push_str(&format!(" -l {}", cfg.block_length));
        }
        cmd.push_str(&format!(" -T {} -Q {}", job.target_hap, job.query_hap));
        cmd.push_str(&format!(" {}", job.target_fasta.display()));
        if job.query_fasta != job.target_fasta {
            cmd.push_str(&format!(" {}", job.query_fasta.display()));
        }
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

    fn pj(tgt: &str, qry: &str, tgt_fa: &str, qry_fa: &str) -> WfmashPansnJob {
        WfmashPansnJob {
            target_hap: tgt.to_string(),
            query_hap: qry.to_string(),
            target_fasta: std::path::PathBuf::from(tgt_fa),
            query_fasta: std::path::PathBuf::from(qry_fa),
        }
    }

    #[test]
    fn wfmash_pansn_emits_T_Q_per_hap_pair_single_fasta() {
        // Same file on both sides → emit with self-map form (one file arg).
        let jobs = vec![
            pj("HG01106#1", "HG00733#2", "/data/pangenome.fa", "/data/pangenome.fa"),
            pj("SGDref#0", "HG01106#1", "/data/pangenome.fa", "/data/pangenome.fa"),
        ];
        let cfg = WfmashPansnEmitConfig {
            output_dir: Path::new("/out"),
            threads: 8,
            block_length: 0,
        };
        let mut buf = Vec::new();
        write_wfmash_pansn_commands(&jobs, &cfg, &mut buf).unwrap();
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
    fn wfmash_pansn_emits_both_files_when_distinct() {
        // Different files on each side → emit both, in wfmash order
        // (target.fa query.fa).
        let jobs = vec![pj("HG01106#1", "HG00733#2", "/data/hg01106.fa", "/data/hg00733.fa")];
        let cfg = WfmashPansnEmitConfig {
            output_dir: Path::new("/out"),
            threads: 4,
            block_length: 0,
        };
        let mut buf = Vec::new();
        write_wfmash_pansn_commands(&jobs, &cfg, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert_eq!(
            s.trim_end(),
            "wfmash -t 4 -T HG01106#1 -Q HG00733#2 /data/hg01106.fa /data/hg00733.fa > /out/HG01106_1_vs_HG00733_2.paf"
        );
    }

    #[test]
    fn wfmash_pansn_includes_block_length_when_set() {
        let jobs = vec![pj("A#0", "B#0", "pg.fa", "pg.fa")];
        let cfg = WfmashPansnEmitConfig {
            output_dir: Path::new("."),
            threads: 4,
            block_length: 20_000,
        };
        let mut buf = Vec::new();
        write_wfmash_pansn_commands(&jobs, &cfg, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("-l 20000"));
        assert!(s.contains("-T A#0 -Q B#0"));
    }
}
