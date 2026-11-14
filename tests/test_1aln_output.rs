//! `--1aln` output holds unchanged records of the input.
use std::collections::HashSet;
use std::path::Path;
use std::process::Command;

/// PAF lines of a .1aln file, from ALNtoPAF.
fn paf_lines(aln: &Path) -> Vec<String> {
    let bin = sweepga::binary_paths::get_embedded_binary_path("ALNtoPAF").unwrap();
    let out = Command::new(bin)
        .arg("-T1")
        .arg(aln.with_extension(""))
        .current_dir(aln.parent().unwrap())
        .output()
        .unwrap();
    assert!(out.status.success(), "ALNtoPAF: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap().lines().map(str::to_string).collect()
}

/// Run sweepga and return its stdout.
fn sweepga(args: &[&str]) -> Vec<u8> {
    let out = Command::new(env!("CARGO_BIN_EXE_sweepga")).args(args).output().unwrap();
    assert!(out.status.success(), "sweepga {args:?}: {}", String::from_utf8_lossy(&out.stderr));
    out.stdout
}

#[test]
fn fasta_to_1aln_output_converts_back() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("in.fa");
    std::fs::copy("data/B-3106.fa", &fasta).unwrap();
    let kept = dir.path().join("kept.1aln");
    std::fs::write(&kept, sweepga(&[fasta.to_str().unwrap(), "--1aln"])).unwrap();
    let lines = paf_lines(&kept);
    assert!(!lines.is_empty());
    for line in &lines {
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 12, "{line}");
        assert_ne!(f[0], f[5], "self mapping kept: {line}");
    }
}

#[test]
fn aln_input_keeps_every_record_when_unfiltered() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("in.fa");
    std::fs::copy("data/B-3106.fa", &fasta).unwrap();
    let produced = fastga_rs::runner::Orchestrator::new_simple(1)
        .align(&fasta, &fasta)
        .unwrap();
    let raw = dir.path().join("raw.1aln");
    std::fs::copy(produced.trim(), &raw).unwrap();
    let raw_paf: HashSet<String> = paf_lines(&raw).iter().map(|l| coords(l)).collect();
    assert!(raw_paf.len() > 10, "raw alignments: {}", raw_paf.len());

    let unfiltered = ["--num-mappings", "many:many", "--self", "--scaffold-jump", "0"];
    let mut args = vec![raw.to_str().unwrap()];
    args.extend(unfiltered);
    let out = String::from_utf8(sweepga(&args)).unwrap();
    let got: HashSet<String> = out.lines().map(coords).collect();
    assert_eq!(got, raw_paf);
}

/// Names, strand and coordinates of a PAF line.
fn coords(line: &str) -> String {
    let f: Vec<&str> = line.split('\t').collect();
    [f[0], f[2], f[3], f[4], f[5], f[7], f[8]].join("\t")
}
