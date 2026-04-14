# SweepGA

Fast all-vs-all genome alignment with plane-sweep filtering and scaffold-aware
chaining. Wraps [FastGA](https://github.com/thegenemyers/FASTGA) (default) and
[wfmash](https://github.com/waveygang/wfmash), then applies a plane sweep to
drop near-duplicate mappings and (optionally) chains nearby alignments into
syntenic scaffolds.

## What it does

- **Align** — one or more FASTA files → PAF or `.1aln`. Uses FastGA or
  wfmash; PanSN-aware; scales k-mer frequency to the cohort size
  automatically.
- **Filter** — take any aligner's PAF/1aln and apply the same plane-sweep
  + scaffold pipeline.
- **Scale** — built-in batching (`--batch-bytes`, `--max-disk`) partitions
  large cohorts into memory/disk-bounded batches.
- **Sparsify** — optionally pre-select which pairs to align (random, tree,
  giant-component, wfmash-density).

Ships two binaries: **`sweepga`** (the aligner/filter) and **`alnstats`**
(a small PAF/1aln statistics tool).

## Install

```bash
# From crates.io
cargo install sweepga

# From source
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo install --force --path .
```

Requires Rust ≥ 1.70.

On mixed-glibc systems (e.g. Debian + Guix) a plain `cargo build` may fail
with `__pthread_barrier_wait@GLIBC_PRIVATE` link errors. Use
`./scripts/build-clean.sh --install` — see
[docs/BUILD-NOTES.md](docs/BUILD-NOTES.md).

## Quick start

```bash
# Self-alignment of one FASTA (FastGA, default)
sweepga genome.fa.gz --output-file aln.paf

# Use wfmash instead
sweepga --wfmash genome.fa.gz --output-file aln.paf

# Pairwise alignment of two FASTAs
sweepga target.fa query.fa --output-file aln.paf

# Filter an existing PAF (any aligner)
sweepga alignments.paf --output-file filtered.paf

# Convert .1aln → PAF
sweepga alignments.1aln --paf > out.paf

# AGC archive input (extract samples, align, filter)
sweepga --agc-samples HG002,HG005 genomes.agc --output-file aln.paf
```

`sweepga` auto-detects input type (FASTA vs. PAF/1aln) and picks a
reasonable pipeline. Run `sweepga --help` for the full flag list.

## Defaults

The defaults are tuned for pangenome alignment (matches impg's historical
settings):

| Flag | Default | Meaning |
|---|---|---|
| `--aligner` | `fastga` | Aligner backend. Use `--wfmash` or `--aligner wfmash` to switch. |
| `--num-mappings` | `many:many` | Pre-scaffold plane-sweep: keep all mappings per query/target. |
| `--scaffold-jump` | `50k` | Scaffolding **is enabled by default**; chains mappings within a 50kb gap. Pass `--scaffold-jump 0` to disable. |
| `--scaffold-mass` | `10k` | Drop scaffold chains shorter than 10kb. |
| `--scaffold-filter` | `many:many` | Keep all non-overlapping scaffolds. Use `1:1` for the aggressive (best-per-chromosome-pair) filter. |
| `--overlap` | `0.95` | Drop mappings whose overlap with a better-scoring one exceeds 95%. |
| `--scoring` | `log-length-ani` | Plane-sweep scoring function. |
| `--threads` | `8` | Parallelism. |

Adaptive behavior: when the input average sequence length is known
(FASTA with `.fai`), `--scaffold-mass` and `--scaffold-jump` are clamped
to the sequence-length scale so short inputs (e.g. pangenome-window
excerpts) don't get filtered empty. Pass `--no-adaptive-scaffolds` to
turn clamping off.

## Scaffolding pipeline

When `--scaffold-jump > 0` (the default), sweepga runs:

1. **Pre-scaffold plane sweep** per (query, target) chromosome pair
   using `--num-mappings` / `--scoring` / `--overlap`.
2. **Scaffold creation** — union-find merges mappings within
   `--scaffold-jump` on both axes; chains shorter than
   `--scaffold-mass` are dropped.
3. **Scaffold plane sweep** using `--scaffold-filter` (default
   `many:many`); `1:1` keeps the single best scaffold per
   chromosome-pair.
4. **Rescue** (`--scaffold-dist > 0`) — recover mappings within the
   Euclidean distance of a kept scaffold anchor. Off by default.

Pass `--no-filter` to skip everything and just emit the raw aligner
output (or the raw input PAF).

## Sparsification

`--sparsify` picks which pairs actually get aligned (for FASTA input)
or how densely wfmash reports mappings:

| Value | Behavior |
|---|---|
| `none` (default) | Align all pairs, report all mappings |
| `auto` | Heuristic pick based on cohort size |
| `<f>` or `random:<f>` | Random `<f>` fraction of pairs |
| `giant:<p>` / `connectivity:<p>` | Giant-component connectivity guarantee |
| `tree:<near>[:<far>[:<random>]]` | Tree-sampled neighbor pairs |
| `wfmash:auto` / `wfmash:<f>` | wfmash mapping density (`-x` flag); auto = `ln(n)/n*10` |

`--sparsify-pairs` is a separate knob that only drives pre-alignment
pair selection (same grammar, no wfmash variant).

## PanSN

Sequence names like `SAMPLE#HAPLOTYPE#CHR` are treated pair-wise
per-chromosome within each genome pair. Prefix delimiter is `#`.

## Output

- `.paf` — standard PAF text.
- `.1aln` — compact binary [ONE](https://github.com/thegenemyers/ONEcode)
  format (50–70% smaller, preserves edit distances). Pass `--paf` to
  convert to PAF on read.

Output path is chosen by `--output-file` (extension auto-detected) or
the explicit `--paf` / `--1aln` flags.

## Alnstats

```bash
alnstats alignments.paf              # summary
alnstats raw.paf filtered.paf        # before/after
alnstats alignments.paf -d           # per-genome-pair breakdown
```

## License

MIT.
