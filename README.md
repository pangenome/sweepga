# SweepGA

Filters highly sensitive whole genome alignments using plane sweep algorithm and optional synteny-based scaffolding with distance-based rescue.

## What it does

By default, SweepGA chains mappings into syntenic scaffolds and applies plane sweep filtering, then rescues nearby mappings. With `-j 0`, it performs only plane sweep filtering without scaffolding. This extracts clean synteny alignments from noisy all-vs-all mappings.

## Installation

Requires Rust 1.70+. Clone and install:

```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo install --force --path .
```

## Basic Usage

SweepGA filters PAF alignments from tools like wfmash, FASTGA, or minimap2. Basic usage:

```bash
# Default: scaffolding + plane sweep + rescue (with -j 10000)
sweepga -i alignments.paf -o filtered.paf

# Plane sweep only (no scaffolding)
sweepga -i alignments.paf -o filtered.paf -j 0

# Custom scaffold parameters
sweepga -i alignments.paf -o filtered.paf -j 50000 -s 5000
```

### Complete example workflow

```bash
# Generate all-vs-all alignments with FASTGA (use -T8 for 8 threads, -pafx for PAF with extended CIGAR)
fastga -T8 -pafx data/scerevisiae8.fa > data/scerevisiae8.raw.paf
# Generates 50,959 raw mappings

# The raw output contains all discovered mappings with CIGAR strings
head -n1 data/scerevisiae8.raw.paf
# SGDref#1#chrI  230218  0  2641  -  SGDref#1#chrIV  1531933  1522805  1525422  2341  2692  255  dv:f:.1135  df:i:351  cg:Z:6=1D2=1X1I6=...

# Apply scaffold-based filtering (default: -j 10000, -s 10000)
sweepga -i data/scerevisiae8.raw.paf -o data/scerevisiae8.filtered.paf
# Reduces to 27,940 mappings: 336 scaffolds + 27,604 rescued

# Check the breakdown
grep -c "st:Z:scaffold" data/scerevisiae8.filtered.paf  # 336 scaffold anchors
grep -c "st:Z:rescued" data/scerevisiae8.filtered.paf   # 27,604 rescued mappings
```

The default parameters work well for most eukaryotic genomes: scaffold mass of 10kb identifies major syntenic blocks, scaffold jump of 10kb allows chaining across gene boundaries, and rescue distance of 100kb captures local rearrangements and smaller homologous features near the main alignments.

## Parameters

### Plane Sweep Filtering (default behavior)

`-n/--num-mappings` controls how many mappings to keep per query position (default "1:1"):
- `"1:1"` - Keep best mapping on both query and target axes
- `"1"` or `"1:∞"` - Keep best on query axis only
- `"many:many"` or `"∞:∞"` - Keep all non-overlapping mappings
- `"M:N"` - Keep top M per query, top N per target

### Scaffolding Parameters

`-j/--scaffold-jump` sets the maximum gap for merging mappings into scaffold chains (default 10k). Set to 0 to disable scaffolding and use plane sweep only.

`-s/--scaffold-mass` sets the minimum length for a chain to be considered a scaffold anchor (default 10k). Only used when scaffolding is enabled (`-j > 0`).

`-D/--scaffold-dist` sets the maximum Euclidean distance for rescue (default 100000). Mappings further than this from any scaffold anchor are discarded.

### Chain Annotations

When scaffolding is enabled, mappings are annotated with:
- `ch:Z:chain_N` - Chain ID for grouped mappings
- `st:Z:{scaffold|rescued|unassigned}` - Status of each mapping

## Example: varying rescue distance

The rescue distance dramatically affects output. Using a yeast chromosome V alignment:

```bash
# Very tight - only mappings within 10kb of scaffolds
./target/release/sweepga -D 10000 < chrV.paf > chrV.strict.paf
# Result: 32 mappings retained

# Default - biological features near main alignments
./target/release/sweepga -D 100000 < chrV.paf > chrV.default.paf
# Result: 71 mappings retained

# Permissive - may include distant paralogs
./target/release/sweepga -D 300000 < chrV.paf > chrV.loose.paf
# Result: 115 mappings retained
```

The scaffold jump parameter has less effect when scaffolds are already well-separated, which is common in finished genomes. The rescue mechanism is strand-agnostic: a reverse strand mapping 50kb from a forward strand scaffold will be rescued if within the distance threshold.

## Algorithm details

The implementation follows wfmash's filterByScaffolds approach. Mappings are first grouped by query-target pair and sorted by query position. Union-find merges mappings where both query and target gaps are below the threshold, allowing small overlaps up to gap/5. The resulting chains are filtered by length to identify scaffolds.

Plane sweep operates independently per query sequence, sweeping across the query axis and keeping the best scoring mappings at each position. The scoring function matches wfmash: identity × log(block_length). For the default "1" mode, this keeps all non-overlapping mappings at each query position; for "1:1" it keeps only the single best mapping per query-target pair. The overlap threshold parameter (default 0.95) filters mappings that are mostly contained within better scoring ones.

The rescue phase calculates Euclidean distance from each non-scaffold mapping to scaffold anchors on the same chromosome pair. Mappings are grouped by (query_chromosome, target_chromosome) and sorted by query position for efficient processing. The distance uses mapping center points in 2D space where axes are query and target positions. Any mapping within the distance threshold of at least one anchor is rescued, regardless of strand orientation.

Output preserves the original PAF records with an added tag indicating filter status: st:Z:scaffold for scaffold anchors, st:Z:rescued for rescued mappings. This allows downstream tools to distinguish high-confidence syntenic anchors from nearby features.

## Citation

SweepGA: PAF filtering using Euclidean distance to scaffold anchors
https://github.com/pangenome/sweepga