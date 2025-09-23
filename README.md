# SweepGA

Filters highly sensitive whole genome alignments to extract maximum synteny alignments using chaining, plane sweep, and distance-based rescue.

## What it does

SweepGA chains mappings into large syntenic scaffolds (>10kb), keeps the best scaffolds per query genome (1:N filtering), then rescues any mappings within 100kb of these anchors. This extracts clean synteny alignments from noisy all-vs-all mappings.

## Installation

Requires Rust 1.70+. Clone and install:

```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo install --force --path .
```

## Usage with FASTGA

FASTGA generates all-vs-all genome alignments in PAF format with full CIGAR strings. These alignments typically need filtering to remove redundant and weak mappings. Basic usage:

```bash
# Generate alignments with FASTGA and filter with SweepGA
fastga -pafm data/scerevisiae8.fa data/scerevisiae8.fa > data/scerevisiae8.raw.paf
sweepga -i data/scerevisiae8.raw.paf -o data/scerevisiae8.filtered.paf
```

### Complete example workflow

```bash
# Generate all-vs-all alignments with FASTGA (use -T8 for 8 threads, -pafx for PAF with extended CIGAR)
fastga -T8 -pafx data/scerevisiae8.fa > data/scerevisiae8.raw.paf
# Generates 50,959 raw mappings

# The raw output contains all discovered mappings with CIGAR strings
head -n1 data/scerevisiae8.raw.paf
# SGDref#1#chrI  230218  0  2641  -  SGDref#1#chrIV  1531933  1522805  1525422  2341  2692  255  dv:f:.1135  df:i:351  cg:Z:6=1D2=1X1I6=...

# Apply scaffold-based filtering (default: scaffolds >10kb, rescue within 100kb)
sweepga -i data/scerevisiae8.raw.paf -o data/scerevisiae8.filtered.paf
# Reduces to 27,940 mappings: 336 scaffolds + 27,604 rescued

# Check the breakdown
grep -c "st:Z:scaffold" data/scerevisiae8.filtered.paf  # 336 scaffold anchors
grep -c "st:Z:rescued" data/scerevisiae8.filtered.paf   # 27,604 rescued mappings
```

The default parameters work well for most eukaryotic genomes: scaffold mass of 10kb identifies major syntenic blocks, scaffold jump of 100kb allows chaining across typical intergenic distances, and rescue distance of 100kb captures local rearrangements and smaller homologous features near the main alignments.

## Parameters

The main parameters control scaffold identification and rescue distance:

`-S/--scaffold-mass` sets the minimum length for a chain to be considered a scaffold anchor (default 10000). Smaller values create more anchors but may include repetitive elements.

`-j/--scaffold-jump` sets the maximum gap for merging mappings into scaffold chains (default 100000). This controls how scattered mappings can be while still being considered part of the same scaffold.

`-D/--scaffold-dist` sets the maximum Euclidean distance for rescue (default 100000). Mappings further than this from any scaffold anchor are discarded.

`-m/--mode` controls the mapping filtering strategy: "1" or "1:∞" (default) keeps all non-overlapping mappings per query position, "1:1" keeps only the best per query-target pair, "N" or "N:N" disables plane sweep filtering entirely.

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