# SweepGA

Fast genome alignment with sophisticated filtering. Wraps FastGA aligner and applies wfmash's plane sweep algorithm with optional synteny-based scaffolding and distance-based rescue.

## What it does

SweepGA can either:
1. **Align FASTA files directly** using integrated FastGA (supports .fa.gz)
2. **Filter existing PAF alignments** from any aligner (wfmash, minimap2, etc.)

By default, it chains mappings into syntenic scaffolds, applies 1:1 filtering to keep best scaffolds per chromosome pair, then rescues nearby mappings. This extracts clean syntenic alignments from noisy all-vs-all mappings.

## Installation

Requires Rust 1.70+. Clone and install:

```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo install --force --path .
```

## Basic Usage

### Direct alignment (FASTA input)

```bash
# Self-alignment
sweepga genome.fa.gz -o output.paf

# Pairwise alignment
sweepga target.fa query.fa -o output.paf

# With custom parameters
sweepga genome.fa.gz -o output.paf -j 20k -d 50k
```

### PAF filtering (existing alignments)

```bash
# Default: scaffolding + 1:1 filter + rescue
sweepga -i alignments.paf -o filtered.paf

# Plane sweep only (no scaffolding)
sweepga -i alignments.paf -o filtered.paf -j 0

# Custom scaffold parameters
sweepga -i alignments.paf -o filtered.paf -j 50k -s 20k
```

### Complete example workflow

```bash
# Method 1: Direct alignment with integrated FastGA
sweepga data/scerevisiae8.fa.gz -o data/scerevisiae8.filtered.paf
# Runs FastGA alignment then applies filtering
# Result: ~17,800 mappings from ~30K raw alignments

# Method 2: Filter existing PAF
fastga -T8 -pafx data/scerevisiae8.fa > data/scerevisiae8.raw.paf
sweepga -i data/scerevisiae8.raw.paf -o data/scerevisiae8.filtered.paf

# Check the breakdown
grep -c "st:Z:scaffold" data/scerevisiae8.filtered.paf  # Scaffold anchors
grep -c "st:Z:rescued" data/scerevisiae8.filtered.paf   # Rescued mappings
```

### Default parameters

- **No pre-filtering** (`-n many`): All mappings participate in scaffold building
- **Scaffold creation** (`-j 10k`): Merge mappings within 10kb gaps
- **Min scaffold** (`-s 10k`): Scaffolds must be ≥10kb
- **1:1 scaffold filter** (`-m 1:1`): Keep best scaffold per chromosome pair
- **Rescue distance** (`-d 20k`): Rescue mappings within 20kb of scaffolds

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

`-d/--scaffold-dist` sets the maximum Euclidean distance for rescue (default 20k). Mappings further than this from any scaffold anchor are discarded.

`-m/--scaffold-filter` controls scaffold filtering (default "1:1"). Options: "1:1", "1" (1:many), "many" (no filter).

### Chain Annotations

When scaffolding is enabled, mappings are annotated with:
- `ch:Z:chain_N` - Chain ID for grouped mappings
- `st:Z:{scaffold|rescued|unassigned}` - Status of each mapping

## Example: varying rescue distance

The rescue distance dramatically affects output. Using a yeast chromosome V alignment:

```bash
# Very tight - only mappings within 10kb of scaffolds
sweepga -d 10k < chrV.paf > chrV.strict.paf
# Result: 32 mappings retained

# Default (20kb) - biological features near main alignments
sweepga < chrV.paf > chrV.default.paf
# Result: 45 mappings retained

# Permissive - may include distant paralogs
sweepga -d 100k < chrV.paf > chrV.loose.paf
# Result: 71 mappings retained
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