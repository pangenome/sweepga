# SweepGA

Fast genome alignment with sophisticated filtering. Wraps FastGA aligner and applies plane sweep filtering to keep the best non-overlapping alignments.

## What it does

SweepGA can either:
1. **Align FASTA files directly** using integrated FastGA (supports .fa.gz)
2. **Filter existing PAF alignments** from any aligner (wfmash, minimap2, etc.)

By default, it applies 1:1 plane sweep filtering to keep the single best mapping per query-target chromosome pair. Optionally, you can enable scaffolding (`-j > 0`) to chain nearby mappings and rescue distant features.

## Tools Included

This package includes two binaries:

- **`sweepga`** - Genome alignment and filtering tool
- **`alnstats`** - Alignment statistics and validation tool

Use `alnstats` to verify filtering results:
```bash
# Show statistics for a PAF file
alnstats alignments.paf

# Compare before/after filtering
alnstats raw.paf filtered.paf

# Detailed per-genome-pair breakdown
alnstats alignments.paf -d
```

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
# Self-alignment with 1:1 filtering
sweepga genome.fa.gz > output.paf

# Pairwise alignment (target, query order)
sweepga target.fa query.fa > output.paf

# With 2 threads
sweepga genome.fa.gz -t 2 > output.paf
```

### PAF filtering (pipe from stdin)

```bash
# Default: 1:1 plane sweep filtering
cat alignments.paf | sweepga > filtered.paf

# Keep best mapping per query only (1:∞)
cat alignments.paf | sweepga -n 1 > filtered.paf

# No filtering, just pass through
cat alignments.paf | sweepga -n many > output.paf
```

### Reading PAF files directly

```bash
# Read from file instead of stdin
sweepga alignments.paf > filtered.paf
```

### Example: 8 yeast genomes

```bash
# Direct alignment and filtering in one step
sweepga data/scerevisiae8.fa.gz > scerevisiae8.paf

# Result: ~26K mappings (1:1 filtered)
# - Each genome pair gets best alignment per chromosome pair
# - Self-mappings excluded by default (use --self to include)
```

## Parameters

### Core Filtering Parameters

**`-n/--num-mappings`** - Plane sweep filtering (default: `1:1`)
- `"1:1"` - Keep best mapping on both query and target axes (default)
- `"1"` - Keep best mapping per query position only
- `"many"` - No filtering, keep all mappings
- `"M:N"` - Keep top M per query, top N per target

**`-o/--overlap`** - Maximum overlap ratio (default: 0.95)
- Mappings with >95% overlap with a better-scoring mapping are removed

**`-t/--threads`** - Number of threads (default: 8)

**`--self`** - Include self-mappings (excluded by default)

### Optional: Scaffolding Mode

Scaffolding is **disabled by default** (`-j 0`). When enabled, it chains nearby mappings and rescues distant features.

**`-j/--scaffold-jump`** - Enable scaffolding by setting gap distance (default: 0 = disabled)
- Example: `-j 10k` merges mappings within 10kb gaps into chains

**`-s/--scaffold-mass`** - Minimum scaffold length (default: 10k)
- Only chains ≥ this length become scaffold anchors

**`-d/--scaffold-dist`** - Maximum rescue distance (default: 20k)
- Mappings within this distance of scaffolds are rescued

**`-m/--scaffold-filter`** - Scaffold chain filtering (default: `many`)
- `"1:1"` - Keep best scaffold per chromosome pair
- `"many"` - Keep all scaffolds (default)

When scaffolding is enabled, output includes:
- `ch:Z:chain_N` - Chain ID for grouped mappings
- `st:Z:scaffold` - Member of a scaffold chain
- `st:Z:rescued` - Rescued due to proximity to scaffold

## Example: Scaffolding Mode

Scaffolding groups nearby mappings into chains and rescues distant features:

```bash
# Enable scaffolding with 10kb gap merging
cat alignments.paf | sweepga -j 10k > scaffolded.paf

# Tighter scaffolding with 1:1 filtering
cat alignments.paf | sweepga -j 10k -m 1:1 > strict.paf

# Adjust rescue distance
cat alignments.paf | sweepga -j 10k -d 50k > permissive.paf
```

Scaffolding is useful for:
- Grouping fragmented alignments from draft assemblies
- Rescuing small features near main syntenic blocks
- Creating cleaner visualizations by chaining related mappings

The rescue mechanism is strand-agnostic: reverse-strand mappings can be rescued by forward-strand scaffolds if within the distance threshold.

## How It Works

### Plane Sweep Filtering (Default)

The plane sweep algorithm operates per query-target chromosome pair:

1. **Sort** mappings by query position
2. **Score** each mapping: `identity × log(block_length)` (matches wfmash)
3. **Sweep** left-to-right, keeping best mappings based on `-n` setting:
   - `1:1`: Keep single best mapping per position on both query and target
   - `1`: Keep best mapping per query position (multiple targets allowed)
   - `many`: Keep all non-overlapping mappings
4. **Filter** mappings with >95% overlap (configurable with `-o`)

### Optional: Scaffolding Mode (`-j > 0`)

When scaffolding is enabled:

1. **Chain** nearby mappings using union-find:
   - Merge if gap < `-j` on both query and target
   - Allow small overlaps (up to gap/5)
2. **Filter** chains by minimum length (`-s`)
3. **Apply** scaffold filtering (`-m`) to chains
4. **Rescue** mappings within distance (`-d`) of scaffold anchors
   - Distance = Euclidean distance between mapping centers
   - Checked per chromosome pair only

Output includes PAF tags:
- `ch:Z:chain_N` - Chain identifier (when scaffolding enabled)
- `st:Z:scaffold` - Member of scaffold chain
- `st:Z:rescued` - Rescued by proximity to scaffold

## Citation

SweepGA: PAF filtering using Euclidean distance to scaffold anchors
https://github.com/pangenome/sweepga