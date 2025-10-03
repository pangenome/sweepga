# SweepGA

Fast genome alignment with plane sweep filtering. Wraps FastGA aligner and applies plane sweep filtering to keep the best non-overlapping alignments.

## What it does

SweepGA can either:
1. **Align FASTA files directly** using integrated FastGA (supports .fa.gz)
2. **Filter existing PAF alignments** from any aligner (wfmash, minimap2, etc.)

By default, it applies 1:1 plane sweep filtering to keep the single best mapping per query-target chromosome pair.

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

**`-n/--num-mappings`** - Plane sweep filtering (default: `1:1`)
- `"1:1"` - Keep best mapping on both query and target axes (default)
- `"1"` - Keep best mapping per query position only
- `"many"` - No filtering, keep all mappings
- `"M:N"` - Keep top M per query, top N per target

**`-o/--overlap`** - Maximum overlap ratio (default: 0.95)
- Mappings with >95% overlap with a better-scoring mapping are removed

**`-l/--min-block-length`** - Minimum alignment block length (default: 0)

**`-i/--min-identity`** - Minimum identity threshold (0-1 fraction, 1-100%, or "aniN")

**`-t/--threads`** - Number of threads (default: 8)

**`--self`** - Include self-mappings (excluded by default)

**`-f/--no-filter`** - Disable all filtering

## How It Works

The plane sweep algorithm operates per query-target chromosome pair:

1. **Sort** mappings by query position
2. **Score** each mapping: `identity × log(block_length)` (matches wfmash)
3. **Sweep** left-to-right, keeping best mappings based on `-n` setting:
   - `1:1`: Keep single best mapping per position on both query and target
   - `1`: Keep best mapping per query position (multiple targets allowed)
   - `many`: Keep all non-overlapping mappings
4. **Filter** mappings with >95% overlap (configurable with `-o`)

## Citation

SweepGA: Fast plane sweep filtering for genome alignments
https://github.com/pangenome/sweepga