# SweepGA

Fast genome alignment filtering using plane sweep algorithm with scaffold-based rescue.

## Overview

SweepGA is a Rust implementation of sophisticated PAF (Pairwise mApping Format) filtering algorithms inspired by wfmash. It applies plane sweep filtering to genome alignments and uses a scaffold-based rescue mechanism to retain biologically relevant mappings near high-confidence anchors.

## Key Features

- **Plane sweep filtering**: Efficiently removes overlapping/weaker alignments
- **Scaffold-based rescue**: Retains mappings within a specified distance of high-confidence scaffold anchors
- **Cross-strand rescue**: Mappings on opposite strands can be rescued if near scaffold anchors
- **Union-find chaining**: Groups nearby mappings into coherent chains
- **Multiple filtering modes**: 1:1, 1:∞ (default), and N:N filtering strategies
- **Stream processing**: Efficient handling of large PAF files

## Installation

### Prerequisites
- Rust 1.70 or later
- Cargo

### Build from source
```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo build --release
```

The binary will be at `target/release/sweepga`.

## Quick Demo

Using FASTGA with S. cerevisiae data:

```bash
# First, run FASTGA to generate initial mappings
cd deps/FASTGA
cargo build --release
./target/release/fastga -t 8 ../../data/scerevisiae8.fa > ../../scerevisiae8.paf

# Apply SweepGA filtering with default parameters
cd ../..
./target/release/sweepga < scerevisiae8.paf > scerevisiae8.filtered.paf

# Check the results
echo "Input mappings: $(wc -l < scerevisiae8.paf)"
echo "Filtered mappings: $(wc -l < scerevisiae8.filtered.paf)"
echo "Scaffold anchors: $(grep -c 'st:Z:scaffold' scerevisiae8.filtered.paf)"
echo "Rescued mappings: $(grep -c 'st:Z:rescued' scerevisiae8.filtered.paf)"
```

## Usage

```bash
sweepga [OPTIONS] < input.paf > output.paf
```

### Key Options

#### Scaffold Parameters
- `-S, --scaffold-mass <N>`: Minimum scaffold length [default: 10000]
- `-j, --scaffold-jump <N>`: Maximum gap for chaining into scaffolds [default: 100000]
- `-D, --scaffold-dist <N>`: Maximum distance from scaffold for rescue [default: 100000]

#### Filtering Parameters
- `-m, --mode <MODE>`: Filtering mode: "1:1", "1:∞"/"map" (default), or "N:N"
- `-n, --mappings <N>`: Maximum mappings per segment (for 1:N mode)
- `-l, --block-length <N>`: Minimum block length [default: 0]
- `-O, --overlap <RATIO>`: Maximum overlap ratio [default: 0.95]

#### Other Options
- `-c, --chain-jump <N>`: Chain gap distance [default: 2000]
- `--no-merge`: Keep fragment mappings separate (don't merge)
- `--self`: Keep self-mappings (excluded by default)
- `--no-filter`: Disable all filtering
- `--quiet`: Suppress progress output

### Example Commands

```bash
# Basic filtering with default parameters (1:∞ mode)
sweepga < input.paf > output.paf

# Strict 1:1 filtering
sweepga -m 1:1 < input.paf > output.paf

# Adjust scaffold parameters for tighter clustering
sweepga -S 5000 -j 50000 -D 50000 < input.paf > output.paf

# Keep all mappings (no filtering) but still annotate scaffolds
sweepga --no-filter < input.paf > output.paf

# Debug scaffolds only (outputs synthetic PAF for scaffold chains)
sweepga --scaffolds-only < input.paf > scaffolds.paf
```

## Algorithm Details

### 1. Scaffold Identification
- Merges nearby mappings into chains using union-find
- Identifies chains exceeding minimum scaffold length
- Applies plane sweep to remove overlapping scaffolds

### 2. Rescue Mechanism
- Calculates Euclidean distance from non-scaffold mappings to scaffold anchors
- Rescues mappings within threshold distance
- Works across strands (forward/reverse mappings can be rescued by any scaffold)

### 3. Output Annotation
Each mapping in the output includes a status tag:
- `st:Z:scaffold` - Part of a scaffold anchor
- `st:Z:rescued` - Rescued due to proximity to scaffold
- `st:Z:unassigned` - Passed filters but not scaffolded/rescued

## Comparison with wfmash

SweepGA implements the same core filtering algorithms as wfmash's `filterByScaffolds`:
- Compatible chain merging with union-find
- Same distance calculations for rescue
- Equivalent plane sweep implementation

Key differences:
- Rust implementation (vs C++ for wfmash)
- Stream processing of PAF records
- Simplified interface focused on filtering

## Development

### Running tests
```bash
cargo test
```

### Building with debug symbols
```bash
cargo build
```

## Citation

If you use SweepGA in your research, please cite:

```
SweepGA: Fast genome alignment filtering with scaffold-based rescue
https://github.com/pangenome/sweepga
```

## License

MIT

## Authors

- Erik Garrison
- Implementations based on algorithms from wfmash by Andrea Guarracino and the Pangenome Consortium