# SweepGA

PAF filtering using Euclidean distance to scaffold anchors.

## What it does

SweepGA filters genome alignment PAF files by identifying large scaffold alignments and rescuing nearby mappings based on Euclidean distance. The tool first finds high-confidence scaffold chains through merging and plane sweep filtering, then retains any mapping within a specified Euclidean distance threshold of these anchors.

The algorithm works in three phases. First, it merges nearby mappings into chains using union-find, where mappings are chained if they're within a gap threshold in both query and target coordinates. Second, it identifies scaffold chains that exceed a minimum length threshold and removes overlapping scaffolds through plane sweep. Third, it rescues non-scaffold mappings by calculating their Euclidean distance to scaffold anchors in query-target coordinate space, keeping those within the distance threshold.

The Euclidean distance is calculated between mapping center points: for a mapping M and scaffold anchor S, the distance is sqrt((M_query_center - S_query_center)² + (M_target_center - S_target_center)²). This allows rescue across strands, so reverse complement mappings near forward strand scaffolds are retained.

## Installation

Requires Rust 1.70+. Clone and build:

```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo build --release
```

## Usage with FASTGA

FASTGA generates all-vs-all genome alignments in PAF format with full CIGAR strings (pafxm style). These alignments typically need filtering to remove redundant and weak mappings. Here's a complete workflow:

```bash
# Generate all-vs-all alignments with FASTGA
cd deps/FASTGA
cargo build --release
./target/release/fastga -t 8 ../../data/scerevisiae8.fa > ../../scerevisiae8.raw.paf

# The raw output contains all discovered mappings with CIGAR strings
head -n1 ../../scerevisiae8.raw.paf
# S288C#1#chrI    230218  1748    2170    +    SK1#1#chrI    228864  1742    2164    422    422    60    NM:i:0    cg:Z:422M

# Apply scaffold-based filtering (default: scaffolds >10kb, rescue within 100kb)
cd ../..
./target/release/sweepga < scerevisiae8.raw.paf > scerevisiae8.filtered.paf

# Each output mapping is annotated with its filter status
grep "st:Z:" scerevisiae8.filtered.paf | head -n3
# Shows mappings tagged as st:Z:scaffold or st:Z:rescued
```

The default parameters work well for most eukaryotic genomes: scaffold mass of 10kb identifies major syntenic blocks, scaffold jump of 100kb allows chaining across typical intergenic distances, and rescue distance of 100kb captures local rearrangements and smaller homologous features near the main alignments.

## Parameters

The main parameters control scaffold identification and rescue distance:

`-S/--scaffold-mass` sets the minimum length for a chain to be considered a scaffold anchor (default 10000). Smaller values create more anchors but may include repetitive elements.

`-j/--scaffold-jump` sets the maximum gap for merging mappings into scaffold chains (default 100000). This controls how scattered mappings can be while still being considered part of the same scaffold.

`-D/--scaffold-dist` sets the maximum Euclidean distance for rescue (default 100000). Mappings further than this from any scaffold anchor are discarded.

`-m/--mode` controls the filtering strategy: "1:∞" (default) keeps all non-overlapping mappings per query position, "1:1" keeps only the best per query-target pair, "N:N" or equivalently "∞:∞" disables plane sweep filtering.

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

Plane sweep removes overlapping scaffolds by maintaining a list of non-overlapping chains per query sequence, rejecting new chains that overlap existing ones by more than the overlap threshold (default 95%). This ensures scaffold anchors don't overlap significantly.

The rescue phase calculates Euclidean distance from each non-scaffold mapping to all scaffold anchors on the same target sequence. The distance uses mapping center points in 2D space where axes are query and target positions. Any mapping within the distance threshold of at least one anchor is rescued, regardless of strand orientation.

Output preserves the original PAF records with an added tag indicating filter status: st:Z:scaffold for scaffold anchors, st:Z:rescued for rescued mappings. This allows downstream tools to distinguish high-confidence syntenic anchors from nearby features.

## Citation

SweepGA: PAF filtering using Euclidean distance to scaffold anchors
https://github.com/pangenome/sweepga