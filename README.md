# SweepGA

Fast genome alignment with plane sweep filtering and scaffolding. Wraps FastGA aligner and applies plane sweep and other filtering methods to keep the best non-overlapping alignments.

## What it does

SweepGA can:
1. **Align FASTA files directly** using integrated FastGA (supports .fa.gz)
2. **Filter existing alignments** from any aligner (wfmash, minimap2, etc.)
3. **Apply scaffolding/chaining** to merge nearby alignments into syntenic regions
4. **Output multiple formats**: PAF (text) or .1aln (binary ONE format)

By default, it applies a two-stage filtering approach:
- **Stage 1**: Removes duplicate mappings (1:1 plane sweep)
- **Stage 2**: Creates scaffolds, keeps all non-overlapping scaffold chains (N:N), and rescues nearby mappings

## Tools Included

This package includes two binaries:

- **`sweepga`** - Genome alignment, filtering, and scaffolding tool
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

### From crates.io (recommended)

```bash
cargo install sweepga
```

This installs both `sweepga` and `alnstats` binaries from the published crate.

### From source

Requires Rust 1.70+. Clone and install:

```bash
git clone https://github.com/pangenome/sweepga.git
cd sweepga
cargo install --force --path .
```

### Build Issues on Mixed glibc Systems

**Symptoms:** Build fails with linker errors like:
```
ld: /usr/lib/x86_64-linux-gnu/librt.so: undefined reference to '__pthread_barrier_wait@GLIBC_PRIVATE'
```

This occurs on systems with multiple package managers (e.g., Debian + Guix) providing different glibc versions.

**Fix:** Use the clean build script to isolate from environment conflicts:

```bash
./scripts/build-clean.sh --install
```

See [docs/BUILD-NOTES.md](docs/BUILD-NOTES.md) for details.

### With Guix

Adapted from https://issues.genenetwork.org/topics/rust/guix-rust-bootstrap:

```shell
# Update Guix
mkdir -p $HOME/opt
guix pull -p $HOME/opt/guix-pull-20251012 --url=https://codeberg.org/guix/guix

# Be sure to use the updated Guix
alias guix=$HOME/opt/guix-pull-20251012/bin/guix

# Update Rust and Cargo
mkdir -p ~/.cargo ~/.rustup # to prevent rebuilds
guix shell --share=$HOME/.cargo  --share=$HOME/.rustup -C -N -D -F -v 3 guix gcc-toolchain make libdeflate pkg-config xz coreutils sed zstd zlib nss-certs openssl curl
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. ~/.cargo/env
rustup default stable
exit

# Clone the repository
git clone https://github.com/pangenome/sweepga.git
cd sweepga

guix shell --share=$HOME/.cargo  --share=$HOME/.rustup -C -N -D -F -v 3 guix gcc-toolchain make libdeflate pkg-config xz coreutils sed zstd zlib nss-certs openssl curl cmake clang # we need cmake and clang too for building
. ~/.cargo/env
export LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib
cargo build --release

# Check the lib path and put it into your ~/.bashrc or ~/.zshrc
echo $GUIX_ENVIRONMENT/
   #/gnu/store/whgjblccmr4kdmsi4vg8h0p53m5f7sch-profile/
exit
echo "export GUIX_ENVIRONMENT=/gnu/store/whgjblccmr4kdmsi4vg8h0p53m5f7sch-profile/" >> ~/.bashrc # or ~/.zshrc
source ~/.bashrc # or ~/.zshrc

# Use the executable in sweepga/target/release
env LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib ./target/release/sweepga --help
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

# Output in .1aln format (binary, more compact)
sweepga genome.fa.gz --output-file output.1aln
```

### Filtering existing alignments

```bash
# Filter PAF from stdin
cat alignments.paf | sweepga > filtered.paf

# Read PAF file directly
sweepga alignments.paf > filtered.paf

# Filter .1aln format (preserves format)
sweepga alignments.1aln --output-file filtered.1aln

# Convert .1aln to PAF
sweepga alignments.1aln --paf > output.paf
```

### Scaffolding and chaining

Scaffolding merges nearby alignments into syntenic chains, then filters and rescues alignments:

```bash
# Enable scaffolding with 10kb gap distance
sweepga genome.fa.gz -j 10000 > output.paf

# Aggressive scaffolding with 1:1 filtering and rescue
sweepga alignments.paf -j 10k -s 10k -m 1:1 -d 20k > filtered.paf

# Permissive: keep all scaffolds without filtering
sweepga alignments.paf -j 50k -m N:N > filtered.paf
```

## Parameters

### Basic Filtering

**`-n/--num-mappings`** - Pre-scaffold filter: n:m-best mappings in query:target dimensions (default: `1:1`)
- `"1:1"` - Orthogonal: keep best mapping on both query and target axes
- `"1"` or `"1:∞"` - Keep best mapping per query position only
- `"many"` or `"N:N"` - No pre-filtering, pass all to scaffolding
- `"n:m"` - Keep top n per query, top m per target

**`-o/--overlap`** - Maximum overlap ratio for mappings (default: 0.95)
- Mappings with >95% overlap with a better-scoring mapping are removed

**`-b/--block-length`** - Minimum alignment block length (default: 0)

**`-i/--min-identity`** - Minimum identity threshold (0-1 fraction, 1-100%, or "aniN")

**`--scoring`** - Scoring function for plane sweep (default: `log-length-ani`)
- `ani` - Sort by alignment identity %
- `length` - Sort by alignment length
- `length-ani` - Sort by length × identity
- `log-length-ani` - Sort by log(length) × identity
- `matches` - Sort by number of matching bases

**`--self`** - Include self-mappings (excluded by default)

**`-N/--no-filter`** - Disable all filtering

### Scaffolding and Chaining

**`-j/--scaffold-jump`** - Gap distance for merging alignments into scaffolds (default: `10k`)
- `0` - Scaffolding disabled
- `10000` or `10k` - Merge alignments within 10kb gaps (default, moderate)
- Higher values create longer scaffold chains
- Accepts k/m/g suffix (e.g., `50k`, `1m`)

**`-s/--scaffold-mass`** - Minimum scaffold chain length (default: `10k`)
- Filters out short scaffold chains
- Accepts k/m/g suffix

**`-m/--scaffold-filter`** - Scaffold filter mode (default: `many:many` = no filtering)
- `"1:1"` - Keep best scaffold per chromosome pair (**90-99% reduction**)
- `"1"` or `"1:∞"` - One scaffold per query, many per target
- `"N:N"` or `"many:many"` - Keep all non-overlapping scaffolds (**30-70% reduction**)
- `"n:m"` - Keep top n per query, top m per target

**`-O/--scaffold-overlap`** - Overlap threshold for scaffold filtering (default: 0.5)

**`-d/--scaffold-dist`** - Maximum distance for rescuing alignments near scaffolds (default: `20k`)
- `0` - No rescue (most aggressive)
- `20000` or `20k` - Rescue alignments within 20kb of scaffold anchors
- Higher values rescue more alignments
- Distance is Euclidean: `sqrt((q_dist)² + (t_dist)²)`

**`-Y/--min-scaffold-identity`** - Minimum scaffold identity threshold (defaults to `-i` value)

**`--scaffolds-only`** - Output only scaffold chains for debugging

### Output Options

**`--output-file <FILE>`** - Write output to file (auto-detects format from extension)
- `.paf` - PAF text format
- `.1aln` - Binary ONE format (more compact)

**`--paf`** - Force PAF output format (overrides .1aln default for .1aln inputs)

### Advanced Filtering

**`-x/--sparsify`** - Tree sparsification pattern (default: `1.0` = keep all)
- `"0.5"` - Keep 50% of alignments
- `"tree:3"` - Tree pattern with depth 3
- `"tree:3,2,0.1"` - Complex tree pattern

**`--ani-method`** - ANI calculation method (default: `n100`)
- `all` - Use all bases
- `orthogonal` - Use orthogonal alignments only
- `nX[-sort]` - Use top X% (e.g., `n50`, `n90-identity`, `n100-score`)

### General Options

**`-t/--threads`** - Number of threads for parallel processing (default: 8)

**`--quiet`** - Quiet mode (no progress output)

**`--tempdir <DIR>`** - Temporary directory for intermediate files

**`--check-fastga`** - Check FastGA binary locations and exit (diagnostic)

**`--aligner <ALIGNER>`** - Aligner for FASTA input (default: `fastga`)

**`-f/--frequency <N>`** - FastGA k-mer frequency threshold

**`--all-pairs`** - Align all genome pairs separately (for many genomes)

## How Scaffolding Works

When scaffolding is enabled (default: `-j 10k`), the filtering process follows these steps:

1. **Input Processing** - Filter by minimum block length, exclude self-mappings
2. **Pre-scaffold Filter** - Apply plane sweep filter to individual mappings (`-n`, default: `1:1`)
   - Removes duplicate/overlapping mappings
   - Keeps best mapping per query-target pair
3. **Scaffold Creation** - Merge nearby alignments into chains using union-find algorithm
   - Alignments within `-j` gap distance on both query and target are merged
   - Chains shorter than `-s` minimum length are discarded
4. **Scaffold Filter** - Apply plane sweep filter to scaffold chains (`-m`, default: `N:N`)
   - `1:1` mode: Keep single best scaffold per chromosome pair (aggressive)
   - `N:N` mode: Keep all non-overlapping scaffolds (moderate, default)
5. **Rescue Phase** - Recover alignments near kept scaffolds (`-d`, default: `20k`)
   - Alignments within Euclidean distance of scaffold anchors are rescued
   - Works per chromosome pair only

**Example effects** (514k input alignments):
- **Default** (`-n 1:1 -j 10k -m N:N -d 20k`): 180k alignments (65% reduction, moderate)
- **Aggressive** (`-n 1:1 -j 10k -m 1:1 -d 20k`): 13k alignments (97% reduction)
- **No scaffolding** (`-n 1:1 -j 0`): 476k alignments (7% reduction from 1:1 plane sweep only)

## How Plane Sweep Works

The plane sweep algorithm operates per query-target chromosome pair:

1. **Group** by chromosome pairs
2. **Sort** alignments by query position
3. **Score** each alignment using `--scoring` function (default: `log(length) × identity`)
4. **Sweep** left-to-right, keeping best alignments based on multiplicity:
   - `1:1`: Keep single best alignment per position on both query and target
   - `1:∞`: Keep best alignment per query position (multiple targets allowed)
   - `N:N`: Keep all non-overlapping alignments
5. **Filter** alignments with overlap exceeding threshold (`-o`)

## PanSN Naming Convention

When using PanSN-style sequence names (e.g., `genome#haplotype#chromosome`):

- Filtering operates per chromosome pair within genome pairs
- Example: `SGDref#1#chrI` paired with `DBVPG6765#1#chrI`
- Ensures best alignments are kept for each chromosome pair within each genome pair
- Maintains ~100% coverage for highly similar genomes

## Output Formats

### PAF Format (default for FASTA input)
```
query_name  query_len  query_start  query_end  strand  target_name  target_len  target_start  target_end  matches  block_len  mapq  [tags...]
```

### .1aln Format (binary ONE format)
- More compact than PAF (typically 50-70% smaller)
- Preserves all alignment information including X records (edit distances)
- Faster to read/write for large files
- Use `--paf` flag to convert to PAF for visualization

## Citation

SweepGA: Fast plane sweep filtering for genome alignments
https://github.com/pangenome/sweepga

## License

MIT License - see LICENSE file for details
