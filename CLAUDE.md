# SweepGA Filtering Algorithm

## Tools Available

### alnstats - Alignment Statistics Tool

**Purpose:** Calculate and compare PAF alignment statistics, especially useful for validating filtering operations.

**Usage:**
```bash
# Show statistics for a PAF file
alnstats alignments.paf

# Compare before/after filtering
alnstats raw.paf filtered.paf

# Detailed per-genome-pair breakdown
alnstats alignments.paf -d
```

**Key Metrics:**
- Total mappings, bases, average identity
- Per-genome-pair coverage statistics
- Genome pairs with >95% coverage
- Inter-chromosomal and inter-genome mapping counts

**When to Use:**
1. **After filtering operations** - Verify genome pairs aren't lost
2. **Before/after comparisons** - Quantify filtering effects
3. **Debugging coverage issues** - Use `-d` to see per-pair breakdown
4. **Validating changes** - Ensure modifications preserve expected coverage

**Example validation workflow:**
```bash
# Before making changes
alnstats z.paf > before_stats.txt

# Make changes, rerun filtering
sweepga z.paf > filtered.paf

# Compare results
alnstats z.paf filtered.paf

# Verify all genome pairs preserved
alnstats filtered.paf -d | grep "cov"
```

## CRITICAL: Git Commit Rules - NEVER VIOLATE THESE

**ABSOLUTELY FORBIDDEN - NO EXCEPTIONS:**
- **NEVER EVER use `git add -A`, `git add .`, or `git add --all`**
- **NEVER add test data, FASTA files (except data/scerevisiae8.fa.gz), PAF files, or intermediate results**
- **NEVER add files with extensions: .fa, .fasta, .paf, .gdb, .ktab, .bps, .log (unless specifically code logs)**

**REQUIRED PROCEDURE:**
1. Always list files to be added explicitly: `git add src/specific_file.rs`
2. Always check with `git status` before committing
3. If you accidentally stage wrong files, immediately `git reset HEAD <file>`
4. All test data must be generated from data/scerevisiae8.fa.gz only

**BEFORE EVERY COMMIT:**
```bash
git status  # Check what's staged
git diff --cached --name-only  # List all staged files
# Only proceed if list contains ONLY source code files
```

## CRITICAL FIX: 1:1 Filtering with ~100% Coverage (SOLVED)

### The Problem
When applying 1:1 filtering to 99% identical yeast genomes, we were only getting ~23% coverage instead of the expected ~99%. This was unacceptable for genome alignment where we need complete coverage between highly similar genomes.

### The Solution
The key insight is that **1:1 filtering must operate at the chromosome pair level, not the genome pair level**.

#### What Was Wrong:
- We were grouping by genome prefix pairs (e.g., "SGDref#1" → "DBVPG6765#1")
- Within each genome pair group (~460 alignments across all chromosomes), we only kept 1 alignment total
- This meant only 1 chromosome pair got an alignment per genome pair = terrible coverage

#### The Fix:
- Group by full chromosome names (e.g., "SGDref#1#chrI" → "DBVPG6765#1#chrI")
- Apply 1:1 filtering within each chromosome pair group
- Use `plane_sweep_both()` for true 1:1 (respecting both query and target constraints)
- Result: Keep the best alignment for EACH chromosome pair = ~100% coverage

### Implementation Details
In `paf_filter.rs`:
```rust
// For 1:1 filtering: group by full chromosome names (includes genome prefix)
let query_group = meta.query_name.clone();
let target_group = meta.target_name.clone();
```

In `plane_sweep_exact.rs`:
```rust
// True 1:1 filtering applies constraints on BOTH axes
let kept_in_group = if mappings_to_keep == 1 {
    plane_sweep_both(&mut group_mappings, 1, 1, overlap_threshold)
} else {
    plane_sweep_query(&mut group_mappings, mappings_to_keep, overlap_threshold)
};
```

### Results
- Before fix: 1,738 alignments, 23.1% coverage, 0% genome pairs >95% coverage
- After fix: 26,272 alignments, 100.1% coverage, 100% genome pairs >95% coverage

This maintains the expected property that 99% identical genomes should have ~100% reciprocal coverage while still removing redundant/overlapping alignments within each chromosome pair.

## Core Algorithm (Corrected Implementation)

The filtering process follows this exact sequence:

1. **Input Processing**
   - Take input mappings/alignments from any source (PAF format)
   - Filter by minimum block length and optionally exclude self-mappings
   - Store all original mappings for potential rescue later

2. **Primary Mapping Filter** (default: 1:1)
   - Apply plane sweep filtering to raw mappings BEFORE scaffold creation
   - Respects PanSN prefix grouping when `-Y` is set
   - Controlled by `-n/--num-mappings` (default: "1:1")
   - Options: "1:1", "1" (same as "1:∞"), "N" (no filtering)

3. **Scaffold Creation** (if `-s` > 0)
   - Create scaffolds from the (optionally pre-filtered) mappings
   - Merge nearby mappings into scaffold chains (using `-j/--scaffold-jump` parameter, default: 10kb)
   - Filter chains by minimum length (`-s/--scaffold-mass`, default: 10kb)
   - Scaffolds define high-confidence syntenic regions

4. **Scaffold Filter** (default: 1:1)
   - Apply plane sweep filtering to the SCAFFOLD chains
   - Respects PanSN prefix grouping when `-Y` is set
   - Controlled by `-m/--scaffold-filter` (default: "1:1")
   - Options: "1:1", "1" (same as "1:∞"), "N" (no filtering)
   - Keeps best scaffolds per query/target pair within each prefix group

5. **Rescue Phase**
   - Identify "anchors": mappings that are members of kept scaffold chains
   - For ALL original mappings (before any filtering):
     - Keep if it's an anchor
     - Otherwise, calculate Euclidean distance to nearest anchor on SAME chromosome pair
     - Rescue if within `-d/--scaffold-dist` distance (default: 20kb)

## Key Parameters

- `-n/--num-mappings`: Primary mapping filter before scaffolds (default: "many" = no filtering)
- `-m/--scaffold-filter`: Filter for scaffold chains (default: "1:1")
- `-s/--scaffold-mass`: Minimum scaffold length (default: 10kb)
- `-j/--scaffold-jump`: Maximum gap to merge mappings into scaffolds (default: 10kb)
- `-d/--scaffold-dist`: Maximum distance for rescue (default: 20kb)
- `--aligner`: Alignment method for FASTA input (default: "fastga")

## Implementation Notes

- Both pre-scaffold and scaffold filtering respect PanSN prefix grouping (`-Y`)
- Small mappings can be rescued regardless of size if within rescue distance of a scaffold
- Rescue happens per chromosome pair (both query and target chromosomes must match)
- When `-S` is 0, only pre-scaffold filtering is applied (no scaffolding/rescue)
- The rescue phase checks ALL original mappings, not just filtered ones

## Future: Filtering DSL Concept

The current parameter space has grown complex with multiple filtering stages:
- Pre-scaffold plane sweep filtering (`-n`)
- Scaffold creation with merge distance (`-j`)
- Post-scaffold filtering (`--scaffold-filter`)
- Rescue operations (`-d`)
- Various multiplicities at each stage

A potential future improvement would be a filtering DSL/pipeline language:

### Example DSL Syntax:
- `"raw | scaffold(100k) | filter(1:1) | rescue(100k)"` - Current wfmash default
- `"filter(1:1) | scaffold(100k) | filter(1:1)"` - Aggressive filtering
- `"scaffold(100k) | filter(many:1)"`  - Keep many queries per target

### Operations:
- `filter(M:N)` - Plane sweep filtering with multiplicity
- `scaffold(distance)` - Union-find merging with gap distance
- `rescue(distance)` - Rescue mappings near anchors
- `length(min)` - Filter by minimum length

### Preset Definitions:
- `--pipeline wfmash` = No pre-filter, scaffold, 1:1 filter, rescue
- `--pipeline aggressive` = 1:1 pre-filter, scaffold, 1:1 post-filter
- `--pipeline permissive` = Scaffold only, no filtering

This would make the pipeline stages explicit and composable, reducing confusion about what operations happen when.

## Common Usage Patterns

1. **Default (wfmash-like)**: No pre-filtering, 1:1 scaffold filtering
   ```
   sweepga  # Creates scaffolds, keeps best per chromosome pair
   ```

2. **Aggressive filtering**: 1:1 pre-filtering, 1:1 scaffold filtering
   ```
   sweepga -n 1:1 -m 1:1
   ```

3. **Pre-filter only**: No scaffolding
   ```
   sweepga -n 1:1 -j 0  # Apply 1:1 filter to mappings, no scaffolds
   ```

4. **Keep more mappings**: 1:∞ for both filters
   ```
   sweepga -n 1 -m 1   # "1" is shorthand for "1:∞"
   ```

## Key Differences from Original Description

1. **Pre-scaffold filtering is optional** (default: N:N = no filtering)
2. **Scaffold filtering has separate control** (default: 1:1)
3. **Both filters respect prefix grouping** when `-Y` is set
4. **Rescue uses ALL original mappings**, not post-filter mappings

## Testing Coverage - Yeast 8 Genomes

### Expected Behavior for 1:1 Scaffold Filtering
When testing with highly similar genomes (e.g., 8 yeast genomes with ~99% identity):

1. **Input**: ~30K alignments from self-alignment
2. **After plane sweep**: ~29K alignments (no pre-filtering by default)
3. **Scaffold creation**: Groups into ~2-3K scaffold chains
4. **After 1:1 scaffold filter**: Keeps best non-overlapping scaffolds per chromosome pair
5. **Final output**: ~17-18K alignments (scaffold members + rescued)

### Key Testing Commands
```bash
# Test with yeast 8 genomes
./target/release/sweepga data/scerevisiae8.fa.gz -o output.paf

# Check coverage statistics
grep "^SGDref#1#chrI" output.paf | grep "DBVPG6044#1#chrI" | wc -l
# Expected: Multiple alignments per chromosome pair (scaffold members)

# Count scaffold chains per chromosome pair
grep "^SGDref#1#chrI" output.paf | grep "DBVPG6044#1#chrI" | grep -o "ch:Z:chain_[0-9]*" | sort -u | wc -l
# Expected: 3-7 scaffold chains (non-overlapping regions)

# Verify all genome pairs have alignments
awk -F'\t' '{
    split($1, q, "#"); split($6, t, "#");
    if (q[1] != t[1]) pairs[q[1]"->"t[1]] = 1
} END {
    for (p in pairs) count++;
    print count" genome pairs covered"
}' output.paf
# Expected: 56 pairs for 8 genomes (complete coverage)
```

### Understanding Scaffold Output
- Each scaffold is a **chain of nearby alignments** (merged within `-j` distance)
- Scaffold members all share the same `ch:Z:chain_N` tag
- 1:1 filtering applies to **scaffolds**, not individual alignments
- Result: Multiple alignments per chromosome pair, organized into non-overlapping scaffold chains