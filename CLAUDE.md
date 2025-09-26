# SweepGA Filtering Algorithm

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

3. **Scaffold Creation** (if `-S` > 0)
   - Create scaffolds from the (optionally pre-filtered) mappings
   - Merge nearby mappings into scaffold chains (using `-j/--scaffold-jump` parameter)
   - Filter chains by minimum length (`-S/--scaffold-mass`, default: 10kb)
   - Scaffolds define high-confidence syntenic regions

4. **Scaffold Filter** (default: 1:1)
   - Apply plane sweep filtering to the SCAFFOLD chains
   - Respects PanSN prefix grouping when `-Y` is set
   - Controlled by `--scaffold-filter` (default: "1:1")
   - Options: "1:1", "1" (same as "1:∞"), "N" (no filtering)
   - Keeps best scaffolds per query/target pair within each prefix group

5. **Rescue Phase**
   - Identify "anchors": mappings that are members of kept scaffold chains
   - For ALL original mappings (before any filtering):
     - Keep if it's an anchor
     - Otherwise, calculate Euclidean distance to nearest anchor on SAME chromosome pair
     - Rescue if within `-D/--scaffold-dist` distance (default: 100kb)

## Key Parameters

- `-n/--num-mappings`: Primary mapping filter before scaffolds (default: "1:1")
- `--scaffold-filter`: Filter for scaffold chains (default: "1:1")
- `-S/--scaffold-mass`: Minimum scaffold length (default: 10kb)
- `-j/--scaffold-jump`: Maximum gap to merge mappings into scaffolds (default: 100kb)
- `-D/--scaffold-dist`: Maximum distance for rescue (default: 100kb)
- `-Y/--group-prefix`: Delimiter for prefix grouping (default: "#" for PanSN format)

## Implementation Notes

- Both pre-scaffold and scaffold filtering respect PanSN prefix grouping (`-Y`)
- Small mappings can be rescued regardless of size if within rescue distance of a scaffold
- Rescue happens per chromosome pair (both query and target chromosomes must match)
- When `-S` is 0, only pre-scaffold filtering is applied (no scaffolding/rescue)
- The rescue phase checks ALL original mappings, not just filtered ones

## Common Usage Patterns

1. **Default (wfmash-like)**: No pre-filtering, 1:1 scaffold filtering
   ```
   sweepga -S 10000  # Creates scaffolds, keeps best per genome pair
   ```

2. **Aggressive filtering**: 1:1 pre-filtering, 1:1 scaffold filtering
   ```
   sweepga -S 10000 -n 1:1 --scaffold-filter 1:1
   ```

3. **Pre-filter only**: No scaffolding
   ```
   sweepga -n 1:1 -j 0  # Apply 1:1 filter to mappings, no scaffolds
   ```

4. **Keep more mappings**: 1:∞ for both filters
   ```
   sweepga -S 10000 -n 1 --scaffold-filter 1   # "1" is shorthand for "1:∞"
   ```

## Key Differences from Original Description

1. **Pre-scaffold filtering is optional** (default: N:N = no filtering)
2. **Scaffold filtering has separate control** (default: 1:1)
3. **Both filters respect prefix grouping** when `-Y` is set
4. **Rescue uses ALL original mappings**, not post-filter mappings