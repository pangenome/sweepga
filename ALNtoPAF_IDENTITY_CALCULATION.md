# ALNtoPAF Identity Calculation

## Summary

ALNtoPAF calculates sequence identity using a specific formula that differs from naive approaches. Understanding this is critical for achieving equivalence between PAF and .1aln filtering workflows.

## The Calculation

### Variables from ALNtoPAF.c

```c
int  del = 0;                              // Accumulated deletions (gaps in query)
int  query_span = path->aepos - path->abpos;  // Query coordinates span
int  blocksum = query_span + del;          // Total alignment length
int  iid = blocksum - path->diffs;         // Identity matches
```

Where:
- `del`: Total number of deletions (D operations in CIGAR) - gaps in query relative to target
- `path->diffs`: From .1aln 'D' record - total edit distance including ALL differences
- `query_span`: Length of aligned region in query coordinates
- `blocksum`: Total alignment length (query_span + deletions)

### PAF Output

```c
// Column 10: Number of matching bases
itoa(iid, buf, out);

// Column 11: Alignment block length
itoa(blocksum, buf, out);

// dv:f: tag - Divergence fraction
int x = 10000 + (10000ll * (query_span - iid)) / query_span;
// Output as 0.XXXX
```

### Divergence Formula (dv:f: tag)

```
divergence = (query_span - iid) / query_span
```

Expanding `iid`:
```
iid = blocksum - path->diffs
    = (query_span + del) - path->diffs
```

Therefore:
```
query_span - iid = query_span - (query_span + del - path->diffs)
                 = path->diffs - del
```

So:
```
divergence = (path->diffs - del) / query_span
           = (total_diffs - deletions) / query_span
           = (substitutions + insertions) / query_span
```

**Key Insight**: The divergence excludes deletions from the numerator but uses query_span (which doesn't include deletions) as the denominator!

### Identity Calculation

```
identity = 1 - divergence
         = 1 - (path->diffs - del) / query_span
         = (query_span - path->diffs + del) / query_span
```

Or equivalently:
```
identity = iid / blocksum
         = (blocksum - path->diffs) / blocksum
         = (query_span + del - path->diffs) / (query_span + del)
```

## What This Means

### The 'D' Record

The .1aln 'D' record contains:
```
path->diffs = substitutions + insertions + deletions
```

This is the **total edit distance** including all types of differences.

### The Problem

When we read .1aln files, we only have access to:
- `path->diffs` (from 'D' record) = total edit distance
- `query_span` = query_end - query_start
- `target_span` = target_end - target_start

But we need `del` (deletions) to calculate divergence correctly!

### The Solution

**Option 1: Calculate `del` from coordinates**
```rust
let query_span = aln.query_end - aln.query_start;
let target_span = aln.target_end - aln.target_start;
let del = target_span.saturating_sub(query_span); // Deletions in query

let diffs = aln.mismatches; // From 'D' record
let divergence_numerator = diffs - del;
let divergence = divergence_numerator as f64 / query_span as f64;
let identity = 1.0 - divergence;
```

**Option 2: Use blocksum**
```rust
let query_span = aln.query_end - aln.query_start;
let target_span = aln.target_end - aln.target_start;
let del = target_span.saturating_sub(query_span);
let blocksum = query_span + del;

let diffs = aln.mismatches; // From 'D' record
let matches = blocksum - diffs;
let identity = matches as f64 / blocksum as f64;
```

Both formulas are equivalent and match ALNtoPAF's calculation!

## Example from Test Data

Record 2 from B-3106.fa:
```
Query:  2-3343 (span = 3341)
Target: 0-3339 (span = 3339)
D record: diffs = 101

del = target_span - query_span = 3339 - 3341 = -2 (WAIT, THIS IS WRONG!)
```

Actually, we need to think about strand:
- For forward strand: `del = max(0, target_span - query_span)`
- For reverse strand: need to be careful about coordinate systems

Let me look at the actual coordinates more carefully...

Actually, the issue is more subtle. The `del` variable in ALNtoPAF is calculated from the CIGAR string during trace point reconstruction, not from simple coordinate arithmetic.

## The Real Solution: Read X Records

The .1aln format has:
- 'T' record: Trace points (target coordinates)
- 'X' record: Per-tracepoint edit distances (INT_LIST)

The sum of X values gives the actual edit distance used by ALNtoPAF, which may differ from the 'D' record!

To achieve perfect equivalence:
1. Read X records using `onecode-rs::OneFile::int_list()`
2. Sum the X values to get the actual divergence
3. Use the formula: `identity = (query_span - sum_X) / query_span`

This avoids needing to calculate `del` separately.
