# Solution: Correct Identity Calculation from .1aln Format

## Problem Statement

When reading .1aln files in sweepga, the calculated identity values did not match those produced by ALNtoPAF. This prevented format-equivalent filtering between .1aln and PAF workflows.

## Root Cause

The issue had two parts:

1. **Wrong field used**: The 'D' record in .1aln format is trace-point related data, NOT the actual edit distance
2. **Wrong formula**: The actual edit distances are stored in the 'X' record (INT_LIST), and require a specific formula

## The Correct Formula

ALNtoPAF calculates divergence as:

```
del = target_span - query_span
divergence = (sum(X) - del) / query_span / 2.0
identity = 1.0 - divergence
```

Where:
- **sum(X)** = Sum of all values in the X record INT_LIST (per-tracepoint edit distances)
- **del** = Deletions in query (target_span - query_span)
- **query_span** = query_end - query_start
- **target_span** = target_end - target_start

The key insight is that `sum(X)` represents a "symmetric" divergence that must be divided by 2.

## Implementation

### Changes to fastga-rs (~/fastga-rs/src/onelib.rs)

Added X record reading in `AlnReader::read_alignment()`:

```rust
'X' => {
    // X record contains INT_LIST of per-tracepoint edit distances
    if let Some(x_values) = self.file.int_list() {
        x_list_sum = x_values.iter().map(|&v| v as u64).sum();
        has_x_records = true;
    }
}
```

Updated identity calculation:

```rust
let query_span = (a_end - a_beg) as i64;
let target_span = (b_end - b_beg) as i64;
let diffs = if has_x_records { x_list_sum } else { diffs_d_record };

let del = target_span - query_span;
let divergence = if query_span > 0 {
    ((diffs as i64 - del) as f64 / query_span as f64) / 2.0
} else {
    0.0
};

let identity = (1.0 - divergence).max(0.0).min(1.0);

// Calculate matches from identity
let calculated_matches = (identity * query_span as f64) as u64;
let final_matches = if matches > 0 { matches } else { calculated_matches };
```

### Validation Results

Tested on 100 records from yeast genome alignments:
- **Match rate: 100%** (all identities within 0.0001 tolerance)
- **Maximum difference: 0.000098**

Example:
```
Rec  FastGA_ID   PAF_ID     Diff
1    0.975936    0.976000   0.000064  ✓
2    0.985553    0.985600   0.000047  ✓
3    0.951141    0.951200   0.000059  ✓
...
```

## Files Modified

1. **~/fastga-rs/src/onelib.rs**: Added X record reading and correct identity calculation
2. **~/sweepga/Cargo.toml**: Changed to use local fastga-rs path for testing

## Testing

Run verification:
```bash
cargo build --release --example verify_identity
./target/release/examples/verify_identity
```

Expected output:
```
✅ SUCCESS: All identities match between .1aln and PAF!
The X record reading and divergence calculation is working correctly.
```

## Next Steps

1. ✅ Verify unified_filter.rs works with corrected identity values
2. ✅ Create format equivalence integration test
3. Push changes to fastga-rs repository
4. Update sweepga to use updated fastga-rs commit

## References

- Original analysis: `ALNtoPAF_IDENTITY_CALCULATION.md`
- Test data: `test_output.1aln` (yeast 8-genome alignment)
- Comparison: `test_output_converted.paf` (via ALNtoPAF)
