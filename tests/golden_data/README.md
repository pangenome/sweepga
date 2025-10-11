# Golden Test Data

This directory contains reference outputs used for regression testing.
Tests verify that current outputs match these golden files via SHA256 checksums.

## Purpose

Golden files lock down exact behavior to prevent unintended changes:
- Coordinate conversion stability
- Filtering consistency
- Output format preservation

## Files

- `checksums.txt` - SHA256 hashes of all golden outputs
- `*.1aln` - Binary .1aln reference files
- `*.paf` - Text PAF reference files
- `generate_golden.sh` - Script to regenerate golden files (use with caution!)

## Regenerating Golden Files

**WARNING**: Only regenerate when intentionally changing behavior!

```bash
cd tests/golden_data
./generate_golden.sh
```

This updates all golden files and checksums. Commit the changes with a clear
explanation of why output changed.

## Test Failure

If golden tests fail:

1. **Expected failure**: You changed filtering/coordinates intentionally
   - Review the diff carefully
   - Regenerate golden files
   - Update CHANGELOG.md explaining the change

2. **Unexpected failure**: You introduced a regression
   - Revert your changes
   - Fix the bug
   - Tests should pass again

## Checksum Format

```
<sha256> <filename>
```

Example:
```
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  test_output.paf
```
