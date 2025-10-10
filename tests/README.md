# SweepGA Test Suite

## Test Categories

### Unit Tests & Working Integration Tests
Fast tests that run automatically in CI using the library API:

- **Unit tests** - Test individual modules and functions
- **Working integration tests** (`test_working_integration.rs`) - Test end-to-end workflows using:
  - Small test data (`data/B-3106.fa` - 32KB)
  - Synthetic generated sequences
  - Library API directly (no binary required)

```bash
cargo test
```

### Binary Integration Tests (Ignored)
These tests require the compiled binary and are marked `#[ignore]`:

- `tests/fastga_integration.rs` - Tests that invoke sweepga binary with FastGA
- `tests/test_centromere_plane_sweep.rs` - Centromere inversion tests
- `tests/test_chaining_stability.rs` - Chaining behavior tests with yeast data
- `tests/test_chain_monotonicity.rs` - Chain monotonicity tests
- `tests/test_genome_pair_grouping.rs` - Requires z.paf file (not in repo)

Run these manually after building:
```bash
cargo build --release
cargo test --release -- --ignored
```

## Test Strategy

### Working Tests (Not Ignored)
Use the sweepga library API directly:
- ✅ Run in CI automatically
- ✅ Fast execution
- ✅ Test core functionality
- ✅ Gracefully skip if FastGA binaries unavailable

### Binary Tests (Ignored)
Invoke the compiled `sweepga` binary:
- ⚠️ Require binary to be built first
- ⚠️ May need data files not in repo
- ⚠️ Slower, better for manual testing
- Use for regression testing and end-to-end validation

## Running Tests in CI

CI runs:
1. `cargo test` - Runs all non-ignored unit tests
2. `cargo test --release` - Runs non-ignored tests in release mode

To run ignored tests locally:
```bash
# After ensuring sweepga is built
cargo test --release -- --ignored --test-threads=1
```

## Test Data

- `data/scerevisiae8.fa.gz` - 8 yeast genomes for integration testing (committed to repo)
- `z.paf` - Large PAF file used by some tests (NOT in repo, generate locally)
