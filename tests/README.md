# SweepGA Test Suite

## Test Categories

### Unit Tests
Fast tests that don't require external dependencies. Run automatically in CI.
```bash
cargo test
```

### Integration Tests (Ignored)
These tests require:
- The `sweepga` binary to be compiled
- FastGA binaries (built automatically by fastga-rs)
- Test data files (e.g., `data/scerevisiae8.fa.gz`)

Run these manually after building:
```bash
cargo build --release
cargo test --release -- --ignored
```

### Test Files Marked as Ignored

- `tests/fastga_integration.rs` - Tests that use FastGA binaries (FAtoGDB, etc.)
- `tests/test_centromere_plane_sweep.rs` - Binary-dependent centromere tests
- `tests/test_chaining_stability.rs` - Chaining behavior tests with yeast data
- `tests/test_chain_monotonicity.rs` - Chain monotonicity tests
- `tests/test_genome_pair_grouping.rs` - Requires z.paf file

## Why Are Tests Ignored?

Integration tests that invoke the compiled `sweepga` binary are marked `#[ignore]` because:

1. **Build Order**: Tests run before binaries are fully built in debug mode
2. **CI Environment**: Some tests need data files not in the repository  
3. **Performance**: Some tests are slow and better suited for manual runs

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
