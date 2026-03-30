# Research: Disk Space Monitoring & Temp File Lifecycle in SweepGA

## 1. Temporary Files and Directories Created During a Run

### 1.1 FastGA Index Files (largest disk consumers)

Created by `FastGAIntegration::create_gdb_only()` and `create_index_only()` (`src/fastga_integration.rs:250-298`):

| File Pattern | Created By | Purpose | Relative Size |
|---|---|---|---|
| `<base>.1gdb` | `prepare_gdb()` | Genome database (sequence storage) | ~1x input FASTA |
| `<base>.gix` | `create_index()` | K-mer index header | Small |
| `.<base>.ktab.<N>` | `create_index()` | K-mer frequency tables (per-thread) | **Largest** (~3-8x input) |
| `.<base>.ktab.<N>.zst` | `compress_index()` | Zstd-compressed ktab (optional, `--zstd`) | ~0.5x uncompressed |
| `.<base>.bps` | `prepare_gdb()` | Base-pair sequence data | ~1x input FASTA |

The **ktab files are the dominant disk consumer** -- one per thread, each storing k-mer frequency tables for the entire genome database. With 8 threads and a 100MB genome, expect ~800MB of ktab files.

### 1.2 Batch Alignment Temp Directories

Created by `run_batch_alignment_generic()` (`src/batch_align.rs:557`) and `process_agc_batched()` (`src/main.rs:2158`):

| Directory Pattern | Purpose |
|---|---|
| `<tempdir>/sweepga_batch_<PID>/batch_<N>/genomes.fa` | Per-batch FASTA files |
| `<tempdir>/sweepga_agc_batches_<PID>/batch_<N>.fa` | AGC-extracted batch FASTAs |
| `<tempdir>/sweepga_pairs_<PID>/<sample>.fa` | Per-pair extracted FASTAs |

### 1.3 Intermediate Alignment Files (NamedTempFile)

Created via Rust's `tempfile::NamedTempFile` throughout the pipeline:

| Location | File Pattern | Purpose |
|---|---|---|
| `src/main.rs:691` | `NamedTempFile::new()` | Temp filtered PAF for tree mode |
| `src/main.rs:1038` | `NamedTempFile::with_suffix(".paf")` | .1aln-to-PAF conversion output |
| `src/main.rs:2855` | `NamedTempFile::with_suffix(".1aln")` | Tree-filtered .1aln |
| `src/main.rs:2894` | `NamedTempFile::with_suffix(".1aln")` | Temp .1aln for filter pipeline |
| `src/main.rs:3461` | `NamedTempFile::with_suffix(".paf")` | PAF output temp |
| `src/fastga_integration.rs:391-398` | `_tmp_<PID>_<NANOS>.1aln` | FastGA alignment intermediate |
| `src/fastga_integration.rs:632` | `NamedTempFile::new()` | Temp .1aln holder |
| `src/fastga_integration.rs:751` | `NamedTempFile::new()` | Temp PAF from align_to_temp_paf |
| `src/agc.rs:163-172` | `sweepga_agc_*.fa` | AGC extraction temp FASTA |
| `src/paf_merge.rs:107` | `NamedTempFile::new()` | Merged PAF output |

`NamedTempFile` auto-deletes when dropped (Rust RAII), so these are cleaned up when the owning variable goes out of scope.

### 1.4 FastGA Direct Alignment Intermediates

In `align_direct_paf()` (`src/fastga_integration.rs:355-489`):
- Creates `_tmp_<PID>_<NANOS>.1aln` in the target's working directory
- This file is explicitly deleted after ALNtoPAF conversion (line 443)
- On FastGA failure, it's also cleaned up (line 430)

## 2. Existing Disk Space Monitoring

### 2.1 `disk_usage` Module (`src/disk_usage.rs`)

A **complete tracking infrastructure** already exists:

**Global atomic counters:**
- `CURRENT_USAGE` -- bytes currently on disk from tracked files
- `PEAK_USAGE` -- high-water mark
- `CUMULATIVE_WRITTEN` -- total bytes ever written

**Key functions:**
- `track_file_created(path)` -- record a new file's size (line 17)
- `track_file_deleted(path)` -- subtract a deleted file's size (line 51)
- `add_bytes(n)` / `remove_bytes(n)` -- manual byte adjustments (lines 69, 89)
- `scan_fastga_index_files(dir)` -- scan for `.1gdb`, `.1gix`, `.1bps`, `.ktab.*` (line 185)
- `track_directory(path)` -- recursive directory tracking (line 167)
- `summary()` / `log_summary()` -- report current/peak/cumulative (lines 126, 135)
- `format_bytes(n)` -- human-readable formatting (line 109)

### 2.2 Integration Points

The disk_usage module is actively used in:

1. **Timing log output** (`src/main.rs:87-95`): Every `timing.log()` call reports `disk:<current> peak_disk:<peak>` on stderr.
2. **Batch alignment** (`src/batch_align.rs:97-127`): `prepare_target()` tracks index file size with `add_bytes()`, `cleanup_target()` removes with `remove_bytes()`.
3. **Single-run alignment** (`src/fastga_integration.rs:716-738`): Background thread polls index directory every 1s, tracking size deltas.
4. **Temp PAF files** (`src/fastga_integration.rs:762`): `track_file_created()` called on temp PAF output.

### 2.3 What's Missing: Available Disk Space Query

**There is NO existing mechanism to query available disk space at runtime.** The module only tracks what SweepGA itself has written. There is no:
- `statvfs` / `statfs` call
- Check against available filesystem space before starting
- Pre-flight validation that enough space exists for the planned operation
- Warning when approaching filesystem limits

## 3. Temp File Cleanup Lifecycle

### 3.1 Batch Alignment Lifecycle (`src/batch_align.rs`)

```
Phase 0: Create batch_dir (line 557-558)
Phase 1: Write batch FASTAs to batch_dir/batch_N/genomes.fa (lines 590-603)
         Aligner.prepare_all() -> create GDBs for all batches (line 606)

For each target batch i:
    Phase 2a: Aligner.prepare_target(i) -> create index (lines 90-110)
              [DISK PEAK: GDB for all batches + index for batch i]
    Phase 2b: For each query batch j:
                Aligner.align(j, i) -> produce PAF bytes (line 628)
    Phase 2c: Aligner.cleanup_target(i) -> delete index files (lines 117-129)
              [Index freed]

Phase 3: Aligner.cleanup_all() -> delete all GDB files (lines 132-138)
         remove_dir_all(batch_dir) (line 645)
         [All temp files freed]
```

### 3.2 AGC Batched Lifecycle (`src/main.rs:2107-2264`)

```
Phase 0: Create batch_dir (line 2158-2159)
Phase 1: Extract all batches to batch_N.fa (lines 2162-2173)
Phase 2: Create GDB for all batches (lines 2189-2193)

For each target batch i:
    Phase 3a: Build index (line 2211)
              [DISK PEAK: all batch FASTAs + all GDBs + target index]
    Phase 3b: Align all query batches against target i (lines 2218-2240)
    Phase 3c: Cleanup index (line 2248-2252)

Phase 4: remove_dir_all(batch_dir) (line 2257)
```

### 3.3 Pair-by-Pair Lifecycle (`src/main.rs:1948-2097`)

```
For each pair:
    Extract query + target FASTA to work_dir (lines 2025-2026)
    Create GDB + index for target (lines 2046-2047)
    Align (line 2056)
    cleanup_all(target_gdb) (line 2062)
    remove_file(query/target FASTA) (lines 2038, 2071-2072)

After all pairs: remove_dir_all(work_dir) (line 2090)
```

### 3.4 Single-Run (Non-Batch) Lifecycle

When running without `--batch-bytes`, FastGA creates and cleans up internally:
- `align_to_temp_paf()` (`src/fastga_integration.rs:668-765`): Uses `NamedTempFile` (auto-cleanup on drop). Background thread monitors index files during alignment; indexes are cleaned up by FastGA internally after `align_with_existing_indices()`.
- `align_direct_paf()` (`src/fastga_integration.rs:355-489`): Creates `_tmp_*.1aln`, explicitly deleted at line 443.

### 3.5 NamedTempFile Auto-Cleanup

All `NamedTempFile` instances are cleaned up when their owning variable goes out of scope (Rust RAII). This means:
- Temp PAF files from alignment are cleaned up when the function returns or the variable is reassigned
- No explicit cleanup needed, but files persist until the owning scope ends
- **Risk**: If the process is killed (SIGKILL), temp files in the OS temp directory are NOT cleaned up

## 4. Relationship Between Batch Size, Number of Genomes, and Disk Usage

### 4.1 Disk Usage Formula

For a FASTA input alignment run, the peak disk usage comes from:

```
Peak_disk = batch_FASTAs + GDBs_all_batches + index_one_batch + temp_PAF_output

Where:
  batch_FASTAs      = total_sequence_bp (all genomes, written to batch dirs)
  GDBs_all_batches  = N_batches * avg_batch_bp * GDB_MULTIPLIER
  index_one_batch   = avg_batch_bp * INDEX_MULTIPLIER * N_threads
  temp_PAF_output   = N_alignment_pairs * avg_alignments_per_pair * ~150 bytes/line

Multipliers (empirical, from FastGA):
  GDB_MULTIPLIER    ~= 2.0  (1gdb + .bps files, ~1x each)
  INDEX_MULTIPLIER  ~= 1.0 per thread  (ktab files, one per thread)
```

### 4.2 Concrete Example: 8 Yeast Genomes (~12MB each)

Single batch (no `--batch-bytes`):
```
Input FASTA:       ~96 MB total
GDB files:         ~192 MB (2x input)
Index (8 threads): ~768 MB (8x input, ktab per thread)
Temp PAF:          ~5 MB
Peak total:        ~1,061 MB
```

With `--batch-bytes 50M` (2 genomes per batch = 4 batches):
```
Batch FASTAs:      ~96 MB (copied to batch dirs)
GDB files:         ~96 MB (4 batches * 24MB per batch GDB)
Index (1 batch):   ~192 MB (24MB * 8 threads)
Temp PAF:          ~5 MB
Peak total:        ~389 MB
```

### 4.3 Scaling Laws

| Parameter | Effect on Peak Disk |
|---|---|
| N_genomes (no batching) | Linear: GDB + index scale with total sequence |
| N_genomes (batched) | Sub-linear: only one batch's index at a time |
| batch_size | Controls peak: larger batches = larger index per target |
| N_threads | Linear: one ktab file per thread |
| --zstd | ~2x reduction in ktab files on disk |
| Genome size | Linear: all intermediate files scale with sequence size |
| N_batches^2 | Total alignments = N_batches^2 (all-pairs), cumulative I/O |

### 4.4 Key Insight: Batching Reduces Peak But Not Total I/O

Batching reduces peak disk usage (only one target's index at a time) but **increases total I/O** because:
- Each genome's GDB is created once but the index is built N_batches times (once as target for each batch round)
- More batch FASTA files are written

## 5. Can We Estimate Disk Usage Before Running?

### 5.1 What We Can Compute Pre-Flight

From the input FASTA and CLI parameters, we can estimate:

```rust
fn estimate_peak_disk(
    total_sequence_bp: u64,  // from parse_genome_sizes()
    n_genomes: usize,
    batch_max_bp: Option<u64>,
    n_threads: usize,
    zstd_compress: bool,
) -> u64 {
    let (n_batches, max_batch_bp) = if let Some(max_bp) = batch_max_bp {
        let n = (total_sequence_bp + max_bp - 1) / max_bp;
        (n as usize, max_bp)
    } else {
        (1, total_sequence_bp)
    };

    let gdb_factor = 2.0;   // .1gdb + .bps
    let ktab_per_thread = 1.0;  // each ktab ~= sequence size
    let zstd_factor = if zstd_compress { 0.5 } else { 1.0 };

    let batch_fastas = total_sequence_bp;  // copies of input
    let all_gdbs = (total_sequence_bp as f64 * gdb_factor) as u64;
    let one_index = (max_batch_bp as f64 * ktab_per_thread * n_threads as f64 * zstd_factor) as u64;
    let paf_estimate = total_sequence_bp / 10;  // rough upper bound

    batch_fastas + all_gdbs + one_index + paf_estimate
}
```

### 5.2 Available Functions for Size Estimation

Already available:
- `batch_align::parse_genome_sizes()` (`src/batch_align.rs:274`) -- parses FASTA to get per-genome bp counts
- `agc.get_sample_sizes()` (`src/agc.rs:62`) -- gets bp per sample from AGC without extraction
- `batch_align::partition_into_batches_by_bp()` (`src/batch_align.rs:378`) -- computes batch layout

Not yet available:
- FastGA GDB/index size multiplier constants (would need empirical calibration)
- Filesystem free space query

### 5.3 Empirical Calibration Needed

The GDB_MULTIPLIER and INDEX_MULTIPLIER values above are estimates. To get accurate numbers:
1. Run a test alignment and measure actual GDB, ktab, and gix file sizes relative to input
2. These ratios may vary with k-mer frequency (`-f`) and sequence complexity
3. Could be measured once at startup on a small subsample

## 6. What's Missing for Disk-Budget Enforcement

### 6.1 No Available Disk Space Query

The `disk_usage` module tracks **what SweepGA has written** but never checks **how much space is available**. On Linux, this would require:

```rust
use std::os::unix::fs::MetadataExt;
use nix::sys::statvfs::statvfs;

fn available_disk_bytes(path: &str) -> Result<u64> {
    let stat = statvfs(path)?;
    Ok(stat.blocks_available() * stat.fragment_size())
}
```

Or without nix dependency, using libc directly:
```rust
fn available_bytes(path: &str) -> std::io::Result<u64> {
    use std::ffi::CString;
    let c_path = CString::new(path)?;
    let mut stat: libc::statvfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::statvfs(c_path.as_ptr(), &mut stat) } != 0 {
        return Err(std::io::Error::last_os_error());
    }
    Ok(stat.f_bavail as u64 * stat.f_frsize as u64)
}
```

### 6.2 No Pre-Flight Validation

There is no check that says "this run will need ~X GB of temp space, you only have ~Y GB available." The batch alignment code will simply fail with I/O errors when disk fills up.

### 6.3 No Runtime Budget Enforcement

During a run, there is no mechanism to:
- Pause alignment if disk usage exceeds a threshold
- Abort early if free space drops below a safety margin
- Dynamically reduce batch size if disk is running low

### 6.4 No Cleanup on Failure

If an alignment fails mid-batch (e.g., disk full), the cleanup code at `batch_align.rs:645` and `main.rs:2257` may not execute because the error propagates before reaching cleanup. The `remove_dir_all` calls are after the main loop, not in a `Drop` guard.

**Exception**: `NamedTempFile` instances do auto-delete on drop, even on panic (unless SIGKILL).

### 6.5 Missing Signal Handler for Cleanup

If the process receives SIGINT/SIGTERM, batch directories may be left behind. There is no signal handler to clean up `sweepga_batch_*` or `sweepga_agc_*` directories.

## Summary

| Aspect | Status |
|---|---|
| Temp file creation | Well-structured, uses NamedTempFile and explicit cleanup |
| Disk usage tracking | **Exists** (`disk_usage.rs`), reports current/peak/cumulative |
| Disk usage in logs | **Active** -- every `timing.log()` shows disk stats |
| Available disk query | **Missing** -- no `statvfs` or equivalent |
| Pre-flight estimation | **Missing** -- could be built from `parse_genome_sizes()` + empirical multipliers |
| Budget enforcement | **Missing** -- no abort/pause on disk limits |
| Cleanup on error | **Partial** -- NamedTempFile OK, batch dirs may leak on failure |
| Signal-safe cleanup | **Missing** -- no SIGINT/SIGTERM handler |
| Batch size vs. disk | **Well understood** -- `--batch-bytes` directly controls peak index size |
