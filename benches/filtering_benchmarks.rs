/// Performance benchmarks for filtering operations
///
/// Run with: cargo bench
///
/// These benchmarks track performance over time to detect regressions.
/// Note: These are skeleton benchmarks - the actual filtering API doesn't
/// expose enough public methods for detailed internal benchmarking yet.
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::fs;
use std::process::Command;
use tempfile::TempDir;

/// Generate synthetic PAF data for benchmarking
fn generate_synthetic_paf(num_alignments: usize) -> String {
    let mut paf_lines = Vec::with_capacity(num_alignments);

    for i in 0..num_alignments {
        let query_name = format!("seq{}", i % 100);
        let target_name = format!("seq{}", (i + 1) % 100);
        let query_start = (i * 1000) % 50000;
        let query_end = query_start + 1000;
        let target_start = (i * 800) % 40000;
        let target_end = target_start + 1000;

        paf_lines.push(format!(
            "{query_name}\t100000\t{query_start}\t{query_end}\t+\t{target_name}\t80000\t{target_start}\t{target_end}\t950\t1000\t60"
        ));
    }

    paf_lines.join("\n")
}

/// Benchmark: End-to-end PAF filtering pipeline
fn bench_paf_filtering_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("paf_filtering");

    for size in [100, 1000, 10000].iter() {
        group.throughput(Throughput::Elements(*size as u64));
        group.sample_size(10); // Reduce sample size for faster benchmarks

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            let paf_data = generate_synthetic_paf(size);

            b.iter_with_setup(
                || {
                    // Setup: create temp files
                    let temp_dir = TempDir::new().unwrap();
                    let input_paf = temp_dir.path().join("input.paf");
                    let output_paf = temp_dir.path().join("output.paf");
                    fs::write(&input_paf, &paf_data).unwrap();
                    (temp_dir, input_paf, output_paf)
                },
                |(temp_dir, input_paf, output_paf)| {
                    // Benchmark: run the filtering binary
                    Command::new("cargo")
                        .args([
                            "run",
                            "--release",
                            "--quiet",
                            "--bin",
                            "sweepga",
                            "--",
                            input_paf.to_str().unwrap(),
                        ])
                        .stdout(fs::File::create(black_box(&output_paf)).unwrap())
                        .status()
                        .unwrap();

                    drop(temp_dir); // Clean up
                },
            );
        });
    }

    group.finish();
}

/// Benchmark: 1:1 filtering
fn bench_one_to_one_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("filtering_1to1");

    for size in [100, 1000, 10000].iter() {
        group.throughput(Throughput::Elements(*size as u64));
        group.sample_size(10);

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            let paf_data = generate_synthetic_paf(size);

            b.iter_with_setup(
                || {
                    let temp_dir = TempDir::new().unwrap();
                    let input_paf = temp_dir.path().join("input.paf");
                    let output_paf = temp_dir.path().join("output.paf");
                    fs::write(&input_paf, &paf_data).unwrap();
                    (temp_dir, input_paf, output_paf)
                },
                |(temp_dir, input_paf, output_paf)| {
                    Command::new("cargo")
                        .args([
                            "run",
                            "--release",
                            "--quiet",
                            "--bin",
                            "sweepga",
                            "--",
                            input_paf.to_str().unwrap(),
                            "-n",
                            "1:1",
                        ])
                        .stdout(fs::File::create(black_box(&output_paf)).unwrap())
                        .status()
                        .unwrap();

                    drop(temp_dir);
                },
            );
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_paf_filtering_pipeline,
    bench_one_to_one_filtering
);

criterion_main!(benches);
