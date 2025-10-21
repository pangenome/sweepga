/// alnstats - Statistics for alignment files (PAF, 1aln)
///
/// Calculates per-genome-pair coverage, mapping counts, and other statistics
/// for alignment files. Future support for 1aln format planned.
use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashMap;
use std::io::BufRead;

#[derive(Parser)]
#[clap(
    name = "alnstats",
    about = "Statistics for alignment files (PAF, 1aln)"
)]
struct Args {
    /// First alignment file (PAF format)
    file1: String,

    /// Optional second file for comparison
    file2: Option<String>,

    /// Show detailed per-genome-pair statistics
    #[clap(short = 'd', long)]
    detailed: bool,
}

#[derive(Debug, Default)]
struct AlignmentStats {
    total_mappings: usize,
    total_bases: u64,
    total_matches: u64,
    self_mappings: usize,
    inter_chromosomal: usize,
    inter_genome: usize,
    chr_pair_count: usize,
    genome_pair_bases: HashMap<(String, String), u64>,
    genome_pair_matches: HashMap<(String, String), u64>,
    genome_sizes: HashMap<String, u64>,
}

impl AlignmentStats {
    fn calculate_coverage_stats(&self) -> CoverageStats {
        // Calculate total genome sizes per genome prefix
        let mut genome_totals: HashMap<String, u64> = HashMap::new();
        for (seq_name, &size) in &self.genome_sizes {
            let genome = extract_genome_prefix(seq_name);
            *genome_totals.entry(genome).or_insert(0) += size;
        }

        // Calculate coverage for each genome pair
        let mut coverages = Vec::new();
        for ((q_genome, t_genome), &bases) in &self.genome_pair_bases {
            if let Some(&genome_size) = genome_totals.get(q_genome) {
                let coverage = 100.0 * bases as f64 / genome_size as f64;
                coverages.push((q_genome.clone(), t_genome.clone(), coverage, bases));
            }
        }

        let avg_coverage = if !coverages.is_empty() {
            coverages.iter().map(|(_, _, c, _)| c).sum::<f64>() / coverages.len() as f64
        } else {
            0.0
        };

        let above_95_pct = coverages.iter().filter(|(_, _, c, _)| *c > 95.0).count();

        CoverageStats {
            avg_coverage,
            genome_pairs: coverages.len(),
            above_95_pct,
            per_pair: coverages,
        }
    }

    fn calculate_avg_identity(&self) -> f64 {
        if self.total_bases > 0 {
            self.total_matches as f64 / self.total_bases as f64
        } else {
            0.0
        }
    }
}

#[derive(Debug)]
struct CoverageStats {
    avg_coverage: f64,
    genome_pairs: usize,
    above_95_pct: usize,
    per_pair: Vec<(String, String, f64, u64)>, // (query_genome, target_genome, coverage%, bases)
}

/// Extract genome prefix from sequence name
/// "SGDref#1#chrI" -> "SGDref#1#"
fn extract_genome_prefix(seq_name: &str) -> String {
    if let Some(last_pos) = seq_name.rfind('#') {
        seq_name[..=last_pos].to_string()
    } else {
        seq_name.to_string()
    }
}

/// Parse a PAF file and collect statistics
fn parse_paf(path: &str) -> Result<AlignmentStats> {
    let reader = sweepga::paf::open_paf_input(path).context(format!("Failed to open {path}"))?;

    let mut stats = AlignmentStats::default();
    let mut chr_pairs = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 11 {
            continue; // Skip invalid lines
        }

        let query = fields[0];
        let query_len: u64 = fields[1].parse().context("Invalid query length")?;
        let query_start: u64 = fields[2].parse().context("Invalid query start")?;
        let query_end: u64 = fields[3].parse().context("Invalid query end")?;
        let target = fields[5];
        let target_len: u64 = fields[6].parse().context("Invalid target length")?;
        let matches: u64 = fields[9].parse().context("Invalid match count")?;
        let _block_len: u64 = fields[10].parse().context("Invalid block length")?;

        stats.total_mappings += 1;
        let mapping_len = query_end - query_start;
        stats.total_bases += mapping_len;
        stats.total_matches += matches;

        // Track genome sizes
        stats.genome_sizes.insert(query.to_string(), query_len);
        stats.genome_sizes.insert(target.to_string(), target_len);

        // Extract genome prefixes
        let q_genome = extract_genome_prefix(query);
        let t_genome = extract_genome_prefix(target);

        // Count mapping types
        if query == target {
            stats.self_mappings += 1;
        } else if q_genome != t_genome {
            stats.inter_genome += 1;

            // Track genome pair bases and matches
            let pair = (q_genome, t_genome);
            *stats.genome_pair_bases.entry(pair.clone()).or_insert(0) += mapping_len;
            *stats.genome_pair_matches.entry(pair).or_insert(0) += matches;
        } else {
            // Same genome, different chromosomes
            stats.inter_chromosomal += 1;
        }

        // Track chromosome pairs
        chr_pairs.insert((query.to_string(), target.to_string()));
    }

    stats.chr_pair_count = chr_pairs.len();

    Ok(stats)
}

fn print_stats(path: &str, stats: &AlignmentStats, detailed: bool) {
    let coverage = stats.calculate_coverage_stats();
    let avg_identity = stats.calculate_avg_identity();

    println!("\nStatistics for {path}:");
    println!("{}", "=".repeat(60));
    println!(
        "Total mappings:        {:>12}",
        format_number(stats.total_mappings)
    );
    println!(
        "Total bases:           {:>12}",
        format_number(stats.total_bases as usize)
    );
    println!("Average identity:      {:>11.1}%", avg_identity * 100.0);
    println!(
        "Self mappings:         {:>12}",
        format_number(stats.self_mappings)
    );
    println!(
        "Inter-chromosomal:     {:>12}",
        format_number(stats.inter_chromosomal)
    );
    println!(
        "Inter-genome:          {:>12}",
        format_number(stats.inter_genome)
    );
    println!(
        "Chromosome pairs:      {:>12}",
        format_number(stats.chr_pair_count)
    );
    println!("Genome pairs:          {:>12}", coverage.genome_pairs);
    println!("Average coverage:      {:>11.1}%", coverage.avg_coverage);
    println!(
        "Pairs >95% coverage:   {:>12}",
        format!("{}/{}", coverage.above_95_pct, coverage.genome_pairs)
    );

    if detailed && !coverage.per_pair.is_empty() {
        println!("\nPer-genome-pair statistics:");
        println!("{}", "-".repeat(60));
        let mut pairs = coverage.per_pair.clone();
        pairs.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap()); // Sort by coverage descending

        for (q, t, cov, bases) in pairs {
            // Calculate identity for this pair
            let pair_matches = stats
                .genome_pair_matches
                .get(&(q.clone(), t.clone()))
                .unwrap_or(&0);
            let identity = if bases > 0 {
                *pair_matches as f64 / bases as f64 * 100.0
            } else {
                0.0
            };

            println!(
                "{:20} -> {:20} {:6.1}% cov, {:6.1}% id, {:>10} bp",
                q.trim_end_matches('#'),
                t.trim_end_matches('#'),
                cov,
                identity,
                format_number(bases as usize)
            );
        }
    }
}

fn compare_stats(file1: &str, file2: &str, stats1: &AlignmentStats, stats2: &AlignmentStats) {
    let coverage1 = stats1.calculate_coverage_stats();
    let coverage2 = stats2.calculate_coverage_stats();
    let identity1 = stats1.calculate_avg_identity();
    let identity2 = stats2.calculate_avg_identity();

    println!("\nComparison: {file1} vs {file2}");
    println!("{}", "=".repeat(60));

    print_comparison("Mappings", stats1.total_mappings, stats2.total_mappings);
    print_comparison(
        "Total bases",
        stats1.total_bases as usize,
        stats2.total_bases as usize,
    );

    println!("\nAverage identity:");
    println!("  {:30} {:>11.1}%", file1, identity1 * 100.0);
    println!("  {:30} {:>11.1}%", file2, identity2 * 100.0);
    println!(
        "  {:30} {:>+10.1}%",
        "Change",
        (identity2 - identity1) * 100.0
    );

    print_comparison(
        "Inter-chromosomal",
        stats1.inter_chromosomal,
        stats2.inter_chromosomal,
    );
    print_comparison(
        "Chromosome pairs",
        stats1.chr_pair_count,
        stats2.chr_pair_count,
    );

    println!("\nAverage genome pair coverage:");
    println!("  {:30} {:>11.1}%", file1, coverage1.avg_coverage);
    println!("  {:30} {:>11.1}%", file2, coverage2.avg_coverage);
    println!(
        "  {:30} {:>+10.1}%",
        "Change",
        coverage2.avg_coverage - coverage1.avg_coverage
    );

    println!("\nGenome pairs with >95% coverage:");
    println!(
        "  {:30} {:>12}",
        file1,
        format!("{}/{}", coverage1.above_95_pct, coverage1.genome_pairs)
    );
    println!(
        "  {:30} {:>12}",
        file2,
        format!("{}/{}", coverage2.above_95_pct, coverage2.genome_pairs)
    );
}

fn print_comparison(label: &str, val1: usize, val2: usize) {
    println!("\n{label}:");
    println!("  {:30} {:>12}", "Before", format_number(val1));
    println!("  {:30} {:>12}", "After", format_number(val2));

    let diff = val2 as i64 - val1 as i64;
    let pct = if val1 > 0 {
        100.0 * diff as f64 / val1 as f64
    } else {
        0.0
    };
    println!(
        "  {:30} {:>12} ({:+.1}%)",
        "Change",
        format_signed(diff),
        pct
    );
}

fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

fn format_signed(n: i64) -> String {
    if n >= 0 {
        format!("+{}", format_number(n as usize))
    } else {
        format!("-{}", format_number((-n) as usize))
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let stats1 = parse_paf(&args.file1)?;

    if let Some(file2) = args.file2 {
        let stats2 = parse_paf(&file2)?;
        compare_stats(&args.file1, &file2, &stats1, &stats2);
    } else {
        print_stats(&args.file1, &stats1, args.detailed);
    }

    Ok(())
}
