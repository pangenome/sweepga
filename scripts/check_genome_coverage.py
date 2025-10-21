#!/usr/bin/env python3
"""Check actual genome-level coverage in PAF files."""

import sys
from collections import defaultdict

def check_genome_coverage(paf_file):
    """Calculate total coverage between genome pairs."""

    # Track total aligned bases per genome pair
    genome_pair_bases = defaultdict(int)
    genome_sizes = {}

    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            target = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])

            # Extract genome names
            q_genome = '#'.join(query.split('#')[:2]) if '#' in query else query
            t_genome = '#'.join(target.split('#')[:2]) if '#' in target else target

            # Skip self mappings
            if q_genome == t_genome:
                continue

            # Track genome sizes (use max seen for each chromosome)
            genome_sizes[query] = query_len
            genome_sizes[target] = target_len

            # Track aligned bases
            pair = f"{q_genome} -> {t_genome}"
            genome_pair_bases[pair] += (query_end - query_start)

    # Calculate total genome sizes
    genome_total_size = defaultdict(int)
    for seq_name, size in genome_sizes.items():
        genome = '#'.join(seq_name.split('#')[:2]) if '#' in seq_name else seq_name
        genome_total_size[genome] += size

    print(f"\nAnalyzing {paf_file}")
    print("=" * 80)
    print(f"Found {len(genome_pair_bases)} genome pairs")
    print(f"Genome sizes: {len(genome_total_size)} genomes")

    # Calculate coverage for each genome pair
    coverage_stats = []
    for pair, aligned_bases in genome_pair_bases.items():
        q_genome, t_genome = pair.split(' -> ')
        q_size = genome_total_size[q_genome]
        t_size = genome_total_size[t_genome]

        if q_size > 0:
            coverage = 100.0 * aligned_bases / q_size
            coverage_stats.append((coverage, pair, aligned_bases, q_size))

    # Sort by coverage
    coverage_stats.sort(reverse=True)

    print(f"\nTop genome pair coverages:")
    for cov, pair, bases, total in coverage_stats[:20]:
        print(f"  {pair}: {cov:.1f}% ({bases:,} / {total:,} bases)")

    # Summary statistics
    coverages = [c for c, _, _, _ in coverage_stats]
    if coverages:
        avg_cov = sum(coverages) / len(coverages)
        above_90 = sum(1 for c in coverages if c > 90)
        above_95 = sum(1 for c in coverages if c > 95)
        above_99 = sum(1 for c in coverages if c > 99)

        print(f"\nCoverage summary:")
        print(f"  Average coverage: {avg_cov:.1f}%")
        print(f"  Pairs with >90% coverage: {above_90}/{len(coverages)} ({100*above_90/len(coverages):.1f}%)")
        print(f"  Pairs with >95% coverage: {above_95}/{len(coverages)} ({100*above_95/len(coverages):.1f}%)")
        print(f"  Pairs with >99% coverage: {above_99}/{len(coverages)} ({100*above_99/len(coverages):.1f}%)")

if __name__ == "__main__":
    for paf_file in sys.argv[1:]:
        check_genome_coverage(paf_file)