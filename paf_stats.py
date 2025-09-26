#!/usr/bin/env python3
"""Quick PAF statistics for comparing filtering effects"""

import sys
from collections import defaultdict

def get_paf_stats(filename):
    """Compute statistics for a PAF file"""

    total_mappings = 0
    total_bases = 0
    inter_chromosomal = 0
    inter_genome = 0
    self_mappings = 0

    # Track coverage by genome pair
    genome_pair_bases = defaultdict(int)
    chr_pair_mappings = defaultdict(int)
    genome_sizes = {}

    with open(filename, 'r') as f:
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

            total_mappings += 1
            mapping_len = query_end - query_start
            total_bases += mapping_len

            # Extract genome and chromosome names
            if '#' in query:
                q_parts = query.split('#')
                q_genome = '#'.join(q_parts[:2]) if len(q_parts) >= 2 else query
                q_chr = q_parts[2] if len(q_parts) >= 3 else query
            else:
                q_genome = query
                q_chr = query

            if '#' in target:
                t_parts = target.split('#')
                t_genome = '#'.join(t_parts[:2]) if len(t_parts) >= 2 else target
                t_chr = t_parts[2] if len(t_parts) >= 3 else target
            else:
                t_genome = target
                t_chr = target

            # Track genome sizes
            genome_sizes[query] = query_len
            genome_sizes[target] = target_len

            # Count mapping types
            if query == target:
                self_mappings += 1
            elif q_genome != t_genome:
                inter_genome += 1
                genome_pair_bases[(q_genome, t_genome)] += mapping_len
            elif q_chr != t_chr:
                inter_chromosomal += 1

            # Track chromosome pairs
            chr_pair_mappings[(query, target)] += 1

    # Calculate genome sizes
    genome_totals = defaultdict(int)
    for seq_name, size in genome_sizes.items():
        if '#' in seq_name:
            genome = '#'.join(seq_name.split('#')[:2])
        else:
            genome = seq_name
        genome_totals[genome] += size

    # Calculate coverage for genome pairs
    genome_coverages = []
    for (q_genome, t_genome), bases in genome_pair_bases.items():
        if q_genome in genome_totals:
            coverage = 100.0 * bases / genome_totals[q_genome]
            genome_coverages.append(coverage)

    avg_coverage = sum(genome_coverages) / len(genome_coverages) if genome_coverages else 0

    return {
        'total_mappings': total_mappings,
        'total_bases': total_bases,
        'self_mappings': self_mappings,
        'inter_chromosomal': inter_chromosomal,
        'inter_genome': inter_genome,
        'chr_pairs': len(chr_pair_mappings),
        'genome_pairs': len(genome_pair_bases),
        'avg_coverage': avg_coverage,
        'coverages_above_95': sum(1 for c in genome_coverages if c > 95),
        'total_genome_pairs': len(genome_coverages),
    }

def compare_paf_files(file1, file2):
    """Compare two PAF files and show differences"""
    stats1 = get_paf_stats(file1)
    stats2 = get_paf_stats(file2)

    print(f"\nComparison: {file1} vs {file2}")
    print("=" * 60)

    print(f"\nMappings:")
    print(f"  {file1:30s}: {stats1['total_mappings']:,}")
    print(f"  {file2:30s}: {stats2['total_mappings']:,}")
    diff = stats2['total_mappings'] - stats1['total_mappings']
    pct = 100.0 * diff / stats1['total_mappings'] if stats1['total_mappings'] > 0 else 0
    print(f"  {'Change':30s}: {diff:+,} ({pct:+.1f}%)")

    print(f"\nInter-chromosomal mappings:")
    print(f"  {file1:30s}: {stats1['inter_chromosomal']:,}")
    print(f"  {file2:30s}: {stats2['inter_chromosomal']:,}")
    diff = stats2['inter_chromosomal'] - stats1['inter_chromosomal']
    pct = 100.0 * diff / stats1['inter_chromosomal'] if stats1['inter_chromosomal'] > 0 else 0
    print(f"  {'Change':30s}: {diff:+,} ({pct:+.1f}%)")

    print(f"\nChromosome pairs:")
    print(f"  {file1:30s}: {stats1['chr_pairs']:,}")
    print(f"  {file2:30s}: {stats2['chr_pairs']:,}")
    diff = stats2['chr_pairs'] - stats1['chr_pairs']
    pct = 100.0 * diff / stats1['chr_pairs'] if stats1['chr_pairs'] > 0 else 0
    print(f"  {'Change':30s}: {diff:+,} ({pct:+.1f}%)")

    print(f"\nAverage genome pair coverage:")
    print(f"  {file1:30s}: {stats1['avg_coverage']:.1f}%")
    print(f"  {file2:30s}: {stats2['avg_coverage']:.1f}%")
    diff = stats2['avg_coverage'] - stats1['avg_coverage']
    print(f"  {'Change':30s}: {diff:+.1f}%")

    print(f"\nGenome pairs with >95% coverage:")
    print(f"  {file1:30s}: {stats1['coverages_above_95']}/{stats1['total_genome_pairs']}")
    print(f"  {file2:30s}: {stats2['coverages_above_95']}/{stats2['total_genome_pairs']}")

if __name__ == "__main__":
    if len(sys.argv) == 2:
        # Single file stats
        stats = get_paf_stats(sys.argv[1])
        print(f"\nStatistics for {sys.argv[1]}:")
        print("=" * 40)
        print(f"Total mappings:        {stats['total_mappings']:,}")
        print(f"Total bases:           {stats['total_bases']:,}")
        print(f"Self mappings:         {stats['self_mappings']:,}")
        print(f"Inter-chromosomal:     {stats['inter_chromosomal']:,}")
        print(f"Inter-genome:          {stats['inter_genome']:,}")
        print(f"Chromosome pairs:      {stats['chr_pairs']:,}")
        print(f"Average coverage:      {stats['avg_coverage']:.1f}%")
        print(f"Pairs >95% coverage:   {stats['coverages_above_95']}/{stats['total_genome_pairs']}")
    elif len(sys.argv) == 3:
        # Compare two files
        compare_paf_files(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python3 paf_stats.py <paf_file> [<paf_file2>]")
        print("  Single file: show statistics")
        print("  Two files: compare statistics")
        sys.exit(1)