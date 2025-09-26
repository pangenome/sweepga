#!/usr/bin/env python3
"""
Analyze chromosome pair coverage in PAF files.
For yeast genomes, we expect most homologous chromosomes to have 1:1 mappings.
"""

import sys
from collections import defaultdict

def parse_paf_line(line):
    """Parse a PAF line and extract query, target, and lengths."""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None

    query = fields[0]
    query_len = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    target = fields[5]
    target_len = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])

    return {
        'query': query,
        'target': target,
        'query_len': query_len,
        'target_len': target_len,
        'query_cov': (query_end - query_start) / query_len,
        'target_cov': (target_end - target_start) / target_len,
    }

def extract_genome_and_chr(seq_name):
    """Extract genome and chromosome from PanSN format: genome#haplotype#chromosome"""
    parts = seq_name.split('#')
    if len(parts) >= 3:
        genome = f"{parts[0]}#{parts[1]}"
        chromosome = parts[2]
    else:
        genome = seq_name
        chromosome = seq_name
    return genome, chromosome

def analyze_coverage(paf_file):
    """Analyze chromosome pair coverage in a PAF file."""

    # Track alignments between chromosome pairs
    chr_pairs = defaultdict(list)  # (query_chr, target_chr) -> list of alignments
    genome_pairs = defaultdict(set)  # (query_genome, target_genome) -> set of chr pairs

    with open(paf_file, 'r') as f:
        for line in f:
            aln = parse_paf_line(line)
            if not aln:
                continue

            q_genome, q_chr = extract_genome_and_chr(aln['query'])
            t_genome, t_chr = extract_genome_and_chr(aln['target'])

            # Skip self-mappings at genome level
            if q_genome == t_genome:
                continue

            chr_pair = (aln['query'], aln['target'])
            chr_pairs[chr_pair].append(aln)

            genome_pair = (q_genome, t_genome)
            genome_pairs[genome_pair].add((q_chr, t_chr))

    # Analyze coverage for each chromosome pair
    print(f"\nAnalyzing {paf_file}")
    print("=" * 80)

    # Count chromosome pairs by coverage level
    coverage_bins = defaultdict(int)
    chr_coverage = {}

    for (q_chr, t_chr), alignments in chr_pairs.items():
        # Calculate total coverage
        q_genome, q_chr_name = extract_genome_and_chr(q_chr)
        t_genome, t_chr_name = extract_genome_and_chr(t_chr)

        # Merge overlapping alignments to get true coverage
        total_q_cov = sum(a['query_cov'] for a in alignments)
        total_t_cov = sum(a['target_cov'] for a in alignments)
        avg_cov = (total_q_cov + total_t_cov) / 2

        chr_coverage[(q_chr, t_chr)] = avg_cov

        # Bin the coverage
        if avg_cov < 0.1:
            coverage_bins['<10%'] += 1
        elif avg_cov < 0.5:
            coverage_bins['10-50%'] += 1
        elif avg_cov < 0.8:
            coverage_bins['50-80%'] += 1
        elif avg_cov < 0.95:
            coverage_bins['80-95%'] += 1
        elif avg_cov <= 1.05:
            coverage_bins['95-105% (1:1)'] += 1
        else:
            coverage_bins['>105%'] += 1

    # Print summary statistics
    print(f"\nTotal chromosome pairs with alignments: {len(chr_pairs)}")
    print(f"Total genome pairs: {len(genome_pairs)}")

    print("\nChromosome pair coverage distribution:")
    for bin_name in ['<10%', '10-50%', '50-80%', '80-95%', '95-105% (1:1)', '>105%']:
        count = coverage_bins[bin_name]
        pct = 100 * count / max(1, len(chr_pairs))
        print(f"  {bin_name:15s}: {count:4d} pairs ({pct:5.1f}%)")

    # Find expected homologous pairs (same chromosome name)
    homologous_pairs = {}
    for (q_chr, t_chr), cov in chr_coverage.items():
        q_genome, q_chr_name = extract_genome_and_chr(q_chr)
        t_genome, t_chr_name = extract_genome_and_chr(t_chr)

        if q_chr_name == t_chr_name:  # Same chromosome (e.g., chrI to chrI)
            homologous_pairs[(q_chr, t_chr)] = cov

    print(f"\nHomologous chromosome pairs (same chr name): {len(homologous_pairs)}")

    # Check coverage of homologous pairs
    good_coverage = sum(1 for cov in homologous_pairs.values() if 0.95 <= cov <= 1.05)
    print(f"  With ~100% coverage (0.95-1.05): {good_coverage} ({100*good_coverage/max(1,len(homologous_pairs)):.1f}%)")

    # List problematic homologous pairs
    problematic = [(pair, cov) for pair, cov in homologous_pairs.items() if cov < 0.95 or cov > 1.05]
    if problematic:
        print("\n  Problematic homologous pairs (coverage != ~1.0):")
        for (q_chr, t_chr), cov in sorted(problematic, key=lambda x: x[1])[:10]:
            q_genome, q_chr_name = extract_genome_and_chr(q_chr)
            t_genome, t_chr_name = extract_genome_and_chr(t_chr)
            print(f"    {q_genome} {q_chr_name} -> {t_genome} {t_chr_name}: {cov:.1%} coverage")

    # Check for missing homologous pairs
    all_queries = set()
    all_targets = set()
    for q_chr, t_chr in chr_pairs:
        all_queries.add(q_chr)
        all_targets.add(t_chr)

    print(f"\nUnique query sequences: {len(all_queries)}")
    print(f"Unique target sequences: {len(all_targets)}")

    # For each query chromosome, check if it has alignments to the same chromosome in other genomes
    missing_homologs = []
    for q_chr in all_queries:
        q_genome, q_chr_name = extract_genome_and_chr(q_chr)
        has_homolog = False

        for t_chr in all_targets:
            t_genome, t_chr_name = extract_genome_and_chr(t_chr)

            if q_genome != t_genome and q_chr_name == t_chr_name:
                if (q_chr, t_chr) in chr_pairs:
                    has_homolog = True
                    break

        if not has_homolog:
            missing_homologs.append(q_chr)

    if missing_homologs:
        print(f"\nChromosomes missing homologous alignments: {len(missing_homologs)}")
        for chr_name in missing_homologs[:10]:
            genome, chr = extract_genome_and_chr(chr_name)
            print(f"  {genome} {chr}")

    return chr_pairs, genome_pairs, chr_coverage

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_coverage.py <paf_file> [<paf_file2> ...]")
        sys.exit(1)

    for paf_file in sys.argv[1:]:
        analyze_coverage(paf_file)