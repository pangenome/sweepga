# SweepGA Usage Guide

## Overview
SweepGA integrates FastGA alignment with sophisticated filtering algorithms from wfmash. It can either:
1. Run FastGA alignment on FASTA files and filter the results
2. Filter existing PAF alignments

## Usage Modes

### FastGA Alignment Mode (with FASTA files)

#### Self-alignment (one file)
```bash
sweepga genome.fa -o output.paf
```
Aligns a genome to itself, useful for finding repeats and duplications.

#### Pairwise alignment (two files)
```bash
sweepga target.fa query.fa -o output.paf
```
Aligns query sequences to target sequences. Note the order: **target first, query second**.

### PAF Filtering Mode (existing alignments)

#### From file
```bash
sweepga -i alignments.paf -o filtered.paf
```

#### From stdin
```bash
cat alignments.paf | sweepga -o filtered.paf
```

## Key Parameters

### Filtering Options
- `-b, --block-length <N>`: Minimum alignment block length (supports k/m/g suffixes)
- `-p, --overlap <FLOAT>`: Maximum overlap ratio for plane sweep (default: 0.95)
- `-n, --num-mappings <MODE>`: Mapping filter mode (e.g., "1:1", "1:many", "N")
- `-j, --scaffold-jump <N>`: Gap distance for scaffolding (0 = disable)
- `-s, --scaffold-mass <N>`: Minimum scaffold length

### Performance Options
- `-t, --threads <N>`: Number of threads (default: 8)
- `--quiet`: Suppress progress output

## Examples

### Find all repeats in a genome
```bash
sweepga genome.fa -j 10k -s 50k -o repeats.paf
```

### Align two assemblies with strict filtering
```bash
sweepga reference.fa assembly.fa -n 1:1 -b 10k -o filtered.paf
```

### Filter existing alignments with scaffolding
```bash
sweepga -i raw_alignments.paf -j 100k -s 10k -o scaffolded.paf
```

## Notes
- When using FastGA mode, the alignment is run through embedded FFI bindings (no external FastGA needed)
- Temporary PAF files are automatically cleaned up after processing
- The tool respects PanSN naming conventions when grouping sequences