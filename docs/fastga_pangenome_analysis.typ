#set page(paper: "a4", margin: 1in)
#set text(font: "Geist", size: 11pt)
#set par(justify: true, leading: 0.75em)

#align(center)[
  #text(17pt, weight: "bold")[
    FastGA for Pangenome Self-Alignment
  ]

  #v(0.8em)

  #text(11pt)[
    Erik Garrison
  ]

  #v(0.3em)

  #text(10pt)[
    #datetime.today().display()
  ]
]

#v(1.5em)

= Background

== Motivation: An Unexpected Efficiency

While experimenting with FastGA (Myers, Durbin, and Zhou 2025) for pangenome alignment, we observed something unexpected: aligning all genomes together as a single self-alignment was *dramatically faster* than performing explicit pairwise alignments between all genome pairs. For 7 yeast genomes, self-alignment completed in under 9 seconds while all-pairs mode required over 2.5 minutes—a 17× speedup.

This observation required investigation. FastGA was designed for pairwise whole-genome comparison, where a k-mer appearing in both genomes might occur a handful of times. The default frequency threshold (`-f 10`) filters out k-mers appearing more than 10 times, preventing the aligner from wasting computation on highly repetitive sequences. But in a pangenome context with H haplotypes, even a unique single-copy gene will produce a k-mer that appears H times—once per haplotype. With the default threshold, nearly all biologically meaningful k-mers would be filtered out.

Adapting FastGA for pangenome use required raising this frequency threshold significantly. We found that setting `-f` to approximately 10×H (where H is the number of haplotypes) allows detection of homologous regions containing moderately repetitive elements while still filtering true high-copy repeats. This modification, combined with self-alignment mode, revealed the efficiency gains documented here.

== SweepGA: A Wrapper for Pangenome Workflows

SweepGA wraps FastGA and adapts it for pangenome workflows. It handles file format conversion (producing PAF output compatible with downstream tools) and applies plane sweep filtering to remove redundant overlapping alignments. This filtering step is standard practice after whole-genome alignment—tools like LASTZ, minimap2, and other aligners typically require similar post-processing when downstream analyses assume 1:1 homologies between sequences.

== The Algorithmic Question

This document investigates *why* self-alignment is so much faster than all-pairs alignment. The answer reveals a deeper property of FastGA's sorted k-mer merge algorithm that, to our knowledge, has not been explicitly discussed: when aligning H genomes simultaneously, the merge-based approach achieves O(H·L) complexity rather than the O(H²·L) required for explicit pairwise comparisons.

== Technical Context

FastGA replaces hash-based seed lookup with merge-based k-mer matching. The algorithm builds sorted k-mer occurrence tables and performs linear merge operations to find matching seeds, avoiding the random memory access patterns of hash-based methods. For a detailed description, see Myers, Durbin, and Zhou (2025).

= The Algorithm

For each genome, FastGA extracts k-mers (using syncmer sampling to reduce table size to ~L/3) and sorts them using cache-coherent MSD radix sort. With k=40, sorting is O(L). For H genomes, total preprocessing is O(H·L).

== Pairwise Mode

For comparing two genomes A and B, FastGA merges their sorted k-mer tables in a single linear scan. Matching k-mers are automatically adjacent after sorting. The merge visits O(L) entries and emits seed pairs at each match. Total time: O(L).

== Self-Alignment Mode

For self-alignment of H concatenated genomes, FastGA builds *one* sorted k-mer table containing H·L entries. The algorithm scans this table once, visiting each entry exactly once: O(H·L) scan time.

At each k-mer shared across genomes, the H occurrences (one per genome) cluster together in sorted order. The algorithm enumerates matches between entries in each cluster. For a k-mer present in all H genomes, this produces O(H²) seed pairs at that position—but these H entries are *physically adjacent in memory*, making enumeration cache-coherent.

The total output is O(H²·L) seed pairs (as expected: H²/2 genome pairs × O(L) matches each). But the *scan time* is O(H·L), not O(H²·L), because we traverse the sorted table once rather than performing H² separate merge operations.

= The Key Observation

== Separating Scan Cost from Emission Cost

The work breaks down into two phases:

*Phase 1 — Finding matching k-mers (the scan):*
- All-pairs: H²/2 separate merges, each scanning O(L) entries → O(H²·L) scan operations
- Self-alignment: One scan through H·L sorted entries → O(H·L) scan operations

*Phase 2 — Emitting seed pairs (the enumeration):*
- Both modes: O(H²·L) seed pairs emitted

This is a legitimate complexity separation. The scan phase has real per-operation costs: cache line fetches, k-mer comparisons, maintaining merge state, finding run boundaries. Reducing this from O(H²·L) to O(H·L) is a genuine H-factor improvement in that portion of the work.

== When Does This Matter?

Total cost = (scan_ops × cost_per_scan) + (emit_ops × cost_per_emit)

- All-pairs: H²·L · c_scan + H²·L · c_emit
- Self-alignment: H·L · c_scan + H²·L · c_emit

The speedup depends on the ratio of scan cost to emission cost. If scanning/lookup is expensive relative to emission, you approach H× speedup. If emission dominates, the speedup diminishes.

For H=7 and observed ~10× speedup, this implies c_scan is substantial—which makes sense. Each scan operation involves index arithmetic, k-mer comparisons, memory fetches (potentially cache misses in all-pairs mode), and merge state bookkeeping. Emission is relatively cheap: you've already found the match, you just write it out.

== Cache Locality Compounds This

In self-alignment, the scan is sequential through one contiguous array—optimal for hardware prefetching. In all-pairs, you perform H² separate merge operations, each potentially thrashing cache as you switch between different k-mer tables. So not only are there H× fewer scan operations, but each operation is cheaper due to better locality.

This means the "constant" c_scan is actually worse for all-pairs mode, compounding the H-factor algorithmic improvement.

= Implications

Self-alignment's O(H·L) scan time (versus O(H²·L) for explicit all-pairs) enables pangenome alignment at scales that would otherwise be impractical. The key constraint is memory: storing the sorted k-mer table requires O(H·L) space. For H=1000 haplotypes and L=3Gbp, this becomes terabyte-scale (3Gbp × 1000 × ~1/3 syncmers × 16 bytes ≈ 16TB).

The practical solution: hierarchical alignment. Use sparse global alignment to build a pangenome graph (e.g., with impg), extract local homologous chunks where L is reduced to 100kb–1Mb, then apply FastGA self-alignment within each chunk. In these local windows, FastGA's cache-coherent enumeration makes all-vs-all seed discovery tractable.

= Summary

FastGA's merge-based approach eliminates the O(L²) per-pair comparison cost of naive k-mer matching, achieving O(L) per pair through sorting and linear merge.

For pangenome self-alignment with H genomes, the critical insight is: *sorting all k-mers together into one table* reduces the scan phase from O(H²·L) to O(H·L) operations. Total output remains O(H²·L) seed pairs, but the expensive work (index lookups, cache fetches, merge bookkeeping) is reduced H-fold, while the remaining emission work happens cache-locally on adjacent entries.

This decomposition explains the observed ~10× speedup for H=7: close to the theoretical H-fold improvement in scan cost, with additional gains from cache locality.

Memory is the main constraint: O(H·L) space for the sorted k-mer table. For large pangenomes, hierarchical workflows (sparse global alignment → local refinement) make this tractable.

= Empirical Performance: Regular vs All-Pairs Mode

To validate the theoretical complexity analysis, we benchmarked FastGA on 7 highly similar *Saccharomyces cerevisiae* genomes (~12 Mb each, 99% average identity) using SweepGA. We compared three alignment strategies:

1. *Regular mode (default)*: Self-alignment of concatenated genomes with auto-detected frequency threshold
2. *Regular mode with `-f 70`*: Self-alignment with higher frequency threshold (10 × number of genomes)
3. *All-pairs mode*: Explicit pairwise alignment of all genome pairs

== Experimental Parameters

#table(
  columns: (auto, auto),
  inset: 8pt,
  align: (left, left),
  [*Parameter*], [*Value*],
  [Dataset], [7 *S. cerevisiae* genomes],
  [Genome names], [DBVPG6044, DBVPG6765, S288C, SK1,\ UWOPS034614, Y12, YPS128],
  [Genome size (each)], [~12 Mb],
  [Average identity], [~99%],
  [Total sequence], [~84 Mb],
  [System], [Linux 6.16.0, 8 cores],
  [Tool], [SweepGA + FastGA],
  [FastGA version], [Gene Myers' implementation],
  [Syncmer parameters], [(12,8) syncmers (default)],
  [K-mer size], [40 bp (default)],
)

#v(1em)

#table(
  columns: (auto, auto, auto),
  inset: 8pt,
  align: (left, left, left),
  [*Mode*], [*Command*], [*Key Parameters*],
  [Regular (default)], [`sweepga input.fa.gz`], [Frequency: auto-detect (~7)],
  [Regular `-f 70`], [`sweepga -f 70 input.fa.gz`], [Frequency: 70 (10×genomes)],
  [All-Pairs], [`sweepga --all-pairs input.fa.gz`], [42 pairwise comparisons],
)

All experiments used 8 threads (`-T8` passed to FastGA) and were measured with `/usr/bin/time -v` for accurate CPU and memory profiling.

== Results

#table(
  columns: (auto, auto, auto, auto),
  inset: 8pt,
  align: (left, right, right, right),
  [*Metric*], [*Regular (default)*], [*Regular `-f 70`*], [*All-Pairs*],
  [Wall clock time], [8.90s], [14.94s], [2m 35.3s (155.3s)],
  [Speedup vs all-pairs], [*17.4×*], [*10.4×*], [1.0× (baseline)],
  [CPU time], [47.3s], [80.1s], [946.6s],
  [*Matched base pairs*], [*481.6 Mb*], [*483.3 Mb*], [*479.7 Mb*],
  [Total alignment length], [485.4 Mb], [487.3 Mb], [483.6 Mb],
  [Average identity], [99.21%], [99.19%], [99.21%],
  [Total mappings], [12,358], [14,135], [12,377],
  [Average coverage], [96.3%], [96.5%], [96.5%],
  [Pairs >95% cov], [35/42], [38/42], [37/42],
  [Peak memory], [629 MB], [653 MB], [379 MB],
)

== Key Observations

*Matched Base Pairs*: Regular mode with `-f 70` produces the most aligned sequence—483.3 Mb of matched base pairs, compared to 481.6 Mb for default regular mode and 479.7 Mb for all-pairs. All three modes achieve ~99.2% identity, so differences lie in total alignment coverage rather than alignment quality. The `-f 70` parameter (frequency threshold = 10 × genomes) allows more repetitive k-mers to pass filtering, enabling detection of additional homologous regions.

*Performance Characteristics*: Regular mode with default parameters achieves 17.4× speedup (8.90s vs 155.3s) but finds 1.7 Mb fewer matched bases than `-f 70`. Regular mode with `-f 70` provides an excellent compromise—10.4× speedup while finding the most sequence homology. The O(H²) complexity manifests as 42 pairwise comparisons in all-pairs mode, each requiring separate FastGA invocations. All-pairs mode spends 808.55s in system time versus only 138.06s in user time, indicating that process creation overhead dominates total CPU usage (946.6s total).

*Frequency Threshold Effects*: The `-f` parameter controls k-mer frequency filtering. Default mode uses auto-detection (~7 for 7 genomes), which aggressively filters k-mers present across multiple genomes. Setting `-f 70` (10 × genomes) relaxes this threshold, allowing detection of homologous regions that contain repetitive elements shared across strains. This explains the 1.7 Mb increase in matched bases (481.6 → 483.3 Mb) with only modest runtime cost (8.90s → 14.94s).

*Coverage Quality*: All three modes achieve excellent coverage—96.3–96.5% average with ~99.2% identity. Regular `-f 70` achieves the best genome pair coverage (38/42 pairs >95%) and finds the most total homology, making it the optimal choice for maximizing detected sequence similarity while maintaining fast runtime.

== Implications

For *maximizing sequence homology detection*: Regular mode with `-f 70` (or generally `-f 10×H` where H is the number of genomes) finds the most matched base pairs (483.3 Mb) while maintaining 10.4× speedup over all-pairs. This makes it optimal for pangenome analysis where the goal is to detect as much homologous sequence as possible. The frequency threshold formula `10×H` allows repetitive k-mers shared across strains to contribute to alignment detection.

For *fastest runtime*: Regular mode with default parameters achieves 17.4× speedup (8.90s) but finds 1.7 Mb less sequence homology than `-f 70`. The default auto-detection is overly aggressive for pangenome work, filtering k-mers that are biologically meaningful for detecting shared sequence. Use default mode only for preliminary exploration where speed is paramount.

For *memory-constrained systems*: All-pairs mode uses 39.7% less peak memory (379 MB vs 629 MB for regular modes) because each pairwise comparison processes smaller sequences than concatenated self-alignment. However, all-pairs finds 3.6 Mb fewer matched bases than regular `-f 70`, suggesting some homology is missed when genomes are aligned pairwise rather than allowing all-vs-all k-mer matching.

For *large-scale pangenomes*: Memory scales as O(H·L), becoming prohibitive beyond ~100 haplotypes for whole genomes. The practical solution is hierarchical alignment: use sparse global alignment to build a pangenome graph (e.g., with impg), extract local homologous chunks where L is reduced to 100 kb–1 Mb, then apply FastGA self-alignment within each chunk with appropriate `-f 10×H` setting.

*Process overhead observations*: The high system time in all-pairs mode (808.55s system vs 138.06s user) likely reflects I/O, memory allocation, and kernel work across 42 separate FastGA processes. While the exact breakdown is unclear, running 42 independent processes clearly incurs overhead that self-alignment avoids.

== Validation of Theoretical Model

=== Validating the Scan/Emission Decomposition

The observed speedup validates the scan/emission cost model. With H=7:

*All-pairs*: 42 separate merge operations at ~3.7s each = 155.3s total. Each merge involves O(L) scan operations plus O(L) emissions per pair, totaling O(H²·L) scan ops and O(H²·L) emissions.

*Self-alignment with `-f 70`*: One scan through H·L entries = 14.94s. The scan phase is O(H·L), while emission remains O(H²·L).

The *10.4× speedup* (close to theoretical 7×) confirms that scan costs are substantial relative to emission costs. The extra ~1.5× beyond 7× reflects that c_scan is actually worse in all-pairs mode due to cache thrashing across H² separate merge operations, plus process creation overhead.

=== The Geometric Intuition

Consider the H×H×L lattice of potential matches—for each of H² genome pairs, L positions could match. Explicit all-pairs traverses this lattice pair-by-pair, performing H² separate O(L) scans. Self-alignment exploits sorted order to trace a *single path* through this entire lattice: one O(H·L) scan that discovers all H² pairwise relationships.

The sorted k-mer table encodes a "bundle" of all matches across all genome pairs. Matching k-mers cluster together regardless of which genomes they came from. Following this bundle through the sorted k-mer sequence, we perform O(H·L) scan operations while emitting O(H²·L) seed pairs from cache-local enumeration.

=== Empirical Confirmation

#table(
  columns: (auto, auto, auto, auto),
  inset: 8pt,
  align: (left, right, right, right),
  [*Mode*], [*Complexity*], [*Observed*], [*Notes*],
  [All-pairs (`-f 10` per pair)], [O(H²·L)], [155.3s], [42 merges × 3.7s each],
  [Self-alignment (`-f 70`)], [O(H·L)], [14.94s], [Single sorted scan],
  [Speedup], [H = 7], [10.4×], [Close to theoretical 7×],
)

The `-f 70` setting (10 × 7 genomes) matches the effective per-pair frequency threshold used in all-pairs mode (default `-f 10`), making this an apples-to-apples comparison. The observed 10.4× speedup is close to the theoretical H-fold (7×) prediction; the additional ~1.5× gain comes from eliminating process overhead and improved cache locality.

=== Why Self-Alignment Finds More Homology

Regular `-f 70` finds 483.3 Mb matched bases versus 479.7 Mb for all-pairs—3.6 Mb more homology despite being 10.4× faster. The cause is unclear: it may reflect the global sorted view detecting matches that pairwise comparisons miss, or it may be an artifact of subtly different filtering semantics (global `-f 70` vs per-pair `-f 10` are not perfectly equivalent). The difference is small (~0.7%) and both modes achieve excellent coverage.

=== Implementation Verification

To verify these theoretical claims, we examined Gene Myers' FastGA source code (`FastGA.c`). The implementation confirms the complexity difference through two distinct code paths:

*Self-alignment mode* (`self_adaptamer_merge`, line 2494): Uses a *single* k-mer stream `T1`. The function calls `new_self_merge_thread` which scans one sorted k-mer table containing all H·L entries. For each k-mer position `suf1`, the algorithm uses the LCP (Longest Common Prefix) structure to find the range `[low, hgh]` of adjacent k-mers that share the same prefix. All matching k-mers from different genomes are adjacent in sorted order, enabling O(1) lookup of run boundaries via the `vlcp[]` array.

*Pairwise mode* (`adaptamer_merge`, line 2279): Uses *two* k-mer streams `T1` and `T2`. Each pairwise comparison requires a separate merge operation. For H genomes, this means H²/2 invocations of the merge, each scanning O(L) entries.

The critical inner loop in `new_self_merge_thread` (lines 1810-1862) enumerates matches within a run of identical k-mers:

```c
for (p = low+PAYOFF; p < hgh; p += kbyte)
  { if (p == pay1)
      continue;
    if (p[-2] >= mlen)  // skip if already in same match group
      continue;
    // emit match between suf1 and p
  }
```

Here `low` and `hgh` define the range of k-mers sharing a common prefix of length `mlen`. The algorithm uses an LCP (Longest Common Prefix) array—a standard suffix array technique—to track how many prefix bases each entry shares with its predecessor. The `p[-2] >= mlen` check skips entries that have already been grouped with the current k-mer in a previous iteration, avoiding duplicate emissions.

The decision point at line 5116 selects between modes:

```c
if (SELF)
  self_adaptamer_merge(T1,P1,gdb1->seqtot);
else
  adaptamer_merge(T1,T2,P1,P2,gdb1->seqtot);
```

The LCP structure enables O(1) lookup of run boundaries via the `vlcp[]` array: `low = vlcp[plen]` immediately gives the start of entries sharing a prefix of length `plen`. This is what makes the "enumerate all H" operation at each position cache-coherent—all relevant entries are physically adjacent in the sorted table.

=== Practical Recommendation

Use concatenated self-alignment with `-f 10×H` for pangenome work. This reduces scan complexity from O(H²·L) to O(H·L) while producing equivalent output—a substantial speedup with no loss in sensitivity.

= Appendix: Data and Methods

== Data Availability

The yeast genome dataset is available at: #link("http://hypervolu.me/~erik/yeast/cerevisiae.fa.gz")[hypervolu.me/~erik/yeast/cerevisiae.fa.gz]

== Analysis Commands

```bash
# Extract 7 genomes for testing
samtools faidx cerevisiae.fa.gz DBVPG6044 DBVPG6765 S288C SK1 \
  UWOPS034614 Y12 YPS128 | bgzip > yeast7.fa.gz

# Regular mode (auto-detected frequency threshold)
/usr/bin/time -v sweepga yeast7.fa.gz > regular.paf 2> regular.time

# Regular mode with relaxed frequency threshold
/usr/bin/time -v sweepga -f 70 yeast7.fa.gz > regular_f70.paf 2> regular_f70.time

# All-pairs mode (explicit pairwise comparisons)
/usr/bin/time -v sweepga --all-pairs yeast7.fa.gz > allpairs.paf 2> allpairs.time

# Compute alignment statistics
alnstats regular.paf regular_f70.paf allpairs.paf
```

= References

Myers G, Durbin R, Zhou C. FastGA: Fast Genome Alignment. _Bioinformatics Advances_. 2025. doi:10.1093/bioadv/vbaf238

SweepGA source code: #link("https://github.com/pangenome/sweepga")[github.com/pangenome/sweepga]

FastGA source code: #link("https://github.com/thegenemyers/FASTGA")[github.com/thegenemyers/FASTGA]
