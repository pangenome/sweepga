use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use crate::mapping::ChainStatus;
use crate::paf::open_paf_input;
use crate::plane_sweep_exact::PlaneSweepMapping;
use crate::sequence_index::SequenceIndex;

/// Scoring function for plane sweep
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScoringFunction {
    Identity,          // Identity only
    Length,            // Length only
    LengthIdentity,    // Length * Identity
    LogLengthIdentity, // log(Length) * Identity (default)
    Matches,           // Total matches only (gap-neutral)
}

/// Filtering mode
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterMode {
    OneToOne,   // 1:1 - best mapping per query AND per target
    OneToMany,  // 1:N - best mapping per query, N per target
    ManyToMany, // N:N - N mappings per query and per target
}

/// Filter configuration
#[derive(Clone)]
#[allow(dead_code)]
pub struct FilterConfig {
    pub chain_gap: u64,        // -c/--chain-jump
    pub min_block_length: u64, // -l/--block-length

    // Primary mapping filter (applied to raw mappings before scaffold creation)
    pub mapping_filter_mode: FilterMode, // Default: N:N (no filtering)
    pub mapping_max_per_query: Option<usize>,
    pub mapping_max_per_target: Option<usize>,
    pub plane_sweep_secondaries: usize, // -n parameter: number of secondaries to keep in plane sweep

    // Scaffold filter (applied to scaffold chains)
    pub scaffold_filter_mode: FilterMode, // Default: 1:1
    pub scaffold_max_per_query: Option<usize>,
    pub scaffold_max_per_target: Option<usize>,

    pub overlap_threshold: f64,   // -O/--overlap
    pub sparsity: f64,            // -x/--sparsify
    pub no_merge: bool,           // -M/--no-merge
    pub scaffold_gap: u64,        // -j/--scaffold-jump
    pub min_scaffold_length: u64, // -S/--scaffold-mass
    pub scaffold_overlap_threshold: f64,
    pub scaffold_max_deviation: u64, // -D/--scaffold-dist
    pub prefix_delimiter: char,
    pub skip_prefix: bool,

    // Scoring and identity filtering
    pub scoring_function: ScoringFunction,
    pub min_identity: f64, // Minimum block identity threshold (0.0-1.0)
    pub min_scaffold_identity: f64, // Minimum scaffold identity threshold (0.0-1.0)
}

/// Record metadata for filtering without modifying records
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RecordMeta {
    pub rank: usize, // 0-based index in original file
    pub query_name: String,
    pub target_name: String,
    pub query_start: u64,
    pub query_end: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub block_length: u64,
    pub identity: f64,         // Block identity: matches / alignment_length
    pub matches: u64,          // Number of matching bases
    pub alignment_length: u64, // Total alignment length (including gaps)
    pub strand: char,
    pub chain_id: Option<String>,
    pub chain_status: ChainStatus,
    pub discard: bool,
    pub overlapped: bool,
}

/// Compact record metadata using sequence IDs instead of strings
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct CompactRecordMeta {
    rank: usize, // 0-based index in original file
    query_id: u32,
    target_id: u32,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
    block_length: u64,
    identity: f64,
    strand: char,
    chain_id: Option<String>, // Still a string for now
    chain_status: ChainStatus,
    discard: bool,
    overlapped: bool,
}

impl CompactRecordMeta {
    /// Convert from regular RecordMeta using a sequence index
    #[allow(dead_code)]
    fn from_record_meta(meta: &RecordMeta, seq_index: &mut SequenceIndex) -> Self {
        Self {
            rank: meta.rank,
            query_id: seq_index.get_or_insert(&meta.query_name),
            target_id: seq_index.get_or_insert(&meta.target_name),
            query_start: meta.query_start,
            query_end: meta.query_end,
            target_start: meta.target_start,
            target_end: meta.target_end,
            block_length: meta.block_length,
            identity: meta.identity,
            strand: meta.strand,
            chain_id: meta.chain_id.clone(),
            chain_status: meta.chain_status.clone(),
            discard: meta.discard,
            overlapped: meta.overlapped,
        }
    }

    /// Convert back to RecordMeta for output
    #[allow(dead_code)]
    fn to_record_meta(&self, seq_index: &SequenceIndex) -> RecordMeta {
        RecordMeta {
            rank: self.rank,
            query_name: seq_index.name(self.query_id).to_string(),
            target_name: seq_index.name(self.target_id).to_string(),
            query_start: self.query_start,
            query_end: self.query_end,
            target_start: self.target_start,
            target_end: self.target_end,
            block_length: self.block_length,
            identity: self.identity,
            matches: (self.identity * self.block_length as f64) as u64, // Estimate matches from identity
            alignment_length: self.block_length, // Use block_length as alignment length
            strand: self.strand,
            chain_id: self.chain_id.clone(),
            chain_status: self.chain_status.clone(),
            discard: self.discard,
            overlapped: self.overlapped,
        }
    }
}

/// Represents a merged chain for scaffold filtering
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct MergedChain {
    query_name: String,
    target_name: String,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
    strand: char,
    total_length: u64,
    weighted_identity: f64, // Average identity of mapped regions (matches/mapped_length)
    sum_matches: u64,       // Sum of matches from all mappings
    sum_block_lengths: u64, // Sum of actual mapped lengths
    member_indices: Vec<usize>, // Indices of original mappings in this chain
}

impl MergedChain {
    /// Calculate score for this chain based on the scoring function
    fn score(&self, scoring: ScoringFunction) -> f64 {
        match scoring {
            ScoringFunction::Identity => {
                // Use average identity of mapped regions
                self.weighted_identity
            }
            ScoringFunction::Length => {
                // Use sum of mapped lengths (not span, to be gap-neutral)
                self.sum_block_lengths as f64
            }
            ScoringFunction::LengthIdentity => {
                // Mapped length * weighted identity
                self.sum_block_lengths as f64 * self.weighted_identity
            }
            ScoringFunction::LogLengthIdentity => {
                // log(mapped length) * weighted identity
                if self.sum_block_lengths > 0 {
                    (self.sum_block_lengths as f64).ln() * self.weighted_identity
                } else {
                    0.0
                }
            }
            ScoringFunction::Matches => {
                // Total matches only - more matches = better
                self.sum_matches as f64
            }
        }
    }
}

/// Compact merged chain using sequence IDs
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct CompactMergedChain {
    query_id: u32,
    target_id: u32,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
    strand: char,
    total_length: u64,
    member_indices: Vec<usize>,
}

/// PAF filter that preserves original records
pub struct PafFilter {
    config: FilterConfig,
    #[allow(dead_code)]
    temp_dir: Option<String>,
    keep_self: bool,
    scaffolds_only: bool,
}

#[allow(dead_code)]
impl PafFilter {
    pub fn new(config: FilterConfig) -> Self {
        PafFilter {
            config,
            temp_dir: std::env::var("TMPDIR")
                .ok()
                .or_else(|| Some("/tmp".to_string())),
            keep_self: false, // Exclude self-mappings by default
            scaffolds_only: false,
        }
    }

    /// Extract prefix from sequence name for grouping
    fn extract_prefix(&self, name: &str) -> String {
        if self.config.skip_prefix {
            // Skip prefix extraction, use full name
            return name.to_string();
        }

        // Extract prefix: take everything up to and including the last delimiter
        // e.g., "r#1#2" -> "r#1#", "Rabacal-1#Chr1" -> "Rabacal-1#", "r#1" -> "r#"
        if let Some(pos) = name.rfind(self.config.prefix_delimiter) {
            name[..=pos].to_string() // Include the delimiter itself
        } else {
            // No delimiter found, use full name
            name.to_string()
        }
    }

    pub fn with_keep_self(mut self, keep_self: bool) -> Self {
        self.keep_self = keep_self;
        self
    }

    pub fn with_scaffolds_only(mut self, scaffolds_only: bool) -> Self {
        self.scaffolds_only = scaffolds_only;
        self
    }

    /// Main filtering pipeline using record ranks
    pub fn filter_paf<P: AsRef<Path>>(&self, input_path: P, output_path: P) -> Result<()> {
        // First pass: extract metadata for all records
        let metadata = self.extract_metadata(&input_path)?;

        // Apply filters to get passing record ranks
        let passing_ranks = self.apply_filters(metadata)?;

        // Second pass: write passing records with annotations
        self.write_filtered_output(&input_path, &output_path, passing_ranks)?;

        Ok(())
    }

    /// Extract metadata from PAF without modifying records (private implementation)
    fn extract_metadata<P: AsRef<Path>>(&self, path: P) -> Result<Vec<RecordMeta>> {
        let reader = open_paf_input(path.as_ref())?;
        let mut metadata = Vec::new();
        let mut has_cigar = false;
        let mut checked_cigar = false;

        for (rank, line) in reader.lines().enumerate() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() < 11 {
                continue;
            }

            // Parse essential fields
            let query_name = fields[0].to_string();
            let _query_len = fields[1].parse::<u64>().unwrap_or(0);
            let query_start = fields[2].parse::<u64>().unwrap_or(0);
            let query_end = fields[3].parse::<u64>().unwrap_or(0);
            let strand = if fields[4] == "+" { '+' } else { '-' };
            let target_name = fields[5].to_string();
            let _target_len = fields[6].parse::<u64>().unwrap_or(0);
            let target_start = fields[7].parse::<u64>().unwrap_or(0);
            let target_end = fields[8].parse::<u64>().unwrap_or(0);
            let matches = fields[9].parse::<u64>().unwrap_or(0);
            let block_length = fields[10].parse::<u64>().unwrap_or(1);

            // Calculate block identity: matches / alignment_length
            // alignment_length is the block_length (denominator in PAF format)
            let alignment_length = block_length;
            let mut identity = matches as f64 / alignment_length.max(1) as f64;
            let mut exact_matches = matches; // Default to PAF matches field

            // Look for tags: divergence (dv:f:) and CIGAR (cg:Z:)
            for field in &fields[11..] {
                if let Some(div_str) = field.strip_prefix("dv:f:") {
                    if let Ok(div) = div_str.parse::<f64>() {
                        identity = 1.0 - div;
                    }
                } else if let Some(cigar_str) = field.strip_prefix("cg:Z:") {
                    // Parse CIGAR to get exact match count
                    if let Ok((cigar_matches, _, _, _)) = crate::paf::parse_cigar_counts(cigar_str)
                    {
                        if cigar_matches > 0 {
                            exact_matches = cigar_matches;
                            // Also update identity based on CIGAR matches
                            identity = cigar_matches as f64 / alignment_length.max(1) as f64;
                            has_cigar = true;
                        }
                    }
                }
            }

            // Check CIGAR availability on first record when using matches scoring
            if !checked_cigar && self.config.scoring_function == ScoringFunction::Matches {
                checked_cigar = true;
                if !has_cigar {
                    eprintln!("[sweepga] WARNING: Using 'matches' scoring but input lacks CIGAR strings (cg:Z: tags).");
                    eprintln!("[sweepga] WARNING: Match counts will be estimated from PAF matches field (column 10).");
                    eprintln!("[sweepga] WARNING: For exact match counts, use aligners that output extended CIGAR (e.g., wfmash).");
                }
            }

            metadata.push(RecordMeta {
                rank,
                query_name,
                target_name,
                query_start,
                query_end,
                target_start,
                target_end,
                block_length,
                identity,
                matches: exact_matches, // Use exact matches from CIGAR if available
                alignment_length,
                strand,
                chain_id: None,
                chain_status: ChainStatus::Unassigned,
                discard: false,
                overlapped: false,
            });
        }

        Ok(metadata)
    }

    /// Apply filtering pipeline following wfmash's algorithm
    pub fn apply_filters(
        &self,
        mut metadata: Vec<RecordMeta>,
    ) -> Result<HashMap<usize, RecordMeta>> {
        // 1. Filter by minimum block length, self-mappings, and minimum identity
        metadata.retain(|m| {
            m.block_length >= self.config.min_block_length
                && (self.keep_self || m.query_name != m.target_name)
                && m.identity >= self.config.min_identity
        });

        // Keep all original mappings for rescue phase (before any plane sweep)
        let all_original_mappings = metadata.clone();

        // 2. Apply plane sweep as the default filtering method
        // IMPORTANT: Plane sweep must be applied PER QUERY SEQUENCE, not per query-target pair!
        // This is how wfmash does it - all mappings from a query compete regardless of target

        // The plane sweep is already correctly implemented in apply_plane_sweep_to_mappings
        // which groups by query sequence internally. We just need to pass ALL mappings to it.
        let before_plane_sweep = metadata.len();
        metadata = self.apply_plane_sweep_to_mappings(&metadata)?;
        let after_plane_sweep = metadata.len();

        // Report plane sweep if it filtered anything
        if before_plane_sweep != after_plane_sweep {
            eprintln!("[sweepga] Plane sweep: {before_plane_sweep} → {after_plane_sweep} mappings");
        }

        // If no scaffolding (scaffold_gap == 0), we're done - return the plane-swept mappings
        if self.config.scaffold_gap == 0 {
            // Calculate statistics for plane sweep output
            let total_kept_bases: u64 = metadata.iter().map(|m| m.block_length).sum();
            let total_kept_matches: u64 = metadata.iter().map(|m| m.matches).sum();
            let avg_kept_identity = if total_kept_bases > 0 {
                total_kept_matches as f64 / total_kept_bases as f64
            } else {
                0.0
            };

            eprintln!("[sweepga] Plane sweep filtering (no scaffolding)");
            eprintln!(
                "[sweepga] Summary: {} → {} mappings ({:.1}% kept)",
                all_original_mappings.len(),
                metadata.len(),
                (metadata.len() as f64 / all_original_mappings.len().max(1) as f64) * 100.0
            );
            eprintln!(
                "[sweepga]   Output: {:.1} Mb total, {:.1}% avg identity",
                total_kept_bases as f64 / 1_000_000.0,
                avg_kept_identity * 100.0
            );

            let mut result = HashMap::new();
            for m in metadata {
                result.insert(m.rank, m);
            }
            return Ok(result);
        }

        // 3. Apply scaffold filtering (wfmash's filterByScaffolds)
        // According to CLAUDE.md: create scaffolds from the PLANE-SWEPT mappings,
        // then rescue from ALL ORIGINAL mappings

        // Calculate input statistics
        let total_mapped_bases: u64 = metadata.iter().map(|m| m.block_length).sum();
        let avg_identity = if !metadata.is_empty() {
            metadata.iter().map(|m| m.identity).sum::<f64>() / metadata.len() as f64
        } else {
            0.0
        };

        eprintln!("[sweepga] Scaffold creation");
        eprintln!(
            "[sweepga]   Input: {} mappings, {:.1} Mb total, {:.1}% avg identity",
            metadata.len(),
            total_mapped_bases as f64 / 1_000_000.0,
            avg_identity * 100.0
        );

        // Use scaffold_gap for merging into scaffolds
        let merged_chains = self.merge_mappings_into_chains(&metadata, self.config.scaffold_gap)?;
        eprintln!(
            "[sweepga]   Merged into {} chains (gap ≤ {})",
            merged_chains.len(),
            self.config.scaffold_gap
        );

        // Step 2: Filter chains by minimum scaffold length and identity
        let before_filter = merged_chains.len();
        let mut filtered_chains: Vec<MergedChain> = merged_chains
            .into_iter()
            .filter(|chain| {
                chain.total_length >= self.config.min_scaffold_length
                    && chain.weighted_identity >= self.config.min_scaffold_identity
            })
            .collect();
        let filtered_out = before_filter - filtered_chains.len();
        eprintln!(
            "[sweepga]   Length/identity filter: {} chains kept, {} removed",
            filtered_chains.len(),
            filtered_out
        );
        if self.config.min_scaffold_length > 0 || self.config.min_scaffold_identity > 0.0 {
            eprintln!(
                "[sweepga]   Min length: {}, min identity: {:.1}%",
                self.config.min_scaffold_length,
                self.config.min_scaffold_identity * 100.0
            );
        }

        // Step 3: Apply plane sweep to scaffolds
        eprintln!("[sweepga] Scaffold sweep");
        let before_sweep = filtered_chains.len();

        // Track which mappings are in scaffolds BEFORE plane sweep
        let mut pre_sweep_scaffold_members: HashSet<usize> = HashSet::new();
        for chain in &filtered_chains {
            for &member_rank in &chain.member_indices {
                pre_sweep_scaffold_members.insert(member_rank);
            }
        }

        filtered_chains = self.apply_scaffold_plane_sweep(filtered_chains)?;
        eprintln!(
            "[sweepga]   Scaffold sweep: {} → {} scaffolds",
            before_sweep,
            filtered_chains.len()
        );

        // If scaffolds_only mode, return the actual mappings that form scaffolds
        if self.scaffolds_only {
            let mut scaffold_mappings = HashMap::new();

            // Build a map from rank to original metadata for quick lookup
            // Use all_original_mappings since member_indices refers to original ranks
            let mut rank_to_meta: HashMap<usize, &RecordMeta> = HashMap::new();
            for meta in &all_original_mappings {
                rank_to_meta.insert(meta.rank, meta);
            }

            // Collect all member mappings from scaffold chains
            for (chain_idx, chain) in filtered_chains.iter().enumerate() {
                // Create a unique chain ID
                let chain_id = format!("chain_{}", chain_idx + 1);

                for &member_rank in &chain.member_indices {
                    // member_indices contains ranks of original mappings
                    if let Some(meta) = rank_to_meta.get(&member_rank) {
                        let mut scaffold_meta = (*meta).clone();
                        scaffold_meta.chain_status = ChainStatus::Scaffold;
                        scaffold_meta.chain_id = Some(chain_id.clone());
                        scaffold_mappings.insert(member_rank, scaffold_meta);
                    }
                }
            }

            return Ok(scaffold_mappings);
        }

        // Step 4: Identify anchors - use the actual member mappings of scaffold chains
        // NOTE: member_indices actually contains ranks, not array indices!
        let mut anchor_ranks = HashSet::new();
        let mut rank_to_chain_id: HashMap<usize, String> = HashMap::new();

        for (chain_idx, chain) in filtered_chains.iter().enumerate() {
            let chain_id = format!("chain_{}", chain_idx + 1);

            // The member_indices field contains the ranks of mappings in this chain
            for &member_rank in &chain.member_indices {
                anchor_ranks.insert(member_rank);
                rank_to_chain_id.insert(member_rank, chain_id.clone());
            }
        }

        // Identify mappings that were in filtered-out scaffolds
        // These should NOT be rescued even if they're near surviving anchors
        let filtered_scaffold_members: HashSet<usize> = pre_sweep_scaffold_members
            .difference(&anchor_ranks)
            .copied()
            .collect();

        eprintln!("[sweepga] Rescue phase");
        eprintln!(
            "[sweepga]   Anchors: {} mappings in {} scaffolds",
            anchor_ranks.len(),
            filtered_chains.len()
        );

        // Map ranks to indices in all_original_mappings
        let mut rank_to_idx = HashMap::new();
        for (idx, meta) in all_original_mappings.iter().enumerate() {
            rank_to_idx.insert(meta.rank, idx);
        }

        // Step 5: Rescue mappings within scaffold_max_deviation of anchors
        // OPTIMIZED: Group by chromosome pair and sort for efficient rescue

        // First, group all mappings by (query_chr, target_chr) pair
        let mut mappings_by_chr_pair: HashMap<(String, String), Vec<usize>> = HashMap::new();
        for (idx, mapping) in all_original_mappings.iter().enumerate() {
            let key = (mapping.query_name.clone(), mapping.target_name.clone());
            mappings_by_chr_pair.entry(key).or_default().push(idx);
        }

        // Sort each chromosome pair's mappings by query position for binary search
        for indices in mappings_by_chr_pair.values_mut() {
            indices.sort_by_key(|&idx| all_original_mappings[idx].query_start);
        }

        // Collect anchors by chromosome pair for efficient lookup
        let mut anchors_by_chr_pair: HashMap<(String, String), Vec<usize>> = HashMap::new();
        for &anchor_rank in &anchor_ranks {
            if let Some(&anchor_idx) = rank_to_idx.get(&anchor_rank) {
                let anchor = &all_original_mappings[anchor_idx];
                let key = (anchor.query_name.clone(), anchor.target_name.clone());
                anchors_by_chr_pair.entry(key).or_default().push(anchor_idx);
            }
        }

        let mut kept_mappings = Vec::new();
        let mut kept_status = HashMap::new();
        let max_deviation = self.config.scaffold_max_deviation;

        // Process each chromosome pair independently (can be parallelized)
        for ((query_chr, target_chr), mapping_indices) in &mappings_by_chr_pair {
            let chr_key = (query_chr.clone(), target_chr.clone());
            let chr_anchors = anchors_by_chr_pair
                .get(&chr_key)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            if chr_anchors.is_empty() {
                continue; // No anchors on this chromosome pair, skip
            }

            // For each mapping on this chromosome pair
            for &mapping_idx in mapping_indices {
                let mapping = &all_original_mappings[mapping_idx];

                if anchor_ranks.contains(&mapping.rank) {
                    // This is an anchor (member of a scaffold chain) - keep it with its chain ID
                    let mut anchor_mapping = mapping.clone();
                    if let Some(chain_id) = rank_to_chain_id.get(&mapping.rank) {
                        anchor_mapping.chain_id = Some(chain_id.clone());
                    }
                    kept_mappings.push(anchor_mapping);
                    // All anchors are scaffold members by definition
                    kept_status.insert(mapping.rank, ChainStatus::Scaffold);
                } else if filtered_scaffold_members.contains(&mapping.rank) {
                    // This mapping was part of a plane-sweep-filtered scaffold
                    // Do NOT rescue it, even if it's near an anchor
                    continue;
                } else if max_deviation > 0 {
                    // Only attempt rescue if max_deviation > 0
                    // Check if within deviation distance of any anchor
                    // Use binary search to find anchors within query range
                    let mapping_q_center = (mapping.query_start + mapping.query_end) / 2;
                    let mapping_t_center = (mapping.target_start + mapping.target_end) / 2;

                    let mut min_distance = u64::MAX;
                    let mut closest_anchor_rank = None;

                    // Only check anchors that are within max_deviation in query space
                    for &anchor_idx in chr_anchors {
                        let anchor = &all_original_mappings[anchor_idx];
                        let anchor_q_center = (anchor.query_start + anchor.query_end) / 2;

                        // Early exit if anchor is too far in query space
                        let q_diff =
                            (mapping_q_center as i64 - anchor_q_center as i64).unsigned_abs();
                        if q_diff > max_deviation {
                            continue; // Too far in query dimension alone
                        }

                        let anchor_t_center = (anchor.target_start + anchor.target_end) / 2;
                        let t_diff =
                            (mapping_t_center as i64 - anchor_t_center as i64).unsigned_abs();

                        // Euclidean distance
                        let distance = ((q_diff * q_diff + t_diff * t_diff) as f64).sqrt() as u64;
                        if distance < min_distance {
                            min_distance = distance;
                            closest_anchor_rank = Some(anchor.rank);
                        }

                        // Early exit if we found a close enough anchor
                        if min_distance <= max_deviation {
                            break;
                        }
                    }

                    if min_distance <= max_deviation {
                        // This mapping is rescued - assign it to the same chain as its closest anchor
                        let mut rescued_mapping = mapping.clone();
                        if let Some(anchor_rank) = closest_anchor_rank {
                            if let Some(chain_id) = rank_to_chain_id.get(&anchor_rank) {
                                rescued_mapping.chain_id = Some(chain_id.clone());
                            }
                        }
                        kept_mappings.push(rescued_mapping);
                        kept_status.insert(mapping.rank, ChainStatus::Rescued);
                    }
                }
                // If max_deviation == 0 and not an anchor, mapping is not kept
            }
        }

        // Count anchors vs rescued
        let anchor_count = anchor_ranks.len();
        let rescued_count = kept_mappings.len() - anchor_count;
        let total_kept_bases: u64 = kept_mappings.iter().map(|m| m.block_length).sum();
        let avg_kept_identity = if !kept_mappings.is_empty() {
            kept_mappings.iter().map(|m| m.identity).sum::<f64>() / kept_mappings.len() as f64
        } else {
            0.0
        };

        eprintln!(
            "[sweepga]   Rescued: {} mappings within {}bp of anchors",
            rescued_count, self.config.scaffold_max_deviation
        );
        eprintln!(
            "[sweepga] Summary: {} → {} mappings ({:.1}% kept)",
            all_original_mappings.len(),
            kept_mappings.len(),
            (kept_mappings.len() as f64 / all_original_mappings.len().max(1) as f64) * 100.0
        );
        eprintln!(
            "[sweepga]   Output: {:.1} Mb total, {:.1}% avg identity",
            total_kept_bases as f64 / 1_000_000.0,
            avg_kept_identity * 100.0
        );

        // Build result map directly from kept mappings
        let mut passing = HashMap::new();
        for meta in kept_mappings {
            let mut result = meta.clone();
            // Set the chain status based on what we determined earlier
            result.chain_status = kept_status
                .get(&meta.rank)
                .cloned()
                .unwrap_or(ChainStatus::Scaffold);
            passing.insert(meta.rank, result);
        }

        Ok(passing)
    }

    /// Merge mappings into chains using wfmash's union-find approach
    fn merge_mappings_into_chains(
        &self,
        metadata: &[RecordMeta],
        max_gap: u64,
    ) -> Result<Vec<MergedChain>> {
        use crate::union_find::UnionFind;

        // Group by (query, target, strand) - this is like wfmash's refSeqId grouping
        // Store the original rank, not the position in metadata array
        let mut groups: HashMap<(String, String, char), Vec<(usize, usize)>> = HashMap::new();

        for (idx, meta) in metadata.iter().enumerate() {
            let key = (
                meta.query_name.clone(),
                meta.target_name.clone(),
                meta.strand,
            );
            groups.entry(key).or_default().push((meta.rank, idx));
        }

        let mut all_chains = Vec::new();

        for ((query, target, strand), indices) in groups {
            // Sort by query start position (like wfmash's sort)
            let mut sorted_indices = indices.clone();
            sorted_indices.sort_by_key(|&(_rank, idx)| metadata[idx].query_start);

            // Best-buddy chaining: track best predecessor for each mapping
            let mut best_pred_score: Vec<u64> = vec![u64::MAX; sorted_indices.len()];
            let mut best_pred_idx: Vec<Option<usize>> = vec![None; sorted_indices.len()];

            // Phase 1: Find all best-buddy relationships
            for i in 0..sorted_indices.len() {
                let (_rank_i, idx_i) = sorted_indices[i];
                let search_bound = metadata[idx_i].query_end + max_gap;

                let mut best_j = None;
                let mut best_score = u64::MAX;

                for j in (i + 1)..sorted_indices.len() {
                    let (_rank_j, idx_j) = sorted_indices[j];

                    if metadata[idx_j].query_start > search_bound {
                        break;
                    }

                    // Calculate query distance/overlap (using actual overlap as distance, not 0)
                    let q_gap = if metadata[idx_j].query_start >= metadata[idx_i].query_end {
                        // Normal gap
                        metadata[idx_j].query_start - metadata[idx_i].query_end
                    } else {
                        // Overlap: use overlap amount as distance, penalize if too large
                        let overlap = metadata[idx_i].query_end - metadata[idx_j].query_start;
                        if overlap <= max_gap / 5 {
                            overlap // Small overlap: use actual overlap distance
                        } else {
                            max_gap + 1 // Large overlap: reject
                        }
                    };

                    // Calculate target distance/overlap
                    let r_gap = if strand == '+' {
                        if metadata[idx_j].target_start >= metadata[idx_i].target_end {
                            metadata[idx_j].target_start - metadata[idx_i].target_end
                        } else {
                            let overlap = metadata[idx_i].target_end - metadata[idx_j].target_start;
                            if overlap <= max_gap / 5 {
                                overlap
                            } else {
                                max_gap + 1
                            }
                        }
                    } else if metadata[idx_i].target_start >= metadata[idx_j].target_end {
                        metadata[idx_i].target_start - metadata[idx_j].target_end
                    } else {
                        let overlap = metadata[idx_j].target_end - metadata[idx_i].target_start;
                        if overlap <= max_gap / 5 {
                            overlap
                        } else {
                            max_gap + 1
                        }
                    };

                    if q_gap <= max_gap && r_gap <= max_gap {
                        let dist_sq = q_gap * q_gap + r_gap * r_gap;

                        // Best-buddy: only link if i is the best predecessor for j
                        if dist_sq < best_score && dist_sq < best_pred_score[j] {
                            best_score = dist_sq;
                            best_j = Some(j);
                        }
                    }
                }

                // Record best-buddy relationship for j
                if let Some(j) = best_j {
                    best_pred_score[j] = best_score;
                    best_pred_idx[j] = Some(i);
                }
            }

            // Phase 2: Create chains using union-find on best-buddy pairs only
            let mut uf = UnionFind::new(sorted_indices.len());
            for (j, pred) in best_pred_idx.iter().enumerate() {
                if let Some(i) = pred {
                    uf.union(*i, j);
                }
            }

            // Extract chains from union-find sets
            let sets = uf.get_sets();
            let chains: Vec<Vec<(usize, usize)>> = sets
                .into_iter()
                .map(|set_indices| set_indices.into_iter().map(|i| sorted_indices[i]).collect())
                .collect();

            // Create merged chains from groups of mappings
            for chain_indices in chains {
                if chain_indices.is_empty() {
                    continue;
                }

                // Calculate merged chain boundaries and weighted identity
                let mut q_min = u64::MAX;
                let mut q_max = 0;
                let mut t_min = u64::MAX;
                let mut t_max = 0;
                let mut member_ranks = Vec::new();
                let mut sum_matches = 0u64;
                let mut sum_block_lengths = 0u64;

                for &(rank, idx) in &chain_indices {
                    let meta = &metadata[idx];
                    q_min = q_min.min(meta.query_start);
                    q_max = q_max.max(meta.query_end);
                    t_min = t_min.min(meta.target_start);
                    t_max = t_max.max(meta.target_end);
                    member_ranks.push(rank); // Store the original rank

                    // Sum up matches and block lengths for weighted identity
                    sum_matches += meta.matches;
                    sum_block_lengths += meta.block_length;
                }

                let total_length = q_max - q_min;

                // Calculate identity with log-compressed gaps
                // Gaps between alignments are penalized logarithmically
                // This gives a middle ground: gaps hurt but not linearly
                let gap_length = total_length.saturating_sub(sum_block_lengths);
                let log_compressed_gap = if gap_length > 0 {
                    (gap_length as f64).ln().max(0.0)
                } else {
                    0.0
                };
                let effective_length = sum_block_lengths as f64 + log_compressed_gap;

                let weighted_identity = if effective_length > 0.0 {
                    sum_matches as f64 / effective_length
                } else {
                    0.0
                };

                all_chains.push(MergedChain {
                    query_name: query.clone(),
                    target_name: target.clone(),
                    query_start: q_min,
                    query_end: q_max,
                    target_start: t_min,
                    target_end: t_max,
                    strand,
                    total_length,
                    weighted_identity,
                    sum_matches,
                    sum_block_lengths,
                    member_indices: member_ranks, // Now storing ranks, not indices
                });
            }
        }

        Ok(all_chains)
    }

    /// Apply plane sweep to filter overlapping chains
    /// Calculate overlap fraction between two ranges
    fn calculate_overlap(start1: u32, end1: u32, start2: u32, end2: u32) -> f64 {
        let overlap_start = start1.max(start2);
        let overlap_end = end1.min(end2);

        if overlap_start >= overlap_end {
            return 0.0; // No overlap
        }

        let overlap_len = overlap_end - overlap_start;
        let len1 = end1 - start1;
        let len2 = end2 - start2;
        let min_len = len1.min(len2);

        if min_len == 0 {
            return 0.0;
        }

        overlap_len as f64 / min_len as f64
    }

    /// Apply 1:1 filtering within a prefix group
    fn filter_one_to_one_within_group(&self, mappings: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        if mappings.is_empty() {
            return Ok(Vec::new());
        }

        // Find best mapping (longest)
        let best = mappings.into_iter().max_by_key(|m| m.block_length).unwrap();

        Ok(vec![best])
    }

    /// Apply plane sweep to raw mappings (before scaffold filtering)
    /// This implements the wfmash plane sweep algorithm exactly
    /// CRITICAL: Plane sweep must run PER QUERY SEQUENCE, not globally
    fn apply_plane_sweep_to_mappings(&self, mappings: &[RecordMeta]) -> Result<Vec<RecordMeta>> {
        if mappings.is_empty() || mappings.len() <= 1 {
            return Ok(mappings.to_vec());
        }

        // Convert RecordMeta to PlaneSweepMapping with grouping keys
        let plane_sweep_mappings: Vec<(PlaneSweepMapping, String, String)> = mappings
            .iter()
            .enumerate()
            .map(|(idx, meta)| {
                let mapping = PlaneSweepMapping {
                    idx, // Store original index
                    query_start: meta.query_start,
                    query_end: meta.query_end,
                    target_start: meta.target_start,
                    target_end: meta.target_end,
                    identity: meta.identity, // Use actual identity from PAF
                    flags: 0,
                };
                // For 1:1 filtering with prefix grouping:
                // We want to keep the best alignment for each chromosome WITHIN each genome pair
                // So we group by full chromosome names (which include the genome prefix)
                let query_group = meta.query_name.clone();
                let target_group = meta.target_name.clone();
                (mapping, query_group, target_group)
            })
            .collect();

        let overlap_threshold = self.config.overlap_threshold;

        use crate::plane_sweep_exact::{plane_sweep_query, plane_sweep_target};

        let query_limit = match self.config.mapping_filter_mode {
            FilterMode::OneToOne => 1,
            FilterMode::OneToMany => self.config.mapping_max_per_query.unwrap_or(1),
            FilterMode::ManyToMany => self.config.mapping_max_per_query.unwrap_or(usize::MAX),
        };

        let target_limit = match self.config.mapping_filter_mode {
            FilterMode::OneToOne => 1,
            FilterMode::OneToMany => self.config.mapping_max_per_target.unwrap_or(usize::MAX),
            FilterMode::ManyToMany => self.config.mapping_max_per_target.unwrap_or(usize::MAX),
        };

        // Helper function to extract genome prefix from sequence name
        // For "SGDref#1#chrI", returns "SGDref#1#"
        let extract_genome_prefix = |seq_name: &str| -> String {
            if let Some(last_pos) = seq_name.rfind('#') {
                seq_name[..=last_pos].to_string()
            } else {
                seq_name.to_string()
            }
        };

        // CRITICAL: Group by (query_genome_prefix, target_genome_prefix) pairs FIRST
        // This ensures plane sweep runs independently for each genome pair
        let mut genome_pair_groups: HashMap<(String, String), Vec<usize>> = HashMap::new();

        for (i, (_, q, t)) in plane_sweep_mappings.iter().enumerate() {
            let query_genome = extract_genome_prefix(q);
            let target_genome = extract_genome_prefix(t);
            genome_pair_groups
                .entry((query_genome, target_genome))
                .or_default()
                .push(i);
        }

        // Process each genome pair independently
        let mut all_kept_indices = Vec::new();

        for (_genome_pair, genome_pair_indices) in genome_pair_groups {
            // Within this genome pair, apply query and target sweeps with intersection

            // Query axis sweep: group by query chr within this genome pair
            let mut query_kept_set = HashSet::new();
            let mut by_query: HashMap<String, Vec<usize>> = HashMap::new();

            for &idx in &genome_pair_indices {
                let (_, q, _t) = &plane_sweep_mappings[idx];
                by_query.entry(q.clone()).or_default().push(idx);
            }

            for (_q_chr, indices) in by_query {
                let mut query_mappings: Vec<_> =
                    indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                let kept = plane_sweep_query(
                    &mut query_mappings,
                    query_limit,
                    overlap_threshold,
                    self.config.scoring_function,
                );
                for k in kept {
                    query_kept_set.insert(indices[k]);
                }
            }

            // Target axis sweep: group by target chr within this genome pair
            let mut target_kept_set = HashSet::new();
            let mut by_target: HashMap<String, Vec<usize>> = HashMap::new();

            for &idx in &genome_pair_indices {
                let (_, _q, t) = &plane_sweep_mappings[idx];
                by_target.entry(t.clone()).or_default().push(idx);
            }

            for (_t_chr, indices) in by_target {
                let mut target_mappings: Vec<_> =
                    indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                let kept = plane_sweep_target(
                    &mut target_mappings,
                    target_limit,
                    overlap_threshold,
                    self.config.scoring_function,
                );
                for k in kept {
                    target_kept_set.insert(indices[k]);
                }
            }

            // Intersection: keep only indices that passed both sweeps within this genome pair
            for idx in query_kept_set.intersection(&target_kept_set) {
                all_kept_indices.push(*idx);
            }
        }

        let kept_indices = all_kept_indices;

        // Convert back to RecordMeta
        let result: Vec<RecordMeta> = kept_indices
            .iter()
            .map(|&idx| mappings[idx].clone())
            .collect();

        Ok(result)
    }

    /// Apply scaffold plane sweep - SAME ALGORITHM as regular mappings, just different params
    fn apply_scaffold_plane_sweep(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        if chains.is_empty() || chains.len() <= 1 {
            return Ok(chains);
        }

        // Convert MergedChain to PlaneSweepMapping - same as regular mappings
        let plane_sweep_mappings: Vec<(PlaneSweepMapping, String, String)> = chains
            .iter()
            .enumerate()
            .map(|(idx, chain)| {
                let mapping = PlaneSweepMapping {
                    idx,
                    query_start: chain.query_start,
                    query_end: chain.query_end,
                    target_start: chain.target_start,
                    target_end: chain.target_end,
                    identity: chain.weighted_identity, // Use weighted identity for scoring
                    flags: 0,
                };
                let query_group = chain.query_name.clone();
                let target_group = chain.target_name.clone();
                (mapping, query_group, target_group)
            })
            .collect();

        // Get the scaffold filter params (vs mapping filter params)
        let overlap_threshold = self.config.scaffold_overlap_threshold;

        // Apply plane sweep - EXACT SAME FUNCTION as regular mappings
        // Always use independent query and target sweeps
        use crate::plane_sweep_exact::{plane_sweep_query, plane_sweep_target};

        let kept_indices: Vec<usize> = match self.config.scaffold_filter_mode {
            FilterMode::OneToOne => {
                // True 1:1: Apply plane sweep independently on query and target axes
                // 1. Sweep across each query chromosome (scaffolds to different targets compete)
                // 2. Sweep across each target chromosome (scaffolds from different queries compete)
                // 3. Keep only scaffolds that pass BOTH sweeps

                // Query axis sweep: group by query chr
                let mut query_kept_set = HashSet::new();
                let mut by_query: std::collections::HashMap<String, Vec<usize>> =
                    std::collections::HashMap::new();

                for (i, (_, q, _t)) in plane_sweep_mappings.iter().enumerate() {
                    by_query.entry(q.clone()).or_default().push(i);
                }

                for (_q_chr, indices) in by_query {
                    let mut query_mappings: Vec<_> =
                        indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                    let kept = plane_sweep_query(
                        &mut query_mappings,
                        1,
                        overlap_threshold,
                        self.config.scoring_function,
                    );
                    for k in kept {
                        query_kept_set.insert(indices[k]);
                    }
                }

                // Target axis sweep: group by target chr
                let mut target_kept_set = HashSet::new();
                let mut by_target: std::collections::HashMap<String, Vec<usize>> =
                    std::collections::HashMap::new();

                for (i, (_, _q, t)) in plane_sweep_mappings.iter().enumerate() {
                    by_target.entry(t.clone()).or_default().push(i);
                }

                for (_t_chr, indices) in by_target {
                    let mut target_mappings: Vec<_> =
                        indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                    let kept = plane_sweep_target(
                        &mut target_mappings,
                        1,
                        overlap_threshold,
                        self.config.scoring_function,
                    );
                    for k in kept {
                        target_kept_set.insert(indices[k]);
                    }
                }

                // Intersection: keep only indices that passed both sweeps
                query_kept_set
                    .intersection(&target_kept_set)
                    .copied()
                    .collect()
            }
            FilterMode::OneToMany | FilterMode::ManyToMany => {
                // Same independent sweep pattern as 1:1, just different limits
                let query_limit = match self.config.scaffold_filter_mode {
                    FilterMode::OneToMany => self.config.scaffold_max_per_query.unwrap_or(1),
                    FilterMode::ManyToMany => {
                        self.config.scaffold_max_per_query.unwrap_or(usize::MAX)
                    }
                    _ => unreachable!(),
                };

                let target_limit = match self.config.scaffold_filter_mode {
                    FilterMode::OneToMany => {
                        self.config.scaffold_max_per_target.unwrap_or(usize::MAX)
                    }
                    FilterMode::ManyToMany => {
                        self.config.scaffold_max_per_target.unwrap_or(usize::MAX)
                    }
                    _ => unreachable!(),
                };

                // Query axis sweep: group by query chr
                let mut query_kept_set = HashSet::new();
                let mut by_query: std::collections::HashMap<String, Vec<usize>> =
                    std::collections::HashMap::new();

                for (i, (_, q, _t)) in plane_sweep_mappings.iter().enumerate() {
                    by_query.entry(q.clone()).or_default().push(i);
                }

                for (_q_chr, indices) in by_query {
                    let mut query_mappings: Vec<_> =
                        indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                    let kept = plane_sweep_query(
                        &mut query_mappings,
                        query_limit,
                        overlap_threshold,
                        self.config.scoring_function,
                    );
                    for k in kept {
                        query_kept_set.insert(indices[k]);
                    }
                }

                // Target axis sweep: group by target chr
                let mut target_kept_set = HashSet::new();
                let mut by_target: std::collections::HashMap<String, Vec<usize>> =
                    std::collections::HashMap::new();

                for (i, (_, _q, t)) in plane_sweep_mappings.iter().enumerate() {
                    by_target.entry(t.clone()).or_default().push(i);
                }

                for (_t_chr, indices) in by_target {
                    let mut target_mappings: Vec<_> =
                        indices.iter().map(|&i| plane_sweep_mappings[i].0).collect();

                    let kept = plane_sweep_target(
                        &mut target_mappings,
                        target_limit,
                        overlap_threshold,
                        self.config.scoring_function,
                    );
                    for k in kept {
                        target_kept_set.insert(indices[k]);
                    }
                }

                // Intersection: keep only indices that passed both sweeps
                query_kept_set
                    .intersection(&target_kept_set)
                    .copied()
                    .collect()
            }
        };

        // Convert back to MergedChain
        let result: Vec<MergedChain> = kept_indices
            .iter()
            .map(|&idx| chains[idx].clone())
            .collect();

        Ok(result)
    }

    /// DEPRECATED - Now using unified apply_scaffold_plane_sweep
    #[allow(dead_code)]
    fn scaffold_one_to_one_filter(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        use crate::plane_sweep_exact::{plane_sweep_both, PlaneSweepMapping};

        // Group by (query_chr, target_chr) pairs
        let mut by_chr_pair: HashMap<(String, String), Vec<MergedChain>> = HashMap::new();

        for chain in chains {
            let key = (chain.query_name.clone(), chain.target_name.clone());
            by_chr_pair.entry(key).or_default().push(chain);
        }

        let mut filtered_chains = Vec::new();

        // Apply plane sweep within each chromosome pair
        for (_pair, pair_chains) in by_chr_pair {
            if pair_chains.is_empty() {
                continue;
            }

            // Convert chains to mappings for plane sweep
            let mut mappings: Vec<PlaneSweepMapping> = pair_chains
                .iter()
                .enumerate()
                .map(|(idx, chain)| PlaneSweepMapping {
                    idx,
                    query_start: chain.query_start,
                    query_end: chain.query_end,
                    target_start: chain.target_start,
                    target_end: chain.target_end,
                    identity: chain.weighted_identity, // Use weighted identity for scoring
                    flags: 0,
                })
                .collect();

            // Apply plane sweep with 1:1 filtering on both query and target axes
            // This keeps all non-overlapping scaffolds plus best overlapping ones
            let kept_indices = plane_sweep_both(
                &mut mappings,
                1, // max per query position
                1, // max per target position
                self.config.scaffold_overlap_threshold,
                self.config.scoring_function,
            );

            // Collect the kept chains
            for idx in kept_indices {
                filtered_chains.push(pair_chains[idx].clone());
            }
        }

        Ok(filtered_chains)
    }

    /// DEPRECATED - Now using unified apply_scaffold_plane_sweep
    #[allow(dead_code)]
    fn scaffold_one_to_many_filter(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        use crate::plane_sweep_exact::{plane_sweep_query, PlaneSweepMapping};

        // Group by (query_chr, target_chr) pairs
        let mut by_chr_pair: HashMap<(String, String), Vec<MergedChain>> = HashMap::new();

        for chain in chains {
            let key = (chain.query_name.clone(), chain.target_name.clone());
            by_chr_pair.entry(key).or_default().push(chain);
        }

        let mut filtered_chains = Vec::new();

        // Apply plane sweep within each chromosome pair
        for (_pair, pair_chains) in by_chr_pair {
            if pair_chains.is_empty() {
                continue;
            }

            // Convert chains to mappings for plane sweep
            let mut mappings: Vec<PlaneSweepMapping> = pair_chains
                .iter()
                .enumerate()
                .map(|(idx, chain)| PlaneSweepMapping {
                    idx,
                    query_start: chain.query_start,
                    query_end: chain.query_end,
                    target_start: chain.target_start,
                    target_end: chain.target_end,
                    identity: chain.weighted_identity, // Use weighted identity for scoring
                    flags: 0,
                })
                .collect();

            // Apply plane sweep with 1:∞ filtering (limit on query axis only)
            let max_per_query = self.config.scaffold_max_per_query.unwrap_or(1);
            let kept_indices = plane_sweep_query(
                &mut mappings,
                max_per_query,
                self.config.scaffold_overlap_threshold,
                self.config.scoring_function,
            );

            // Collect the kept chains
            for idx in kept_indices {
                filtered_chains.push(pair_chains[idx].clone());
            }
        }

        Ok(filtered_chains)
    }

    /// Keep all non-overlapping chains (for scaffold filtering with -n inf behavior)
    fn plane_sweep_keep_all_non_overlapping(
        &self,
        chains: Vec<MergedChain>,
    ) -> Result<Vec<MergedChain>> {
        // Group chains by query
        let mut by_query: HashMap<String, Vec<MergedChain>> = HashMap::new();
        for chain in chains {
            by_query
                .entry(chain.query_name.clone())
                .or_default()
                .push(chain);
        }

        let mut filtered = Vec::new();

        for (_query, mut query_chains) in by_query {
            // Sort by chain score (descending) to prioritize better chains
            let scoring = self.config.scoring_function;
            query_chains.sort_by(|a, b| {
                let score_a = a.score(scoring);
                let score_b = b.score(scoring);
                score_b
                    .partial_cmp(&score_a)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            let mut kept: Vec<MergedChain> = Vec::new();

            for chain in query_chains {
                let mut has_significant_overlap = false;

                // Check overlap with already kept chains
                for existing in &kept {
                    let overlap_start = chain.query_start.max(existing.query_start);
                    let overlap_end = chain.query_end.min(existing.query_end);

                    if overlap_start < overlap_end {
                        let overlap_len = overlap_end - overlap_start;
                        let chain_len = chain.query_end - chain.query_start;
                        let existing_len = existing.query_end - existing.query_start;
                        let min_len = chain_len.min(existing_len);

                        // Use scaffold_overlap_threshold to determine significant overlap
                        let overlap_frac = overlap_len as f64 / min_len as f64;
                        if overlap_frac > self.config.scaffold_overlap_threshold {
                            has_significant_overlap = true;
                            break;
                        }
                    }
                }

                // Keep chain if it doesn't significantly overlap with any existing
                if !has_significant_overlap {
                    kept.push(chain);
                }
            }

            filtered.extend(kept);
        }

        Ok(filtered)
    }

    fn plane_sweep_filter_chains(&self, chains: Vec<MergedChain>) -> Result<Vec<MergedChain>> {
        // Group chains by query
        let mut by_query: HashMap<String, Vec<MergedChain>> = HashMap::new();
        for chain in chains {
            by_query
                .entry(chain.query_name.clone())
                .or_default()
                .push(chain);
        }

        let mut filtered = Vec::new();

        for (_query, mut query_chains) in by_query {
            // Sort by query start
            query_chains.sort_by_key(|c| c.query_start);

            // Plane sweep to remove overlaps
            let mut kept: Vec<MergedChain> = Vec::new();
            let scoring_func = self.config.scoring_function;

            for chain in query_chains {
                // Check overlap with already kept chains
                let mut has_significant_overlap = false;

                kept.retain(|existing| {
                    let overlap_start = chain.query_start.max(existing.query_start);
                    let overlap_end = chain.query_end.min(existing.query_end);

                    if overlap_start < overlap_end {
                        let overlap_len = overlap_end - overlap_start;
                        let chain_len = chain.query_end - chain.query_start;
                        let existing_len = existing.query_end - existing.query_start;
                        let min_len = chain_len.min(existing_len);

                        // If >50% overlap (wfmash uses scaffold_overlap_threshold)
                        if overlap_len as f64 / min_len as f64
                            > self.config.scaffold_overlap_threshold
                        {
                            // Keep the better scoring chain
                            if existing.score(scoring_func) > chain.score(scoring_func) {
                                has_significant_overlap = true;
                                true // Keep existing
                            } else {
                                false // Remove existing, will add current
                            }
                        } else {
                            true // No significant overlap, keep both
                        }
                    } else {
                        true // No overlap, keep
                    }
                });

                if !has_significant_overlap {
                    kept.push(chain);
                }
            }

            filtered.extend(kept);
        }

        Ok(filtered)
    }

    /// Apply final filtering mode (1:1, 1:N, N:N)
    fn apply_filtering_mode(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        match self.config.mapping_filter_mode {
            FilterMode::OneToOne => self.one_to_one_filter(metadata),
            FilterMode::OneToMany => self.one_to_many_filter(metadata),
            FilterMode::ManyToMany => self.many_to_many_filter(metadata),
        }
    }

    /// One-to-one filtering
    fn one_to_one_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        let mut best_per_query: HashMap<String, RecordMeta> = HashMap::new();

        // Find best mapping per query
        for meta in metadata {
            best_per_query
                .entry(meta.query_name.clone())
                .and_modify(|e| {
                    if meta.block_length > e.block_length {
                        *e = meta.clone();
                    }
                })
                .or_insert(meta);
        }

        // From those, keep best per target
        let mut best_per_target: HashMap<String, RecordMeta> = HashMap::new();
        for meta in best_per_query.into_values() {
            best_per_target
                .entry(meta.target_name.clone())
                .and_modify(|e| {
                    if meta.block_length > e.block_length {
                        *e = meta.clone();
                    }
                })
                .or_insert(meta);
        }

        Ok(best_per_target.into_values().collect())
    }

    /// One-to-many filtering
    fn one_to_many_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        let mut result = Vec::new();

        // Group by query
        let mut by_query: HashMap<String, Vec<RecordMeta>> = HashMap::new();
        for meta in metadata {
            by_query
                .entry(meta.query_name.clone())
                .or_default()
                .push(meta);
        }

        // Keep best per query, optionally limit per target
        for (_query, mut metas) in by_query {
            metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));

            if let Some(limit) = self.config.mapping_max_per_target {
                metas.truncate(limit);
            }

            result.extend(metas);
        }

        Ok(result)
    }

    /// Many-to-many filtering
    fn many_to_many_filter(&self, metadata: Vec<RecordMeta>) -> Result<Vec<RecordMeta>> {
        // For scaffold filtering, we typically don't need additional limiting
        // but we can apply per-query/target limits if specified

        if self.config.mapping_max_per_query.is_none()
            && self.config.mapping_max_per_target.is_none()
        {
            return Ok(metadata);
        }

        let mut result = metadata;

        // Apply per-query limit
        if let Some(limit) = self.config.mapping_max_per_query {
            let mut by_query: HashMap<String, Vec<RecordMeta>> = HashMap::new();
            for meta in result {
                by_query
                    .entry(meta.query_name.clone())
                    .or_default()
                    .push(meta);
            }

            result = Vec::new();
            for (_query, mut metas) in by_query {
                metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));
                metas.truncate(limit);
                result.extend(metas);
            }
        }

        // Apply per-target limit
        if let Some(limit) = self.config.mapping_max_per_target {
            let mut by_target: HashMap<String, Vec<RecordMeta>> = HashMap::new();
            for meta in result {
                by_target
                    .entry(meta.target_name.clone())
                    .or_default()
                    .push(meta);
            }

            result = Vec::new();
            for (_target, mut metas) in by_target {
                metas.sort_by_key(|m| std::cmp::Reverse(m.block_length));
                metas.truncate(limit);
                result.extend(metas);
            }
        }

        Ok(result)
    }

    /// Write filtered output with annotations
    fn write_filtered_output<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: P,
        passing: HashMap<usize, RecordMeta>,
    ) -> Result<()> {
        let output_file = File::create(output_path)?;
        let mut writer = BufWriter::new(output_file);

        // Scaffolds_only mode is now handled the same as normal mode
        // since we're returning actual mappings, not synthetic records

        // Normal mode - read input and filter
        let reader = open_paf_input(input_path)?;

        for (rank, line) in reader.lines().enumerate() {
            if let Some(meta) = passing.get(&rank) {
                let mut line = line?;

                // Add our annotations as tags
                if let Some(ref chain_id) = meta.chain_id {
                    line.push_str(&format!("\tch:Z:{chain_id}"));
                }
                // Output status tag like wfmash does
                let status_str = match meta.chain_status {
                    ChainStatus::Scaffold => "scaffold",
                    ChainStatus::Rescued => "rescued",
                    ChainStatus::Unassigned => "unassigned",
                };
                line.push_str(&format!("\tst:Z:{status_str}"));

                writeln!(writer, "{line}")?;
            }
        }

        writer.flush()?;
        Ok(())
    }
}

/// Public function to extract PAF metadata without filtering (for testing/debugging)
#[allow(dead_code)]
pub fn extract_metadata<P: AsRef<Path>>(path: P) -> Result<(Vec<RecordMeta>, ())> {
    let config = FilterConfig {
        chain_gap: 0,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::ManyToMany,
        scaffold_max_per_query: None,
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: 0,
        min_scaffold_length: 0,
        scaffold_overlap_threshold: 0.95,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    let metadata = filter.extract_metadata(path)?;
    Ok((metadata, ()))
}
