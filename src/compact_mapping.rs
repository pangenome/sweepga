#![allow(dead_code)]
/// Compact mapping structures using sequence IDs instead of strings
use crate::mapping::ChainStatus;

/// Compact record metadata using sequence IDs instead of strings
#[derive(Debug, Clone)]
pub struct CompactRecordMeta {
    pub rank: usize, // 0-based index in original file
    pub query_id: u32,
    pub target_id: u32,
    pub query_start: u64,
    pub query_end: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub block_length: u64,
    pub identity: f64,
    pub strand: char,
    pub chain_id: Option<String>, // Still a string for now
    pub chain_status: ChainStatus,
    pub discard: bool,
    pub overlapped: bool,
}

/// Compact merged chain using sequence IDs
#[derive(Debug, Clone)]
pub struct CompactMergedChain {
    pub query_id: u32,
    pub target_id: u32,
    pub query_start: u64,
    pub query_end: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub strand: char,
    pub total_length: u64,
    pub member_indices: Vec<usize>,
}
