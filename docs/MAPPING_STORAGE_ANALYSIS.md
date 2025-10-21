# Mapping Storage Analysis & Optimization Proposal

## Current Storage Situation

### Data Structures

1. **RecordMeta** (main storage in paf_filter.rs):
   - Uses u32 (32-bit) for positions - **LIMITS TO 4GB SEQUENCES**
   - Stores full String for query_name and target_name - **INEFFICIENT**
   - Size: ~100+ bytes per mapping due to String overhead

2. **PlaneSweepMapping** (temporary during filtering):
   - Uses u32 for positions - **SAME 4GB LIMIT**
   - Compact: only 28 bytes per mapping
   - But loses sequence name information

3. **Grouping Information**:
   - Currently computed on-the-fly by cloning full sequence names
   - Creates temporary tuples like `(mapping, "SGDref#1#chrI", "Y12#1#chrI")`
   - Very inefficient with string cloning for every mapping

## Issues

1. **32-bit position limits**: Can't handle chromosomes/contigs > 4GB
2. **String storage overhead**: Each mapping stores 2 full sequence names
3. **Repeated string cloning**: Grouping creates new string copies
4. **No sequence ID mapping**: Can't use compact integer IDs

## Proposed Optimization

### Phase 1: Use 64-bit positions
- Change all position fields from u32 to u64
- Allows handling of arbitrarily large sequences
- Small memory increase (16 bytes per mapping)

### Phase 2: Implement sequence ID mapping
```rust
struct SequenceIndex {
    names: Vec<String>,           // All unique sequence names
    name_to_id: HashMap<String, u32>, // Map name -> ID
}

struct CompactRecordMeta {
    rank: usize,
    query_id: u32,      // 4 bytes instead of String
    target_id: u32,     // 4 bytes instead of String
    query_start: u64,   // 8 bytes for large sequences
    query_end: u64,     // 8 bytes
    target_start: u64,  // 8 bytes
    target_end: u64,    // 8 bytes
    // ... other fields
}
```

### Phase 3: Pre-compute groups and store by group
```rust
struct GroupedMappings {
    // Group key -> vector of mappings in that group
    groups: HashMap<(u32, u32), Vec<CompactRecordMeta>>,
    sequence_index: SequenceIndex,
}
```

Benefits:
- Direct iteration over groups without re-grouping
- No string operations during filtering
- Cache-friendly access patterns

## Memory Impact

Current (with 30K mappings):
- ~100 bytes/mapping × 30K = ~3MB (mostly String overhead)

Proposed:
- 48 bytes/mapping × 30K = ~1.4MB
- Plus one-time sequence index: ~1KB for typical datasets
- **Net savings: >50%**

## Implementation Priority

1. **URGENT**: Change positions to u64 (fixes 4GB limit)
2. **HIGH**: Implement sequence ID mapping (major memory savings)
3. **MEDIUM**: Pre-grouped storage (performance optimization)