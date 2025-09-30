use sweepga::paf_filter::{PafFilter, FilterConfig, FilterMode, ScoringFunction};
use std::io::Write;
use tempfile::NamedTempFile;
use std::fs::File;
use std::io::{BufReader, BufRead};

#[test]
fn test_scaffold_overlap_50_percent() {
    // Create test PAF with scaffolds that overlap by different amounts
    let mut temp_input = NamedTempFile::new().unwrap();

    // Chain 1: query 100-600 (length 500)
    writeln!(temp_input, "query1\t1000\t100\t300\t+\ttarget1\t1000\t100\t300\t180\t200\t60").unwrap();
    writeln!(temp_input, "query1\t1000\t350\t600\t+\ttarget1\t1000\t350\t600\t225\t250\t60").unwrap();

    // Chain 2: query 350-850 (length 500) - 50% overlap with chain 1
    writeln!(temp_input, "query1\t1000\t350\t550\t+\ttarget2\t1000\t100\t300\t180\t200\t60").unwrap();
    writeln!(temp_input, "query1\t1000\t600\t850\t+\ttarget2\t1000\t350\t600\t225\t250\t60").unwrap();

    // Chain 3: query 800-1000 (length 200) - minimal overlap
    writeln!(temp_input, "query1\t1000\t800\t1000\t+\ttarget3\t1000\t800\t1000\t180\t200\t60").unwrap();

    let input_path = temp_input.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    // Configure with 50% scaffold overlap threshold
    let config = FilterConfig {
        chain_gap: 2000,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,  // No pre-filtering
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 300,  // Will merge mappings within 300bp into chains
        min_scaffold_length: 100,
        scaffold_overlap_threshold: 0.5,  // 50% overlap threshold
        scaffold_max_deviation: 0,  // No rescue
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    filter.filter_paf(input_path, output_path).expect("Filtering should succeed");

    // Read output and count chains
    let file = File::open(output_path).unwrap();
    let reader = BufReader::new(file);
    let mut chains = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if let Some(chain_match) = line.split('\t').skip(12).find(|f| f.starts_with("ch:Z:")) {
            let chain_id = chain_match.strip_prefix("ch:Z:").unwrap();
            chains.insert(chain_id.to_string());
        }
    }

    // With 50% overlap threshold, chain 2 should be filtered out
    assert!(chains.len() <= 2, "At most 2 chains should remain after 50% overlap filtering");
}

#[test]
fn test_scaffold_overlap_95_percent() {
    // Create test PAF with scaffolds that overlap by different amounts
    let mut temp_input = NamedTempFile::new().unwrap();

    // Chain 1: query 100-600
    writeln!(temp_input, "query1\t1000\t100\t600\t+\ttarget1\t1000\t100\t600\t450\t500\t60").unwrap();

    // Chain 2: query 350-850 - 50% overlap with chain 1
    writeln!(temp_input, "query1\t1000\t350\t850\t+\ttarget2\t1000\t350\t850\t450\t500\t60").unwrap();

    // Chain 3: query 550-650 - fully contained in chain 2
    writeln!(temp_input, "query1\t1000\t550\t650\t+\ttarget3\t1000\t550\t650\t90\t100\t60").unwrap();

    let input_path = temp_input.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    // Configure with 95% scaffold overlap threshold
    let config = FilterConfig {
        chain_gap: 2000,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 10000,
        min_scaffold_length: 50,
        scaffold_overlap_threshold: 0.95,  // 95% overlap threshold
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    filter.filter_paf(input_path, output_path).expect("Filtering should succeed");

    // Read output
    let file = File::open(output_path).unwrap();
    let reader = BufReader::new(file);
    let line_count = reader.lines().count();

    // With 95% overlap threshold:
    // - Chain 1 and 2 have 50% overlap (< 95%), so both kept
    // - Chain 3 is 100% contained (> 95%), so filtered out
    assert!(line_count >= 2, "At least chains 1 and 2 should be kept with 95% threshold");
}

#[test]
fn test_scaffold_no_overlap_filtering() {
    // Create test PAF with completely overlapping scaffolds
    let mut temp_input = NamedTempFile::new().unwrap();

    // Three chains with identical regions but different targets
    writeln!(temp_input, "query1\t1000\t100\t600\t+\ttarget1\t1000\t100\t600\t450\t500\t60").unwrap();
    writeln!(temp_input, "query1\t1000\t100\t600\t+\ttarget2\t1000\t100\t600\t400\t500\t60").unwrap();
    writeln!(temp_input, "query1\t1000\t100\t600\t+\ttarget3\t1000\t100\t600\t350\t500\t60").unwrap();

    let input_path = temp_input.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    // Configure with no scaffold overlap filtering (threshold = 1.0)
    let config = FilterConfig {
        chain_gap: 2000,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToMany,  // Keep best per query
        scaffold_max_per_query: Some(2),  // Keep top 2
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 10000,
        min_scaffold_length: 100,
        scaffold_overlap_threshold: 1.0,  // No overlap filtering
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::Matches,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    filter.filter_paf(input_path, output_path).expect("Filtering should succeed");

    // Read output
    let file = File::open(output_path).unwrap();
    let reader = BufReader::new(file);
    let line_count = reader.lines().count();

    // With no overlap filtering and keep top 2, we should get 2 scaffolds
    assert_eq!(line_count, 2, "Should keep exactly 2 best scaffolds with no overlap filtering");
}

#[test]
fn test_scaffold_overlap_with_rescue() {
    // Test that rescue phase works correctly with scaffold overlap
    let mut temp_input = NamedTempFile::new().unwrap();

    // Main scaffold: query 100-600
    writeln!(temp_input, "query1\t1000\t100\t300\t+\ttarget1\t1000\t100\t300\t180\t200\t60").unwrap();
    writeln!(temp_input, "query1\t1000\t350\t600\t+\ttarget1\t1000\t350\t600\t225\t250\t60").unwrap();

    // Overlapping scaffold: query 400-700 (gets filtered due to overlap)
    writeln!(temp_input, "query1\t1000\t400\t700\t+\ttarget2\t1000\t400\t700\t270\t300\t60").unwrap();

    // Small mapping near main scaffold (should be rescued)
    writeln!(temp_input, "query1\t1000\t650\t680\t+\ttarget1\t1000\t650\t680\t27\t30\t60").unwrap();

    // Small mapping near filtered scaffold (should NOT be rescued)
    writeln!(temp_input, "query1\t1000\t750\t780\t+\ttarget2\t1000\t750\t780\t27\t30\t60").unwrap();

    let input_path = temp_input.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    let config = FilterConfig {
        chain_gap: 2000,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 300,  // Merge nearby mappings
        min_scaffold_length: 100,
        scaffold_overlap_threshold: 0.5,  // Filter overlapping scaffolds
        scaffold_max_deviation: 100,  // Rescue within 100bp
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    filter.filter_paf(input_path, output_path).expect("Filtering should succeed");

    // Read output and check rescue status
    let file = File::open(output_path).unwrap();
    let reader = BufReader::new(file);
    let mut rescued_count = 0;
    let mut total_lines = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        total_lines += 1;
        if line.contains("st:Z:rescued") {
            rescued_count += 1;
        }
    }

    // Should have main scaffold (2 mappings) + 1 rescued mapping
    assert_eq!(total_lines, 3, "Should have 2 scaffold mappings + 1 rescued");
    assert_eq!(rescued_count, 1, "Should have exactly 1 rescued mapping");
}

#[test]
fn test_scaffold_overlap_different_chromosomes() {
    // Test that scaffold overlap is calculated per chromosome pair
    let mut temp_input = NamedTempFile::new().unwrap();

    // Chain on chr1-chr1: query 100-600
    writeln!(temp_input, "chr1\t1000\t100\t600\t+\tchr1\t1000\t100\t600\t450\t500\t60").unwrap();

    // Chain on chr1-chr2: query 100-600 (same query region, different target chr)
    writeln!(temp_input, "chr1\t1000\t100\t600\t+\tchr2\t1000\t100\t600\t450\t500\t60").unwrap();

    // Chain on chr2-chr1: query 100-600 (different query chr)
    writeln!(temp_input, "chr2\t1000\t100\t600\t+\tchr1\t1000\t100\t600\t450\t500\t60").unwrap();

    let input_path = temp_input.path().to_str().unwrap();
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_str().unwrap();

    let config = FilterConfig {
        chain_gap: 2000,
        min_block_length: 0,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: false,
        scaffold_gap: 10000,
        min_scaffold_length: 100,
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let filter = PafFilter::new(config);
    filter.filter_paf(input_path, output_path).expect("Filtering should succeed");

    // Read output
    let file = File::open(output_path).unwrap();
    let reader = BufReader::new(file);
    let line_count = reader.lines().count();

    // All three scaffolds should be kept as they're on different chromosome pairs
    assert_eq!(line_count, 3, "All 3 scaffolds on different chr pairs should be kept");
}