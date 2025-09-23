use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;

/// Helper to create a temporary PAF file
fn create_temp_paf(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    writeln!(file, "{}", content).unwrap();
    file
}

/// Helper to count lines matching a pattern in a file
fn count_lines_matching(path: &str, pattern: &str) -> usize {
    let content = fs::read_to_string(path).unwrap();
    content.lines()
        .filter(|line| line.contains(pattern))
        .count()
}

/// Helper to check if specific inter-chromosomal rescues exist
fn check_inter_chromosomal_rescues(path: &str) -> usize {
    let content = fs::read_to_string(path).unwrap();
    let mut inter_chrom_count = 0;

    for line in content.lines() {
        if line.contains("st:Z:rescued") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 5 {
                let query = fields[0];
                let target = fields[5];

                // Extract chromosome names
                let query_chr = query.split('#').last().unwrap_or("");
                let target_chr = target.split('#').last().unwrap_or("");

                if query_chr != target_chr {
                    inter_chrom_count += 1;
                }
            }
        }
    }

    inter_chrom_count
}

#[test]
fn test_no_inter_chromosomal_rescue_bug() {
    // Create test PAF with mappings that should NOT be rescued across chromosomes
    let paf_content = "\
S288C#1#chrI\t230218\t9\t12602\t+\tY12#1#chrI\t219929\t9\t12627\t11986\t12619\t255\tdv:f:0.05\n\
S288C#1#chrI\t230218\t12530\t12746\t+\tY12#1#chrVIII\t547529\t12530\t12746\t210\t216\t255\tdv:f:0.03\n\
S288C#1#chrVIII\t543240\t100\t10100\t+\tY12#1#chrVIII\t547529\t100\t10100\t9950\t10000\t255\tdv:f:0.005";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run sweepga filter
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    // Check that chrI->chrVIII mapping is NOT rescued
    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let rescued_inter_chrom = output_content.lines()
        .find(|line| {
            line.contains("S288C#1#chrI") &&
            line.contains("Y12#1#chrVIII") &&
            line.contains("st:Z:rescued")
        });

    assert!(rescued_inter_chrom.is_none(),
            "Bug: S288C#chrI -> Y12#chrVIII should NOT be rescued by chrVIII->chrVIII anchor!");

    // Count total inter-chromosomal rescues (should be 0)
    let inter_chrom_rescues = check_inter_chromosomal_rescues(output_file.path().to_str().unwrap());
    assert_eq!(inter_chrom_rescues, 0,
               "Found {} inter-chromosomal rescues, expected 0", inter_chrom_rescues);
}

#[test]
fn test_scaffold_filtering_basic() {
    // Test basic scaffold creation and filtering
    let paf_content = "\
query1\t100000\t0\t11000\t+\ttarget1\t100000\t0\t11000\t10900\t11000\t255\tdv:f:0.01\n\
query1\t100000\t11100\t21000\t+\ttarget1\t100000\t11100\t21000\t9800\t9900\t255\tdv:f:0.01\n\
query1\t100000\t21100\t31000\t+\ttarget1\t100000\t21100\t31000\t9800\t9900\t255\tdv:f:0.01\n\
query1\t100000\t50000\t51000\t+\ttarget1\t100000\t50000\t51000\t990\t1000\t255\tdv:f:0.01";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run with scaffold filtering (default: 10kb minimum)
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-S")
        .arg("10000")  // 10kb minimum scaffold
        .arg("-D")
        .arg("10000")  // 10kb max deviation (instead of default 100kb)
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    // Check scaffolds were created
    let scaffold_count = count_lines_matching(output_file.path().to_str().unwrap(), "st:Z:scaffold");
    assert!(scaffold_count > 0, "Should have created scaffolds");

    // Small isolated mapping should not be kept
    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let small_mapping = output_content.lines()
        .find(|line| line.contains("50000\t51000"));
    assert!(small_mapping.is_none(), "Small isolated mapping should be filtered out");
}

#[test]
fn test_one_to_one_filtering() {
    // Test 1:1 filtering mode
    let paf_content = "\
query1\t100000\t0\t5000\t+\ttarget1\t100000\t0\t5000\t4900\t5000\t255\tdv:f:0.02\n\
query1\t100000\t10000\t20000\t+\ttarget2\t100000\t0\t10000\t9800\t10000\t255\tdv:f:0.02\n\
query2\t100000\t0\t15000\t+\ttarget1\t100000\t10000\t25000\t14700\t15000\t255\tdv:f:0.02\n\
query2\t100000\t20000\t25000\t+\ttarget2\t100000\t20000\t25000\t4900\t5000\t255\tdv:f:0.02";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run with 1:1 filtering
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-m")
        .arg("1:1")
        .arg("-S")
        .arg("0")  // Disable scaffold filtering for this test
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    // In 1:1 mode, each query and target should have at most 1 mapping
    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();

    // Count mappings per query
    let query1_count = lines.iter().filter(|l| l.starts_with("query1")).count();
    let query2_count = lines.iter().filter(|l| l.starts_with("query2")).count();

    assert!(query1_count <= 1, "query1 should have at most 1 mapping in 1:1 mode");
    assert!(query2_count <= 1, "query2 should have at most 1 mapping in 1:1 mode");
}

#[test]
fn test_plane_sweep_filtering() {
    // Test plane sweep (1:∞) filtering - keeps best overlapping mapping per query position
    let paf_content = "\
query1\t100000\t0\t5000\t+\ttarget1\t100000\t0\t5000\t4900\t5000\t255\tdv:f:0.02\n\
query1\t100000\t2000\t8000\t+\ttarget2\t100000\t0\t6000\t5880\t6000\t255\tdv:f:0.02\n\
query1\t100000\t1000\t3000\t+\ttarget3\t100000\t0\t2000\t1960\t2000\t255\tdv:f:0.02\n\
query1\t100000\t20000\t30000\t+\ttarget1\t100000\t20000\t30000\t9800\t10000\t255\tdv:f:0.02";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run with plane sweep (1:∞)
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-m")
        .arg("1:∞")
        .arg("-S")
        .arg("0")  // Disable scaffold filtering for this test
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();

    // Should keep the longest overlapping mapping in each query region
    // The 2000-8000 mapping (length 6000) should win over 0-5000 (length 5000)
    // and 1000-3000 (length 2000) in the overlapping region
    let has_long_mapping = lines.iter().any(|l| l.contains("2000\t8000"));
    assert!(has_long_mapping, "Should keep the 2000-8000 mapping (longest in overlap)");

    // The 20000-30000 mapping should be kept (no overlap)
    let has_distant_mapping = lines.iter().any(|l| l.contains("20000\t30000"));
    assert!(has_distant_mapping, "Should keep the non-overlapping 20000-30000 mapping");
}

#[test]
fn test_rescue_within_chromosome_pair() {
    // Test that rescue works correctly WITHIN the same chromosome pair
    let paf_content = "\
S288C#chrI\t230218\t0\t12000\t+\tY12#chrI\t219929\t0\t12000\t11800\t12000\t255\tdv:f:0.02\n\
S288C#chrI\t230218\t12100\t12500\t+\tY12#chrI\t219929\t12100\t12500\t390\t400\t255\tdv:f:0.02\n\
S288C#chrII\t813597\t0\t11000\t+\tY12#chrII\t800000\t0\t11000\t10800\t11000\t255\tdv:f:0.02\n\
S288C#chrII\t813597\t50000\t50400\t+\tY12#chrII\t800000\t50000\t50400\t390\t400\t255\tdv:f:0.02";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run with default scaffold filtering
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-S")
        .arg("10000")  // 10kb minimum scaffold
        .arg("-D")
        .arg("1000")   // 1kb max deviation for rescue
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content = fs::read_to_string(output_file.path()).unwrap();

    // The 12000bp mapping should be a scaffold
    let has_scaffold = output_content.lines()
        .any(|l| l.contains("0\t12000") && l.contains("st:Z:scaffold"));
    assert!(has_scaffold, "Large mapping should be marked as scaffold");

    // The nearby small mapping (12100-12500) should be rescued
    let has_rescued_nearby = output_content.lines()
        .any(|l| l.contains("12100\t12500") && l.contains("st:Z:rescued"));
    assert!(has_rescued_nearby, "Small mapping near scaffold should be rescued");

    // The distant small mapping (50000-50400) should NOT be rescued (too far)
    let has_distant = output_content.lines()
        .any(|l| l.contains("50000\t50400"));
    assert!(!has_distant, "Distant small mapping should NOT be rescued (too far from scaffold)");
}

#[test]
fn test_self_mapping_filtering() {
    // Test that self-mappings are filtered by default
    let paf_content = "\
chr1\t100000\t0\t5000\t+\tchr1\t100000\t10000\t15000\t4900\t5000\t255\tdv:f:0.02\n\
chr1\t100000\t0\t5000\t+\tchr2\t100000\t0\t5000\t4900\t5000\t255\tdv:f:0.02";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run without --self flag (default: filter self-mappings)
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-S")
        .arg("0")  // Disable scaffold filtering
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();

    // Should only have chr1->chr2 mapping, not chr1->chr1
    assert_eq!(lines.len(), 1, "Should have exactly 1 mapping (self-mapping filtered)");
    assert!(lines[0].contains("chr2"), "Should keep chr1->chr2 mapping");

    // Now test WITH --self flag
    let output_file2 = NamedTempFile::new().unwrap();
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file2.path())
        .arg("--self")
        .arg("-S")
        .arg("0")  // Disable scaffold filtering
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content2 = fs::read_to_string(output_file2.path()).unwrap();
    let lines2: Vec<&str> = output_content2.lines().collect();

    assert_eq!(lines2.len(), 2, "Should have 2 mappings with --self flag");
}

#[test]
fn test_prefix_grouping() {
    // Test that filtering modes work on prefix pairs, not chromosome pairs
    // genome1#chr1 -> genome2#chr1 should compete with genome1#chr2 -> genome2#chr2
    // in 1:1 mode when grouping by prefix
    let paf_content = "\
genome1#chr1\t100000\t0\t10000\t+\tgenome2#chr1\t100000\t0\t10000\t9800\t10000\t255\tdv:f:0.02\n\
genome1#chr2\t100000\t0\t15000\t+\tgenome2#chr2\t100000\t0\t15000\t14700\t15000\t255\tdv:f:0.02\n\
genome1#chr1\t100000\t20000\t25000\t+\tgenome2#chr3\t100000\t0\t5000\t4900\t5000\t255\tdv:f:0.02";

    let input_file = create_temp_paf(paf_content);
    let output_file = NamedTempFile::new().unwrap();

    // Run with 1:1 filtering and prefix grouping enabled (default)
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .arg("-m")
        .arg("1:1")
        .arg("-Y")
        .arg("#")  // Group by prefix before #
        .arg("-S")
        .arg("0")  // Disable scaffold filtering
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content = fs::read_to_string(output_file.path()).unwrap();
    let lines: Vec<&str> = output_content.lines().collect();

    // With prefix grouping, genome1->genome2 is ONE pair
    // In 1:1 mode, should keep only the longest mapping (15kb chr2->chr2)
    assert_eq!(lines.len(), 1, "Should keep only 1 mapping in 1:1 mode with prefix grouping");
    assert!(lines[0].contains("chr2"), "Should keep the longest mapping (chr2->chr2)");

    // Now test WITHOUT prefix grouping
    let output_file2 = NamedTempFile::new().unwrap();
    let status = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file2.path())
        .arg("-m")
        .arg("1:1")
        .arg("-Y")
        .arg("")  // Disable prefix grouping
        .arg("-S")
        .arg("0")
        .status()
        .expect("Failed to run sweepga");

    assert!(status.success());

    let output_content2 = fs::read_to_string(output_file2.path()).unwrap();
    let lines2: Vec<&str> = output_content2.lines().collect();

    // Without prefix grouping, each chromosome is independent
    // But in 1:1 mode, each query chromosome can only map to one target
    // genome1#chr1 has two mappings (->chr1 and ->chr3), only best is kept (10kb)
    // genome1#chr2 has one mapping (->chr2), it's kept (15kb)
    assert_eq!(lines2.len(), 2, "Should keep 2 mappings in 1:1 mode without prefix grouping");
}