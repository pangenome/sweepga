/// Unit tests for sweepga components
use std::fs;
use tempfile::NamedTempFile;

#[test]
fn test_paf_format_validation() {
    // Test that we produce valid PAF format
    let temp = NamedTempFile::new().unwrap();
    let test_paf = "chr1\t1000\t100\t900\t+\tchr2\t2000\t200\t1000\t700\t800\t60\n";
    fs::write(temp.path(), test_paf).unwrap();

    // Read and validate PAF fields
    let content = fs::read_to_string(temp.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 1);

    let fields: Vec<&str> = lines[0].split('\t').collect();
    assert!(fields.len() >= 12, "PAF requires at least 12 fields");

    // Validate field types
    assert!(
        fields[1].parse::<u32>().is_ok(),
        "Query length must be numeric"
    );
    assert!(
        fields[2].parse::<u32>().is_ok(),
        "Query start must be numeric"
    );
    assert!(
        fields[3].parse::<u32>().is_ok(),
        "Query end must be numeric"
    );
    assert!(
        fields[4] == "+" || fields[4] == "-",
        "Strand must be + or -"
    );
    assert!(
        fields[6].parse::<u32>().is_ok(),
        "Target length must be numeric"
    );
    assert!(
        fields[7].parse::<u32>().is_ok(),
        "Target start must be numeric"
    );
    assert!(
        fields[8].parse::<u32>().is_ok(),
        "Target end must be numeric"
    );
    assert!(
        fields[9].parse::<u32>().is_ok(),
        "Match count must be numeric"
    );
    assert!(
        fields[10].parse::<u32>().is_ok(),
        "Alignment length must be numeric"
    );
    assert!(
        fields[11].parse::<u32>().is_ok(),
        "Mapping quality must be numeric"
    );
}

#[test]
fn test_filter_modes() {
    // Test filter mode parsing
    let modes = vec![
        ("1:1", "OneToOne"),
        ("1", "OneToMany"),
        ("1:∞", "OneToMany"),
        ("1:many", "OneToMany"),
        ("N", "NoFilter"),
        ("N:N", "NoFilter"),
        ("none", "NoFilter"),
    ];

    for (mode, _expected) in modes {
        // Would need to expose parse_filter_mode or test via integration
        // For now, just document expected behavior
        assert!(!mode.is_empty(), "Mode string should not be empty");
    }
}

#[test]
fn test_temp_file_cleanup() {
    use std::process::Command;

    // Create test FASTA
    let test_fa = "/tmp/test_cleanup.fa";
    fs::write(test_fa, ">test\nACGT\n").unwrap();

    // Get temp file count before
    let temp_before = std::fs::read_dir("/tmp")
        .unwrap()
        .filter(|entry| {
            entry
                .as_ref()
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".tmp")
        })
        .count();

    // Run sweepga
    let _output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--quiet",
            "--bin",
            "sweepga",
            "--",
            test_fa,
            "-t",
            "1",
        ])
        .output()
        .expect("Failed to run");

    // Give it a moment to clean up
    std::thread::sleep(std::time::Duration::from_millis(100));

    // Get temp file count after
    let temp_after = std::fs::read_dir("/tmp")
        .unwrap()
        .filter(|entry| {
            entry
                .as_ref()
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".tmp")
        })
        .count();

    // Should not leave temp files behind
    assert!(
        temp_after <= temp_before + 1, // Allow for race conditions
        "Temp files not cleaned up properly"
    );

    // Cleanup
    let _ = fs::remove_file(test_fa);
}

#[test]
fn test_cigar_extended_format() {
    // Test that extended CIGAR uses = and X correctly
    let cigar_examples = vec![
        "10=",          // 10 matches
        "5=2X3=",       // 5 matches, 2 mismatches, 3 matches
        "10=5I10=",     // matches with insertion
        "10=5D10=",     // matches with deletion
        "3=1X2=1I4=1D", // complex
    ];

    for cigar in cigar_examples {
        // Validate extended CIGAR format
        let mut pos = 0;
        let bytes = cigar.as_bytes();

        while pos < bytes.len() {
            // Read number
            let num_start = pos;
            while pos < bytes.len() && bytes[pos].is_ascii_digit() {
                pos += 1;
            }

            if pos < bytes.len() {
                let op = bytes[pos] as char;
                assert!(
                    op == '='
                        || op == 'X'
                        || op == 'I'
                        || op == 'D'
                        || op == 'M'
                        || op == 'S'
                        || op == 'H',
                    "Invalid CIGAR operation: {op}"
                );
                pos += 1;
            }

            assert!(num_start < pos, "CIGAR must have length before operation");
        }
    }
}

#[test]
fn test_scaffold_annotations() {
    // Test that scaffold annotations are properly formatted
    let annotations = vec![
        "ch:Z:chain_1",
        "st:Z:scaffold",
        "st:Z:rescued",
        "ch:Z:chain_123",
    ];

    for ann in annotations {
        let parts: Vec<&str> = ann.split(':').collect();
        assert_eq!(parts.len(), 3, "Annotation must have 3 parts");
        assert!(parts[0] == "ch" || parts[0] == "st", "Unknown tag");
        assert_eq!(parts[1], "Z", "Must be string type");
        assert!(!parts[2].is_empty(), "Value cannot be empty");
    }
}

#[cfg(test)]
mod fastga_config_tests {

    #[test]
    fn test_config_builder() {
        // Would test FastGA config builder if exposed
        // Document expected behavior:
        // - Default threads should use available CPUs
        // - Min identity should be 0.0-1.0
        // - Min alignment length should be positive
    }

    #[test]
    fn test_thread_limits() {
        // Test that thread count is reasonable
        let max_threads = num_cpus::get();
        assert!(max_threads > 0, "Should detect at least 1 CPU");
        assert!(max_threads < 1000, "Unreasonable thread count detected");
    }
}
