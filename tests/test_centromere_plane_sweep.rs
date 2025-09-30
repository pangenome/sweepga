/// Test that large reverse-strand scaffolds aren't filtered out by plane sweep
/// when they should win over smaller forward-strand scaffolds
///
/// This reproduces a bug where an 8Mb reverse-strand centromere inversion
/// was being filtered out in favor of smaller forward-strand alignments
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_reverse_strand_scaffold_plane_sweep() {
    // Create synthetic mappings based on real HG002 chr1 centromere data
    // Forward strand: 129-133Mb region, ~3.8Mb, 76% identity
    // Reverse strand: 129-137Mb region, ~8Mb, 76% identity
    // The reverse strand should WIN because it's larger!
    
    let mut paf = NamedTempFile::new().unwrap();
    
    // Forward-strand alignment spanning 129-133Mb
    writeln!(paf, "query\t250000000\t129142789\t132986703\t+\ttarget\t250000000\t129142789\t132986703\t2938926\t3843914\t60\tNM:i:904988\tcg:Z:2938926=904988X").unwrap();
    
    // Reverse-strand alignment spanning 129-137Mb (LARGER, should win)
    writeln!(paf, "query\t250000000\t129213003\t137240549\t-\ttarget\t250000000\t131937578\t139967018\t6372479\t8027546\t60\tNM:i:1655067\tcg:Z:6372479=1655067X").unwrap();
    
    paf.flush().unwrap();
    
    // With Y=0 (no identity filter) and scaffolding enabled
    let output = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-Y")
        .arg("0")
        .arg("-j")
        .arg("100000")  // Enable scaffolding
        .output()
        .expect("Failed to run sweepga");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    
    eprintln!("Output:\n{}", stderr);
    eprintln!("Mappings:\n{}", stdout);
    
    // Check which alignment(s) survived
    let has_forward = stdout.contains("\t+\t");
    let has_reverse = stdout.contains("\t-\t");
    
    // Count lines
    let forward_count = stdout.lines().filter(|l| l.contains("\t+\t") && !l.starts_with('[')).count();
    let reverse_count = stdout.lines().filter(|l| l.contains("\t-\t") && !l.starts_with('[')).count();
    
    eprintln!("Forward-strand mappings: {}", forward_count);
    eprintln!("Reverse-strand mappings: {}", reverse_count);
    
    // The 8Mb reverse-strand alignment should be kept
    // It's twice the size of the forward-strand alignment!
    assert!(has_reverse, 
            "8Mb reverse-strand alignment should NOT be filtered out by plane sweep");
    
    // In fact, the reverse-strand should win because it's larger
    assert!(reverse_count > 0, 
            "Reverse-strand alignment (8Mb) should beat forward-strand alignment (3.8Mb)");
}

#[test]
fn test_reverse_vs_forward_scaffold_scoring() {
    // Simpler test: two overlapping scaffolds, reverse is larger
    let mut paf = NamedTempFile::new().unwrap();
    
    // Forward: 1Mb at 95% identity
    writeln!(paf, "query\t100000000\t10000000\t11000000\t+\ttarget\t100000000\t10000000\t11000000\t950000\t1000000\t60\tNM:i:50000\tcg:Z:950000=50000X").unwrap();
    
    // Reverse: 2Mb at 95% identity, overlaps the same query region
    // Should WIN because it's 2x larger with same identity
    writeln!(paf, "query\t100000000\t10000000\t12000000\t-\ttarget\t100000000\t20000000\t22000000\t1900000\t2000000\t60\tNM:i:100000\tcg:Z:1900000=100000X").unwrap();
    
    paf.flush().unwrap();
    
    let output = std::process::Command::new("./target/release/sweepga")
        .arg("-i")
        .arg(paf.path())
        .arg("-Y")
        .arg("0")
        .arg("-j")
        .arg("100000")
        .output()
        .expect("Failed to run sweepga");
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    
    eprintln!("Test 2 Output:\n{}", stderr);
    eprintln!("Test 2 Mappings:\n{}", stdout);
    
    let has_forward = stdout.contains("\t+\t");
    let has_reverse = stdout.contains("\t-\t");
    
    // The 2Mb reverse-strand should win over the 1Mb forward-strand
    assert!(has_reverse, "Larger reverse-strand alignment should be kept");
    
    // In an ideal world, only the reverse would be kept
    // But at minimum, the reverse MUST be present
    eprintln!("Has forward: {}, Has reverse: {}", has_forward, has_reverse);
}
