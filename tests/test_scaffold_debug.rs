use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn debug_scaffold_output() {
    // Create test PAF
    let mut file = NamedTempFile::new().unwrap();
    writeln!(
        file,
        "chrA\t500000\t0\t10000\t+\tchr1\t600000\t0\t10000\t10000\t10000\t60\tcg:Z:10000M"
    ).unwrap();
    writeln!(
        file,
        "chrA\t500000\t10000\t20000\t+\tchr1\t600000\t10000\t20000\t10000\t10000\t60\tcg:Z:10000M"
    ).unwrap();
    writeln!(
        file,
        "chrA\t500000\t20000\t30000\t+\tchr1\t600000\t20000\t30000\t10000\t10000\t60\tcg:Z:10000M"
    ).unwrap();
    file.flush().unwrap();

    // Run with scaffolding
    let cmd = std::process::Command::new("./target/release/sweepga")
        .arg(file.path().to_str().unwrap())
        .arg("-j").arg("20k")  // Enable scaffolding
        .arg("-s").arg("15k")  // Min scaffold size
        .arg("-d").arg("0")    // No rescue
        .output()
        .expect("Failed to run sweepga");

    println!("STDERR:\n{}", String::from_utf8_lossy(&cmd.stderr));
    println!("STDOUT:\n{}", String::from_utf8_lossy(&cmd.stdout));

    let output = String::from_utf8_lossy(&cmd.stdout);
    println!("OUTPUT:\n{}", output);

    assert!(!output.is_empty(), "Output should not be empty");
}