/// Test that all required FastGA helper binaries are installed and accessible
///
/// This test verifies that the build.rs script correctly installed all helper
/// binaries needed by FastGA to function properly.
use sweepga::binary_paths;

#[test]
fn test_all_fastga_binaries_installed() {
    // List of all required FastGA binaries
    // FAtoGDB, GIXmake, and GIXrm are called by FastGA itself via system()
    let required_binaries = vec![
        "FastGA",   // Main aligner
        "FAtoGDB",  // FASTA to GDB converter (called by FastGA.c:4733-4739)
        "GIXmake",  // Genome index creator (called by FastGA.c:4743-4750)
        "GIXrm",    // Index cleanup (called by FastGA.c:170-184)
        "ALNtoPAF", // Alignment format converter
        "PAFtoALN", // PAF to alignment converter
    ];

    let mut all_found = true;
    let mut missing_binaries = Vec::new();

    for binary in &required_binaries {
        match binary_paths::get_embedded_binary_path(binary) {
            Ok(path) => {
                // Verify the binary actually exists and is executable
                assert!(
                    path.exists(),
                    "{} binary path returned but file doesn't exist: {}",
                    binary,
                    path.display()
                );

                #[cfg(unix)]
                {
                    use std::os::unix::fs::PermissionsExt;
                    let metadata = std::fs::metadata(&path)
                        .unwrap_or_else(|_| panic!("Could not read metadata for {binary}"));
                    let is_executable = metadata.permissions().mode() & 0o111 != 0;
                    assert!(
                        is_executable,
                        "{} is not executable: {}",
                        binary,
                        path.display()
                    );
                }

                eprintln!("✓ {} found at: {}", binary, path.display());
            }
            Err(e) => {
                eprintln!("✗ {binary} NOT FOUND: {e}");
                missing_binaries.push(*binary);
                all_found = false;
            }
        }
    }

    if !all_found {
        panic!(
            "Missing required FastGA binaries: {missing_binaries:?}\n\
             This means the build.rs script didn't install all required binaries.\n\
             FastGA will fail at runtime when it tries to call these via system()."
        );
    }
}

#[test]
fn test_oneview_optional_binary() {
    // ONEview is optional - the test should pass whether it's present or not
    match binary_paths::get_embedded_binary_path("ONEview") {
        Ok(path) => {
            eprintln!("✓ ONEview found (optional): {}", path.display());
            assert!(path.exists(), "ONEview path returned but doesn't exist");
        }
        Err(_) => {
            eprintln!("  ONEview not found (this is OK, it's optional)");
        }
    }
}

#[test]
fn test_binary_directory_accessible() {
    // Test that we can determine a binary directory
    // This is what fastga-rs uses to set PATH
    let binaries_to_check = ["FastGA"];

    for binary in &binaries_to_check {
        let path =
            binary_paths::get_embedded_binary_path(binary).expect("FastGA binary should be found");

        let dir = path
            .parent()
            .expect("Binary should have a parent directory");

        assert!(
            dir.is_dir(),
            "Binary directory should exist: {}",
            dir.display()
        );

        eprintln!("Binary directory: {}", dir.display());

        // Verify that FAtoGDB and GIXmake are in the same directory
        // (required for FastGA's system() calls to work)
        for helper in &["FAtoGDB", "GIXmake", "GIXrm"] {
            let helper_path = dir.join(helper);
            assert!(
                helper_path.exists(),
                "{} should be in the same directory as FastGA: {}",
                helper,
                dir.display()
            );
        }
    }
}
