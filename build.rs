/// Build script for sweepga - copies FastGA helper binaries for cargo install
///
/// This script ensures that when users do `cargo install sweepga`, all the FastGA
/// helper binaries (FAtoGDB, GIXmake, GIXrm, etc.) are available at runtime.
///
/// The binaries are built by the fastga-rs dependency and placed in its OUT_DIR.
/// We need to copy them to our OUT_DIR so they get installed alongside sweepga.
use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Find the fastga-rs build directory
    // During build, dependencies' OUT_DIRs are in target/<profile>/build/<dep-name>-<hash>/out/
    let target_dir = out_dir
        .ancestors()
        .find(|p| p.ends_with("build"))
        .expect("Could not find build directory");

    // List of FastGA utilities we need to copy
    let utilities = [
        "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN", "ONEview",
    ];

    println!("cargo:warning=Searching for FastGA binaries...");

    // Search for fastga-rs build output directory
    if let Ok(entries) = std::fs::read_dir(target_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir()
                && path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .map(|n| n.starts_with("fastga-rs-"))
                    .unwrap_or(false)
            {
                let fastga_out_dir = path.join("out");
                if fastga_out_dir.exists() {
                    println!(
                        "cargo:warning=Found FastGA binaries in: {}",
                        fastga_out_dir.display()
                    );

                    // Copy each utility to our OUT_DIR
                    for utility in &utilities {
                        let src = fastga_out_dir.join(utility);
                        if src.exists() {
                            let dst = out_dir.join(utility);
                            if let Err(e) = std::fs::copy(&src, &dst) {
                                println!("cargo:warning=Failed to copy {utility}: {e}");
                            } else {
                                println!("cargo:warning=Copied {utility}");

                                // Make executable on Unix
                                #[cfg(unix)]
                                {
                                    use std::os::unix::fs::PermissionsExt;
                                    if let Ok(metadata) = std::fs::metadata(&dst) {
                                        let mut perms = metadata.permissions();
                                        perms.set_mode(0o755);
                                        let _ = std::fs::set_permissions(&dst, perms);
                                    }
                                }
                            }
                        } else if utility != &"ONEview" {
                            // ONEview is optional
                            println!(
                                "cargo:warning={} not found at: {}",
                                utility,
                                src.display()
                            );
                        }
                    }

                    break;
                }
            }
        }
    }

    // Also copy binaries to $CARGO_HOME/lib/sweepga/ for cargo install
    // This happens during every build, ensuring binaries are available after install
    if let Ok(cargo_home) =
        env::var("CARGO_HOME").or_else(|_| env::var("HOME").map(|h| format!("{h}/.cargo")))
    {
        let lib_dir = PathBuf::from(cargo_home).join("lib").join("sweepga");

        // Create the directory
        if std::fs::create_dir_all(&lib_dir).is_ok() {
            println!(
                "cargo:warning=Installing FastGA helper binaries to: {}",
                lib_dir.display()
            );

            // Copy binaries to lib directory
            for utility in &utilities {
                let src = out_dir.join(utility);
                if src.exists() {
                    let dst = lib_dir.join(utility);
                    if std::fs::copy(&src, &dst).is_ok() {
                        #[cfg(unix)]
                        {
                            use std::os::unix::fs::PermissionsExt;
                            if let Ok(metadata) = std::fs::metadata(&dst) {
                                let mut perms = metadata.permissions();
                                perms.set_mode(0o755);
                                let _ = std::fs::set_permissions(&dst, perms);
                            }
                        }
                    }
                }
            }
        }
    }

    // Set environment variable pointing to our OUT_DIR
    // This will be baked into the binary at compile time
    println!("cargo:rustc-env=SWEEPGA_BINARY_DIR={}", out_dir.display());
}
