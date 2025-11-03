// ! Binary path resolution for embedded FastGA tools
//!
//! This module provides robust binary path resolution that:
//! 1. Prioritizes embedded binaries built by cargo in target/build directories
//! 2. Falls back to PATH only if embedded binaries aren't found
//! 3. Provides clear error messages when binaries are missing

use anyhow::{anyhow, Result};
use std::env;
use std::path::PathBuf;

/// Get the path to an embedded FastGA binary.
///
/// Search order:
/// 1. OUT_DIR (set during build time only)
/// 2. Installed location: $CARGO_HOME/lib/sweepga/ or ~/.cargo/lib/sweepga/
/// 3. Development build: target/{debug,release}/build/fastga-rs-*/out/{binary_name}
/// 4. PATH (system fallback - least preferred)
///
/// This ensures we use the FastGA binaries built as part of our dependency chain,
/// not system-installed versions which may be incompatible.
pub fn get_embedded_binary_path(binary_name: &str) -> Result<PathBuf> {
    // Try OUT_DIR first (set during build)
    if let Ok(out_dir) = env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join(binary_name);
        if path.exists() {
            return Ok(path);
        }
    }

    // Check installed location (for cargo install)
    // Try $CARGO_HOME/lib/sweepga/ or ~/.cargo/lib/sweepga/
    let cargo_home = env::var("CARGO_HOME").ok().map(PathBuf::from).or_else(|| {
        env::var("HOME")
            .ok()
            .map(|h| PathBuf::from(h).join(".cargo"))
    });

    if let Some(cargo_home) = cargo_home {
        let installed_path = cargo_home.join("lib").join("sweepga").join(binary_name);
        if installed_path.exists() {
            return Ok(installed_path);
        }
    }

    // Get the executable's directory to find target/
    let exe_dir = env::current_exe()
        .ok()
        .and_then(|exe| exe.parent().map(|p| p.to_path_buf()));

    if let Some(mut target_dir) = exe_dir {
        // Navigate up to find the target/ directory
        // Executable might be in target/release/ or target/release/deps/
        while target_dir.file_name().is_some_and(|n| n != "target") {
            if !target_dir.pop() {
                break;
            }
        }

        // Now target_dir should point to target/
        if target_dir.ends_with("target") {
            // Search in target/{debug,release}/build/fastga-rs-*/out/
            for profile in &["release", "debug"] {
                let build_dir = target_dir.join(profile).join("build");
                if build_dir.exists() {
                    if let Ok(entries) = std::fs::read_dir(&build_dir) {
                        for entry in entries.flatten() {
                            let path = entry.path();
                            if path.is_dir()
                                && path
                                    .file_name()
                                    .and_then(|n| n.to_str())
                                    .is_some_and(|n| n.starts_with("fastga-rs-"))
                            {
                                let binary_path = path.join("out").join(binary_name);
                                if binary_path.exists() {
                                    return Ok(binary_path);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Last resort: check PATH (system binaries)
    if let Ok(output) = std::process::Command::new("which")
        .arg(binary_name)
        .output()
    {
        if output.status.success() {
            let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
            let path = PathBuf::from(path);

            // Warn if using system binary
            if !path.to_string_lossy().contains("/build/fastga-rs-") {
                eprintln!(
                    "WARNING: Using system binary for {}: {}",
                    binary_name,
                    path.display()
                );
                eprintln!("         This may cause compatibility issues. Consider running 'cargo build' again.");
            }

            return Ok(path);
        }
    }

    Err(anyhow!(
        "FastGA binary '{}' not found in embedded build directories or PATH.\n\
         Try running 'cargo build' or 'cargo build --release' first.",
        binary_name
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_embedded_binaries() {
        // This test verifies we can find the embedded binaries
        // Note: FAtoGDB, GIXmake, and GIXrm are required by FastGA itself
        let bins = [
            "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN",
        ];

        for bin in &bins {
            match get_embedded_binary_path(bin) {
                Ok(path) => {
                    assert!(
                        path.exists(),
                        "{} binary not found at: {}",
                        bin,
                        path.display()
                    );

                    // Verify it's an embedded binary (in build dir)
                    let path_str = path.to_string_lossy();
                    let is_embedded = path_str.contains("/build/fastga-rs-");

                    if !is_embedded {
                        println!(
                            "WARNING: {} is not an embedded binary: {}",
                            bin,
                            path.display()
                        );
                    }
                }
                Err(e) => {
                    panic!("{} binary not found: {}", bin, e);
                }
            }
        }
    }
}
