//! Binary path resolution for embedded FastGA tools
//!
//! Scans target/build/fastga-rs-*/out/ for binaries built by the fastga-rs
//! dependency. Falls back to PATH if not found.

use anyhow::{anyhow, Result};
use std::env;
use std::path::PathBuf;

/// Get the path to an embedded FastGA binary.
///
/// Search order:
/// 1. target/{debug,release}/build/fastga-rs-*/out/{binary_name}
/// 2. PATH (system fallback)
pub fn get_embedded_binary_path(binary_name: &str) -> Result<PathBuf> {
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
            return Ok(path);
        }
    }

    Err(anyhow!(
        "FastGA binary '{binary_name}' not found in embedded build directories or PATH.\n\
         Try running 'cargo build' or 'cargo build --release' first."
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_embedded_binaries() {
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
                }
                Err(e) => {
                    panic!("{bin} binary not found: {e}");
                }
            }
        }
    }
}
