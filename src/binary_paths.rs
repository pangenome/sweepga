//! Binary path resolution for FastGA and wfmash tools.
//!
//! Search order:
//! 1. ~/.cache/sweepga/{cache_key}/  (installed by build.rs)
//! 2. target/{profile}/build/{dep}-*/out/  (development fallback)
//! 3. PATH  (system fallback)

use anyhow::{anyhow, Result};
use std::env;
use std::path::PathBuf;

/// Cache key baked in at compile time by build.rs.
/// `None` only if build.rs didn't run (shouldn't happen in normal builds).
const CACHE_KEY: Option<&str> = option_env!("SWEEPGA_CACHE_KEY");

/// Return the versioned cache directory if it exists on disk.
pub fn cache_dir() -> Option<PathBuf> {
    let key = CACHE_KEY?;
    let base = if let Ok(xdg) = env::var("XDG_CACHE_HOME") {
        PathBuf::from(xdg)
    } else {
        PathBuf::from(env::var("HOME").ok()?).join(".cache")
    };
    let dir = base.join("sweepga").join(key);
    if dir.is_dir() {
        Some(dir)
    } else {
        None
    }
}

/// Locate a binary by name.
pub fn get_embedded_binary_path(binary_name: &str) -> Result<PathBuf> {
    // 1. Cache directory
    if let Some(dir) = cache_dir() {
        let p = dir.join(binary_name);
        if p.exists() {
            return Ok(p);
        }
    }

    // 2. Cargo build tree (development)
    if let Some(path) = scan_build_tree(binary_name) {
        return Ok(path);
    }

    // 3. System PATH
    if let Ok(output) = std::process::Command::new("which")
        .arg(binary_name)
        .output()
    {
        if output.status.success() {
            let p = PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
            if p.exists() {
                return Ok(p);
            }
        }
    }

    Err(anyhow!(
        "Binary '{binary_name}' not found.\n\
         Searched: cache (~/.cache/sweepga/), build tree, PATH.\n\
         Try running 'cargo build --release' first."
    ))
}

/// Prepend the binary cache/build directory to PATH and set WFMASH_BIN_DIR.
///
/// Call once at startup so that:
/// - FastGA's internal `system()` calls find GIXmake, FAtoGDB, etc.
/// - wfmash-rs finds the wfmash binary via WFMASH_BIN_DIR.
pub fn setup_binary_env() {
    let dir = cache_dir().or_else(|| {
        // Fall back: use whatever directory we find FastGA in
        get_embedded_binary_path("FastGA")
            .ok()
            .and_then(|p| p.parent().map(|d| d.to_path_buf()))
    });

    if let Some(dir) = dir {
        let dir_str = dir.display().to_string();

        // Prepend to PATH for FastGA's system() calls
        let current_path = env::var("PATH").unwrap_or_default();
        env::set_var("PATH", format!("{dir_str}:{current_path}"));

        // Tell wfmash-rs where to find its binary
        env::set_var("WFMASH_BIN_DIR", &dir_str);
    }
}

/// Walk the cargo build tree to find a binary in any dependency's out/ directory.
fn scan_build_tree(binary_name: &str) -> Option<PathBuf> {
    let exe_dir = env::current_exe()
        .ok()
        .and_then(|p| p.canonicalize().ok())
        .and_then(|e| e.parent().map(|p| p.to_path_buf()))?;

    // Walk up to target/
    let mut target_dir = exe_dir;
    while target_dir.file_name().is_some_and(|n| n != "target") {
        if !target_dir.pop() {
            return None;
        }
    }
    if !target_dir.ends_with("target") {
        return None;
    }

    for profile in &["release", "debug"] {
        let build_dir = target_dir.join(profile).join("build");
        if let Ok(entries) = std::fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    let binary_path = path.join("out").join(binary_name);
                    if binary_path.exists() {
                        return Some(binary_path.canonicalize().unwrap_or(binary_path));
                    }
                }
            }
        }
    }
    None
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
