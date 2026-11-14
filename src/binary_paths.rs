//! Binary path resolution for FastGA and wfmash tools.
//!
//! Search order:
//! 1. ~/.cache/sweepga/{cache_key}/  (installed by build.rs)
//! 2. target/{profile}/build/{dep}-{hash}/out/ or .../build/{dep}/{hash}/out/
//!    (development fallback, both cargo build-directory layouts)
//! 3. <exe_dir>/../libexec/<exe_stem>/  (packaged installs, e.g. conda)
//! 4. PATH  (system fallback)

use anyhow::{anyhow, Result};
use std::env;
use std::path::{Path, PathBuf};

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

    // 3. Private install directory (packaged installs, e.g. conda)
    if let Some(dir) = libexec_dir() {
        let p = dir.join(binary_name);
        if p.exists() {
            return Ok(p);
        }
    }

    // 4. System PATH
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
         Searched: cache (~/.cache/sweepga/), build tree, libexec, PATH.\n\
         Try running 'cargo build --release' first."
    ))
}

/// Private install directory next to the running executable:
/// `<exe_dir>/../libexec/<exe_stem>/`.
///
/// Package managers use this to keep bundled wfmash and FastGA binaries out of
/// `bin/`, where they would overwrite the separately packaged versions.
pub fn libexec_dir() -> Option<PathBuf> {
    let exe = env::current_exe().ok()?;
    let dir = exe
        .parent()?
        .parent()?
        .join("libexec")
        .join(exe.file_stem()?);
    if dir.is_dir() {
        Some(dir)
    } else {
        None
    }
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
        if let Some(p) = newest_in_build(&target_dir.join(profile).join("build"), binary_name) {
            return Some(p.canonicalize().unwrap_or(p));
        }
    }
    None
}

/// Newest `binary_name` in any `out` dir of `build_dir`.
/// Layouts: `{dep}-{hash}/out` and `{dep}/{hash}/out`.
fn newest_in_build(build_dir: &Path, binary_name: &str) -> Option<PathBuf> {
    let dirs = |d: &Path| {
        std::fs::read_dir(d)
            .into_iter()
            .flatten()
            .flatten()
            .map(|e| e.path())
            .collect::<Vec<_>>()
    };
    dirs(build_dir)
        .into_iter()
        .flat_map(|d| {
            let mut v = dirs(&d);
            v.push(d);
            v
        })
        .map(|d| d.join("out").join(binary_name))
        .filter_map(|p| Some((p.metadata().ok()?.modified().ok()?, p)))
        .max_by_key(|(t, _)| *t)
        .map(|(_, p)| p)
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

    #[test]
    fn test_newest_in_build_both_layouts() {
        let tmp = tempfile::tempdir().unwrap();
        let build = tmp.path().join("build");
        let put = |rel: &str, age: u64| {
            let p = build.join(rel).join("out").join("FastGA");
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            let f = std::fs::File::create(&p).unwrap();
            let t = std::time::SystemTime::now() - std::time::Duration::from_secs(age);
            f.set_modified(t).unwrap();
            p
        };
        assert_eq!(newest_in_build(&build, "FastGA"), None);
        let old = put("fastga-rs-aaaa", 300);
        assert_eq!(newest_in_build(&build, "FastGA"), Some(old.clone()));
        let new = put("fastga-rs/bbbb", 200);
        assert_eq!(newest_in_build(&build, "FastGA"), Some(new));
        std::fs::remove_dir_all(build.join("fastga-rs-aaaa")).unwrap();
        let stale = put("fastga-rs/cccc", 400);
        assert_ne!(newest_in_build(&build, "FastGA"), Some(stale));
        let old2 = put("fastga-rs-dddd", 100);
        assert_eq!(newest_in_build(&build, "FastGA"), Some(old2));
        assert_eq!(newest_in_build(&build, "wfmash"), None);
    }
}
