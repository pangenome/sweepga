/// Build script for sweepga — caches FastGA and wfmash binaries so that
/// `cargo install sweepga` produces a self-contained installation.
///
/// At build time the dependency crates (fastga-rs, wfmash-rs) compile their
/// helper binaries into target/{profile}/build/{dep}-*/out/.  We copy those
/// into  ~/.cache/sweepga/{cache_key}/  and remove any stale version dirs.
///
/// At runtime binary_paths.rs checks the cache first, so the binaries are
/// always available regardless of whether a cargo target/ tree still exists.
use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=Cargo.lock");

    // OUT_DIR = target/{profile}/build/sweepga-{hash}/out
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    let build_dir = out_dir
        .ancestors()
        .find(|p| p.file_name().map(|n| n == "build").unwrap_or(false))
        .map(|p| p.to_path_buf());

    let build_dir = match build_dir {
        Some(d) => d,
        None => {
            println!("cargo:warning=Could not locate build directory");
            return;
        }
    };

    // ── cache key ──────────────────────────────────────────────────────
    let cache_key = cache_key();
    println!("cargo:rustc-env=SWEEPGA_CACHE_KEY={cache_key}");

    // ── cache directory ────────────────────────────────────────────────
    let cache_dir = cache_base().join(&cache_key);

    if let Err(e) = std::fs::create_dir_all(&cache_dir) {
        println!("cargo:warning=Failed to create cache dir: {e}");
        return;
    }

    // ── copy binaries ──────────────────────────────────────────────────
    let mut installed = 0usize;

    // FastGA utilities
    let fastga_bins = [
        "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN", "ONEview",
    ];
    if let Some(dir) = find_dep_out(&build_dir, "fastga-rs", "FastGA") {
        installed += copy_binaries(&dir, &cache_dir, &fastga_bins);
    } else {
        println!("cargo:warning=fastga-rs build output not found");
    }

    // wfmash
    if let Some(dir) = find_dep_out(&build_dir, "wfmash-rs", "wfmash") {
        installed += copy_binaries(&dir, &cache_dir, &["wfmash"]);
    } else {
        println!("cargo:warning=wfmash-rs build output not found (may not be built yet)");
    }

    if installed > 0 {
        println!(
            "cargo:warning=Cached {installed} binaries in {}",
            cache_dir.display()
        );
        cleanup_old_versions(&cache_dir, &cache_key);
    }
}

// ── helpers ────────────────────────────────────────────────────────────

/// Build a cache key: `git describe --always --dirty` or CARGO_PKG_VERSION.
fn cache_key() -> String {
    if let Ok(out) = Command::new("git")
        .args(["describe", "--always", "--dirty"])
        .output()
    {
        if out.status.success() {
            let desc = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if !desc.is_empty() {
                return desc;
            }
        }
    }
    env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "unknown".into())
}

/// $XDG_CACHE_HOME/sweepga  or  ~/.cache/sweepga
fn cache_base() -> PathBuf {
    env::var("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            PathBuf::from(env::var("HOME").expect("HOME not set")).join(".cache")
        })
        .join("sweepga")
}

/// Scan `build_dir` for `{prefix}-*/out/` containing `marker_bin`.
fn find_dep_out(build_dir: &Path, prefix: &str, marker_bin: &str) -> Option<PathBuf> {
    let pat = format!("{prefix}-");
    for entry in std::fs::read_dir(build_dir).ok()?.flatten() {
        let path = entry.path();
        if path.is_dir() {
            if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                if name.starts_with(&pat) {
                    let out = path.join("out");
                    if out.join(marker_bin).exists() {
                        return Some(out);
                    }
                }
            }
        }
    }
    None
}

/// Copy listed binaries from `src` to `dst`.  Returns number of successes.
///
/// Uses atomic rename to avoid ETXTBSY ("Text file busy") when the
/// destination binary is currently being executed by another process.
fn copy_binaries(src: &Path, dst: &Path, names: &[&str]) -> usize {
    let mut n = 0;
    for name in names {
        let s = src.join(name);
        if !s.exists() {
            continue;
        }
        let d = dst.join(name);
        let tmp = dst.join(format!(".{name}.tmp"));
        match std::fs::copy(&s, &tmp) {
            Ok(_) => {
                make_executable(&tmp);
                match std::fs::rename(&tmp, &d) {
                    Ok(_) => n += 1,
                    Err(e) => {
                        println!("cargo:warning=Failed to rename {name}: {e}");
                        let _ = std::fs::remove_file(&tmp);
                    }
                }
            }
            Err(e) => println!("cargo:warning=Failed to copy {name}: {e}"),
        }
    }
    n
}

#[cfg(unix)]
fn make_executable(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    if let Ok(m) = std::fs::metadata(path) {
        let mut p = m.permissions();
        p.set_mode(0o755);
        let _ = std::fs::set_permissions(path, p);
    }
}

#[cfg(not(unix))]
fn make_executable(_path: &Path) {}

/// Delete sibling directories under the sweepga cache that aren't `current_key`.
fn cleanup_old_versions(current_dir: &Path, current_key: &str) {
    let parent = match current_dir.parent() {
        Some(p) => p,
        None => return,
    };
    let entries = match std::fs::read_dir(parent) {
        Ok(e) => e,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        if entry.path().is_dir() {
            if let Some(name) = entry.file_name().to_str() {
                if name != current_key {
                    println!(
                        "cargo:warning=Removing old cache: {}",
                        entry.path().display()
                    );
                    let _ = std::fs::remove_dir_all(entry.path());
                }
            }
        }
    }
}
