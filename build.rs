/// Build script for sweepga — caches FastGA and wfmash binaries so that
/// `cargo install sweepga` produces a self-contained installation.
///
/// At build time the dependency crates (fastga-rs, wfmash-rs) compile their
/// helper binaries into target/{profile}/build/{dep}-{hash}/out/ (cargo up to
/// 1.99) or target/{profile}/build/{dep}/{hash}/out/ (cargo 1.100).  We copy those
/// into  ~/.cache/sweepga/{cache_key}/  and remove any stale version dirs.
///
/// At runtime binary_paths.rs checks the cache first, so the binaries are
/// always available regardless of whether a cargo target/ tree still exists.
use std::env;
use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=Cargo.lock");

    // OUT_DIR is below target/{profile}/build
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
        println!(
            "cargo:warning=FastGA not found: no fastga-rs out dir under {}",
            build_dir.display()
        );
    }

    // wfmash
    if let Some(dir) = find_dep_out(&build_dir, "wfmash-rs", "wfmash") {
        installed += copy_binaries(&dir, &cache_dir, &["wfmash"]);
    } else {
        println!(
            "cargo:warning=wfmash not found: no wfmash-rs out dir under {}",
            build_dir.display()
        );
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

/// Build a cache key from package version + dependency git revisions.
///
/// The key changes when `Cargo.lock` changes (which triggers a build.rs
/// rerun), ensuring cached binaries always match the compiled deps.
/// We do NOT use `git describe` because git state changes without any
/// file changing, which means Cargo won't re-run build.rs and the
/// cache key goes stale.
fn cache_key() -> String {
    let version = env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "0".into());

    // Extract short git revisions from Cargo.lock for our binary deps
    let fastga_rev = dep_rev_from_lockfile("fastga-rs").unwrap_or_default();
    let wfmash_rev = dep_rev_from_lockfile("wfmash-rs").unwrap_or_default();

    if fastga_rev.is_empty() && wfmash_rev.is_empty() {
        version
    } else {
        format!("{version}_f{fastga_rev}_w{wfmash_rev}")
    }
}

/// Extract the short git revision for a dependency from Cargo.lock.
/// Looks for lines like: `source = "git+...#01eb1ffe5df54b45..."`
fn dep_rev_from_lockfile(dep_name: &str) -> Option<String> {
    let lock = std::fs::read_to_string("Cargo.lock").ok()?;
    let mut in_dep = false;
    for line in lock.lines() {
        if line.starts_with("name = ") && line.contains(dep_name) {
            in_dep = true;
        } else if in_dep && line.starts_with("source = ") {
            // source = "git+https://...#01eb1ffe5df54b45123fee40cf69a8beb04fd1fc"
            if let Some(hash) = line.rsplit('#').next() {
                let hash = hash.trim_end_matches('"');
                return Some(hash[..8.min(hash.len())].to_string());
            }
        } else if line.starts_with("[[") {
            in_dep = false;
        }
    }
    None
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

/// Newest dependency `out` dir holding `marker_bin`.
/// Layouts: `{prefix}-{hash}/out` and `{prefix}/{hash}/out`.
fn find_dep_out(build_dir: &Path, prefix: &str, marker_bin: &str) -> Option<PathBuf> {
    let pat = format!("{prefix}-");
    let mut outs = Vec::new();
    let entries = std::fs::read_dir(build_dir).ok()?;
    for path in entries.flatten().map(|e| e.path()) {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with(&pat) {
            outs.push(path.join("out"));
        } else if name == prefix {
            let subs = std::fs::read_dir(&path).into_iter().flatten();
            outs.extend(subs.flatten().map(|e| e.path().join("out")));
        }
    }
    outs.into_iter()
        .filter_map(|o| Some((o.join(marker_bin).metadata().ok()?.modified().ok()?, o)))
        .max_by_key(|(t, _)| *t)
        .map(|(_, o)| o)
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
            println!("cargo:warning={name} not found in {}", src.display());
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
