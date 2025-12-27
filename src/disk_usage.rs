//! Disk usage tracking for sweepga
//!
//! Tracks disk space used by temporary files and alignment outputs.
//! Reports cumulative and peak usage.

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// Global disk usage tracker
static CURRENT_USAGE: AtomicU64 = AtomicU64::new(0);
static PEAK_USAGE: AtomicU64 = AtomicU64::new(0);
static CUMULATIVE_WRITTEN: AtomicU64 = AtomicU64::new(0);
static TRACKED_FILES: Mutex<Vec<(String, u64)>> = Mutex::new(Vec::new());

/// Record that a file was created/written
pub fn track_file_created<P: AsRef<Path>>(path: P) {
    if let Ok(metadata) = std::fs::metadata(path.as_ref()) {
        let size = metadata.len();
        let path_str = path.as_ref().to_string_lossy().to_string();

        // Update cumulative written
        CUMULATIVE_WRITTEN.fetch_add(size, Ordering::Relaxed);

        // Update current usage
        let new_current = CURRENT_USAGE.fetch_add(size, Ordering::Relaxed) + size;

        // Update peak if necessary
        let mut peak = PEAK_USAGE.load(Ordering::Relaxed);
        while new_current > peak {
            match PEAK_USAGE.compare_exchange_weak(
                peak,
                new_current,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(p) => peak = p,
            }
        }

        // Track the file
        if let Ok(mut files) = TRACKED_FILES.lock() {
            files.push((path_str, size));
        }
    }
}

/// Record that a file was deleted
#[allow(dead_code)]
pub fn track_file_deleted<P: AsRef<Path>>(path: P) {
    let path_str = path.as_ref().to_string_lossy().to_string();

    if let Ok(mut files) = TRACKED_FILES.lock() {
        if let Some(pos) = files.iter().position(|(p, _)| p == &path_str) {
            let (_, size) = files.remove(pos);
            CURRENT_USAGE.fetch_sub(size, Ordering::Relaxed);
        }
    }
}

/// Track a file by checking its size (for files we didn't create)
#[allow(dead_code)]
pub fn track_existing_file<P: AsRef<Path>>(path: P) {
    track_file_created(path);
}

/// Add bytes to current usage (for index files tracked by size delta)
pub fn add_bytes(bytes: u64) {
    CUMULATIVE_WRITTEN.fetch_add(bytes, Ordering::Relaxed);
    let new_current = CURRENT_USAGE.fetch_add(bytes, Ordering::Relaxed) + bytes;

    // Update peak if necessary
    let mut peak = PEAK_USAGE.load(Ordering::Relaxed);
    while new_current > peak {
        match PEAK_USAGE.compare_exchange_weak(
            peak,
            new_current,
            Ordering::Relaxed,
            Ordering::Relaxed,
        ) {
            Ok(_) => break,
            Err(p) => peak = p,
        }
    }
}

/// Remove bytes from current usage (for index files deleted)
pub fn remove_bytes(bytes: u64) {
    CURRENT_USAGE.fetch_sub(bytes, Ordering::Relaxed);
}

/// Get current disk usage in bytes
pub fn current_usage() -> u64 {
    CURRENT_USAGE.load(Ordering::Relaxed)
}

/// Get peak disk usage in bytes
pub fn peak_usage() -> u64 {
    PEAK_USAGE.load(Ordering::Relaxed)
}

/// Get cumulative bytes written
pub fn cumulative_written() -> u64 {
    CUMULATIVE_WRITTEN.load(Ordering::Relaxed)
}

/// Format bytes in human-readable form
pub fn format_bytes(bytes: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = KB * 1024;
    const GB: u64 = MB * 1024;

    if bytes >= GB {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    } else {
        format!("{} bytes", bytes)
    }
}

/// Get a summary of disk usage
pub fn summary() -> DiskUsageSummary {
    DiskUsageSummary {
        current: current_usage(),
        peak: peak_usage(),
        cumulative: cumulative_written(),
    }
}

/// Log disk usage summary to stderr
pub fn log_summary() {
    let s = summary();
    eprintln!("[sweepga::disk] Disk usage summary:");
    eprintln!("[sweepga::disk]   Current:    {}", format_bytes(s.current));
    eprintln!("[sweepga::disk]   Peak:       {}", format_bytes(s.peak));
    eprintln!("[sweepga::disk]   Cumulative: {}", format_bytes(s.cumulative));
}

/// Disk usage summary
#[derive(Debug, Clone)]
pub struct DiskUsageSummary {
    pub current: u64,
    pub peak: u64,
    pub cumulative: u64,
}

/// Reset all counters (useful for testing)
#[allow(dead_code)]
pub fn reset() {
    CURRENT_USAGE.store(0, Ordering::Relaxed);
    PEAK_USAGE.store(0, Ordering::Relaxed);
    CUMULATIVE_WRITTEN.store(0, Ordering::Relaxed);
    if let Ok(mut files) = TRACKED_FILES.lock() {
        files.clear();
    }
}

/// Track disk usage for a directory (recursive)
#[allow(dead_code)]
pub fn track_directory<P: AsRef<Path>>(path: P) -> std::io::Result<u64> {
    let mut total = 0u64;
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if metadata.is_file() {
            let size = metadata.len();
            total += size;
            track_file_created(entry.path());
        } else if metadata.is_dir() {
            total += track_directory(entry.path())?;
        }
    }
    Ok(total)
}

/// Scan a directory for FastGA index files and return their total size
/// Index patterns: *.1gdb, .*.ktab.*, *.1gix, *.1bps
pub fn scan_fastga_index_files<P: AsRef<Path>>(dir: P) -> std::io::Result<u64> {
    let dir = dir.as_ref();
    let mut total = 0u64;

    if !dir.is_dir() {
        return Ok(0);
    }

    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if !metadata.is_file() {
            continue;
        }

        let name = entry.file_name();
        let name_str = name.to_string_lossy();

        // Check for FastGA index patterns
        let is_index = name_str.ends_with(".1gdb")
            || name_str.ends_with(".1gix")
            || name_str.ends_with(".1bps")
            || name_str.contains(".ktab.");

        if is_index {
            total += metadata.len();
        }
    }

    Ok(total)
}

/// Track FastGA index files in a directory (records them as created)
/// Returns the total size of index files found
#[allow(dead_code)]
pub fn track_fastga_indexes<P: AsRef<Path>>(dir: P) -> std::io::Result<u64> {
    let dir = dir.as_ref();
    let mut total = 0u64;

    if !dir.is_dir() {
        return Ok(0);
    }

    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if !metadata.is_file() {
            continue;
        }

        let name = entry.file_name();
        let name_str = name.to_string_lossy();

        // Check for FastGA index patterns
        let is_index = name_str.ends_with(".1gdb")
            || name_str.ends_with(".1gix")
            || name_str.ends_with(".1bps")
            || name_str.contains(".ktab.");

        if is_index {
            let size = metadata.len();
            total += size;
            track_file_created(entry.path());
        }
    }

    Ok(total)
}

/// Untrack FastGA index files (call after alignment when they're cleaned up)
#[allow(dead_code)]
pub fn untrack_fastga_indexes<P: AsRef<Path>>(dir: P) {
    let dir = dir.as_ref();

    if let Ok(mut files) = TRACKED_FILES.lock() {
        // Remove entries that look like FastGA indexes in this directory
        let dir_str = dir.to_string_lossy();
        files.retain(|(path, size)| {
            let is_in_dir = path.starts_with(dir_str.as_ref());
            let name = Path::new(path)
                .file_name()
                .map(|n| n.to_string_lossy())
                .unwrap_or_default();
            let is_index = name.ends_with(".1gdb")
                || name.ends_with(".1gix")
                || name.ends_with(".1bps")
                || name.contains(".ktab.");

            if is_in_dir && is_index {
                // Subtract from current usage
                CURRENT_USAGE.fetch_sub(*size, Ordering::Relaxed);
                false // remove from list
            } else {
                true // keep in list
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_disk_tracking() {
        reset();

        // Create a temp file
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "Hello, world! This is test data.").unwrap();
        temp.flush().unwrap();

        // Track it
        track_file_created(temp.path());

        assert!(current_usage() > 0);
        assert!(peak_usage() > 0);
        assert!(cumulative_written() > 0);

        let size_before = current_usage();

        // Track deletion
        track_file_deleted(temp.path());

        assert_eq!(current_usage(), 0);
        assert_eq!(peak_usage(), size_before);
    }
}
