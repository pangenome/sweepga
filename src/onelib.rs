/// ONElib FFI bindings for reading/writing .1aln files
///
/// This provides safe Rust wrappers around the ONElib C library
/// for working with .1aln alignment files directly with integer IDs.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::path::Path;
use std::ptr;
use anyhow::{Result, bail, Context};

use crate::mapping::{Mapping, MappingAux, ChainStatus};

// ONElib C types and functions
#[repr(C)]
struct OneFile {
    _private: [u8; 0], // Opaque type
}

#[repr(C)]
struct OneSchema {
    _private: [u8; 0], // Opaque type
}

extern "C" {
    // Schema creation
    fn make_Aln_Schema() -> *mut OneSchema;
    fn oneSchemaDestroy(schema: *mut OneSchema);

    // File operations
    fn oneFileOpenRead(
        path: *const c_char,
        schema: *mut OneSchema,
        file_type: *const c_char,
        nthreads: c_int,
    ) -> *mut OneFile;

    fn oneFileOpenWriteNew(
        path: *const c_char,
        schema: *mut OneSchema,
        file_type: *const c_char,
        is_binary: bool,
        nthreads: c_int,
    ) -> *mut OneFile;

    fn oneFileClose(file: *mut OneFile);

    // Reading
    fn oneReadLine(file: *mut OneFile) -> c_char;

    // Field accessors (from rust_onelib.c)
    fn one_int(file: *mut OneFile, index: c_int) -> i64;
    fn one_real(file: *mut OneFile, index: c_int) -> f64;
    fn one_char(file: *mut OneFile, index: c_int) -> c_char;
    fn one_line_type(file: *mut OneFile) -> c_char;
    fn one_line_count(file: *mut OneFile) -> i64;

    // Field setters (from rust_onelib.c)
    fn one_int_set(file: *mut OneFile, index: c_int, value: i64);
    fn one_real_set(file: *mut OneFile, index: c_int, value: f64);
    fn one_char_set(file: *mut OneFile, index: c_int, value: c_char);

    // Writing
    fn oneWriteLine(file: *mut OneFile, line_type: c_char, list_len: i64, list_buf: *const u8);
}

/// Reader for .1aln files that outputs Mapping structs
pub struct AlnReader {
    file: *mut OneFile,
    schema: *mut OneSchema,
}

impl AlnReader {
    /// Open a .1aln file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;
        let path_cstr = CString::new(path_str)?;

        unsafe {
            let schema = make_Aln_Schema();
            if schema.is_null() {
                bail!("Failed to create .1aln schema");
            }

            let type_cstr = CString::new("aln")?;
            let file = oneFileOpenRead(
                path_cstr.as_ptr(),
                schema,
                type_cstr.as_ptr(),
                1, // single threaded for now
            );

            if file.is_null() {
                oneSchemaDestroy(schema);
                bail!("Failed to open .1aln file: {}", path_str);
            }

            Ok(AlnReader { file, schema })
        }
    }

    /// Read next alignment from the file as Mapping + MappingAux
    /// Returns None when EOF is reached
    pub fn read_mapping(&mut self) -> Result<Option<(Mapping, MappingAux)>> {
        unsafe {
            // Skip to next 'A' record (alignment)
            loop {
                let line_type = one_line_type(self.file);

                // If we're not at an 'A' record, read next line
                if line_type != b'A' as c_char {
                    let next = oneReadLine(self.file);
                    if next == 0 {
                        return Ok(None); // EOF
                    }
                    continue;
                }

                // We're at an 'A' record
                // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
                // Fields: aread, abpos, aepos, bread, bbpos, bepos

                let query_id = one_int(self.file, 0);
                let query_start = one_int(self.file, 1) as u32;
                let query_end = one_int(self.file, 2) as u32;
                let ref_id = one_int(self.file, 3);
                let ref_start = one_int(self.file, 4) as u32;
                let ref_end = one_int(self.file, 5) as u32;

                // Read associated data lines until we hit 'T' (trace) or next 'A'
                let mut matches = 0u64;
                let mut diffs = 0u64;
                let mut is_reverse = false;
                let mut query_len = 0u32;
                let mut ref_len = 0u32;

                loop {
                    let next_type = oneReadLine(self.file);

                    if next_type == 0 {
                        break; // EOF
                    }

                    match next_type as u8 as char {
                        'T' => {
                            // Trace record - marks end of this alignment's data
                            // Read next line to position for next alignment
                            oneReadLine(self.file);
                            break;
                        }
                        'M' => matches = one_int(self.file, 0) as u64,
                        'D' => diffs = one_int(self.file, 0) as u64,
                        'R' => is_reverse = true,
                        'L' => {
                            query_len = one_int(self.file, 0) as u32;
                            ref_len = one_int(self.file, 1) as u32;
                        }
                        'A' => {
                            // Hit next alignment without seeing 'T'
                            // This is valid - not all alignments have traces
                            break;
                        }
                        _ => {
                            // Skip other records (X, C, Q, etc.)
                            continue;
                        }
                    }
                }

                // Calculate identity: matches / (matches + diffs)
                let identity = if matches + diffs > 0 {
                    matches as f64 / (matches + diffs) as f64
                } else {
                    0.0
                };

                // Create mapping structure
                let mut mapping = Mapping {
                    ref_seq_id: ref_id as u32,
                    ref_start_pos: ref_start,
                    query_start_pos: query_start,
                    block_length: ref_end - ref_start,
                    n_merged: 1,
                    conserved_sketches: 0,
                    nuc_identity: 0,
                    flags: 0,
                    kmer_complexity: 100,
                };

                mapping.set_identity(identity);
                if is_reverse {
                    mapping.set_reverse(true);
                }

                // Create aux structure with integer IDs only
                let aux = MappingAux {
                    split_mapping_id: 0,
                    chain_pair_score: f64::MAX,
                    chain_pair_id: i64::MIN,
                    query_seq_id: query_id as i32,
                    query_len,
                    ref_name: String::new(), // Empty - IDs only
                    query_name: String::new(), // Empty - IDs only
                    chain_id: String::new(),
                    chain_status: ChainStatus::Unassigned,
                    ref_len,
                };

                return Ok(Some((mapping, aux)));
            }
        }
    }

    /// Read all mappings from the file
    pub fn read_all(&mut self) -> Result<Vec<(Mapping, MappingAux)>> {
        let mut mappings = Vec::new();

        while let Some(mapping) = self.read_mapping()? {
            mappings.push(mapping);
        }

        Ok(mappings)
    }
}

impl Drop for AlnReader {
    fn drop(&mut self) {
        unsafe {
            if !self.file.is_null() {
                oneFileClose(self.file);
            }
            if !self.schema.is_null() {
                oneSchemaDestroy(self.schema);
            }
        }
    }
}

/// Writer for .1aln files
pub struct AlnWriter {
    file: *mut OneFile,
    schema: *mut OneSchema,
}

impl AlnWriter {
    /// Create a new .1aln file for writing
    pub fn create<P: AsRef<Path>>(path: P, binary: bool) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;
        let path_cstr = CString::new(path_str)?;

        unsafe {
            let schema = make_Aln_Schema();
            if schema.is_null() {
                bail!("Failed to create .1aln schema");
            }

            let type_cstr = CString::new("aln")?;
            let file = oneFileOpenWriteNew(
                path_cstr.as_ptr(),
                schema,
                type_cstr.as_ptr(),
                binary,
                1, // single threaded
            );

            if file.is_null() {
                oneSchemaDestroy(schema);
                bail!("Failed to create .1aln file: {}", path_str);
            }

            Ok(AlnWriter { file, schema })
        }
    }

    /// Write a mapping to the file
    pub fn write_mapping(&mut self, mapping: &Mapping, aux: &MappingAux) -> Result<()> {
        unsafe {
            // Write 'A' record: alignment coordinates
            // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
            // Fields: aread, abpos, aepos, bread, bbpos, bepos
            one_int_set(self.file, 0, aux.query_seq_id as i64);
            one_int_set(self.file, 1, mapping.query_start_pos as i64);
            one_int_set(self.file, 2, (mapping.query_start_pos + mapping.block_length) as i64);
            one_int_set(self.file, 3, mapping.ref_seq_id as i64);
            one_int_set(self.file, 4, mapping.ref_start_pos as i64);
            one_int_set(self.file, 5, (mapping.ref_start_pos + mapping.block_length) as i64);
            oneWriteLine(self.file, b'A' as c_char, 0, std::ptr::null());

            // Write 'L' record: sequence lengths
            one_int_set(self.file, 0, aux.query_len as i64);
            one_int_set(self.file, 1, aux.ref_len as i64);
            oneWriteLine(self.file, b'L' as c_char, 0, std::ptr::null());

            // Write 'R' record if reverse complement
            if mapping.is_reverse() {
                oneWriteLine(self.file, b'R' as c_char, 0, std::ptr::null());
            }

            // Calculate matches and diffs from identity
            // identity = matches / (matches + diffs)
            // So: matches = identity * (matches + diffs)
            // We know block_length, but matches+diffs can be different due to indels
            // For simplicity, assume block_length ≈ matches + diffs
            let identity = mapping.identity() / 100.0; // Convert from percentage
            let block_len = mapping.block_length as u64;
            let matches = (identity * block_len as f64).round() as u64;
            let diffs = block_len - matches;

            // Write 'M' record: matches
            one_int_set(self.file, 0, matches as i64);
            oneWriteLine(self.file, b'M' as c_char, 0, std::ptr::null());

            // Write 'D' record: differences
            one_int_set(self.file, 0, diffs as i64);
            oneWriteLine(self.file, b'D' as c_char, 0, std::ptr::null());

            // Write 'T' record: trace (empty for now - we don't have trace data)
            oneWriteLine(self.file, b'T' as c_char, 0, std::ptr::null());

            Ok(())
        }
    }
}

impl Drop for AlnWriter {
    fn drop(&mut self) {
        unsafe {
            if !self.file.is_null() {
                oneFileClose(self.file);
            }
            if !self.schema.is_null() {
                oneSchemaDestroy(self.schema);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // Only run with actual .1aln files
    fn test_read_1aln() {
        let mut reader = AlnReader::open("test.1aln").unwrap();
        let mappings = reader.read_all().unwrap();
        assert!(mappings.len() > 0);
    }
}
