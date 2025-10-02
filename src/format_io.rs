/// Unified I/O for PAF and .1aln formats
///
/// This module provides format-agnostic reading and writing of alignments.
/// Both formats are converted to/from the common (Mapping, MappingAux) representation.

use anyhow::Result;
use std::path::Path;

use crate::mapping::{Mapping, MappingAux};
use crate::seq_registry::SequenceRegistry;

/// Read alignments from file (PAF or .1aln) into memory
pub fn read_alignments<P: AsRef<Path>>(
    path: P,
    format: &str,
) -> Result<(Vec<(Mapping, MappingAux)>, SequenceRegistry)> {
    match format {
        "paf" => read_paf(path),
        "1aln" => read_1aln(path),
        _ => anyhow::bail!("Unknown format: {}", format),
    }
}

/// Read PAF file into memory
fn read_paf<P: AsRef<Path>>(
    path: P,
) -> Result<(Vec<(Mapping, MappingAux)>, SequenceRegistry)> {
    use crate::paf::PafReader;
    use std::fs::File;

    let file = File::open(path)?;
    let mut reader = PafReader::new(file);
    let mut mappings = Vec::new();

    // Read all records
    while let Some((_paf_record, mapping, aux)) = reader.read_record()? {
        mappings.push((mapping, aux));
    }

    // Extract the registry from the reader
    let registry = reader.into_registry();

    Ok((mappings, registry))
}

/// Read .1aln file into memory
fn read_1aln<P: AsRef<Path>>(
    path: P,
) -> Result<(Vec<(Mapping, MappingAux)>, SequenceRegistry)> {
    use crate::onelib::AlnReader;

    let mut reader = AlnReader::open(path)?;
    let mappings = reader.read_all()?;

    // For .1aln, we don't have string names, just IDs
    // Create an empty registry - names not needed for .1aln→.1aln workflow
    let registry = SequenceRegistry::new();

    Ok((mappings, registry))
}

/// Write alignments to file (PAF or .1aln)
pub fn write_alignments<P: AsRef<Path>>(
    path: P,
    format: &str,
    mappings: &[(Mapping, MappingAux)],
    registry: &SequenceRegistry,
) -> Result<()> {
    match format {
        "paf" => write_paf(path, mappings, registry),
        "1aln" => write_1aln(path, mappings),
        _ => anyhow::bail!("Unknown format: {}", format),
    }
}

/// Write PAF file from memory
fn write_paf<P: AsRef<Path>>(
    path: P,
    mappings: &[(Mapping, MappingAux)],
    registry: &SequenceRegistry,
) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (mapping, aux) in mappings {
        // Look up sequence names from IDs
        let query_name = registry
            .get_query_name(aux.query_seq_id)
            .unwrap_or("unknown");
        let ref_name = registry
            .get_ref_name(mapping.ref_seq_id)
            .unwrap_or("unknown");

        // Calculate identity from mapping
        let identity = mapping.identity() / 100.0; // Convert from percentage

        // Copy packed fields to avoid unaligned references
        let query_start = mapping.query_start_pos;
        let ref_start = mapping.ref_start_pos;
        let block_len = mapping.block_length;
        let matches = (identity * block_len as f64).round() as u64;

        // Write PAF line
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255",
            query_name,
            aux.query_len,
            query_start,
            query_start + block_len,
            if mapping.is_reverse() { '-' } else { '+' },
            ref_name,
            aux.ref_len,
            ref_start,
            ref_start + block_len,
            matches,
            block_len,
        )?;
    }

    Ok(())
}

/// Write .1aln file from memory
fn write_1aln<P: AsRef<Path>>(
    path: P,
    mappings: &[(Mapping, MappingAux)],
) -> Result<()> {
    use crate::onelib::AlnWriter;

    let mut writer = AlnWriter::create(path, true)?; // true = binary format

    for (mapping, aux) in mappings {
        writer.write_mapping(mapping, aux)?;
    }

    Ok(())
}
