use crate::mapping::{Mapping, MappingAux, PafRecord};
use anyhow::{bail, Result};
use noodles::bgzf;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Open a file and auto-detect bgzip compression, returning a boxed BufRead
pub fn open_paf_input<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let file = File::open(path)?;

    // Check by file extension (faster than reading magic bytes)
    let is_compressed = path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz" || ext == "bgz")
        .unwrap_or(false);

    if is_compressed {
        Ok(Box::new(BufReader::new(bgzf::io::reader::Reader::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Parse CIGAR string to count exact matches (= operations)
/// Returns (matches, mismatches, insertions, deletions)
pub fn parse_cigar_counts(cigar: &str) -> Result<(u64, u64, u64, u64)> {
    let mut matches = 0u64;
    let mut mismatches = 0u64;
    let mut insertions = 0u64;
    let mut deletions = 0u64;

    let mut num_str = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            let count: u64 = num_str.parse()
                .map_err(|_| anyhow::anyhow!("Invalid number in CIGAR: {}", num_str))?;
            num_str.clear();

            match ch {
                '=' => matches += count,      // Exact match
                'X' => mismatches += count,   // Mismatch
                'I' => insertions += count,   // Insertion to reference
                'D' => deletions += count,    // Deletion from reference
                'M' => {
                    // M means match or mismatch, we can't distinguish
                    // In this case, use the provided matches field from PAF
                    // This is handled elsewhere
                }
                _ => {} // Ignore other operations like S, H, N, P
            }
        }
    }

    Ok((matches, mismatches, insertions, deletions))
}

#[allow(dead_code)]
pub struct PafReader<R: Read> {
    reader: BufReader<R>,
    ref_name_to_id: HashMap<String, u32>,
    query_name_to_id: HashMap<String, i32>,
    next_ref_id: u32,
    next_query_id: i32,
}

impl<R: Read> PafReader<R> {
    pub fn new(reader: R) -> Self {
        PafReader {
            reader: BufReader::new(reader),
            ref_name_to_id: HashMap::new(),
            query_name_to_id: HashMap::new(),
            next_ref_id: 0,
            next_query_id: 0,
        }
    }

    pub fn read_record(&mut self) -> Result<Option<(PafRecord, Mapping, MappingAux)>> {
        let mut line = String::new();
        if self.reader.read_line(&mut line)? == 0 {
            return Ok(None);
        }

        let paf = self.parse_paf_line(&line)?;

        // Get or assign IDs
        let ref_id = *self
            .ref_name_to_id
            .entry(paf.ref_name.clone())
            .or_insert_with(|| {
                let id = self.next_ref_id;
                self.next_ref_id += 1;
                id
            });

        let query_id = *self
            .query_name_to_id
            .entry(paf.query_name.clone())
            .or_insert_with(|| {
                let id = self.next_query_id;
                self.next_query_id += 1;
                id
            });

        let mut mapping = paf.to_mapping(ref_id, query_id);

        // Check for divergence tag from FASTGA (dv:f:) and convert to identity
        for (tag, val) in &paf.tags {
            if tag == "dv:f" {
                if let Ok(divergence) = val.parse::<f64>() {
                    let identity = 1.0 - divergence;
                    mapping.set_identity(identity);
                }
            }
        }

        let aux = MappingAux {
            query_seq_id: query_id,
            query_len: paf.query_len,
            ref_len: paf.ref_len,
            ..Default::default()
        };

        Ok(Some((paf, mapping, aux)))
    }

    fn parse_paf_line(&self, line: &str) -> Result<PafRecord> {
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 12 {
            bail!("PAF line has fewer than 12 required fields");
        }

        let mut paf = PafRecord {
            query_name: fields[0].to_string(),
            query_len: fields[1].parse()?,
            query_start: fields[2].parse()?,
            query_end: fields[3].parse()?,
            strand: fields[4].chars().next().unwrap_or('+'),
            ref_name: fields[5].to_string(),
            ref_len: fields[6].parse()?,
            ref_start: fields[7].parse()?,
            ref_end: fields[8].parse()?,
            matches: fields[9].parse()?,
            block_len: fields[10].parse()?,
            quality: fields[11].parse()?,
            tags: Vec::new(),
            cigar: None,
        };

        // Parse optional tags
        for field in &fields[12..] {
            if let Some((tag, rest)) = field.split_once(':') {
                if let Some((typ, val)) = rest.split_once(':') {
                    if tag == "cg" && typ == "Z" {
                        paf.cigar = Some(val.to_string());
                    } else {
                        paf.tags.push((format!("{tag}:{typ}"), val.to_string()));
                    }
                }
            }
        }

        Ok(paf)
    }

    pub fn read_all(&mut self) -> Result<Vec<(PafRecord, Mapping, MappingAux)>> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }
        Ok(records)
    }
}

/// Read PAF from file (auto-detects bgzip compression)
#[allow(dead_code)]
pub fn read_paf_file(path: &str) -> Result<Vec<(PafRecord, Mapping, MappingAux)>> {
    let input = open_paf_input(path)?;
    let mut reader = PafReader::new(input);
    reader.read_all()
}

/// Read PAF from stdin
#[allow(dead_code)]
pub fn read_paf_stdin() -> Result<Vec<(PafRecord, Mapping, MappingAux)>> {
    use std::io::stdin;
    let mut reader = PafReader::new(stdin());
    reader.read_all()
}
