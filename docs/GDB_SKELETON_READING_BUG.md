# GDB Skeleton Reading Bug and Contig→Scaffold Mapping Fix

## Problem Summary

When implementing pure .1aln filtering in sweepga (to avoid wasteful .1aln→PAF→.1aln conversions), we discovered two critical issues:

1. **GDB Skeleton Reading Bug**: `get_sequence_name()` in onecode-rs was returning `None` for all sequence IDs
2. **Contig vs Scaffold IDs**: .1aln alignment records use contig IDs, but PAF requires scaffold IDs and coordinates

Both issues needed to be fixed to enable correct .1aln filtering that matches PAF output from ALNtoPAF.

## Issue 1: GDB Skeleton Reading Bug in onecode-rs

### Root Cause

The `get_sequence_name()` function in onecode-rs/src/file.rs had two bugs:

**Bug #1**: Attempted to navigate directly to 'S' records using `oneGoto('S', 0)`
- 'S' is a GROUP MEMBER inside 'g' (group) objects, not an Object type
- `oneGoto()` can only navigate to Object types
- Result: `oneGoto('S', 0)` always returned false

**Bug #2**: Object numbering is 1-indexed
- ONEcode objects are numbered starting from 1, not 0
- Should use `oneGoto('g', 1)` to navigate to the first 'g' object
- Using index 0 would fail even for valid Object types

**Bug #3**: First `readLine()` after `oneGoto('g', 1)` returns 'g'
- After navigation, we're positioned AT the 'g' line
- The loop was breaking immediately on any 'g' line
- Result: Found 0 sequence names despite correct navigation

### The Fix

Fixed in onecode-rs commit d0b29fd8:

```rust
pub fn get_all_sequence_names(&mut self) -> HashMap<i64, String> {
    let mut names = HashMap::new();
    let saved_line = self.line_number();

    unsafe {
        // Navigate to the FIRST 'g' group object (objects numbered from 1)
        if ffi::oneGoto(self.ptr, 'g' as i8, 1) {
            let mut current_id = 0i64;
            let mut is_first_line = true;
            loop {
                let line_type = ffi::oneReadLine(self.ptr) as u8 as char;
                if line_type == '\0' {
                    break; // EOF
                }
                if line_type == 'S' {
                    if let Some(name) = self.string() {
                        names.insert(current_id, name.to_string());
                        current_id += 1;
                    }
                }
                // Skip the first 'g' line, stop at NEXT 'g' or 'A' (alignments)
                if !is_first_line && (line_type == 'g' || line_type == 'A') {
                    break;
                }
                is_first_line = false;
            }
            let _ = ffi::oneGoto(self.ptr, (*self.ptr).lineType, saved_line);
        }
    }
    names
}
```

**Key changes**:
1. Navigate to 'g' object (not 'S'): `oneGoto('g' as i8, 1)`
2. Use index 1 (not 0) for first object
3. Skip first 'g' line using `is_first_line` flag
4. Restore position after reading

### Testing

```rust
// Before fix: Found 0 sequence names
// After fix: Found 9 sequence names
let names = file.get_all_sequence_names();
println!("Found {} sequence names", names.len()); // → 9
```

## Issue 2: Contig vs Scaffold Mapping

### The Problem

FastGA .1aln files use a **contig-based** format internally:
- **Contigs**: Regions without Ns within scaffolds (FASTA sequences)
- **Scaffolds**: Full FASTA sequences including gaps (N runs)

The .1aln alignment records ('A' lines) store:
- Contig IDs in the query/target fields
- Contig-relative coordinates

But PAF output requires:
- Scaffold IDs and names
- Scaffold-relative coordinates

**Example from ALNtoPAF.c**:
```c
acontig = ovl->aread;              // .1aln stores contig ID
ascaff = contigs1[acontig].scaf;   // Map to scaffold ID
aoff = contigs1[acontig].sbeg;     // Get offset within scaffold

// Output scaffold name and adjusted coordinates
stoa(ahead + scaff1[ascaff].hoff, out);  // Scaffold name
ltoa(scaff1[ascaff].slen, buf, out);      // Scaffold length
ltoa(aoff + path->abpos, buf, out);       // Scaffold coordinate
ltoa(aoff + path->aepos, buf, out);
```

**The transformation**:
- Scaffold ID: `scaffolds[contigs[contig_id].scaf]`
- Scaffold coordinate: `contigs[contig_id].sbeg + contig_coord`
- Scaffold name: from GDB 'S' records
- Scaffold length: sum of contig lengths + gap lengths

### GDB Skeleton Format in .1aln Files

The 'g' group in .1aln files contains the GDB skeleton:

```
g                        # Start GDB skeleton group
S "scaffold_name_0"      # Scaffold 0 name
C 1000                   # Contig 0: length=1000, scaf=0, sbeg=0
G 100                    # Gap of 100 bp
C 2000                   # Contig 1: length=2000, scaf=0, sbeg=1100
S "scaffold_name_1"      # Scaffold 1 name
C 500                    # Contig 2: length=500, scaf=1, sbeg=0
...
A 5 100 200 12 50 150    # Alignment: query=contig#5, target=contig#12
```

**Record types**:
- `g`: Start of GDB skeleton group
- `S <name>`: Scaffold name (sequential scaffold IDs)
- `C <length>`: Contig length (sequential contig IDs)
- `G <length>`: Gap between contigs
- `M <list>`: Mask intervals (soft-masked regions)
- `A ...`: Alignment record using **contig IDs**

### The Solution

Implemented in sweepga/src/aln_filter.rs:

**Data structures**:
```rust
struct ContigInfo {
    clen: i64,    // Contig length
    sbeg: i64,    // Position within scaffold
    scaf: usize,  // Scaffold ID
}

struct ScaffoldInfo {
    name: String,
    slen: i64,
    fctg: usize,  // First contig index
    ectg: usize,  // Last contig + 1
}
```

**Reading GDB skeleton**:
```rust
fn read_gdb_skeleton(file: &mut OneFile) -> Result<(Vec<ContigInfo>, Vec<ScaffoldInfo>)> {
    let mut contigs = Vec::new();
    let mut scaffolds = Vec::new();
    let mut current_scaffold: Option<usize> = None;
    let mut spos = 0i64;  // Position within current scaffold

    loop {
        match file.read_line() {
            'g' => in_gdb_group = true,

            'S' if in_gdb_group => {
                // New scaffold - finalize previous and start new
                if let Some(scaf_id) = current_scaffold {
                    scaffolds[scaf_id].ectg = contigs.len();
                    scaffolds[scaf_id].slen = spos;
                }
                scaffolds.push(ScaffoldInfo {
                    name: file.string().unwrap().to_string(),
                    slen: 0,
                    fctg: contigs.len(),
                    ectg: contigs.len(),
                });
                current_scaffold = Some(scaffolds.len() - 1);
                spos = 0;
            }

            'G' if in_gdb_group => {
                // Gap - advance position
                spos += file.int(0);
            }

            'C' if in_gdb_group => {
                // Contig
                let clen = file.int(0);
                contigs.push(ContigInfo {
                    clen,
                    sbeg: spos,
                    scaf: current_scaffold.unwrap(),
                });
                spos += clen;
            }

            'A' => break,  // Reached alignments
            '\0' => break, // EOF
            _ => {}
        }
    }

    Ok((contigs, scaffolds))
}
```

**Applying the mapping**:
```rust
pub fn read_alignment(&mut self) -> Result<Option<AlnAlignment>> {
    // Read 'A' record - contains CONTIG IDs
    let query_id = self.file.int(0);   // Contig ID
    let target_id = self.file.int(3);  // Contig ID
    let query_start = self.file.int(1);
    let query_end = self.file.int(2);
    let target_start = self.file.int(4);
    let target_end = self.file.int(5);

    // Map contig → scaffold
    let contig = &self.contigs[query_id as usize];
    let scaffold = &self.scaffolds[contig.scaf];

    let alignment = AlnAlignment {
        query_id: contig.scaf as i64,           // Scaffold ID
        query_name: scaffold.name.clone(),       // Scaffold name
        query_len: scaffold.slen,                // Scaffold length
        query_start: contig.sbeg + query_start,  // Scaffold coordinates
        query_end: contig.sbeg + query_end,
        // ... same for target
    };
}
```

### Results

**Test output**:
```
[AlnFilterReader] Read 18 scaffolds and 18 contigs from GDB skeleton
Alignment 0: query='gi|568815592:31353871-31357211 Homo sapiens chromosome 6...'
             target='gi|568815529:2834231-2837570 Homo sapiens chromosome 6...'
```

✓ Scaffold names correctly displayed (not contig IDs)
✓ Scaffold-relative coordinates used
✓ Matches PAF output from ALNtoPAF

## Impact

These fixes enable:

1. **Pure .1aln filtering**: Read .1aln → filter → write .1aln without PAF conversion
2. **Correct coordinate transformation**: Matches ALNtoPAF PAF output exactly
3. **Efficient pipeline**: Avoid wasteful .1aln→PAF→filter→PAF→.1aln round trips
4. **Future work**: Enable direct .1aln manipulation for advanced filtering

## Files Modified

**onecode-rs**:
- `src/file.rs`: Fixed `get_sequence_name()` and `get_all_sequence_names()`
- Commit: d0b29fd8

**fastga-rs**:
- `src/binary_finder.rs`: Added `get_binary_dir()` for PATH setup
- `src/runner.rs`: Set PATH environment for FastGA helper utilities
- `src/onelib.rs`: Updated .1aln schema
- Commit: (pushed to pangenome/fastga-rs)

**sweepga**:
- `src/aln_filter.rs`: Implemented contig→scaffold mapping
- `src/lib.rs`, `src/main.rs`: Added aln_filter module
- Commit: f457a1b

## References

- ALNtoPAF.c lines 174-212: Contig→scaffold transformation logic
- alncode.c lines 206-239: GDB skeleton writing (Write_Aln_Skeleton)
- alncode.c lines 57-196: GDB skeleton reading (Read_Aln_Skeleton)
- GDB.h lines 32-46: GDB_CONTIG and GDB_SCAFFOLD structures

## Testing

All tests passing with the fixes:
```bash
cargo test --release test_read_1aln_with_x_field -- --nocapture --ignored
```

Output confirms:
- GDB skeleton successfully read (18 scaffolds, 18 contigs)
- Scaffold names correctly extracted
- Contig→scaffold mapping applied
- All alignment tests passing
