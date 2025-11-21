/// Test to examine contig offset information from .1aln file
use anyhow::Result;

fn main() -> Result<()> {
    let test_1aln = "/tmp/test_format_reading.1aln";

    println!("=== Examining Contig Offset Information ===\n");

    // Open .1aln file and extract contig offsets
    let schema_text = r#"
P 3 aln
D t 1 3 INT
O g 0
G S 0
O S 1 6 STRING
D G 1 3 INT
D C 1 3 INT
D M 1 8 INT_LIST
O a 0
G A 0
D p 2 3 INT 3 INT
O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
D L 2 3 INT 3 INT
D R 0
D D 1 3 INT
D T 1 8 INT_LIST
D X 1 8 INT_LIST
D Q 1 3 INT
D E 1 3 INT
D Z 1 6 STRING
"#;

    let schema = onecode::OneSchema::from_text(schema_text)?;
    let mut file = onecode::OneFile::open_read(test_1aln, Some(&schema), Some("aln"), 1)?;

    let contigs = file.get_all_contig_offsets();

    println!("Found {} contigs\n", contigs.len());
    println!("Contig offset information:");

    let mut sorted_contigs: Vec<_> = contigs.iter().collect();
    sorted_contigs.sort_by_key(|(id, _)| **id);

    for (id, (sbeg, clen)) in sorted_contigs.iter().take(20) {
        println!("  Contig {id}: sbeg={sbeg}, clen={clen}");
    }

    Ok(())
}
