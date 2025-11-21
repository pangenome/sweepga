/// Debug tool to examine raw .1aln coordinates before conversion
use anyhow::Result;

fn main() -> Result<()> {
    let test_1aln = "/tmp/test_format_reading.1aln";

    println!("=== Examining Raw .1aln Coordinates ===\n");

    // Open and read a few alignment records
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

    // Navigate to first A record
    file.goto('A', 0)?;

    let mut count = 0;
    loop {
        let line_type = file.read_line();
        if line_type == '\0' || count >= 15 {
            break;
        }

        if line_type == 'A' {
            count += 1;
            let a_id = file.int(0);
            let a_beg = file.int(1);
            let a_end = file.int(2);
            let b_id = file.int(3);
            let b_beg = file.int(4);
            let b_end = file.int(5);

            // Read associated records to get strand
            let mut is_reverse = false;
            loop {
                let next_type = file.read_line();
                if next_type == 'R' {
                    is_reverse = true;
                } else if next_type == 'T' || next_type == 'A' || next_type == '\0' {
                    break;
                }
            }

            let strand = if is_reverse { '-' } else { '+' };

            println!(
                "Record {count}: a_id={a_id} [{a_beg}-{a_end}] → b_id={b_id} [{b_beg}-{b_end}] {strand}"
            );
        }
    }

    Ok(())
}
