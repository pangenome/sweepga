use fastga_rs::AlnReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Opening test_output.1aln...");
    let mut reader = AlnReader::open("test_output.1aln")?;
    
    println!("✓ File opened successfully!");
    
    println!("Reading first record...");
    if let Some(rec) = reader.read_record()? {
        println!("✓ First record read successfully!");
        println!("  query_id={}, target_id={}", rec.query_id, rec.target_id);
        
        let qname = reader.get_seq_name(rec.query_id, 0)?;
        let tname = reader.get_seq_name(rec.target_id, 1)?;
        println!("  query_name={}", qname);
        println!("  target_name={}", tname);
    }
    
    println!("\n✓ SUCCESS: AlnReader works correctly!");
    Ok(())
}
