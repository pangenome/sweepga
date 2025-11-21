/// Verbose debug program to show what records are being read
use anyhow::Result;

fn main() -> Result<()> {
    let aln_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "test_output.1aln".to_string());

    eprintln!("Reading {aln_path}...\n");

    let mut reader = fastga_rs::AlnReader::open(&aln_path)?;

    // Read just the first alignment and show details
    if let Some(aln) = reader.read_alignment()? {
        eprintln!("First alignment:");
        eprintln!(
            "  Query: {} ({}-{})",
            aln.query_name, aln.query_start, aln.query_end
        );
        eprintln!(
            "  Target: {} ({}-{})",
            aln.target_name, aln.target_start, aln.target_end
        );
        eprintln!("  Matches: {}", aln.matches);
        eprintln!("  Mismatches: {}", aln.mismatches);
        eprintln!("  Block len: {}", aln.block_len);

        let query_span = (aln.query_end - aln.query_start) as f64;
        let identity_from_matches = aln.matches as f64 / query_span;
        eprintln!("  Identity from matches: {identity_from_matches:.6}");

        // Calculate what we expect from PAF
        eprintln!("\nExpected from PAF dv:f:.0240:");
        eprintln!("  Identity should be: {:.6}", 1.0 - 0.024);
    }

    Ok(())
}
