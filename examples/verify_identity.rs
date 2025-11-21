/// Verify that fastga-rs now calculates identity correctly
use anyhow::Result;

fn main() -> Result<()> {
    let aln_path = "test_output.1aln";
    let paf_path = "test_output_converted.paf";

    eprintln!("Comparing identity values from .1aln (via fastga-rs) vs PAF (via ALNtoPAF)...\n");

    // Read PAF divergence values
    let paf_content = std::fs::read_to_string(paf_path)?;
    let paf_dv: Vec<f64> = paf_content
        .lines()
        .take(100)
        .filter_map(|line| {
            line.split('\t')
                .find(|f| f.starts_with("dv:f:"))
                .and_then(|f| f[5..].parse().ok())
        })
        .collect();

    // Read .1aln and extract identity
    let mut reader = fastga_rs::AlnReader::open(aln_path)?;

    let mut matches = 0;
    let mut total = 0;
    let mut max_diff = 0.0f64;

    println!("Rec\tFastGA_ID\tPAF_ID\t\tDiff\t\tMatch?");
    println!("{}", "-".repeat(60));

    for (i, paf_identity) in paf_dv.iter().map(|dv| 1.0 - dv).enumerate() {
        if let Some(aln) = reader.read_alignment()? {
            // Identity is now correctly calculated in fastga-rs!
            let query_span = (aln.query_end - aln.query_start) as f64;
            let fastga_identity = aln.matches as f64 / query_span;

            let diff = (fastga_identity - paf_identity).abs();
            let is_match = diff < 0.0001;

            if is_match {
                matches += 1;
            }
            total += 1;

            max_diff = max_diff.max(diff);

            if i < 20 || !is_match {
                println!(
                    "{}\t{:.6}\t{:.6}\t{:.6}\t{}",
                    i + 1,
                    fastga_identity,
                    paf_identity,
                    diff,
                    if is_match { "✓" } else { "✗" }
                );
            }
        } else {
            break;
        }
    }

    println!("\n{}", "=".repeat(60));
    println!(
        "Result: {matches}/{total} records match (within 0.0001 tolerance)"
    );
    println!("Match rate: {:.1}%", 100.0 * matches as f64 / total as f64);
    println!("Maximum difference: {max_diff:.6}");

    if matches == total {
        println!("\n✅ SUCCESS: All identities match between .1aln and PAF!");
        println!("The X record reading and divergence calculation is working correctly.");
    } else {
        println!("\n⚠️  Some records don't match. This may be due to:");
        println!("   - Precision limits in dv:f: tag storage");
        println!("   - Edge cases in the formula");
    }

    Ok(())
}
