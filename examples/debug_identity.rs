/// Debug program to understand identity calculation from .1aln format
///
/// Compares values from .1aln file with corresponding PAF file from ALNtoPAF
use anyhow::Result;

fn main() -> Result<()> {
    let aln_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "test_output.1aln".to_string());

    eprintln!("Reading {}...\n", aln_path);

    let mut reader = fastga_rs::AlnReader::open(&aln_path)?;
    let id_to_name = reader.get_all_seq_names();

    eprintln!("Found {} sequences\n", id_to_name.len());

    println!("RecordNum,QSpan,TSpan,Matches,Diffs,NaiveID,CalcDV,Strand");

    let mut count = 0;
    while let Some(aln) = reader.read_alignment()? {
        if count >= 50 {
            break;
        }
        count += 1;

        let query_span = (aln.query_end - aln.query_start) as i64;
        let target_span = (aln.target_end - aln.target_start) as i64;
        let matches = aln.matches as i64;
        let diffs = aln.mismatches as i64;

        // Naive: matches / (matches + diffs)
        let alignment_len = matches + diffs;
        let naive_id = if alignment_len > 0 {
            matches as f64 / alignment_len as f64
        } else {
            0.0
        };

        // ALNtoPAF formula: dv = (diffs - del) / query_span
        // where del = target_span - query_span (deletions in query)
        let del = target_span - query_span;
        let calc_dv = if query_span > 0 {
            (diffs - del) as f64 / query_span as f64
        } else {
            0.0
        };

        println!(
            "{},{},{},{},{},{:.6},{:.6},{}",
            count, query_span, target_span, matches, diffs, naive_id, calc_dv, aln.strand
        );
    }

    Ok(())
}
