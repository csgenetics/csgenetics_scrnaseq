// io_count_extract: streaming extractor for the io_count SAM text stage.
//
// Reads SAM records (tab-delimited, no header) on stdin, writes one line per record that
// carries a featureCounts gene assignment in an XT:Z: optional field.
//
// Output per matching record: <barcode>\t<gene>\n
//   barcode = the [A-Z]+ token between the trailing "_..._" of field 1 (read name), else empty
//   gene    = the complete XT:Z: value, with a terminal numeric version suffix removed
//
// Streaming and low-memory (the process is capped at 1 GB): we never hold the whole input.

use std::io::{self, Read, Write};

#[inline]
fn is_upper(b: u8) -> bool {
    b.is_ascii_uppercase()
}

// Match features_names.py: remove only a terminal `.<ASCII digits>` suffix. Gene IDs are
// otherwise opaque bytes. In particular, punctuation and UTF-8 are preserved exactly.
#[inline]
fn normalize_gene_id(gene: &[u8]) -> &[u8] {
    if let Some(dot) = gene.iter().rposition(|&b| b == b'.') {
        let suffix = &gene[dot + 1..];
        if !suffix.is_empty() && suffix.iter().all(u8::is_ascii_digit) {
            return &gene[..dot];
        }
    }
    gene
}

fn extract_gene(line: &[u8]) -> Result<Option<&[u8]>, &'static str> {
    // SAM has eleven mandatory fields. Restrict XT matching to optional fields so an XT:Z:
    // substring in a read name, sequence, or another tag value cannot be mistaken for the tag.
    let mut fields = line.split(|&b| b == b'\t');
    for _ in 0..11 {
        fields
            .next()
            .ok_or("SAM record has fewer than eleven mandatory fields")?;
    }

    let mut gene = None;
    for field in fields {
        // process_line receives records without LF, but tolerate CRLF input on the final field.
        let field = field.strip_suffix(b"\r").unwrap_or(field);
        if !field.starts_with(b"XT:") {
            continue;
        }
        if !field.starts_with(b"XT:Z:") {
            return Err("XT optional field is not a string (expected XT:Z:)");
        }
        if gene.is_some() {
            return Err("SAM record contains more than one XT tag");
        }

        let value = &field[5..];
        if value.is_empty() {
            return Err("XT:Z: gene identifier is empty");
        }
        let value = normalize_gene_id(value);
        if value.is_empty() {
            return Err("XT:Z: gene identifier is empty after version normalization");
        }
        gene = Some(value);
    }

    Ok(gene)
}

fn process_line(line: &[u8], out: &mut Vec<u8>) -> Result<(), &'static str> {
    let Some(gene) = extract_gene(line)? else {
        return Ok(());
    };

    // field 1 = bytes up to first tab (or whole line if no tab)
    let f1_end = line.iter().position(|&b| b == b'\t').unwrap_or(line.len());
    let f1 = &line[..f1_end];

    // barcode: match _[A-Z]+_$ anchored at end of field 1
    // require: f1 ends with '_', preceded by one-or-more [A-Z], preceded by '_'
    let mut barcode: &[u8] = b"";
    if let Some(&last) = f1.last() {
        if last == b'_' {
            let j = f1.len() - 1; // index of trailing '_'
            let mut k = j;
            while k > 0 && is_upper(f1[k - 1]) {
                k -= 1;
            }
            // now f1[k..j] is the [A-Z]+ run (possibly empty); need >=1 and a '_' before it
            if k < j && k > 0 && f1[k - 1] == b'_' {
                barcode = &f1[k..j];
            }
        }
    }

    out.extend_from_slice(barcode);
    out.push(b'\t');
    out.extend_from_slice(gene);
    out.push(b'\n');
    Ok(())
}

fn process_record(line: &[u8], line_number: u64, out: &mut Vec<u8>) -> io::Result<()> {
    process_line(line, out).map_err(|message| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid SAM record at line {line_number}: {message}"),
        )
    })
}

fn main() -> io::Result<()> {
    let mut reader = io::stdin().lock();
    let stdout = io::stdout();
    let mut writer = io::BufWriter::with_capacity(1 << 20, stdout.lock());

    let mut inbuf = vec![0u8; 1 << 20];
    let mut carry: Vec<u8> = Vec::with_capacity(4096);
    let mut out: Vec<u8> = Vec::with_capacity(1 << 21);
    let mut line_number = 0u64;

    loop {
        let n = reader.read(&mut inbuf)?;
        if n == 0 {
            break;
        }
        let mut chunk = &inbuf[..n];
        // if we have a carry, the first line is carry + up-to-first-newline
        if !carry.is_empty() {
            if let Some(nl) = chunk.iter().position(|&b| b == b'\n') {
                carry.extend_from_slice(&chunk[..nl]);
                line_number += 1;
                process_record(&carry, line_number, &mut out)?;
                carry.clear();
                chunk = &chunk[nl + 1..];
            } else {
                carry.extend_from_slice(chunk);
                continue;
            }
        }
        // process whole lines within chunk
        let mut start = 0;
        while let Some(rel) = chunk[start..].iter().position(|&b| b == b'\n') {
            let line = &chunk[start..start + rel];
            line_number += 1;
            process_record(line, line_number, &mut out)?;
            start += rel + 1;
        }
        // leftover partial line -> carry
        if start < chunk.len() {
            carry.extend_from_slice(&chunk[start..]);
        }
        if out.len() > (1 << 20) {
            writer.write_all(&out)?;
            out.clear();
        }
    }
    // final unterminated line (awk processes a final record without trailing newline)
    if !carry.is_empty() {
        line_number += 1;
        process_record(&carry, line_number, &mut out)?;
    }
    writer.write_all(&out)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sam_record(qname: &str, optional_fields: &[&str]) -> Vec<u8> {
        let mut record = format!("{qname}\t0\tchr1\t1\t255\t1M\t*\t0\t0\tA\tI");
        for field in optional_fields {
            record.push('\t');
            record.push_str(field);
        }
        record.into_bytes()
    }

    fn output_for_gene(gene: &str) -> Vec<u8> {
        let record = sam_record("read_ACGTACGTACGTA_", &[&format!("XT:Z:{gene}")]);
        let mut output = Vec::new();
        process_line(&record, &mut output).unwrap();
        output
    }

    #[test]
    fn preserves_complete_punctuated_and_unicode_gene_ids() {
        for gene in [
            "gene-123",
            "gene:alpha",
            "gene.with.words",
            "gene_under_score",
            "12345",
            "gène-δ",
        ] {
            assert_eq!(
                output_for_gene(gene),
                format!("ACGTACGTACGTA\t{gene}\n").as_bytes()
            );
        }
    }

    #[test]
    fn strips_only_a_terminal_numeric_version_suffix() {
        for (input, expected) in [
            ("ENSG000001.42", "ENSG000001"),
            ("gene.0", "gene"),
            ("gene.12a", "gene.12a"),
            ("gene.12.more", "gene.12.more"),
            ("gene.", "gene."),
            ("gene-12", "gene-12"),
        ] {
            assert_eq!(
                output_for_gene(input),
                format!("ACGTACGTACGTA\t{expected}\n").as_bytes()
            );
        }
    }

    #[test]
    fn reads_xt_only_as_a_complete_optional_field() {
        let record = sam_record("read_ACGTACGTACGTA_", &["XX:Z:XT:Z:not-a-tag", "NM:i:0"]);
        let mut output = Vec::new();
        process_line(&record, &mut output).unwrap();
        assert!(output.is_empty());

        let record = sam_record(
            "read_ACGTACGTACGTA_",
            &["XX:Z:value", "XT:Z:gene:with:colons", "NM:i:0"],
        );
        process_line(&record, &mut output).unwrap();
        assert_eq!(output, b"ACGTACGTACGTA\tgene:with:colons\n");
    }

    #[test]
    fn accepts_crlf_and_an_unterminated_final_record() {
        let mut record = sam_record("read_ACGTACGTACGTA_", &["XT:Z:gene-with-cr"]);
        record.push(b'\r');
        let mut output = Vec::new();
        process_line(&record, &mut output).unwrap();
        assert_eq!(output, b"ACGTACGTACGTA\tgene-with-cr\n");
    }

    #[test]
    fn skips_a_well_formed_record_without_an_xt_tag() {
        let record = sam_record("read_ACGTACGTACGTA_", &["NM:i:0"]);
        let mut output = Vec::new();
        process_line(&record, &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn rejects_malformed_xt_tags_and_records() {
        for (fields, expected_error) in [
            (vec!["XT:i:1"], "XT optional field is not a string"),
            (vec!["XT:Z:"], "gene identifier is empty"),
            (vec!["XT:Z:gene-a", "XT:Z:gene-b"], "more than one XT tag"),
            (vec!["XT:Z:.123"], "empty after version normalization"),
        ] {
            let record = sam_record("read_ACGTACGTACGTA_", &fields);
            let mut output = Vec::new();
            let error = process_line(&record, &mut output).unwrap_err();
            assert!(error.contains(expected_error), "unexpected error: {error}");
            assert!(output.is_empty());
        }

        let mut output = Vec::new();
        let error = process_line(b"not\ta\tcomplete\tsam", &mut output).unwrap_err();
        assert!(error.contains("fewer than eleven"));
    }
}
