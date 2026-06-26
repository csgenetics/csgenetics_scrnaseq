// io_count_extract: byte-exact drop-in for the io_count awk text stage.
//
// Reads SAM records (tab-delimited, no header) on stdin, writes one line per record that
// contains the substring "XT:", exactly replicating:
//
//   awk '/XT:/ {match($1, /_[A-Z]+_$/); printf substr($0,RSTART+1,RLENGTH-2);
//               match($0, /XT:Z:[A-Za-z0-9_]+/); print "\t" substr($0,RSTART+5,RLENGTH-5)}'
//
// Output per matching record: <barcode>\t<gene>\n
//   barcode = the [A-Z]+ token between the trailing "_..._" of field 1 (read name), else empty
//   gene    = the [A-Za-z0-9_]+ run after the first "XT:Z:" in the whole line, else empty
//
// Streaming and low-memory (the process is capped at 1 GB): we never hold the whole input.

use std::io::{self, Read, Write};

#[inline]
fn is_upper(b: u8) -> bool {
    b.is_ascii_uppercase()
}

#[inline]
fn is_gene_char(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

// Find substring `needle` in `hay`, return index of first byte after a match start.
#[inline]
fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > hay.len() {
        return None;
    }
    let first = needle[0];
    let mut i = 0;
    let end = hay.len() - needle.len();
    while i <= end {
        if hay[i] == first && &hay[i..i + needle.len()] == needle {
            return Some(i);
        }
        i += 1;
    }
    None
}

fn process_line(line: &[u8], out: &mut Vec<u8>) {
    // /XT:/  -- only emit for records whose raw line contains "XT:"
    if find(line, b"XT:").is_none() {
        return;
    }

    // field 1 = bytes up to first tab (or whole line if no tab)
    let f1_end = line.iter().position(|&b| b == b'\t').unwrap_or(line.len());
    let f1 = &line[..f1_end];

    // barcode: match _[A-Z]+_$ anchored at end of field 1
    // require: f1 ends with '_', preceded by one-or-more [A-Z], preceded by '_'
    let mut barcode: &[u8] = b"";
    if let Some(&last) = f1.last() {
        if last == b'_' {
            let mut j = f1.len() - 1; // index of trailing '_'
            let mut k = j;
            while k > 0 && is_upper(f1[k - 1]) {
                k -= 1;
            }
            // now f1[k..j] is the [A-Z]+ run (possibly empty); need >=1 and a '_' before it
            if k < j && k > 0 && f1[k - 1] == b'_' {
                barcode = &f1[k..j];
            }
            let _ = &mut j;
        }
    }

    // gene: first "XT:Z:" then [A-Za-z0-9_]+ run
    let mut gene: &[u8] = b"";
    if let Some(p) = find(line, b"XT:Z:") {
        let start = p + 5;
        let mut e = start;
        while e < line.len() && is_gene_char(line[e]) {
            e += 1;
        }
        gene = &line[start..e];
    }

    out.extend_from_slice(barcode);
    out.push(b'\t');
    out.extend_from_slice(gene);
    out.push(b'\n');
}

fn main() -> io::Result<()> {
    let mut reader = io::stdin().lock();
    let stdout = io::stdout();
    let mut writer = io::BufWriter::with_capacity(1 << 20, stdout.lock());

    let mut inbuf = vec![0u8; 1 << 20];
    let mut carry: Vec<u8> = Vec::with_capacity(4096);
    let mut out: Vec<u8> = Vec::with_capacity(1 << 21);

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
                process_line(&carry, &mut out);
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
            process_line(line, &mut out);
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
        process_line(&carry, &mut out);
    }
    writer.write_all(&out)?;
    writer.flush()?;
    Ok(())
}
