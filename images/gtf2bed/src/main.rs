use flate2::read::GzDecoder;
use regex::Regex;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::process;

#[derive(Debug, Clone)]
struct Exon {
    chr: String,
    start: u64,
    end: u64,
    strand: char,
    attr: String,
}

struct Transcript {
    exons: Vec<Exon>,
    start_codon: Option<u64>,
    stop_codon: Option<u64>,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: gtf2bed <input.gtf>");
        process::exit(1);
    }

    let filename = &args[1];

    // Open file with potential gzip decompression
    let reader: Box<dyn Read> = if filename.ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(filename).unwrap_or_else(|e| {
            eprintln!("Can't open {}: {}", filename, e);
            process::exit(1);
        })))
    } else {
        Box::new(File::open(filename).unwrap_or_else(|e| {
            eprintln!("Can't open {}: {}", filename, e);
            process::exit(1);
        }))
    };

    let buf_reader = BufReader::new(reader);

    // Regex for extracting transcript_id
    let transcript_re = Regex::new(r#"transcript_id "([^"]+)""#).unwrap();

    // HashMaps to store transcripts and their data
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();
    let mut first_exon: HashMap<String, Exon> = HashMap::new();

    // Parse GTF file
    for line in buf_reader.lines() {
        let line = line.unwrap();

        // Skip comments
        if line.starts_with('#') {
            continue;
        }

        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Extract transcript_id
        let transcript_id = match transcript_re.captures(fields[8]) {
            Some(caps) => caps.get(1).unwrap().as_str().to_string(),
            None => continue,
        };

        if fields[0].is_empty() {
            continue;
        }

        let feature_type = fields[2];

        if feature_type == "exon" {
            let start: u64 = fields[3].parse().unwrap_or_else(|_| {
                eprintln!("no position at exon");
                process::exit(1);
            });
            let end: u64 = fields[4].parse().unwrap();
            let strand = fields[6].chars().next().unwrap();

            let exon = Exon {
                chr: fields[0].to_string(),
                start,
                end,
                strand,
                attr: fields[8].to_string(),
            };

            // Add exon to transcript
            transcripts
                .entry(transcript_id.clone())
                .or_insert_with(|| Transcript {
                    exons: Vec::new(),
                    start_codon: None,
                    stop_codon: None,
                })
                .exons
                .push(exon.clone());

            // Save first exon for transcript if not already saved
            first_exon.entry(transcript_id.clone()).or_insert(exon);

        } else if feature_type == "start_codon" {
            let start: u64 = fields[3].parse().unwrap();
            transcripts
                .entry(transcript_id.clone())
                .or_insert_with(|| Transcript {
                    exons: Vec::new(),
                    start_codon: None,
                    stop_codon: None,
                })
                .start_codon = Some(start);

        } else if feature_type == "stop_codon" {
            let end: u64 = fields[4].parse().unwrap();
            transcripts
                .entry(transcript_id.clone())
                .or_insert_with(|| Transcript {
                    exons: Vec::new(),
                    start_codon: None,
                    stop_codon: None,
                })
                .stop_codon = Some(end);

        } else if feature_type == "miRNA" {
            let start: u64 = fields[3].parse().unwrap();
            let end: u64 = fields[4].parse().unwrap();
            let strand = fields[6].chars().next().unwrap();

            let exon = Exon {
                chr: fields[0].to_string(),
                start,
                end,
                strand,
                attr: fields[8].to_string(),
            };

            transcripts
                .entry(transcript_id.clone())
                .or_insert_with(|| Transcript {
                    exons: Vec::new(),
                    start_codon: None,
                    stop_codon: None,
                })
                .exons
                .push(exon.clone());

            first_exon.entry(transcript_id.clone()).or_insert(exon);
        }
    }

    // Sort transcript IDs by chromosome, position, then transcript ID
    let mut sorted_ids: Vec<String> = first_exon.keys().cloned().collect();
    sorted_ids.sort_by(|a, b| {
        let exon_a = &first_exon[a];
        let exon_b = &first_exon[b];

        if exon_a.chr == exon_b.chr {
            match exon_a.start.cmp(&exon_b.start) {
                std::cmp::Ordering::Equal => {
                    // If same chr and start, sort by transcript ID
                    a.cmp(b)
                }
                other => other,
            }
        } else {
            exon_a.chr.cmp(&exon_b.chr)
        }
    });

    // Output BED12 format
    for transcript_id in sorted_ids {
        let transcript = &transcripts[&transcript_id];
        if transcript.exons.is_empty() {
            continue;
        }

        let first = &first_exon[&transcript_id];
        let chr = &first.chr;
        let strand = first.strand;

        // Sort exons by start position
        let mut exons = transcript.exons.clone();
        exons.sort_by_key(|e| e.start);

        let beg = exons[0].start;
        let end = exons[exons.len() - 1].end;

        // Handle CDS start/end
        let mut cds_start = transcript.start_codon;
        let mut cds_end = transcript.stop_codon;

        if strand == '-' {
            std::mem::swap(&mut cds_start, &mut cds_end);
            if let Some(ref mut cs) = cds_start {
                *cs = cs.saturating_sub(2);
            }
            if let Some(ref mut ce) = cds_end {
                *ce += 2;
            }
        }

        // Default to exon boundaries if not specified
        let cds_start = cds_start.unwrap_or(beg);
        let cds_end = cds_end.unwrap_or(end);

        // Adjust for BED format (0-based)
        let beg = beg - 1;
        let cds_start = cds_start - 1;

        // Calculate exon count, sizes, and starts
        let exon_count = exons.len();
        let exon_sizes: Vec<String> = exons.iter().map(|e| (e.end - e.start + 1).to_string()).collect();
        let exon_starts: Vec<String> = exons.iter().map(|e| (e.start - beg - 1).to_string()).collect();

        // Output BED12 line
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{},\t{},",
            chr,
            beg,
            end,
            transcript_id,
            0,
            strand,
            cds_start,
            cds_end,
            0,
            exon_count,
            exon_sizes.join(","),
            exon_starts.join(",")
        );
    }
}
