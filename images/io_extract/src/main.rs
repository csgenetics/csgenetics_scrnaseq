// Simplified Rust io_extract without external dependencies for arg parsing
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

const BARCODE_LENGTH: usize = 13;
const BUFFER_SIZE: usize = 8 * 1024 * 1024;

struct BarcodeChecker {
    whitelist: HashSet<String>,
    corrections: HashMap<String, String>,
}

impl BarcodeChecker {
    fn load(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut whitelist = HashSet::new();
        let mut corrections = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();

            if parts.is_empty() {
                continue;
            }

            let barcode = parts[0].to_string();
            whitelist.insert(barcode.clone());

            if parts.len() >= 2 {
                for variant in parts[1].split(',') {
                    if !variant.is_empty() && variant != barcode {
                        corrections.insert(variant.to_string(), barcode.clone());
                    }
                }
            }
        }

        eprintln!("Loaded {} barcodes", whitelist.len());
        eprintln!("Loaded {} corrections", corrections.len());

        Ok(Self {
            whitelist,
            corrections,
        })
    }

    fn check(&self, barcode: &str) -> Option<String> {
        if self.whitelist.contains(barcode) {
            Some(barcode.to_string())
        } else {
            self.corrections.get(barcode).cloned()
        }
    }
}

#[derive(Default)]
struct Stats {
    total: usize,
    passed: usize,
    failed: usize,
    direct_match: usize,
    corrected: usize,
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        eprintln!("Usage: {} <r1.fastq.gz> <r2.fastq.gz> <whitelist.tsv> [--sample-id NAME] [--output-dir DIR] [--no-filtered]", args[0]);
        std::process::exit(1);
    }

    let r1_path = &args[1];
    let r2_path = &args[2];
    let whitelist_path = &args[3];

    let mut sample_id = "rust_test".to_string();
    let mut output_dir = ".".to_string();
    let mut no_filtered = false;

    let mut i = 4;
    while i < args.len() {
        match args[i].as_str() {
            "--sample-id" if i + 1 < args.len() => {
                sample_id = args[i + 1].clone();
                i += 2;
            }
            "--output-dir" if i + 1 < args.len() => {
                output_dir = args[i + 1].clone();
                i += 2;
            }
            "--no-filtered" => {
                no_filtered = true;
                i += 1;
            }
            _ => i += 1,
        }
    }

    eprintln!("Rust io_extract - High-performance barcode extraction");
    eprintln!("Buffer size: {} MB", BUFFER_SIZE / (1024 * 1024));

    let checker = BarcodeChecker::load(whitelist_path)?;

    let start = Instant::now();
    let mut stats = Stats::default();

    // Open input files
    let r1_file = File::open(r1_path)?;
    let r2_file = File::open(r2_path)?;

    let r1_gz = MultiGzDecoder::new(BufReader::with_capacity(BUFFER_SIZE, r1_file));
    let r2_gz = MultiGzDecoder::new(BufReader::with_capacity(BUFFER_SIZE, r2_file));

    let mut r1_reader = BufReader::with_capacity(BUFFER_SIZE, r1_gz);
    let mut r2_reader = BufReader::with_capacity(BUFFER_SIZE, r2_gz);

    // Open output files
    let r1_out_path = format!("{}/{}.io_extract.R1.fastq.gz", output_dir, sample_id);
    let r1_out = GzEncoder::new(
        BufWriter::with_capacity(BUFFER_SIZE, File::create(&r1_out_path)?),
        Compression::fast(),
    );
    let mut r1_out = BufWriter::with_capacity(BUFFER_SIZE, r1_out);

    let mut filtered_r1 = if !no_filtered {
        let path = format!("{}/io_extract_{}_filteredOut_R1.fastq.gz", output_dir, sample_id);
        Some(BufWriter::with_capacity(
            BUFFER_SIZE,
            GzEncoder::new(File::create(path)?, Compression::fast()),
        ))
    } else {
        None
    };

    let mut filtered_r2 = if !no_filtered {
        let path = format!("{}/io_extract_{}_filteredOut_R2.fastq.gz", output_dir, sample_id);
        Some(BufWriter::with_capacity(
            BUFFER_SIZE,
            GzEncoder::new(File::create(path)?, Compression::fast()),
        ))
    } else {
        None
    };

    eprintln!("Processing reads...");
    eprintln!("  Filtered output: {}", if no_filtered { "DISABLED" } else { "ENABLED" });

    // Buffers for FASTQ lines
    let mut h1 = String::new();
    let mut s1 = String::new();
    let mut p1 = String::new();
    let mut q1 = String::new();
    let mut h2 = String::new();
    let mut s2 = String::new();
    let mut p2 = String::new();
    let mut q2 = String::new();
    let mut output_buffer = String::new();

    loop {
        // Read R1 record
        h1.clear();
        s1.clear();
        p1.clear();
        q1.clear();

        if r1_reader.read_line(&mut h1)? == 0 {
            break;
        }
        r1_reader.read_line(&mut s1)?;
        r1_reader.read_line(&mut p1)?;
        r1_reader.read_line(&mut q1)?;

        // Read R2 record
        h2.clear();
        s2.clear();
        p2.clear();
        q2.clear();

        r2_reader.read_line(&mut h2)?;
        r2_reader.read_line(&mut s2)?;
        r2_reader.read_line(&mut p2)?;
        r2_reader.read_line(&mut q2)?;

        stats.total += 1;

        // Extract barcode from R2
        let s2_trimmed = s2.trim_end();
        if s2_trimmed.len() < BARCODE_LENGTH {
            stats.failed += 1;
            if let (Some(f1), Some(f2)) = (&mut filtered_r1, &mut filtered_r2) {
                write!(f1, "{}{}{}{}", h1, s1, p1, q1)?;
                write!(f2, "{}{}{}{}", h2, s2, p2, q2)?;
            }
            continue;
        }

        let barcode = &s2_trimmed[..BARCODE_LENGTH];

        match checker.check(barcode) {
            Some(final_barcode) => {
                // Passed - write modified R1 with barcode in header
                output_buffer.clear();

                // Parse header (skip '@')
                let header = h1[1..].trim_end();

                // Split header on first space
                if let Some(space_pos) = header.find(' ') {
                    output_buffer.push('@');
                    output_buffer.push_str(&header[..space_pos]);
                    output_buffer.push('_');
    output_buffer.push_str(&final_barcode);
                    output_buffer.push_str("_ ");
                    output_buffer.push_str(&header[space_pos + 1..]);
                } else {
                    output_buffer.push('@');
                    output_buffer.push_str(header);
                    output_buffer.push('_');
                    output_buffer.push_str(&final_barcode);
                    output_buffer.push('_');
                }
                output_buffer.push('\n');
                output_buffer.push_str(&s1);
                output_buffer.push_str(&p1);
                output_buffer.push_str(&q1);

                r1_out.write_all(output_buffer.as_bytes())?;

                stats.passed += 1;
                if barcode == final_barcode {
                    stats.direct_match += 1;
                } else {
                    stats.corrected += 1;
                }
            }
            None => {
                // Failed
                stats.failed += 1;
                if let (Some(f1), Some(f2)) = (&mut filtered_r1, &mut filtered_r2) {
                    write!(f1, "{}{}{}{}", h1, s1, p1, q1)?;
                    write!(f2, "{}{}{}{}", h2, s2, p2, q2)?;
                }
            }
        }

        // Progress update every 1M reads
        if stats.total % 1_000_000 == 0 {
            let elapsed = start.elapsed().as_secs_f64();
            let rate = stats.total as f64 / elapsed;
            eprintln!("  Processed {} reads ({:.0} reads/sec)", stats.total, rate);
        }
    }

    // Flush and close
    r1_out.flush()?;
    if let Some(f1) = &mut filtered_r1 {
        f1.flush()?;
    }
    if let Some(f2) = &mut filtered_r2 {
        f2.flush()?;
    }

    let duration = start.elapsed().as_secs_f64();

    // Print summary
    eprintln!("\nRust io_extract Complete");
    eprintln!("==================================================");
    eprintln!("Duration: {:.2} seconds", duration);
    eprintln!("Total reads: {}", stats.total);
    eprintln!("Reads passed: {}", stats.passed);
    eprintln!("Reads failed: {}", stats.failed);
    eprintln!("Direct matches: {}", stats.direct_match);
    eprintln!("Corrected: {}", stats.corrected);
    eprintln!("Throughput: {:.0} reads/second", stats.total as f64 / duration);

    // Write log file
    let log_path = format!("{}/{}.io_extract.log", output_dir, sample_id);
    let mut log = File::create(log_path)?;
    writeln!(log, "# Rust io_extract")?;
    writeln!(log, "# Duration: {:.2} seconds", duration)?;
    writeln!(log, "# Throughput: {:.0} reads/second", stats.total as f64 / duration)?;
    writeln!(log)?;
    writeln!(log, "Input Reads: {}", stats.total)?;
    writeln!(log, "Reads output: {}", stats.passed)?;
    writeln!(log, "Filtered cell barcode. Not correctable: {}", stats.failed)?;
    writeln!(log, "False cell barcode. Error-corrected: {}", stats.corrected)?;
    writeln!(log)?;
    writeln!(log, "# Additional Statistics:")?;
    writeln!(log, "# Direct matches: {}", stats.direct_match)?;
    writeln!(log, "# Filtered output: {}", if no_filtered { "No" } else { "Yes" })?;

    Ok(())
}