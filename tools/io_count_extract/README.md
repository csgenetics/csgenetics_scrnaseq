# io_count_extract

A small static binary that replaces the slow `awk` text stage of the `io_count` process.

## What it does

`io_count` turns a deduplicated BAM into `<barcode>\t<gene>` lines for every alignment carrying
a featureCounts gene assignment (`XT:Z:` tag). The original implementation piped `samtools view`
into an `awk` one-liner. The production container (`quay.io/biocontainers/samtools:1.17`) ships
**BusyBox awk**, which is very slow; on a 1.4 GB mixed-species dedup BAM that stage took ~409 s.

`io_count_extract` reads the same `samtools view` text stream on stdin and performs the transform
in compiled code. It is streaming and low-memory (the process is capped at 1 GB RAM). On the
1.4 GB BAM the stage drops from ~409 s to ~22 s.

It originally replaced this command byte-for-byte:
```
awk '/XT:/ {match($1, /_[A-Z]+_$/); printf substr($0,RSTART+1,RLENGTH-2);
            match($0, /XT:Z:[A-Za-z0-9_]+/); print "\t" substr($0,RSTART+5,RLENGTH-5)}'
```

The awk character class silently truncated valid custom-reference identifiers at punctuation.
The extractor now parses `XT:Z:` as a complete SAM optional field and removes only a terminal
numeric version suffix such as `.12`, matching `bin/features_names.py`. Hyphens, colons, other
dots, underscores, spaces, and UTF-8 bytes are preserved. A malformed, empty, non-string, or
duplicate `XT` field fails loudly; a record without `XT` is skipped. Distinct GTF identifiers that
would collide after version normalization are rejected by `features_names.py`.

## Why a committed binary

The binary is statically linked (musl), so it runs inside the existing public samtools container
with no image change. The pipeline's per-process containers are public biocontainers that Wave can
pull anonymously; we cannot push a custom image (private quay is not Wave-pullable), so shipping a
static binary on the pipeline `bin/` PATH is the cleanest way to add a compiled tool. Deployment
targets are x86_64 (AWS Batch / Seqera), matching the binary target.

## Rebuild

```
cd tools/io_count_extract
cargo build --locked --release --target x86_64-unknown-linux-musl
strip target/x86_64-unknown-linux-musl/release/io_count_extract
cp target/x86_64-unknown-linux-musl/release/io_count_extract ../../bin/io_count_extract
```

The source has no third-party dependencies and the lockfile is committed, so the build is
hermetic given a Rust toolchain with the `x86_64-unknown-linux-musl` target installed
(`rustup target add x86_64-unknown-linux-musl`).

CircleCI runs the Rust unit tests and rebuilds with the reviewed Rust 1.93.0 toolchain. After
stripping, the result must match `bin/io_count_extract` byte-for-byte, so source and the
production artifact cannot drift independently.
