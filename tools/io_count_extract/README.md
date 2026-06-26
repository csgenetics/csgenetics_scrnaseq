# io_count_extract

A small static binary that replaces the slow `awk` text stage of the `io_count` process.

## What it does

`io_count` turns a deduplicated BAM into `<barcode>\t<gene>` lines for every alignment carrying
a featureCounts gene assignment (`XT:Z:` tag). The original implementation piped `samtools view`
into an `awk` one-liner. The production container (`quay.io/biocontainers/samtools:1.17`) ships
**BusyBox awk**, which is very slow; on a 1.4 GB mixed-species dedup BAM that stage took ~409 s.

`io_count_extract` reads the same `samtools view` text stream on stdin and does the identical
transform in compiled code, **byte-for-byte** with the awk (verified on human, empty/intergenic
edge cases, and a 1.4 GB barnyard BAM). It is streaming and low-memory (the process is capped at
1 GB RAM). On the 1.4 GB BAM the stage drops from ~409 s to ~22 s.

Replaces exactly:
```
awk '/XT:/ {match($1, /_[A-Z]+_$/); printf substr($0,RSTART+1,RLENGTH-2);
            match($0, /XT:Z:[A-Za-z0-9_]+/); print "\t" substr($0,RSTART+5,RLENGTH-5)}'
```

## Why a committed binary

The binary is statically linked (musl), so it runs inside the existing public samtools container
with no image change. The pipeline's per-process containers are public biocontainers that Wave can
pull anonymously; we cannot push a custom image (private quay is not Wave-pullable), so shipping a
static binary on the pipeline `bin/` PATH is the cleanest way to add a compiled tool. Deployment
targets are x86_64 (AWS Batch / Seqera), matching the binary target.

## Rebuild

```
cd tools/io_count_extract
cargo build --release --target x86_64-unknown-linux-musl
strip target/x86_64-unknown-linux-musl/release/io_count_extract
cp target/x86_64-unknown-linux-musl/release/io_count_extract ../../bin/io_count_extract
```

The source has no third-party dependencies, so the build is hermetic given a Rust toolchain with
the `x86_64-unknown-linux-musl` target installed (`rustup target add x86_64-unknown-linux-musl`).
