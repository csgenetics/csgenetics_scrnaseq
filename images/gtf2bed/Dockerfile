FROM ubuntu:22.04

# Install runtime dependencies including procps for ps command (required by Nextflow)
RUN apt-get update && apt-get install -y \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled static Rust binary
COPY target/x86_64-unknown-linux-musl/release/gtf2bed /usr/local/bin/gtf2bed

# Make it executable
RUN chmod +x /usr/local/bin/gtf2bed

# No ENTRYPOINT - Nextflow needs to run bash scripts, not just the binary
# The gtf2bed binary is in PATH and can be called directly from Nextflow processes
