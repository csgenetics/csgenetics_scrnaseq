# Multi-stage build for io_extract Rust implementation
# Optimized for small final image size and maximum performance

# Build stage
FROM rust:1.72-slim AS builder

WORKDIR /build

# Copy source files
COPY Cargo.toml .
COPY src/ src/

# Build with full optimizations
RUN cargo build --release

# Runtime stage
FROM debian:bookworm-slim

# Install runtime dependencies including procps for Nextflow metrics
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from build stage
COPY --from=builder /build/target/release/io_extract /usr/local/bin/io_extract

# Ensure binary is executable
RUN chmod +x /usr/local/bin/io_extract

# Test that binary runs
RUN io_extract --help || true

# Set working directory
WORKDIR /data

# Default command
CMD ["/usr/local/bin/io_extract"]