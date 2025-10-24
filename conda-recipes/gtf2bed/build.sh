#!/bin/bash

set -ex

# Package pre-built static binary
# The binary should already be built using: cargo build --release --target x86_64-unknown-linux-musl
# Source path points to the release directory containing the 'gtf2bed' binary

# Install the static binary to conda prefix
mkdir -p $PREFIX/bin
cp gtf2bed $PREFIX/bin/gtf2bed
chmod +x $PREFIX/bin/gtf2bed

# Verify the binary exists and is executable
test -x $PREFIX/bin/gtf2bed
echo "gtf2bed static binary successfully installed to $PREFIX/bin/gtf2bed"

# Verify it's actually statically linked
echo "Checking binary dependencies:"
if ldd $PREFIX/bin/gtf2bed 2>&1 | grep -q "not a dynamic executable"; then
    echo "SUCCESS: Static binary confirmed (not a dynamic executable)"
elif ldd $PREFIX/bin/gtf2bed 2>&1 | grep -q "statically linked"; then
    echo "SUCCESS: Static binary confirmed (statically linked)"
else
    echo "WARNING: Binary appears to be dynamically linked:"
    ldd $PREFIX/bin/gtf2bed || true
    exit 1
fi
