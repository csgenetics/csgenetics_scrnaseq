#!/bin/bash

set -ex

# Package pre-built static binary
# The binary should already be built using: cargo build --release --target x86_64-unknown-linux-musl
# Source path points to the release directory containing the 'qc' binary

# Install the static binary to conda prefix
mkdir -p $PREFIX/bin
cp qc $PREFIX/bin/qc
chmod +x $PREFIX/bin/qc

# Verify the binary exists and is executable
test -x $PREFIX/bin/qc
echo "qc static binary successfully installed to $PREFIX/bin/qc"

# Verify it's actually statically linked
echo "Checking binary dependencies:"
if ldd $PREFIX/bin/qc 2>&1 | grep -q "not a dynamic executable"; then
    echo "SUCCESS: Static binary confirmed (not a dynamic executable)"
elif ldd $PREFIX/bin/qc 2>&1 | grep -q "statically linked"; then
    echo "SUCCESS: Static binary confirmed (statically linked)"
else
    echo "WARNING: Binary appears to be dynamically linked:"
    ldd $PREFIX/bin/qc || true
    exit 1
fi
