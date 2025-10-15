#!/bin/bash

set -ex

# Use cargo install which properly handles the binary installation
# Install directly to the conda prefix
cargo install --path . --root $PREFIX

# Verify the binary exists and is executable
test -x $PREFIX/bin/qc
echo "qc binary successfully installed to $PREFIX/bin/qc"
