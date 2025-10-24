#!/bin/bash

set -ex

# Install the Python package
python -m pip install . --no-deps --no-build-isolation -vv
