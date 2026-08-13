"""
Shared pytest configuration for the bin/ Python unit tests.

Adds the repository ``bin/`` directory to ``sys.path`` so the pipeline scripts
can be imported directly as modules (``import categorize_reads`` etc.). The
scripts are designed to run as standalone executables staged into a Nextflow
work dir, so they are not an installable package; pointing ``sys.path`` at
``bin/`` is the cleanest way to exercise their importable functions.
"""

import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")

if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)
