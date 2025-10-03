#!/usr/bin/env python
"""Register unified_qc module as a MultiQC entry point"""

from pathlib import Path
import re

# Direct path to the entry_points.txt file
ep_file = Path('/usr/local/lib/python3.11/site-packages/multiqc-1.14.dist-info/entry_points.txt')

# Read current content
content = ep_file.read_text()

# Check if already registered
if 'unified_qc' not in content:
    # Add unified_qc entry point to the multiqc.modules.v1 section
    # Find the section and add our module
    content = re.sub(
        r'(\[multiqc\.modules\.v1\].*?)(\n\[)',
        r'\1unified_qc = multiqc.modules.unified_qc:MultiqcModule\n\2',
        content,
        flags=re.DOTALL
    )

    # Write back
    ep_file.write_text(content)
    print('✓ Added unified_qc entry point')
else:
    print('✓ unified_qc already registered')
