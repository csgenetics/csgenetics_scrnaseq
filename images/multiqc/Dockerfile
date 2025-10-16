FROM quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0

USER root

# Copy our custom unified_qc module into the MultiQC modules directory
COPY unified_qc /usr/local/lib/python3.11/site-packages/multiqc/modules/unified_qc/

# Add unified_qc pattern to search_patterns.yaml (v1.14 uses utils/search_patterns.yaml)
RUN echo "\nunified_qc:\n  fn_re: '.*\\.(R[12]\\.preQC|R1\\.postQC)\\.fastp\\.json\$'" >> /usr/local/lib/python3.11/site-packages/multiqc/utils/search_patterns.yaml

# Register the unified_qc module as an entry point
COPY register_module.py /tmp/register_module.py
RUN python /tmp/register_module.py && rm /tmp/register_module.py

# Clear Python cache to force reload of entry points
RUN rm -rf /usr/local/lib/python3.11/site-packages/multiqc/__pycache__ && \
    rm -rf /usr/local/lib/python3.11/site-packages/multiqc/modules/__pycache__

# Switch back to default user
USER 1000
