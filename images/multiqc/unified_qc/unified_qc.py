import json
import logging
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """MultiQC module for CS Genetics Unified QC binary"""

    def __init__(self):
        log.info("unified_qc module __init__ called")
        super(MultiqcModule, self).__init__(
            name="Unified QC",
            anchor="unified_qc",
            info="CS Genetics unified QC process (barcode extraction + trimming + filtering)"
        )

        log.info("unified_qc module initialized, starting file search")
        # Parse and group the 3 JSON files per sample
        data_by_sample = self.parse_unified_qc_files()

        if len(data_by_sample) == 0:
            log.debug("No Unified QC files found")
            return

        log.info(f"Found {len(data_by_sample)} Unified QC sample sets")

        # Add to general stats table
        self.unified_qc_general_stats_table(data_by_sample)

        # Add sequence quality plot with Pre-QC / Post-QC toggle
        self.add_section(
            name="Sequence Quality",
            anchor="unified-qc-seq-quality",
            description="Average sequencing quality over each base (Pre-QC vs Post-QC)",
            plot=self.unified_qc_quality_plot(data_by_sample)
        )

        # Add GC content plot with Pre-QC / Post-QC toggle
        self.add_section(
            name="GC Content",
            anchor="unified-qc-gc-content",
            description="Average GC content over each base (Pre-QC vs Post-QC)",
            plot=self.unified_qc_gc_plot(data_by_sample)
        )

        # Add N content plot with Pre-QC / Post-QC toggle
        self.add_section(
            name="N Content",
            anchor="unified-qc-n-content",
            description="Average N content over each base (Pre-QC vs Post-QC)",
            plot=self.unified_qc_n_plot(data_by_sample)
        )

    def parse_unified_qc_files(self):
        """Find and group the 3 JSON files per sample into a nested dict.

        Returns:
            dict: {
                'sample_name': {
                    'r1_preqc': {parsed JSON},
                    'r2_preqc': {parsed JSON},
                    'r1_postqc': {parsed JSON}
                }
            }
        """
        data_by_sample = {}

        # Find all .fastp.json files
        log.debug("Searching for unified_qc files")
        # Important: Don't convert to list - iterate directly so file handles stay open
        for f in self.find_log_files("unified_qc", filehandles=True):
            try:
                parsed_json = json.load(f['f'])

                # Determine file type and sample name from filename
                # Files: {sample}.R1.preQC.fastp.json, {sample}.R2.preQC.fastp.json, {sample}.R1.postQC.fastp.json
                if '.R1.preQC.fastp.json' in f['fn']:
                    file_type = 'r1_preqc'
                    s_name = f['fn'].replace('.R1.preQC.fastp.json', '')
                elif '.R2.preQC.fastp.json' in f['fn']:
                    file_type = 'r2_preqc'
                    s_name = f['fn'].replace('.R2.preQC.fastp.json', '')
                elif '.R1.postQC.fastp.json' in f['fn']:
                    file_type = 'r1_postqc'
                    s_name = f['fn'].replace('.R1.postQC.fastp.json', '')
                else:
                    # Not a recognized file pattern
                    continue

                # Initialize sample dict if needed
                if s_name not in data_by_sample:
                    data_by_sample[s_name] = {}

                # Store parsed JSON by file type
                data_by_sample[s_name][file_type] = parsed_json
                self.add_data_source(f, s_name)

            except json.JSONDecodeError:
                log.warning(f"Could not parse JSON: {f['fn']}")
                continue

        log.debug(f"Parsed {len(data_by_sample)} sample sets from unified_qc files")

        # Filter ignored samples
        data_by_sample = self.ignore_samples(data_by_sample)

        return data_by_sample

    def unified_qc_general_stats_table(self, data_by_sample):
        """Add unified QC stats to general stats table - 3 rows per sample"""
        headers = {
            "total_reads": {
                "title": "M Reads",
                "description": "Total reads (millions)",
                "modify": lambda x: x / 1000000.0,
                "scale": "GnBu",
                "shared_key": "read_count",
                "format": "{:,.2f}"
            },
            "q30_rate": {
                "title": "% > Q30",
                "description": "Percentage of bases > Q30",
                "modify": lambda x: x * 100.0,
                "scale": "GnBu",
                "suffix": "%",
                "format": "{:,.1f}"
            },
            "q30_bases": {
                "title": "Mb Q30",
                "description": "Q30 bases (millions)",
                "modify": lambda x: x / 1000000.0,
                "scale": "GnBu",
                "format": "{:,.1f}",
                "hidden": True
            },
            "gc_content": {
                "title": "GC%",
                "description": "GC content percentage",
                "modify": lambda x: x * 100.0,
                "scale": "GnBu",
                "suffix": "%",
                "format": "{:,.1f}"
            },
        }

        # Create 3 entries per sample: R1 pre-QC, R2 pre-QC, R1 post-QC
        stats_data = {}
        for s_name, data in data_by_sample.items():
            # R1 Pre-QC
            if 'r1_preqc' in data and 'summary' in data['r1_preqc']:
                r1_pre = data['r1_preqc']['summary']['before_filtering']
                stats_data[f"{s_name}: R1 Pre-QC"] = {
                    "total_reads": r1_pre['total_reads'],
                    "q30_rate": r1_pre['q30_rate'],
                    "q30_bases": r1_pre['q30_bases'],
                    "gc_content": r1_pre['gc_content']
                }

            # R2 Pre-QC
            if 'r2_preqc' in data and 'summary' in data['r2_preqc']:
                r2_pre = data['r2_preqc']['summary']['before_filtering']
                stats_data[f"{s_name}: R2 Pre-QC"] = {
                    "total_reads": r2_pre['total_reads'],
                    "q30_rate": r2_pre['q30_rate'],
                    "q30_bases": r2_pre['q30_bases'],
                    "gc_content": r2_pre['gc_content']
                }

            # R1 Post-QC
            if 'r1_postqc' in data and 'summary' in data['r1_postqc']:
                r1_post = data['r1_postqc']['summary']['after_filtering']
                stats_data[f"{s_name}: R1 Post-QC"] = {
                    "total_reads": r1_post['total_reads'],
                    "q30_rate": r1_post['q30_rate'],
                    "q30_bases": r1_post['q30_bases'],
                    "gc_content": r1_post['gc_content']
                }

        self.general_stats_addcols(stats_data, headers)

    def unified_qc_quality_plot(self, data_by_sample):
        """Create quality plot with Pre-QC / Post-QC toggle buttons

        Pre-QC view: Shows R1 and R2 pre-QC data
        Post-QC view: Shows R1 post-QC data only
        """

        # Prepare data for both views
        preqc_data = {}
        postqc_data = {}

        for s_name, data in data_by_sample.items():
            # Pre-QC: Show both R1 and R2
            if 'r1_preqc' in data and 'read1_before_filtering' in data['r1_preqc']:
                r1_before = data['r1_preqc']['read1_before_filtering']
                if 'quality_curves' in r1_before and 'mean' in r1_before['quality_curves']:
                    mean_qual = r1_before['quality_curves']['mean']
                    preqc_data[f"{s_name}: R1 Pre-QC"] = {i+1: q for i, q in enumerate(mean_qual)}

            if 'r2_preqc' in data and 'read1_before_filtering' in data['r2_preqc']:
                r2_before = data['r2_preqc']['read1_before_filtering']
                if 'quality_curves' in r2_before and 'mean' in r2_before['quality_curves']:
                    mean_qual = r2_before['quality_curves']['mean']
                    preqc_data[f"{s_name}: R2 Pre-QC"] = {i+1: q for i, q in enumerate(mean_qual)}

            # Post-QC: Show R1 only
            if 'r1_postqc' in data and 'read1_after_filtering' in data['r1_postqc']:
                r1_after = data['r1_postqc']['read1_after_filtering']
                if 'quality_curves' in r1_after and 'mean' in r1_after['quality_curves']:
                    mean_qual = r1_after['quality_curves']['mean']
                    postqc_data[f"{s_name}: R1 Post-QC"] = {i+1: q for i, q in enumerate(mean_qual)}

        # Configure plot with data_labels for toggle buttons
        pconfig = {
            "id": "unified-qc-seq-quality-plot",
            "title": "Unified QC: Sequence Quality",
            "xlab": "Read Position (bp)",
            "ylab": "Mean Quality Score",
            "ymin": 0,
            "data_labels": [
                {"name": "Pre-QC", "ylab": "Pre-QC: Mean Quality Score"},
                {"name": "Post-QC", "ylab": "Post-QC: Mean Quality Score"}
            ]
        }

        return linegraph.plot([preqc_data, postqc_data], pconfig)

    def unified_qc_gc_plot(self, data_by_sample):
        """Create GC content plot with Pre-QC / Post-QC toggle buttons

        Pre-QC view: Shows R1 and R2 pre-QC data
        Post-QC view: Shows R1 post-QC data only
        """

        preqc_data = {}
        postqc_data = {}

        for s_name, data in data_by_sample.items():
            # Pre-QC: Show both R1 and R2
            if 'r1_preqc' in data and 'read1_before_filtering' in data['r1_preqc']:
                r1_before = data['r1_preqc']['read1_before_filtering']
                if 'content_curves' in r1_before and 'GC' in r1_before['content_curves']:
                    gc_content = r1_before['content_curves']['GC']
                    preqc_data[f"{s_name}: R1 Pre-QC"] = {i+1: gc*100 for i, gc in enumerate(gc_content)}

            if 'r2_preqc' in data and 'read1_before_filtering' in data['r2_preqc']:
                r2_before = data['r2_preqc']['read1_before_filtering']
                if 'content_curves' in r2_before and 'GC' in r2_before['content_curves']:
                    gc_content = r2_before['content_curves']['GC']
                    preqc_data[f"{s_name}: R2 Pre-QC"] = {i+1: gc*100 for i, gc in enumerate(gc_content)}

            # Post-QC: Show R1 only
            if 'r1_postqc' in data and 'read1_after_filtering' in data['r1_postqc']:
                r1_after = data['r1_postqc']['read1_after_filtering']
                if 'content_curves' in r1_after and 'GC' in r1_after['content_curves']:
                    gc_content = r1_after['content_curves']['GC']
                    postqc_data[f"{s_name}: R1 Post-QC"] = {i+1: gc*100 for i, gc in enumerate(gc_content)}

        pconfig = {
            "id": "unified-qc-gc-content-plot",
            "title": "Unified QC: GC Content",
            "xlab": "Read Position (bp)",
            "ylab": "GC Content (%)",
            "ymin": 0,
            "ymax": 100,
            "tt_label": "{point.x}: {point.y:.2f}%",
            "data_labels": [
                {"name": "Pre-QC", "ylab": "Pre-QC: GC Content (%)"},
                {"name": "Post-QC", "ylab": "Post-QC: GC Content (%)"}
            ]
        }

        return linegraph.plot([preqc_data, postqc_data], pconfig)

    def unified_qc_n_plot(self, data_by_sample):
        """Create N content plot with Pre-QC / Post-QC toggle buttons

        Pre-QC view: Shows R1 and R2 pre-QC data
        Post-QC view: Shows R1 post-QC data only
        """

        preqc_data = {}
        postqc_data = {}

        for s_name, data in data_by_sample.items():
            # Pre-QC: Show both R1 and R2
            if 'r1_preqc' in data and 'read1_before_filtering' in data['r1_preqc']:
                r1_before = data['r1_preqc']['read1_before_filtering']
                if 'content_curves' in r1_before and 'N' in r1_before['content_curves']:
                    n_content = r1_before['content_curves']['N']
                    preqc_data[f"{s_name}: R1 Pre-QC"] = {i+1: n*100 for i, n in enumerate(n_content)}

            if 'r2_preqc' in data and 'read1_before_filtering' in data['r2_preqc']:
                r2_before = data['r2_preqc']['read1_before_filtering']
                if 'content_curves' in r2_before and 'N' in r2_before['content_curves']:
                    n_content = r2_before['content_curves']['N']
                    preqc_data[f"{s_name}: R2 Pre-QC"] = {i+1: n*100 for i, n in enumerate(n_content)}

            # Post-QC: Show R1 only
            if 'r1_postqc' in data and 'read1_after_filtering' in data['r1_postqc']:
                r1_after = data['r1_postqc']['read1_after_filtering']
                if 'content_curves' in r1_after and 'N' in r1_after['content_curves']:
                    n_content = r1_after['content_curves']['N']
                    postqc_data[f"{s_name}: R1 Post-QC"] = {i+1: n*100 for i, n in enumerate(n_content)}

        pconfig = {
            "id": "unified-qc-n-content-plot",
            "title": "Unified QC: N Content",
            "xlab": "Read Position (bp)",
            "ylab": "N Content (%)",
            "ymin": 0,
            "ymax": 100,
            "tt_label": "{point.x}: {point.y:.2f}%",
            "data_labels": [
                {"name": "Pre-QC", "ylab": "Pre-QC: N Content (%)"},
                {"name": "Post-QC", "ylab": "Post-QC: N Content (%)"}
            ]
        }

        return linegraph.plot([preqc_data, postqc_data], pconfig)
