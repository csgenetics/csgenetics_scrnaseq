#!/usr/bin/env python

"""
Stats are pulled from the commandline-supplied input files:
    summary_statistics.py $sample_id $h5ad $multiqc_data_json $antisense $dedup $raw_qualimap $filtered_qualimap $annotated_qualimap
"""

import anndata
import sys
import json
import re

class SummaryStatistics:
    def __init__(self):
        self.sample_id = sys.argv[1]
        self.metrics_dict = {}
        with open(sys.argv[3], "r") as json_handle:
            self.multiqc_json_dict = json.load(json_handle)
        # Need to get the order of the various report sources
        # in the report_data_sources section
        # to ensure that we are indexing the
        # correct report item in the report_general_stats_data.
        # We will hold the indices as values to the report name
        # keys
        self.qc_key_to_index_dict = {}
        for i, key in enumerate(self.multiqc_json_dict["report_data_sources"].keys()):
            if key == "featureCounts":
                self.qc_key_to_index_dict["featureCounts"] = i
            elif key == "STAR":
                self.qc_key_to_index_dict["STAR"] = i
            elif key == "fastp":
                self.qc_key_to_index_dict["fastp"] = i
            else:
                raise RuntimeError(f"multiQC key: {key} unrecognised.")

    def generate_metrics(self):
        self.get_sequencing_stats()
        self.anndata = anndata.read_h5ad(sys.argv[2])

        with open(sys.argv[4], "r") as antisense_handle:
            self.antisense_mapped_reads = int(antisense_handle.read().rstrip())
        
    def get_sequencing_stats(self):
        """
        Populate the self.metrics_dict with the stats
        """
        self.get_trimming_qc_stats()
        self.get_mapping_stats()
        self.get_duplication_stats()
        self.get_cell_stats()
    
    def get_allocation_stats(self):
        # TODO update this so that it also includes the filtering of the Intergenic tags
        # Allocated reads (reads with XT tag)
        raise NotImplementedError()

    def get_duplication_stats(self):
        # Reads before deduplication
        # Reads after deduplication
        # Deduplication perc
        with open(sys.argv[5], "r") as dedup_handle:
            for line in dedup_handle:
                if "INFO Reads: Input Reads:" in line:
                    self.metrics_dict["reads_before_deduplication"] = int(line.split()[-1].strip())
                elif "INFO Number of reads out:" in line:
                    self.metrics_dict["reads_after_deduplication"] = int(line.split()[-1].strip())
        self.metrics_dict["duplication_perc"] = self.as_perc(1 - (self.metrics_dict["reads_after_deduplication"]/self.metrics_dict["reads_before_deduplication"]))

    def get_mapping_stats(self):
        # Populate self.metrics_dict with the raw qualimap stats
        self.get_qualimap_stats(sys.argv[6], "raw")
        # Populate self.metrics_dict with the filtered qualimap stats
        self.get_qualimap_stats(sys.argv[7], "filtered")
        # Populate self.metrics_dict with the annotated qualimap stats
        self.get_qualimap_stats(sys.argv[8], "annotated")

    def get_qualimap_stats(self, path, prefix):
        # Read in the qualimap output for the unfiltered mapping
        with open(path, "r") as raw_qualimap_handle:
            lines = [_.strip() for _ in raw_qualimap_handle]
            # Reads aligned to genome
            for i, line in enumerate(lines):
                    # self.metrics_dict[]
                if ">>>>>>> Reads alignment" in line:
                    # When we find the Reads alignment section
                    # we want the following 8 lines (allowing for a blank line)
                    # directly after the header
                    # We capture the following metrics
                    #   reads aligned
                    #   total alignments
                    #   secondary alignments
                    #   non-unique alignments
                    #   aligned to genes
                    #   ambiguous alignments
                    #   no feature assigned
                    #   not aligned
                    for i in range(i+2, i+10):
                        header = re.search("([\w\s\-]*(?=\=))", lines[i]).groups()[0].rstrip().replace(" ", "_")
                        val = int(re.search("((?<=\=\s)[\d,]*\s*$)", lines[i]).groups()[0].rstrip().replace(",",""))
                        self.metrics_dict[f"{prefix}_{header}"] = val
                
                if ">>>>>>> Reads genomic origin" in line:
                    # When we find the Reads genomic origin region
                    # we want the following 4 lines allowing for 
                    # a line directly after the header
                    # We capture the following metrics
                    #   exonic
                    #   intronic
                    #   intergenic
                    #   overlapping exon
                    for i in range(i+2, i+6):
                        header = re.search("([\w\s\-]*(?=\=))", lines[i]).groups()[0].rstrip().replace(" ", "_")
                        val_absolute = int(re.search("((?<=\s)[\d,]*\s+(?=\())", lines[i]).groups()[0].rstrip().replace(",",""))
                        val_percent = re.search(r"((?<=\()[\d]*\.[\d]{2}%(?=\)))", lines[i]).groups()[0].rstrip()
                        self.metrics_dict[f"{prefix}_{header}"] = val_absolute
                        self.metrics_dict[f"{prefix}_{header}_perc"] = val_percent

    @staticmethod
    def as_perc(float_to_convert):
        return f"{float_to_convert:.2f}%"

    def get_trimming_qc_stats(self):
        # Reads pre-QC
        self.metrics_dict["reads_pre_qc"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1"]["before_filtering_total_reads"])
        # Reads containing cellular barcode matching whitelist
        self.metrics_dict["valid_barcode_reads"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.io_extract"]["before_filtering_total_reads"])
        # Proportion of reads containing cellular barcode matching whitelist as percentage of pre-QC reads
        self.metrics_dict["valid_barcode_reads_perc"] = self.as_perc(float(self.metrics_dict["valid_barcode_reads"] / self.metrics_dict["reads_pre_qc"]) * 100)
        # Proportion of barcode bases >= Q30
        self.metrics_dict["barcode_bases_q30_perc"] = self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R2"]["after_filtering_q30_rate"]))
        # Reads after polyX tail and polyA internal trimming
        self.metrics_dict["reads_post_trimming"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_total_reads"])
        # Mean read length after polyX tail and polyA internal trimming
        self.metrics_dict["mean_post_trim_read_length"] = float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_read1_mean_length"])
        # Proportion of bases post trimming >= Q30
        self.metrics_dict["rna_bases_q30_perc"] = self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_q30_rate"]))


if __name__ == "__main__":
    SummaryStatistics().generate_metrics()