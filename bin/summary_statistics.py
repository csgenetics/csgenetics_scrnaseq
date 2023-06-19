#!/usr/bin/env python

"""
Reads in 4 input data files:
1 - _antisense.txt
2 - .multiqc.data.json
3 - _qualimap.txt
4 - .cell_only.count_matrix.h5ad

Input files are provided as command line inputs
summary_statistics.py Sample1 Sample1.cell_only.count_matrix.h5ad Sample1.multiqc.data.json Sample1_antisense.txt 
"""

import anndata
import sys
import json

class SummaryStatistics:
    def __init__(self):
        self.sample_id = sys.argv[1]
        self.metrics_dict = {}
        with open(sys.argv[2], "r") as json_handle:
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
        Get metrics related to sequencing
        # TODO use the comments as tool tips in the html documents.
        # TODO see if we can get rid of the polyA stats file and get
        # everything we need from the fastp outputs.
        
        command line input is: $sample_id $h5ad $multiqc_data_json $antisense $dedup $qualimap

        """
        self.get_trimming_qc_stats()
        self.get_mapping_stats()
        self.get_filtered_mapping_stats()
        # Proportion duplication (reads in to deduplication / reads out of deduplication)
        self.metrics_dict["duplication_proportion"] = self.get_duplication_proportion()
        

    def get_trimming_qc_stats(self):
        # Reads pre-QC
        self.metrics_dict["reads_pre_qc"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1"]["before_filtering_total_reads"])
        # Reads containing cellular barcode matching whitelist
        self.metrics_dict["valid_barcode_reads"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.io_extract"]["before_filtering_total_reads"])
        # Proportion of reads containing cellular barcode matching whitelist as percentage of pre-QC reads
        self.metrics_dict["prop_valid_barcode_reads"] = float(self.metrics_dict["valid_barcode_reads"] / self.metrics_dict["reads_pre_qc"])
        # Proportion of barcode bases >= Q30
        self.metrics_dict["prop_barcode_bases_q30"] = float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R2"]["after_filtering_q30_rate"])
        # Reads after polyX tail and polyA internal trimming
        self.metrics_dict["reads_post_trimming"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_total_reads"])
        # Mean read length after polyX tail and polyA internal trimming
        self.metrics_dict["mean_post_trim_read_length"] = float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_read1_mean_length"])
        # Proportion of bases post trimming >= Q30
        self.metrics_dict["prop_rna_bases_q30"] = float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_q30_rate"])

    def get_duplication_proportion(self):
        with open(sys.argv[5], "r") as dedup_handle:
            dedup_reads_in = None
            dedup_reads_out = None
            for line in dedup_handle:
                if "INFO Reads: Input Reads:" in line:
                    dedup_reads_in = int(line.split()[-1])
                elif "INFO Number of reads out:" in line:
                    dedup_reads_out = int(line.split()[-1])
            
            if dedup_reads_in is None or dedup_reads_out is None:
                raise RuntimeError("Unable to fetch dedup stats")
            
            # If no in reads return 0
            if dedup_reads_in == 0:
                return 0
            else:
                return dedup_reads_out/dedup_reads_in


if __name__ == "__main__":
    SummaryStatistics().generate_metrics()