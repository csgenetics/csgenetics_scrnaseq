#!/usr/bin/env python

"""
Stats are pulled from the commandline-supplied input files:
    summary_statistics.py $sample_id $h5ad $multiqc_data_json $antisense $dedup $raw_qualimap $filtered_qualimap $annotated_qualimap
"""

import anndata
import sys
import json
import re
import numpy as np
import csv

class SummaryStatistics:
    def __init__(self):
        self.sample_id = sys.argv[1]
        self.metrics_dict = {"cell_stats":{},
                             "sequence_stats":{},
                             "alignment_stats":{},
                             "duplication_stats":{}}
        with open(sys.argv[3], "r") as json_handle:
            self.multiqc_json_dict = json.load(json_handle)
        self.qc_key_to_index_dict = self.make_qc_key_to_index_dict()
        
    def make_qc_key_to_index_dict(self):
        """
            Need to get the order of the various report sources
            in multiQC json structure in the report_data_sources section
            to ensure that we are indexing the
            correct report item in the report_general_stats_data.
            We will hold the indices as values to the report name
            keys.
        """
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
        return self.qc_key_to_index_dict

    def generate_metrics(self):
        self.get_sequencing_stats()
        self.get_cell_stats()
        self.write_out_dict_to_csv()
        self.write_out_dict_to_csv()

    def write_out_dict_to_csv(self):
        with open(f"{self.sample_id}.metrics.csv", "w") as csv_handle:
            for dict_id, v in self.metrics_dict.items():
                for dict_id2,v2 in v.items():
                    csv_handle.write(f"{dict_id2},{v2},{dict_id}\n")
    
    # def write_out_dict_to_csv(self):
    #     with open(f"{self.sample_id}.metrics.json","w") as json_handle:
    #         self.json_object = json.dumps(self.metrics_dict)
    #         json_handle.write(self.json_object)

    @staticmethod
    def get_non_zero_sum(np_1d_array):
        return np.sum(np_1d_array[np.nonzero(np_1d_array)])

    def get_cell_stats(self):
        self.anndata = anndata.read_h5ad(sys.argv[2])
        self.anndata.var["is_mito"] = np.where(self.anndata.var['geneSym'].str.startswith("MT-"), True, False)
        # Convert to array for convenience
        anndata_array = self.anndata.X.toarray()
        # Remove genes that have not been detected at all across the sample
        anndata_array_detected_genes = anndata_array[:, ~np.all(anndata_array == 0, axis=0)]
        # Estimated number of cells
        self.metrics_dict["cell_stats"]["num_cells"] = anndata_array.shape[0]

        # Mean counts per cell - all genes (i.e. all counts across all genes considered even those genes not detected in the sample)
        self.metrics_dict["cell_stats"]["mean_counts_per_cell_all_genes"] = np.mean(anndata_array.sum(axis=1))
        # Median counts per cell - all genes (i.e. all counts across all genes considered even those genes not detected in the sample)
        self.metrics_dict["cell_stats"]["median_counts_per_cell_all_genes"] = int(np.median(anndata_array.sum(axis=1)))

        # Mean counts per cell - sample detected genes (i.e. only those genes that were detected in the sample)
        self.metrics_dict["cell_stats"]["mean_counts_per_cell_sample_detected_genes"] = np.mean(anndata_array_detected_genes.sum(axis=1))
        # Median counts per cell - sample detected genes (i.e. only those genes that were detected in the sample)
        self.metrics_dict["cell_stats"]["median_counts_per_cell_sample_detected_genes"] = int(np.median(anndata_array_detected_genes.sum(axis=1)))

        # Mean count per cell -  cell detected genes (i.e. only those genes that were detected in the given cell)
        self.metrics_dict["cell_stats"]["mean_counts_per_cell_cell_detected_genes"] = np.mean(np.apply_along_axis(self.get_non_zero_sum, 1, anndata_array))
        # Median count per cell -  cell detected genes (i.e. only those genes that were detected in the given cell)
        self.metrics_dict["cell_stats"]["median_counts_per_cell_cell_detected_genes"] = int(np.median(np.apply_along_axis(self.get_non_zero_sum, 1, anndata_array)))
        
        # Mean genes detected per cell
        self.metrics_dict["cell_stats"]["mean_genes_detected_per_cell"] = np.mean(anndata_array.astype(bool).sum(axis=1))
        # Median genes detected per cell
        self.metrics_dict["cell_stats"]["median_genes_detected_per_cell"] = int(np.median(anndata_array.astype(bool).sum(axis=1)))

        # Get a subset of the array that doesn't contain the mito genes
        anndata_array_nuc = anndata_array[:,~self.anndata.var["is_mito"]]
        # Mean nuclear genes detected per cell
        self.metrics_dict["cell_stats"]["mean_nuclear_genes_detected_per_cell"] = np.mean(anndata_array_nuc.astype(bool).sum(axis=1))
        # Median nuclear genes detected per cell
        self.metrics_dict["cell_stats"]["median_nuclear_genes_detected_per_cell"] = int(np.median(anndata_array_nuc.astype(bool).sum(axis=1)))

        # Get a subset of the array that contains only the mito genes
        anndata_array_mito = anndata_array[:,self.anndata.var["is_mito"]]
        # Mean mitochondrial genes detected per cell
        self.metrics_dict["cell_stats"]["mean_mito_genes_detected_per_cell"] = np.mean(anndata_array_mito.astype(bool).sum(axis=1))
        # Median mitochondrial genes detected per cell
        self.metrics_dict["cell_stats"]["median_mito_genes_detected_per_cell"] = int(np.median(anndata_array_mito.astype(bool).sum(axis=1)))

        # Unique genes detected across samples (i.e. a gene can only be detected once per sample)
        self.metrics_dict["cell_stats"]["num_unique_genes_detected_across_sample"] = anndata_array_detected_genes.shape[1]

        # Total genes detected across samples (i.e. genes detected per cell summed)
        self.metrics_dict["cell_stats"]["total_genes_detected_across_samples"] = np.count_nonzero(anndata_array_detected_genes)
    
    def get_sequencing_stats(self):
        """
        Populate the self.metrics_dict with the stats
        """
        self.get_trimming_qc_stats()
        self.get_mapping_stats()
        self.get_duplication_stats()


    def get_antisense(self):
        """
        Read in the antisense metric input files and pull out
        the antisense mapped reads and populate into
        self.metrics_dict

        TODO currently we only collect the antisense for the allocated filtered
        bam. Ideally we would want this stat for each level of filtering
        I.e. raw, umr filtered and allocated
        """
        with open(sys.argv[4], "r") as antisense_handle:
            self.metrics_dict["annotated_antisense_mapping_reads"] = int(antisense_handle.read().rstrip())

    def get_duplication_stats(self):
        # Reads before deduplication
        # Reads after deduplication
        # Deduplication perc
        with open(sys.argv[5], "r") as dedup_handle:
            for line in dedup_handle:
                if "INFO Reads: Input Reads:" in line:
                    self.metrics_dict["duplication_stats"]["reads_before_deduplication"] = int(line.split()[-1].strip())
                elif "INFO Number of reads out:" in line:
                    self.metrics_dict["duplication_stats"]["reads_after_deduplication"] = int(line.split()[-1].strip())
        self.metrics_dict["duplication_stats"]["duplication_perc"] = self.as_perc(1 - (self.metrics_dict["duplication_stats"]["reads_after_deduplication"]/self.metrics_dict["duplication_stats"]["reads_before_deduplication"]))

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
                        self.metrics_dict["alignment_stats"][f"{prefix}_{header}"] = val
                
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
                        val_percent = float(re.search(r"((?<=\()[\d]*\.[\d]+)%(?=\))", lines[i]).groups()[0].rstrip())
                        self.metrics_dict["alignment_stats"][f"{prefix}_{header}"] = val_absolute
                        self.metrics_dict["alignment_stats"][f"{prefix}_{header}_perc"] = val_percent

    @staticmethod
    def as_perc(float_to_convert):
        return float_to_convert * 100
        # return f"{float_to_convert:.2f}%"

    def get_trimming_qc_stats(self):
        # Reads pre-QC
        self.metrics_dict["sequence_stats"]["reads_pre_qc"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1"]["before_filtering_total_reads"])
        # Reads containing cellular barcode matching whitelist
        self.metrics_dict["sequence_stats"]["valid_barcode_reads"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.io_extract"]["before_filtering_total_reads"])
        # Proportion of reads containing cellular barcode matching whitelist as percentage of pre-QC reads
        self.metrics_dict["sequence_stats"]["valid_barcode_reads_perc"] = self.as_perc(float(self.metrics_dict["sequence_stats"]["valid_barcode_reads"] / self.metrics_dict["sequence_stats"]["reads_pre_qc"]))
        # Proportion of barcode bases >= Q30
        self.metrics_dict["sequence_stats"]["barcode_bases_q30_perc"] = self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R2"]["after_filtering_q30_rate"]))
        # Reads after polyX tail and polyA internal trimming
        self.metrics_dict["sequence_stats"]["reads_post_trimming"] = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_total_reads"])
        # Mean read length after polyX tail and polyA internal trimming
        self.metrics_dict["sequence_stats"]["mean_post_trim_read_length"] = float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_read1_mean_length"])
        # Proportion of bases post trimming >= Q30
        self.metrics_dict["sequence_stats"]["rna_bases_q30_perc"] = self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_q30_rate"]))


if __name__ == "__main__":
    SummaryStatistics().generate_metrics()