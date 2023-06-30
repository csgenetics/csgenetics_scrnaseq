#!/usr/bin/env python

"""
Stats are pulled from the commandline-supplied input files:
    summary_statistics.py $sample_id $h5ad $multiqc_data_json $antisense $dedup $raw_qualimap $filtered_qualimap $annotated_qualimap

    A csv is written out ({sample_id}.metrics.csv) containing:
        the variable name
        the variable value
        the variable human readable name
        the variable description (used to make tooltips in the html report)
        the wider metric classification (see primary keys below)
"""

import anndata
import sys
import json
import re
import numpy as np
from collections import defaultdict

class SummaryStatistics:
    def __init__(self):
        self.sample_id = sys.argv[1]
        # Primary keys will be:
        #   Read QC
        #   Cell metrics
        #   Deduplication
        #   Post read QC alignment
        #   High confidence read alignment
        #   Annotated reads alignment
        # The primary keys will be used to group the metrics for display
        # in the html reports. They will also be used as headers.
        # For each sub key, the value will be a tuple
        # that is the human readable name of the metric
        # (i.e. what is displayed in the html), the
        # value of the metric, and the tooltip text. They will be accessed as
        # [0] [1], and [2] in the jinja template.
        self.metrics_dict = defaultdict(dict)
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

    def write_out_dict_to_csv(self):
        with open(f"{self.sample_id}.metrics.csv", "w") as csv_handle:
            csv_handle.write("variable_name,value,human_readable_name,description,classification\n")
            for metrtic_classification, inner_stats_dict in self.metrics_dict.items():
                for stat_name, stat_value in inner_stats_dict.items():
                    csv_handle.write(f"{stat_name},{stat_value[1]},{stat_value[0]},{stat_value[2]},{metrtic_classification}\n")
    
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
        self.metrics_dict["Cell metrics"]["num_cells"] = ("Number of cells", anndata_array.shape[0], "Estimated number of cells; Number of barcodes passing the nulcear genes detected threshold.")

        # Mean total counts per cell
        self.metrics_dict["Cell metrics"]["mean_total_counts_per_cell"] = ("Mean total counts per cell", np.mean(anndata_array.sum(axis=1)), "Mean of the sum of counts per cell.")
        # Median total counts per cell
        self.metrics_dict["Cell metrics"]["median_total_reads_per_cell"] = ("Median total counts per cell", np.median(anndata_array.sum(axis=1)), "Median of the sum of counts per cell.")

        # Mean single gene counts per cell - all genes (i.e. all counts across all genes considered even those genes not detected in the sample)
        self.metrics_dict["Cell metrics"]["mean_counts_per_cell_all_genes"] = ("Mean single gene count per cell: all genes", np.mean(anndata_array.mean(axis=1)), "Mean of the mean individual gene count per cell with all genes considered (including those genes not detected in the sample).")
        # Median single gene counts per cell - all genes (i.e. all counts across all genes considered even those genes not detected in the sample)
        self.metrics_dict["Cell metrics"]["median_counts_per_cell_all_genes"] = ("Median single gene count per cell: all genes", np.median(np.median(anndata_array, axis=1)), "Median of the median individual gene count per cell with all genes considered (including those genes not detected in the sample).")

        # Mean single genecounts per cell - sample detected genes (i.e. only those genes that were detected in the sample)
        self.metrics_dict["Cell metrics"]["mean_counts_per_cell_sample_detected_genes"] = ("Mean single gene count per cell: sample-detected genes", np.mean(anndata_array_detected_genes.mean(axis=1)), "Mean of the mean individual gene count per cell with only those genes detected in the sample considered.")
        # Median single gene counts per cell - sample detected genes (i.e. only those genes that were detected in the sample)
        self.metrics_dict["Cell metrics"]["median_counts_per_cell_sample_detected_genes"] = ("Median single gene count per cell: sample-detected genes", np.median(anndata_array_detected_genes.sum(axis=1)), "Median of the median individual gene count per cell with only those genes detected in the sample considered.")

        # Mean single gene count per cell -  cell detected genes (i.e. only those genes that were detected in the given cell)
        self.metrics_dict["Cell metrics"]["mean_counts_per_cell_cell_detected_genes"] = ("Mean single gene count per cell: cell-detected genes", np.mean(np.apply_along_axis(self.get_non_zero_sum, 1, anndata_array)), "Mean of the mean individual gene count per cell with only those genes detected in each cell considered.")
        # Median single gene count per cell -  cell detected genes (i.e. only those genes that were detected in the given cell)
        self.metrics_dict["Cell metrics"]["median_counts_per_cell_cell_detected_genes"] = ("Median single gene count per cell: cell-detected genes", np.median(np.apply_along_axis(self.get_non_zero_sum, 1, anndata_array)), "Median of the median individual gene count per cell with only those genes detected in each cell considered.")
        
        # Mean genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_genes_detected_per_cell"] = ("Mean genes detected per cell", np.mean(anndata_array.astype(bool).sum(axis=1)), "Mean number of genes detected for each cell (including nuclear and mitochondrial genes).")
        # Median genes detected per cell
        self.metrics_dict["Cell metrics"]["median_genes_detected_per_cell"] = ("Median genes detected per cell", int(np.median(anndata_array.astype(bool).sum(axis=1))), "Median number of genes detected for the cells (including nuclear and mitochondrial genes).")

        # Get a subset of the array that doesn't contain the mito genes
        anndata_array_nuc = anndata_array[:,~self.anndata.var["is_mito"]]
        # Mean nuclear genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_nuclear_genes_detected_per_cell"] = ("Mean nuclear (non-mitochondrial) genes detected per cell", np.mean(anndata_array_nuc.astype(bool).sum(axis=1)), "Mean number of nuclear genes detected for each cell.")
        # Median nuclear genes detected per cell
        self.metrics_dict["Cell metrics"]["median_nuclear_genes_detected_per_cell"] = ("Median nuclear (non-mitochondrial) genes detected per cell", np.median(anndata_array_nuc.astype(bool).sum(axis=1)), "Median number of nuclear genes detected for each cell.")

        # Get a subset of the array that contains only the mito genes
        anndata_array_mito = anndata_array[:,self.anndata.var["is_mito"]]
        # Mean mitochondrial genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_mito_genes_detected_per_cell"] = ("Mean mitochondrial genes detected per cell", np.mean(anndata_array_mito.astype(bool).sum(axis=1)), "Mean number of mitochondrial genes detected for each cell.")
        # Median mitochondrial genes detected per cell
        self.metrics_dict["Cell metrics"]["median_mito_genes_detected_per_cell"] = ("Median mitochondrial genes detected per cell", int(np.median(anndata_array_mito.astype(bool).sum(axis=1))), "Median number of mitochondrial genes detected for each cell.")

        # Unique genes detected across samples (i.e. a gene can only be detected once per sample)
        self.metrics_dict["Cell metrics"]["num_unique_genes_detected_across_sample"] = ("Unique genes detected across sample", anndata_array_detected_genes.shape[1], "Number of unique genes detected across the sample (each gene can be counted only once even if found in multiple cells).")

        # Total genes detected across samples (i.e. genes detected per cell summed)
        self.metrics_dict["Cell metrics"]["total_genes_detected_across_samples"] = ("Total genes detected across sample", np.count_nonzero(anndata_array_detected_genes), "Total number of genes detected across the sample (each gene can be counted more than once if detected in more than one cell).")
    
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
            self.metrics_dict["annotated_alignment_stats"]["annotated_antisense_mapping_reads"] = ("Reads aligned antisense", int(antisense_handle.read().rstrip()), "Number of reads aligned antisense.")

    def get_duplication_stats(self):
        # Reads before deduplication
        # Reads after deduplication
        # Deduplication perc
        with open(sys.argv[5], "r") as dedup_handle:
            for line in dedup_handle:
                if "INFO Reads: Input Reads:" in line:
                    try:
                        reads_in = int(line.split()[-1].strip())
                    except ValueError:
                        reads_in = 0
                    self.metrics_dict["Deduplication"]["reads_before_deduplication"] = ("Reads before deduplication", reads_in, "Number of reads before deduplication.")
                elif "INFO Number of reads out:" in line:
                    self.metrics_dict["Deduplication"]["reads_after_deduplication"] = ("Reads after deduplication", int(line.split()[-1].strip()), "Number of reads after deduplication.")
        if reads_in != 0:
            self.metrics_dict["Deduplication"]["duplication_perc"] = ("Deduplication percent", self.as_perc(1 - (self.metrics_dict["Deduplication"]["reads_after_deduplication"][1]/self.metrics_dict["Deduplication"]["reads_before_deduplication"][1])), "Reads before deduplication / reads after deduplication * 100")
        else:
            self.metrics_dict["Deduplication"]["duplication_perc"] = ("Deduplication percent", 0.0, "Reads before deduplication / reads after deduplication * 100")
    
    def get_mapping_stats(self):
        # Populate self.metrics_dict with the raw qualimap stats
        self.get_qualimap_stats(sys.argv[6], "Post read QC alignment")
        # Populate self.metrics_dict with the filtered qualimap stats
        self.get_qualimap_stats(sys.argv[7], "High confidence read alignment")
        # Populate self.metrics_dict with the annotated qualimap stats
        self.get_qualimap_stats(sys.argv[8], "Annotated reads alignment")

    def get_qualimap_stats(self, path, category):
        # Read in the qualimap output for the unfiltered mapping
        with open(path, "r") as raw_qualimap_handle:
            lines = [_.strip() for _ in raw_qualimap_handle]
            # Reads aligned to genome
            for i, line in enumerate(lines):
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
                        self.metrics_dict[category][f"{header}"] = (header.replace("_", " ").capitalize(), val, header.replace("_", " ").capitalize())
                
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
                        self.metrics_dict[category][f"{header}"] = (header.replace("_", " ").capitalize(), val_absolute, header.replace("_", " ").capitalize())
                        self.metrics_dict[category][f"{header}_perc"] = (header.replace("_", " ").capitalize() + " percentage", val_percent, header.replace("_", " ").capitalize() + " percentage")

    @staticmethod
    def as_perc(float_to_convert):
        return float_to_convert * 100
        # return f"{float_to_convert:.2f}%"

    def get_trimming_qc_stats(self):
        # Reads pre-QC
        self.metrics_dict["Read QC"]["reads_pre_qc"] = ("Number of reads pre-QC", int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1"]["before_filtering_total_reads"]), "Number of reads in the input R1 fastq files (after merging if applicable).")
        # Reads containing cellular barcode matching whitelist
        self.metrics_dict["Read QC"]["valid_barcode_reads"] = ("Number of valid barcode-containing reads", int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.io_extract"]["before_filtering_total_reads"]), "Number of reads containing a barcode exactly matching the whitelist.")
        # Percentage of reads containing cellular barcode matching whitelist as percentage of pre-QC reads
        self.metrics_dict["Read QC"]["valid_barcode_reads_perc"] = ("Percentage valid barcode-containing reads", self.as_perc(float(self.metrics_dict["Read QC"]["valid_barcode_reads"][1] / self.metrics_dict["Read QC"]["reads_pre_qc"][1])), "Number of valid barcode-containing reads / Number of reads pre-QC * 100.")
        # Percentage of barcode bases >= Q30
        self.metrics_dict["Read QC"]["barcode_bases_q30_perc"] = ("Barcode bp >= Q30 percentage", self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R2"]["after_filtering_q30_rate"])), "The percentage of the barcode bases with a Phred score >= 30.")
        # Reads after polyX tail and polyA internal trimming
        self.metrics_dict["Read QC"]["reads_post_trimming"] = ("Number of reads post-QC trimming", int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_total_reads"]), "Number of reads after polyX tail and polyA internal trimming.")
        # Mean read length after polyX tail and polyA internal trimming
        self.metrics_dict["Read QC"]["mean_post_trim_read_length"] = ("Mean read length post-QC trimming", float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_read1_mean_length"]), "Mean R1 read length post-QC trimming.")
        # Percentage of bases post trimming >= Q30
        self.metrics_dict["Read QC"]["rna_bases_q30_perc"] = ("R1 bp >= Q30 percentage; post-QC trimming", self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_q30_rate"])), "The percentage of the R1 bases (post-QC trimming) with a Phred score >= 30.")


if __name__ == "__main__":
    SummaryStatistics().generate_metrics()
