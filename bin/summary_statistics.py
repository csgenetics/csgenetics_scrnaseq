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
import pandas as pd

class SummaryStatistics:
    def __init__(self):
        self.sample_id = sys.argv[1]
        if sys.argv[9] == "TRUE":
            self.mixed = True

            # The purity threshold used to classify a barcode as a cell
            # (in addition to the num nuc genes detected threshold)
            # if we are working with a mixed species sample
            self.purity = sys.argv[10]
        else:
            self.mixed = False

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

    def set_cell_stats_to_zero(self):
        """
        Set all of the cell stats to 0
        For the mixed species case we produce the cell stats
        essentially in triplicate, one for total cells (i.e. 
        Hsap + Mmus), one Hsap and one Mmus.
        Because we produce the stats in triplicate we use
        a list of strings that represent each of the base
        stats that will be triplicated.
        """
        base_stats = [
            "num_cells",
            "raw_reads_per_cell",
            "mean_total_counts_per_cell",
            "median_total_counts_per_cell",
            "mean_genes_detected_per_cell",
            "median_genes_detected_per_cell",
            "mean_nuclear_genes_detected_per_cell",
            "median_nuclear_genes_detected_per_cell",
            "mean_mito_genes_detected_per_cell",
            "median_mito_genes_detected_per_cell",
            "percentage_counts_from_mito",
            "num_unique_genes_detected_across_sample",
            "total_genes_detected_across_samples"
        ]
        if self.mixed:
            # Output the full set of mixed species cell_stats
            for base_stat in base_stats:
                for stat_type in ["_total", "_Hsap", "_Mmus"]:
                    setattr(self, f"{base_stat}{stat_type}", 0)
        else:
            for base_stat in base_stats:
                setattr(self, f"{base_stat}", 0)

    def calculate_cell_stats(self):
        # strip the anndata object down to only the is_single_cell cells
        # so that all metrics are only caluclated for single_cells.
        # If single species this means those that meet the total counts threshold
        # If mixed species then the barcodes must also pass the purity threshold.
        self.anndata_sc = self.anndata[self.anndata.obs['is_single_cell']]

        # Subset to get rid of genes that have 0 counts 
        self.anndata_sc = self.anndata_sc[:,~np.all(self.anndata_sc.X.toarray() == 0, axis=0)]
        
        # Convert to array for convenience
        anndata_array_sc = self.anndata_sc.to_df()
        
        if self.mixed:
            # For the Hsap and Mmus single cell df we work with
            # only their associated genes. I.e. the metrics do not
            # include genes / counts from the other species.
            anndata_array_sc_Hsap = self.anndata_sc[self.anndata_sc.obs["species_based_on_purity"] == "Hsap", self.anndata_sc.var["is_hsap"]].to_df()
            anndata_array_sc_Mmus = self.anndata_sc[self.anndata_sc.obs["species_based_on_purity"] == "Mmus", self.anndata_sc.var["is_mmus"]].to_df()
            
            # If mixed then we need to caluculate 3 versions of each of the metrics:
            #   _total, _Hsap, _Mmus
            self.num_cells_total = anndata_array_sc.shape[0]
            self.num_cells_Hsap = anndata_array_sc_Hsap.shape[0]
            self.num_cells_Mmus = anndata_array_sc_Mmus.shape[0]

            self.raw_reads_per_cell_total = self.metrics_dict["Read QC"]["reads_pre_qc"][1] / self.num_cells_total
            self.raw_reads_per_cell_Hsap = self.metrics_dict["Read QC"]["reads_pre_qc"][1] / self.num_cells_Hsap
            self.raw_reads_per_cell_Mmus = self.metrics_dict["Read QC"]["reads_pre_qc"][1] / self.num_cells_Mmus

            # For the calculation of the total values, we make a concat of the two species
            # sum series so that for the _total version of the stat, Hsap cells only count
            # Hsap genes and vice versa for Mmus
            summed_counts_Hsap = anndata_array_sc_Hsap.sum(axis=1)
            summed_counts_Mmus = anndata_array_sc_Mmus.sum(axis=1)
            concat_counts_species_series = pd.concat([summed_counts_Hsap, summed_counts_Mmus])
            
            self.mean_total_counts_per_cell_total = np.mean(concat_counts_species_series)
            self.mean_total_counts_per_cell_Hsap = np.mean(summed_counts_Hsap)
            self.mean_total_counts_per_cell_Mmus = np.mean(summed_counts_Mmus)

            self.median_total_counts_per_cell_total = np.median(concat_counts_species_series)
            self.median_total_counts_per_cell_Hsap = np.median(summed_counts_Hsap)
            self.median_total_counts_per_cell_Mmus = np.median(summed_counts_Mmus)

            # Similar to the counts we create a concatenated mixed species series
            # so that for Hsap cells, only Hsap detected genes are counted
            # and for Mmus cells only Mmus genes are detected
            summed_genes_detected_Hsap = anndata_array_sc_Hsap.astype(bool).sum(axis=1)
            summed_genes_detected_Mmus = anndata_array_sc_Mmus.astype(bool).sum(axis=1)
            concat_genes_detected_species_series = pd.concat([summed_genes_detected_Hsap, summed_genes_detected_Mmus])

            self.mean_genes_detected_per_cell_total = np.mean(concat_genes_detected_species_series)
            self.mean_genes_detected_per_cell_Hsap = np.mean(summed_genes_detected_Hsap)
            self.mean_genes_detected_per_cell_Mmus = np.mean(summed_genes_detected_Mmus)

            self.median_genes_detected_per_cell_total = int(np.median(concat_genes_detected_species_series))
            self.median_genes_detected_per_cell_Hsap = int(np.median(summed_genes_detected_Hsap))
            self.median_genes_detected_per_cell_Mmus = int(np.median(summed_genes_detected_Mmus))

            # Get a subset of the arrays that don't contain the mito genes
            anndata_array_sc_nuc_Hsap = anndata_array_sc_Hsap.loc[:, ~self.anndata.var["is_mito_hsap"]]
            anndata_array_sc_nuc_Mmus = anndata_array_sc_Mmus.loc[:, ~self.anndata.var["is_mito_mmus"]]
            
            # As above we create a concatenated series of the two individual species
            # so that only Hsap genes are detected for Hsap cells and vice versa for Mmus cells
            summed_nuc_genes_detected_Hsap = anndata_array_sc_nuc_Hsap.astype(bool).sum(axis=1)
            summed_nuc_genes_detected_Mmus = anndata_array_sc_nuc_Mmus.astype(bool).sum(axis=1)
            concat_nuc_genes_detected_species_series = pd.concat([summed_nuc_genes_detected_Hsap, summed_nuc_genes_detected_Mmus])

            self.mean_nuclear_genes_detected_per_cell_total = np.mean(concat_nuc_genes_detected_species_series)
            self.mean_nuclear_genes_detected_per_cell_Hsap = np.mean(summed_nuc_genes_detected_Hsap)
            self.mean_nuclear_genes_detected_per_cell_Mmus = np.mean(summed_nuc_genes_detected_Mmus)

            self.median_nuclear_genes_detected_per_cell_total = np.median(concat_nuc_genes_detected_species_series)
            self.median_nuclear_genes_detected_per_cell_Hsap = np.median(summed_nuc_genes_detected_Hsap)
            self.median_nuclear_genes_detected_per_cell_Mmus = np.median(summed_nuc_genes_detected_Mmus)

            # Get a subset of the array that contains only the mito genes
            anndata_array_sc_mito_Hsap = anndata_array_sc_Hsap.loc[:, self.anndata.var["is_mito_hsap"]]
            anndata_array_sc_mito_Mmus = anndata_array_sc_Mmus.loc[:, self.anndata.var["is_mito_mmus"]]

            # As above we create a concatenated series made up of both of the individual species
            # series
            summed_mito_genes_detected_Hsap = anndata_array_sc_mito_Hsap.astype(bool).sum(axis=1)
            summed_mito_genes_detected_Mmus = anndata_array_sc_mito_Mmus.astype(bool).sum(axis=1)
            concat_mito_genes_detected_species_series = pd.concat([summed_mito_genes_detected_Hsap, summed_mito_genes_detected_Mmus])

            self.mean_mito_genes_detected_per_cell_total = np.mean(concat_mito_genes_detected_species_series)
            self.mean_mito_genes_detected_per_cell_Hsap = np.mean(summed_mito_genes_detected_Hsap)
            self.mean_mito_genes_detected_per_cell_Mmus = np.mean(summed_mito_genes_detected_Mmus)

            self.median_mito_genes_detected_per_cell_total = np.median(concat_mito_genes_detected_species_series)
            self.median_mito_genes_detected_per_cell_Hsap = np.median(summed_mito_genes_detected_Hsap)
            self.median_mito_genes_detected_per_cell_Mmus = np.median(summed_mito_genes_detected_Mmus)

            # Calculate percentage of counts of mitochondrial origin
            # Similar to above we only want to count mito counts that
            # come from genes of the species. The concat_counts_species_series.sum()
            # gives us the total mito counts for the two species in this way. Then
            # we get the counts for each individual species from the mito series.
            self.percentage_counts_from_mito_total = self.as_perc((anndata_array_sc_mito_Hsap.sum(axis=1).sum() + anndata_array_sc_mito_Mmus.sum(axis=1).sum()) / concat_counts_species_series.sum())
            self.percentage_counts_from_mito_Hsap = self.as_perc(anndata_array_sc_mito_Hsap.values.sum() / anndata_array_sc_Hsap.values.sum())
            self.percentage_counts_from_mito_Mmus = self.as_perc(anndata_array_sc_mito_Mmus.values.sum() / anndata_array_sc_Mmus.values.sum())

            self.num_unique_genes_detected_across_sample_total = anndata_array_sc.sum(axis=0).astype(bool).sum()
            self.num_unique_genes_detected_across_sample_Hsap = anndata_array_sc_Hsap.sum(axis=0).astype(bool).sum()
            self.num_unique_genes_detected_across_sample_Mmus = anndata_array_sc_Mmus.sum(axis=0).astype(bool).sum()
            
            self.total_genes_detected_across_sample_total = np.count_nonzero(anndata_array_sc)
            self.total_genes_detected_across_sample_Hsap  = np.count_nonzero(anndata_array_sc_Hsap)
            self.total_genes_detected_across_sample_Mmus  = np.count_nonzero(anndata_array_sc_Mmus)
        else:
            # If single species just calculate the single version of each of the metrics.
            self.num_cells = anndata_array_sc.shape[0]
            self.raw_reads_per_cell = self.metrics_dict["Read QC"]["reads_pre_qc"][1] / self.num_cells
            self.mean_total_counts_per_cell = np.mean(anndata_array_sc.sum(axis=1))
            self.median_total_counts_per_cell = np.median(anndata_array_sc.sum(axis=1))
            self.mean_genes_detected_per_cell = np.mean(anndata_array_sc.astype(bool).sum(axis=1))
            self.median_genes_detected_per_cell = int(np.median(anndata_array_sc.astype(bool).sum(axis=1)))

            # Get a subset of the array that doesn't contain the mito genes
            anndata_array_sc_nuc = anndata_array_sc.loc[:, ~self.anndata.var["is_mito"]]

            self.mean_nuclear_genes_detected_per_cell = np.mean(anndata_array_sc_nuc.astype(bool).sum(axis=1))
            self.median_nuclear_genes_detected_per_cell = np.median(anndata_array_sc_nuc.astype(bool).sum(axis=1))

            # Get a subset of the array that contains only the mito genes
            anndata_array_sc_mito = anndata_array_sc.loc[:, self.anndata.var["is_mito"]]

            self.mean_mito_genes_detected_per_cell = np.mean(anndata_array_sc_mito.astype(bool).sum(axis=1))
            self.median_mito_genes_detected_per_cell = int(np.median(anndata_array_sc_mito.astype(bool).sum(axis=1)))

            # Calculate percentage of counts or mitochondrial origin
            self.percentage_counts_from_mito = self.as_perc(anndata_array_sc_mito.values.sum() / anndata_array_sc.values.sum())

            self.num_unique_genes_detected_across_sample = anndata_array_sc.shape[1]

            self.total_genes_detected_across_sample = np.count_nonzero(anndata_array_sc)

    def populate_cell_stats_in_metrics_dict_mixed_species(self):
        # Same as the single_species version below except that the primary key will
        # be the base metric (see the base_stats list above) and then the secondary key will be
        # the the actual metric.
        # The first item in the tuple is the name of the metric that will be displayed
        # in the Html, the second is the value, the third is the tool tip.
        # NOTE the helper text must NOT have a comma in it as it is going to be written out in a csv.
        self.metrics_dict["num_cells"]["num_cells_total"] = ("Number of cells total", self.num_cells_total, "Estimated number of Hsap and Mmus cells combined: Number of barcodes passing the total counts threshold and the purity threshold.")
        self.metrics_dict["num_cells"]["num_cells_Hsap"] = ("Number of Hsap cells", self.num_cells_Hsap, "Estimated number of Hsap cells: Number of barcodes passing the total counts threshold and the purity threshold as Hsap.")
        self.metrics_dict["num_cells"]["num_cells_Mmus"] = ("Number of Mmus cells", self.num_cells_Mmus, "Estimated number of Mmus cells: Number of barcodes passing the total counts threshold and the purity threshold as Mmus.")
        
        self.metrics_dict["raw_reads_per_cell"]["raw_reads_per_cell_total"] = ("Raw reads per cell (Hsap and Mmus)", self.raw_reads_per_cell_total, "Number of reads pre-QC / Number of cells (Hsap and Mmus)")
        self.metrics_dict["raw_reads_per_cell"]["raw_reads_per_cell_Hsap"] = ("Raw reads per Hsap cell", self.raw_reads_per_cell_Hsap, "Number of reads pre-QC / Number of Hsap cells")
        self.metrics_dict["raw_reads_per_cell"]["raw_reads_per_cell_Mmus"] = ("Raw reads per Mmus cell", self.raw_reads_per_cell_Mmus, "Number of reads pre-QC / Number of Mmus cells")
        
        # Mean total counts per cell
        self.metrics_dict["mean_total_counts_per_cell"]["mean_total_counts_per_cell_total"] = ("Mean total counts per cell (Hsap and Mmus)", self.mean_total_counts_per_cell_total, "Mean sum of counts per cell (Hsap and Mmus cells). Only Hsap gene-derived counts are counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["mean_total_counts_per_cell"]["mean_total_counts_per_cell_Hsap"] = ("Mean total counts per Hsap cell", self.mean_total_counts_per_cell_Hsap, "Mean sum of counts per Hsap cell. Mmus gene counts are not included.")
        self.metrics_dict["mean_total_counts_per_cell"]["mean_total_counts_per_cell_Mmus"] = ("Mean total counts per Mmus cell", self.mean_total_counts_per_cell_Mmus, "Mean sum of counts per Mmus cell. Hsap gene counts are not included.")

        # Median total counts per cell
        self.metrics_dict["median_total_reads_per_cell"]["median_total_reads_per_cell_total"] = ("Median total counts per cell (Hsap and Mmus)", self.median_total_counts_per_cell_total, "Median sum of counts per cell (Hsap and Mmus cells). Only Hsap gene-derived counts are counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["median_total_reads_per_cell"]["median_total_reads_per_cell_Hsap"] = ("Median total counts per Hsap cell", self.median_total_counts_per_cell_Hsap, "Median sum of counts per Hsap cell. Mmus gene counts are not included.")
        self.metrics_dict["median_total_reads_per_cell"]["median_total_reads_per_cell_Mmus"] = ("Median total counts per Mmus cell", self.median_total_counts_per_cell_Mmus, "Median sum of counts per Mmus cell. Hsap gene counts are not included.")

        # Mean genes detected per cell
        self.metrics_dict["mean_genes_detected_per_cell"]["mean_genes_detected_per_cell_total"] = ("Mean genes detected per cell (Hsap and Mmus)", self.mean_genes_detected_per_cell_total, "Mean number of genes detected per cell (including nuclear and mitochondrial genes; including Hsap and Mmus cells).  Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["mean_genes_detected_per_cell"]["mean_genes_detected_per_cell_Hsap"] = ("Mean genes detected per Hsap cell", self.mean_genes_detected_per_cell_Hsap, "Mean number of genes detected per Hsap cell (including nuclear and mitochondrial genes). Mmus gene detections are not included.")
        self.metrics_dict["mean_genes_detected_per_cell"]["mean_genes_detected_per_cell_Mmus"] = ("Mean genes detected per Mmus cell", self.mean_genes_detected_per_cell_Mmus, "Mean number of genes detected per Mmus cell (including nuclear and mitochondrial genes). Hsap gene detections are not included.")

        # Median genes detected per cell
        self.metrics_dict["median_genes_detected_per_cell"]["median_genes_detected_per_cell_total"] = ("Median genes detected per cell (Hsap and Mmus)", self.median_genes_detected_per_cell_total, "Median number of genes detected per cell (including nuclear and mitochondrial genes; including Hsap and Mmus cells).  Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["median_genes_detected_per_cell"]["median_genes_detected_per_cell_Hsap"] = ("Median genes detected per Hsap cell", self.median_genes_detected_per_cell_Hsap, "Median number of genes detected per Hsap cell (including nuclear and mitochondrial genes). Mmus gene detections are not included.")
        self.metrics_dict["median_genes_detected_per_cell"]["median_genes_detected_per_cell_Mmus"] = ("Median genes detected per Mmus cell", self.median_genes_detected_per_cell_Mmus, "Median number of genes detected per Mmus cell (including nuclear and mitochondrial genes). Hsap gene detections are not included.")

        # Mean nuclear genes detected per cell
        self.metrics_dict["mean_nuclear_genes_detected_per_cell"]["mean_nuclear_genes_detected_per_cell_total"] = ("Mean nuclear (non-mitochondrial) genes detected per cell (Hsap and Mmus)", self.mean_nuclear_genes_detected_per_cell_total, "Mean number of nuclear genes detected per cell (Hsap and Mmus cells). Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["mean_nuclear_genes_detected_per_cell"]["mean_nuclear_genes_detected_per_cell_Hsap"] = ("Mean nuclear (non-mitochondrial) genes detected per Hsap cell", self.mean_nuclear_genes_detected_per_cell_Hsap, "Mean number of nuclear genes detected per Hsap cell. Mmus gene detections are not included.")
        self.metrics_dict["mean_nuclear_genes_detected_per_cell"]["mean_nuclear_genes_detected_per_cell_Mmus"] = ("Mean nuclear (non-mitochondrial) genes detected per Mmus cell", self.mean_nuclear_genes_detected_per_cell_Mmus, "Mean number of nuclear genes detected per Mmus cell. Hsap gene detections are not included.")

        # Median nuclear genes detected per cell
        self.metrics_dict["median_nuclear_genes_detected_per_cell"]["median_nuclear_genes_detected_per_cell_total"] = ("Median nuclear (non-mitochondrial) genes detected per cell (Hsap and Mmus)", self.median_nuclear_genes_detected_per_cell_total, "Median number of nuclear genes detected per cell (Hsap and Mmus cells). Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["median_nuclear_genes_detected_per_cell"]["median_nuclear_genes_detected_per_cell_Hsap"] = ("Median nuclear (non-mitochondrial) genes detected per Hsap cell", self.median_nuclear_genes_detected_per_cell_Hsap, "Median number of nuclear genes detected per Hsap cell. Mmus gene detections are not included.")
        self.metrics_dict["median_nuclear_genes_detected_per_cell"]["median_nuclear_genes_detected_per_cell_Mmus"] = ("Median nuclear (non-mitochondrial) genes detected per Mmus cell", self.median_nuclear_genes_detected_per_cell_Mmus, "Median number of nuclear genes detected per Mmus cell. Hsap gene detections are not included.")

        # Mean mitochondrial genes detected per cell
        self.metrics_dict["mean_mito_genes_detected_per_cell"]["mean_mito_genes_detected_per_cell_total"] = ("Mean mitochondrial genes detected per cell (Hsap and Mmus)", self.mean_mito_genes_detected_per_cell_total, "Mean number of mitochondrial genes detected per cell (Hsap and Mmus cells). Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["mean_mito_genes_detected_per_cell"]["mean_mito_genes_detected_per_cell_Hsap"] = ("Mean mitochondrial genes detected per Hsap cell", self.mean_mito_genes_detected_per_cell_Hsap, "Mean number of mitochondrial genes detected per Hsap cell. Mmus gene detections are not included.")
        self.metrics_dict["mean_mito_genes_detected_per_cell"]["mean_mito_genes_detected_per_cell_Mmus"] = ("Mean mitochondrial genes detected per Mmus cell", self.mean_mito_genes_detected_per_cell_Mmus, "Mean number of mitochondrial genes detected per Mmus cell. Hsap gene detections are not included.")

        # Median mitochondrial genes detected per cell
        self.metrics_dict["median_mito_genes_detected_per_cell"]["median_mito_genes_detected_per_cell_total"] = ("Median mitochondrial genes detected per cell (Hsap and Mmus)", self.median_mito_genes_detected_per_cell_total, "Median number of mitochondrial genes detected per cell (Hsap and Mmus cells). Mmus gene detections are not counted for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["median_mito_genes_detected_per_cell"]["median_mito_genes_detected_per_cell_Hsap"] = ("Median mitochondrial genes detected per Hsap cell", self.median_mito_genes_detected_per_cell_Hsap, "Median number of mitochondrial genes detected per Hsap cell. Mmus gene detections are not included.")
        self.metrics_dict["median_mito_genes_detected_per_cell"]["median_mito_genes_detected_per_cell_Mmus"] = ("Median mitochondrial genes detected per Mmus cell", self.median_mito_genes_detected_per_cell_Mmus, "Median number of mitochondrial genes detected per Mmus cell. Hsap gene detections are not included.")

        # Percentage of counts from mitochondrial origin
        self.metrics_dict["percentage_counts_from_mito"]["percentage_counts_from_mito_total"] = ("Percentage of counts of mitochondrial origin (Hsap and Mmus) ", self.percentage_counts_from_mito_total, "(Total mitochonrial counts / total counts) * 100; Mmus counts are not included for Hsap cells and vice versa for Mmus cells.")
        self.metrics_dict["percentage_counts_from_mito"]["percentage_counts_from_mito_Hsap"] = ("Percentage of counts of mitochondrial origin from Hsap cells", self.percentage_counts_from_mito_Hsap, "(Hsap mitochonrial counts / all Hsap counts) * 100; Mmus counts are not included.")
        self.metrics_dict["percentage_counts_from_mito"]["percentage_counts_from_mito_Mmus"] = ("Percentage of counts of mitochondrial origin from Mmus cells", self.percentage_counts_from_mito_Mmus, "(Mmus mitochonrial counts / all Mmus counts) * 100; Hsap counts are not included.")

        # Unique genes detected across samples (i.e. a gene can only be detected once per sample)
        self.metrics_dict["num_unique_genes_detected_across_sample"]["num_unique_genes_detected_across_sample_total"] = ("Unique genes detected across sample (Hsap and Mmus) ", self.num_unique_genes_detected_across_sample_total, "Number of unique genes detected across the sample (each gene can be counted only once even if found in multiple cells).")
        self.metrics_dict["num_unique_genes_detected_across_sample"]["num_unique_genes_detected_across_sample_Hsap"] = ("Unique Hsap genes detected in Hsap cells across sample", self.num_unique_genes_detected_across_sample_Hsap, "Number of unique Hsap genes detected in Hsap cells across the sample (each gene can be counted only once even if found in multiple cells).")
        self.metrics_dict["num_unique_genes_detected_across_sample"]["num_unique_genes_detected_across_sample_Mmus"] = ("Unique Mmus genes detected in Mmus cells across sample", self.num_unique_genes_detected_across_sample_Mmus, "Number of unique Mmus genes detected in Mmus cells across the sample (each gene can be counted only once even if found in multiple cells).")

        # Total genes detected across samples (i.e. genes detected per cell summed)
        self.metrics_dict["total_genes_detected_across_samples"]["total_genes_detected_across_samples_total"] = ("Total genes detected across sample (Hsap and Mmus)", self.total_genes_detected_across_sample_total, "Total number of genes detected across the sample (each gene can be counted more than once if detected in more than one cell; Both Hsap and Mmus genes are counted for each single cell).")
        self.metrics_dict["total_genes_detected_across_samples"]["total_genes_detected_across_samples_Hsap"] = ("Total Hsap genes detected across sample", self.total_genes_detected_across_sample_Hsap, "Total number of Hsap genes detected across the sample (each gene can be counted more than once if detected in more than one cell); Mmus gene detections are not included.")
        self.metrics_dict["total_genes_detected_across_samples"]["total_genes_detected_across_samples_Mmus"] = ("Total Mmus genes detected across sample", self.total_genes_detected_across_sample_Mmus, "Total number of Mmus genes detected across the sample (each gene can be counted more than once if detected in more than one cell); Hsap gene detections are not included.")

    def populate_cell_stats_in_metrics_dict_single_species(self):
        # Estimated number of cells
        self.metrics_dict["Cell metrics"]["num_cells"] = ("Number of cells", self.num_cells, "Estimated number of cells; Number of barcodes passing the total counts threshold.")

        # Raw reads per cell
        self.metrics_dict["Cell metrics"]["raw_reads_per_cell"] = ("Raw reads per cell", self.raw_reads_per_cell, "Number of reads pre-QC / Number of cells")

        # Mean total counts per cell
        self.metrics_dict["Cell metrics"]["mean_total_counts_per_cell"] = ("Mean total counts per cell", self.mean_total_counts_per_cell, "Mean sum of counts per cell.")
        # Median total counts per cell
        self.metrics_dict["Cell metrics"]["median_total_reads_per_cell"] = ("Median total counts per cell", self.median_total_counts_per_cell, "Median sum of counts per cell.")

        # Mean genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_genes_detected_per_cell"] = ("Mean genes detected per cell", self.mean_genes_detected_per_cell, "Mean number of genes detected for each cell (including nuclear and mitochondrial genes).")
        # Median genes detected per cell
        self.metrics_dict["Cell metrics"]["median_genes_detected_per_cell"] = ("Median genes detected per cell", self.median_genes_detected_per_cell, "Median number of genes detected for the cells (including nuclear and mitochondrial genes).")

        # Mean nuclear genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_nuclear_genes_detected_per_cell"] = ("Mean nuclear (non-mitochondrial) genes detected per cell", self.mean_nuclear_genes_detected_per_cell, "Mean number of nuclear genes detected for each cell.")
        # Median nuclear genes detected per cell
        self.metrics_dict["Cell metrics"]["median_nuclear_genes_detected_per_cell"] = ("Median nuclear (non-mitochondrial) genes detected per cell", self.median_nuclear_genes_detected_per_cell, "Median number of nuclear genes detected for each cell.")

        # Mean mitochondrial genes detected per cell
        self.metrics_dict["Cell metrics"]["mean_mito_genes_detected_per_cell"] = ("Mean mitochondrial genes detected per cell", self.mean_mito_genes_detected_per_cell, "Mean number of mitochondrial genes detected for each cell.")
        # Median mitochondrial genes detected per cell
        self.metrics_dict["Cell metrics"]["median_mito_genes_detected_per_cell"] = ("Median mitochondrial genes detected per cell", self.median_mito_genes_detected_per_cell, "Median number of mitochondrial genes detected for each cell.")

        # Percentage of counts from mitochondrial origin
        self.metrics_dict["Cell metrics"]["percentage_counts_from_mito"] = ("Percentage of counts of mitochondrial origin", self.percentage_counts_from_mito, "(Mitochonrial counts / all counts) * 100")

        # Unique genes detected across samples (i.e. a gene can only be detected once per sample)
        self.metrics_dict["Cell metrics"]["num_unique_genes_detected_across_sample"] = ("Unique genes detected across sample", self.num_unique_genes_detected_across_sample, "Number of unique genes detected across the sample (each gene can be counted only once even if found in multiple cells).")

        # Total genes detected across samples (i.e. genes detected per cell summed)
        self.metrics_dict["Cell metrics"]["total_genes_detected_across_samples"] = ("Total genes detected across sample", self.total_genes_detected_across_sample, "Total number of genes detected across the sample (each gene can be counted more than once if detected in more than one cell).")

    def get_cell_stats(self):
        # Try to read in the raw h5ad and handle if it is empty
        # by checking for an OSError (empty file) or a 0 barcode count.
        try:
            self.anndata = anndata.read_h5ad(sys.argv[2])
        except OSError: # If the h5ad is an empty file, output empty metrics
            self.set_cell_and_purity_stats_to_zero()
        else:
            if self.anndata.shape[0] == 0:
                self.set_cell_and_purity_stats_to_zero()
            else:
                if self.mixed:
                    # If mixed, we need to check to see if there
                    # are called cells in the count table that
                    # we can generate purity metrics from.
                    # The is_called_cell attribute is annotated in filter_count_matrix.py
                    if sum(self.anndata.obs["is_called_cell"]) > 0:
                        # Then we can populate the cell purity
                        self.calculate_purity_stats()
                    else:
                        # If there are no called cells then we set all stats to 0
                        self.set_purity_metrics_to_zero()
                    self.populate_purity_stats_in_metrics_dict()

                # Finally we need to calculate the cell stats if there
                # are cells to call the stats from
                # else set the cell stats to 0
                # The is_single_cell attribute is annotated in filter_count_matrix.py
                if sum(self.anndata.obs["is_single_cell"]) > 0:
                    # Then we can populate the cell purity
                    self.calculate_cell_stats()
                else:
                    # If there are no called cells then we set all stats to 0
                    self.set_cell_stats_to_zero()
        
        # Once we have either calculated the stats or set them to 0
        # we need to populate the metrics dicts
        # This gets done differently depending on whether we are doing a mixed
        # or single species output.
        # For the mixed species case we populate with subcategory keys for each
        # of the base base_stats.
        # Else for single we use the base_stat directly as the key
        if self.mixed:
            self.populate_cell_stats_in_metrics_dict_mixed_species()
        else:
            self.populate_cell_stats_in_metrics_dict_single_species()

    def calculate_purity_stats(self):
        """
        Calculate num_raw_cells_no_purity
        and num_multiplet_cells
        """
        # The number of cells diregarding purity
        self.num_raw_cells_no_purity = sum(self.anndata.obs["is_called_cell"])

        # The number of cells that are called cells but do not meet the purity theshold
        self.num_multiplet_cells = sum(self.anndata.obs["is_called_cell"]) - sum(self.anndata.obs["is_single_cell"])

    def populate_purity_stats_in_metrics_dict(self):
        """
        The purity stats that are calculated in calculate_purity_stats
        need to be added to the metrics_dict.
        Do something similar to what has been done in populate_cell_stats_in_metrics_dict_mixed_species.
        """
        # NOTE we call the subkey "*_total" so that the metric name
        # matches that triplicate format metrics.
        self.metrics_dict["num_raw_cells"]["num_raw_cells_total"] = (
            "Number of cells disregarding purity", self.num_raw_cells_no_purity,
            "Number of called cells (i.e. meeting the total count threshold irrespective of the purity threshold)."
            )
        
        self.metrics_dict["num_multiplet_cells"]["num_multiplet_cells_total"] = (
            "Number of multiplet cells",
            self.num_multiplet_cells,
            "Total number of multiplet cells (called cells not meeting the purity threshold)."
            )

    def set_purity_metrics_to_zero(self):
        # The number of cells diregarding purity
        self.num_raw_cells_no_purity = 0
        # The number of cells that are called cells but do not meet the purity theshold
        self.num_multiplet_cells = 0

    def set_cell_and_purity_stats_to_zero(self):
        if self.mixed:
            # If mixed we additionally need to set the purity_metrics to 0
            self.set_purity_metrics_to_zero()
        self.set_cell_stats_to_zero()

    def get_sequencing_stats(self):
        """
        Populate the self.metrics_dict with the stats
        """
        self.get_trimming_qc_stats()
        # Uncomment to reenable the mapping stats
        # self.get_mapping_stats()
        self.get_duplication_stats()

    def get_antisense(self):
        """
        Read in the antisense metric input files and pull out
        the antisense mapped reads and populate into
        self.metrics_dict

        NOTE currently we only collect the antisense for the allocated filtered
        bam. Ideally we would want this stat for each level of filtering
        I.e. raw, umr filtered and allocated
        """
        with open(sys.argv[4], "r") as antisense_handle:
            self.metrics_dict["annotated_alignment_stats"]["annotated_antisense_mapping_reads"] = ("Reads aligned antisense", int(antisense_handle.read().rstrip()), "Number of reads aligned antisense.")

    def get_duplication_stats(self):
        # Reads before deduplication
        # Reads after deduplication
        # Sequencing saturation
        with open(sys.argv[5], "r") as dedup_handle:
            for line in dedup_handle:
                if "INFO Reads: Input Reads:" in line:
                    try:
                        reads_in = int(line.split()[-1].strip())
                    except ValueError:
                        reads_in = 0
                    self.metrics_dict["Deduplication"]["reads_before_deduplication"] = ("High confidence gene-annotated reads before deduplication", reads_in, "Number of high confidence (unique alignment; max mismatch <= 3bp) gene-annotated (with XT gene_id annotation) reads before deduplication.")
                elif "INFO Number of reads out:" in line:
                    self.metrics_dict["Deduplication"]["reads_after_deduplication"] = ("High confidence gene-annotated reads after deduplication", int(line.split()[-1].strip()), "Number of high confidence (unique alignment: max mismatch <= 3bp) gene-annotated (with XT gene_id annotation) reads after deduplication.")
        if reads_in != 0:
            self.metrics_dict["Deduplication"]["sequencing_saturation"] = ("Sequencing saturation", self.as_perc(1 - (self.metrics_dict["Deduplication"]["reads_after_deduplication"][1]/self.metrics_dict["Deduplication"]["reads_before_deduplication"][1])), "(1 - (Reads after deduplication / reads before deduplication)) * 100")
        else:
            self.metrics_dict["Deduplication"]["sequencing_saturation"] = ("Sequencing saturation", 0.0, "(1 - (Reads after deduplication / reads before deduplication)) * 100")
    
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
                        val_percent = float(re.search(r"(?<=\()([\d]+\.{0,1}[\d]*)(?=%\))", lines[i]).groups()[0].rstrip())
                        self.metrics_dict[category][f"{header}"] = (header.replace("_", " ").capitalize(), val_absolute, header.replace("_", " ").capitalize())
                        self.metrics_dict[category][f"{header}_perc"] = (header.replace("_", " ").capitalize() + " percentage", val_percent, header.replace("_", " ").capitalize() + " percentage")

    @staticmethod
    def as_perc(float_to_convert):
        return float_to_convert * 100
        # return f"{float_to_convert:.2f}%"

    def get_trimming_qc_stats(self):
        # Reads pre-QC
        reads_pre_qc = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1"]["before_filtering_total_reads"])
        self.metrics_dict["Read QC"]["reads_pre_qc"] = ("Number of reads pre-QC", reads_pre_qc, "Number of reads in the input R1 fastq files (after merging if applicable).")
        # Reads containing cellular barcode matching whitelist
        self.metrics_dict["Read QC"]["valid_barcode_reads"] = ("Number of valid barcode-containing reads", int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.io_extract"]["before_filtering_total_reads"]), "Number of reads containing a barcode exactly matching the whitelist.")
        # Percentage of reads containing cellular barcode matching whitelist as percentage of pre-QC reads
        if reads_pre_qc != 0:
            self.metrics_dict["Read QC"]["valid_barcode_reads_perc"] = ("Percentage valid barcode-containing reads", self.as_perc(float(self.metrics_dict["Read QC"]["valid_barcode_reads"][1] / self.metrics_dict["Read QC"]["reads_pre_qc"][1])), "(Number of valid barcode-containing reads / Number of reads pre-QC) * 100.")
        else:
            self.metrics_dict["Read QC"]["valid_barcode_reads_perc"] = ("Percentage valid barcode-containing reads", 0, "Number of valid barcode-containing reads / Number of reads pre-QC * 100.")
        # Percentage of barcode bases >= Q30
        self.metrics_dict["Read QC"]["barcode_bases_q30_perc"] = ("Barcode bp >= Q30 percentage", self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R2"]["after_filtering_q30_rate"])), "The percentage of the barcode bases with a Phred score >= 30.")
        reads_post_qc = int(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_total_reads"])
        # Reads after polyX tail and polyA internal trimming
        self.metrics_dict["Read QC"]["reads_post_trimming"] = ("Number of reads post-QC trimming", reads_post_qc, "Number of reads after polyX tail and polyA internal trimming.")
        # Reads after polyX tail and polyA internal trimming as percentage of valid barcode reads
        if self.metrics_dict["Read QC"]["valid_barcode_reads"][1] != 0:
            self.metrics_dict["Read QC"]["reads_post_trimming_perc"] = ("Percentage reads post-QC trimming", self.as_perc(float(reads_post_qc/self.metrics_dict["Read QC"]["valid_barcode_reads"][1])), "(Number of reads after polyX tail and polyA internal trimming / Number of valid barcode-containing reads) * 100.")
        else:
            self.metrics_dict["Read QC"]["reads_post_trimming_perc"] = ("Percentage reads post-QC trimming", 0, "(Number of reads after polyX tail and polyA internal trimming / Number of valid barcode-containing reads) * 100.")

        # Mean read length after polyX tail and polyA internal trimming
        self.metrics_dict["Read QC"]["mean_post_trim_read_length"] = ("Mean read length post-QC trimming", float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_read1_mean_length"]), "Mean R1 read length post-QC trimming.")
        # Percentage of bases post trimming >= Q30
        self.metrics_dict["Read QC"]["rna_bases_q30_perc"] = ("R1 bp >= Q30 percentage; post-QC trimming", self.as_perc(float(self.multiqc_json_dict["report_general_stats_data"][self.qc_key_to_index_dict["fastp"]][f"{self.sample_id}_R1.polyA"]["after_filtering_q30_rate"])), "The percentage of the R1 bases (post-QC trimming) with a Phred score >= 30.")


if __name__ == "__main__":
    SummaryStatistics().generate_metrics()
