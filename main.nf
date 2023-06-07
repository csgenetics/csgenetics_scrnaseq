#!/usr/bin/env nextflow

/*
* This script contains the workflow for the daily R&D pipeline.
* The processes used in this workflow live in modules/processes.nf
*/

nextflow.enable.dsl=2
include {
  features_file; merge_lanes; fastqc; barcode; io_extract; fastp;
  trim_extra_polya; star; qualimap; feature_counts; multiqc;
  sort_index_bam; group; dedup; io_count; count_matrix;
  filter_count_matrix; cell_caller; summary_report;experiment_report
  } from './modules/processes.nf'

workflow {
  // Create channels from rows in CSV file
  // tuples are grouped by sample_id so FQs from different lanes may be merged
    Channel
      .fromPath(params.input_csv)
      .ifEmpty { exit 1, "Cannot find input file : ${params.input_csv}" }
      .splitCsv(header:true, sep:',', strip: true)
      .map {row -> [ row.sample_id, file(row.fastq_1), file(row.fastq_2)] }
      .groupTuple(by: 0)
      .set {ch_input}

    // Create path object to the GTF
    gtf = file("${params.gtf_path}")

    // Create the whitelist object
    whitelist = file("${params.whitelist_path}")

    // Create feature file for count_matrix from GTF
    features_file(gtf)
    feature_file_out = features_file.out.modified_gtf

    // This process will merge fastqs split over multiple lanes 
    // and count the number of reads in the merged fastq
    merge_lanes(ch_input)
    ch_merge_lanes_out = merge_lanes.out.merge_lanes_out
    ch_numreads_log = merge_lanes.out.numreads_log

    // Remove fastqs with less than params.depth_min reads
    // and generate a channel where each items is a list of [sample id, filename]

    ch_merge_lanes_out_filtered = ch_merge_lanes_out
          .filter { (it[1].toInteger() >= params.depth_min) }
          .map { [it[0], it[2]] }

    // Generate a fastqc report
    fastqc(ch_merge_lanes_out_filtered)
    ch_fastqc_multiqc = fastqc.out.fastqc_multiqc
    ch_fastqc_out = fastqc.out.fastqc_out

    // Set the barcode_pattern 
    barcode_pattern="CCCCCCCCCCCCC"
    barcode(ch_merge_lanes_out_filtered, barcode_pattern)
    ch_barcode_multiqc = barcode.out.barcode_multiqc

    // Get the whitelist and extract the IOs from the fastqs using umitools
    io_extract(ch_merge_lanes_out_filtered, whitelist, barcode_pattern)
    ch_io_extract_out = io_extract.out.io_extract_out
    ch_io_extract_log = io_extract.out.io_extract_log

    // Trim and remove low quality reads with fastp
    fastp(ch_io_extract_out)
    ch_fastp_out = fastp.out.fastp_out
    ch_fastp_multiqc = fastp.out.fastp_multiqc
    ch_fastp_log = fastp.out.fastp_log

    // Trim extra polyA
    trim_extra_polya(ch_fastp_out)
    ch_trim_extra_polya_out = trim_extra_polya.out.trim_extra_polya_out
    ch_trim_extra_polya_log1 = trim_extra_polya.out.trim_extra_polya_log1
    ch_trim_extra_polya_log2 = trim_extra_polya.out.trim_extra_polya_log2

    // Align with STAR
    ch_star_index = Channel.fromPath("${params.star_index_dir}")
    star(ch_trim_extra_polya_out, ch_star_index.collect())
    ch_star_out = star.out.star_out
    ch_star_multiqc = star.out.star_multiqc

    // Qualimap on STAR output 
    qualimap(ch_star_out, gtf)
    ch_qualimap_txt = qualimap.out.qualimap_txt 


    // Perform featurecount quantification
    feature_counts(ch_star_out, gtf)
    ch_feature_counts_out = feature_counts.out.feature_counts_out
    ch_feature_counts_multiqc = feature_counts.out.feature_counts_multiqc

    // Generate input channel containing all the files needed for multiqc per samples. 
    // The final channel structure is [sample_id, file1, file2, file3, ...]
    ch_fastqc_multiqc
      .mix(ch_barcode_multiqc)
      .mix(ch_fastp_multiqc)
      .mix(ch_star_multiqc)
      .mix(ch_feature_counts_multiqc)
      .groupTuple(by:0, size: 5)
      .map({it.flatten()}).map({[it[0], it.tail()]})
      .set { ch_multiqc_in }

    // Run multiqc  
    multiqc(ch_multiqc_in)
    ch_multiqc_json = multiqc.out.multiqc_json    
    
    // Run Sort and index bam file
    sort_index_bam(ch_feature_counts_out)
    ch_antisense_out = sort_index_bam.out.antisense_out
    ch_sort_index_bam_out = sort_index_bam.out.sort_index_bam_out 

    // Run umi tools group (basically assign reads to cells) and save a lot of output channels
    group(ch_sort_index_bam_out)
    ch_io_group_sam = group.out.io_group_sam
    ch_group_filtered_sam = group.out.io_group_filtered_sam
    ch_group_tsv = group.out.io_group_tsv
    ch_io_group_log = group.out.io_group_log


    // Perform deduplication
    dedup(ch_group_filtered_sam)
    ch_io_dedup_log = dedup.out.io_dedup_log
    ch_io_dedup_sam = dedup.out.io_dedup_sam

    // Generate file for count matrix
    io_count(ch_io_dedup_sam)
    ch_io_count_out = io_count.out.io_count_out
    ch_io_goodumr_count = io_count.out.io_goodumr_count


    // Generate raw count matrix
    count_matrix(ch_io_count_out, whitelist, feature_file_out)
    ch_h5ad = count_matrix.out.h5ad

    // Run cell caller
    cell_caller(ch_h5ad)
    ch_cell_caller_out = cell_caller.out.cell_caller_out //[val(sample_id), int(cell_caller_nuc_gene_threshold)]
    ch_cell_caller_plot = cell_caller.out.cell_caller_plot


    // Sort the groupTuple so that the int is always
    // first and then flatten the tuple list to return a 3mer
    ch_filter_count_matrix_in = ch_cell_caller_out.mix(ch_h5ad)
    .groupTuple(by: 0, size:2, sort:{it.getClass() == sun.nio.fs.UnixPath ? 1 : 0})
    .map{[it[0], it[1][0], it[1][1]]}

    // Output filtered (cells only) count tables
    filter_count_matrix(ch_filter_count_matrix_in)
    ch_filtered_count_matrices = filter_count_matrix.out.cell_only_count_matrix

    // structure of ch_summary_report_in is
    // [sample_id, min_nuc_gene_cutoff, barcodes, features, matrix, multiqc_data, antisense, cell_caller_png, qualimap]
    ch_filtered_count_matrices.map({[it[0], it.tail()]}).transpose()
    .mix(ch_multiqc_json)
    .mix(ch_antisense_out)
    .mix(ch_qualimap_txt)
    .mix(ch_cell_caller_plot)
    .groupTuple(by:0, size: 7, sort:{it.name})
    .mix(ch_cell_caller_out)
    .groupTuple(by:0, size:2).map({it.flatten()})
    .set({ch_summary_report_in})
    
    // Generate Report
    summary_report(ch_summary_report_in)
    // ch_summary_report = summary_report.out.report_html


    ch_metrics_csv = summary_report.out.metrics_csv
 
    // Experiment Report
    experiment_report(ch_metrics_csv.collect())
}
