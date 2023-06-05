#!/usr/bin/env nextflow

/*
* This script contains the workflow for the daily R&D pipeline.
* The processes used in this workflow live in modules/processes.nf
*/

nextflow.enable.dsl=2
include { basespace;  features_file; merge_lanes; fastqc; barcode;io_whitelist; io_extract; fastp; trim_extra_polya; star; qualimap ;feature_counts; multiqc; sort_index_bam; group; dedup; io_count; count_matrix; cell_caller;summary_report } from './modules/processes.nf'


workflow {

    // Create path object to the GTF
    gtf = file("${params.gtf_path}/*.gtf")

    // We can either run the pipeline using fastqs in basespace as input, 
    // or we can use fastqs in s3 as input.
    // Start by initialising an empty input channel
    ch_input_fastqs = Channel.empty()

    // Either fill the channel with fastqs from the provided s3 path, or run
    // the basespace process to obtain the fastqs
    if ( params.fastq ){
      ch_input_fastqs = Channel.fromPath("${params.fastq_path}/*.f*q.gz")
      .ifEmpty { exit 1, "Input *.f*q.gz files are not found in ${params.fastq_path}" }
    }
    else{
      ch_input_fastqs = basespace()
    }

    // Create feature file for count_matrix from GTF
    features_file(gtf)
    feature_file_out = features_file.out.modified_gtf

    // Regardless of where the fastqs came from, exract the sample name from the
    // fastq filenames and add it to the channel, 
    // so each item in the channel is a tuple of [sample name, fastq file].
    ch_input_fastqs
    .flatMap()
    .map { [it.baseName.toString().substring(0, it.baseName.toString().lastIndexOf("_L0")), it] }
    .groupTuple()
    .set { ch_merge_lanes_in }

    // This process will merge fastqs split over multiple lanes 
    // and count the number of reads in the merged fastq
    merge_lanes(ch_merge_lanes_in)
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

    // This process generates an inferred whitelist based on IOs detected in R2
    // NOT run by default. Switch it on by setting params.io_whitelist = true
    io_whitelist(ch_merge_lanes_out_filtered, barcode_pattern)
    ch_io_whitelist_log = io_whitelist.out.io_whitelist_log

    // Get the whitelist and extract the IOs from the fastqs using umitools
    ch_whitelist = Channel.fromPath(params.whitelist_path)
    io_extract(ch_merge_lanes_out_filtered, ch_whitelist.collect(), barcode_pattern)
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
    ch_star_index = Channel.fromPath("${params.genome_path}")
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


    // Generate input channel containing all the files needed for multiqc across all samples. 
    // The final channel structure is just [file1,file2,file3,...]
    ch_fastqc_multiqc
      .mix(ch_barcode_multiqc)
      .mix(ch_fastp_multiqc)
      .mix(ch_star_multiqc)
      .mix(ch_feature_counts_multiqc)
      .transpose()
      .collect()
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
    count_matrix(ch_io_count_out, ch_whitelist.collect(), feature_file_out)
    ch_h5ad = count_matrix.out.h5ad
    ch_raw_matrix = count_matrix.out.raw_matrix
    ch_raw_barcodes = count_matrix.out.raw_barcodes
    ch_raw_features = count_matrix.out.raw_features

    // Run cell caller
    cell_caller(ch_h5ad)
    ch_cell_caller_out = cell_caller.out.cell_caller_out
    ch_cell_caller_plot = cell_caller.out.cell_caller_plot

    // Generate Report
    summary_report(ch_raw_matrix, ch_raw_barcodes, ch_raw_features, ch_multiqc_json.collect(), ch_antisense_out,ch_qualimap_txt.collect(), ch_cell_caller_plot, ch_cell_caller_out)
    ch_summary_report = summary_report.out.report_html

}
