#!/usr/bin/env nextflow

/*
* This script contains the workflow for the daily R&D pipeline.
* The processes used in this workflow live in modules/processes.nf
*/

nextflow.enable.dsl=2
include {
  features_file; merge_lanes; merged_fastp; io_extract; io_extract_fastp;
  trim_extra_polya; post_polyA_fastp; star; raw_qualimap; filter_umr_mismatch; filtered_qualimap;
  feature_counts; filter_for_annotated; annotated_qualimap; multiqc;
  sort_index_bam; dedup; io_count; count_matrix;
  filter_count_matrix; cell_caller; summary_statistics; single_summary_report;
  multi_sample_report
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

    // Create path objects to HTML report templates
    single_sample_report_template = file("templates/single_sample_report_template.html.jinja2")
    multi_sample_report_template = file("templates/multi_sample_report_template.html.jinja2")
    cs_logo = file("templates/csgenetics_logo.png")

    // Create the whitelist object
    whitelist = file("${params.whitelist_path}")

    // Make a file object of the STAR dir
    star_index = file("${params.star_index_dir}")

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
          .map { [it[0], it[2], it[3]] }


    // Set the barcode_pattern 
    barcode_pattern="CCCCCCCCCCCCC"
    merged_fastp(ch_merge_lanes_out_filtered, barcode_pattern)
    ch_merged_fastp_multiqc = merged_fastp.out.merged_fastp_multiqc

    // Get the whitelist and extract the IOs from the fastqs using umitools
    io_extract(ch_merge_lanes_out_filtered, whitelist, barcode_pattern)
    ch_io_extract_out = io_extract.out.io_extract_out
    // Filter out empty fastq files.
    ch_io_extract_out_filtered = ch_io_extract_out
      .filter { it[2].countFastq() > 0}

    ch_io_extract_log = io_extract.out.io_extract_log

    // Trim and remove low quality reads with fastp
    io_extract_fastp(ch_io_extract_out_filtered)
    ch_io_extract_fastp_out = io_extract_fastp.out.fastp_out
    ch_io_extract_fastp_multiqc = io_extract_fastp.out.fastp_multiqc

    // Trim extra polyA
    trim_extra_polya(ch_io_extract_fastp_out)
    ch_trim_extra_polya_out = trim_extra_polya.out.trim_extra_polya_out

    // Send the polyA trimmed reads back through
    // fastp to get total number post QC reads
    // and Q30 percentages.
    post_polyA_fastp(ch_trim_extra_polya_out)
    ch_post_polyA_fastp_out = post_polyA_fastp.out.fastp_out
    ch_post_polyA_fastp_multiqc = post_polyA_fastp.out.fastp_multiqc

    // Align with STAR
    star(ch_post_polyA_fastp_out, star_index)
    ch_star_out_qualimap = star.out.star_out_bams
    ch_star_out_filtering = star.out.star_out_bams
    ch_star_multiqc = star.out.star_multiqc

    // Qualimap on STAR output
    raw_qualimap(ch_star_out_qualimap, gtf)

    // Filter the mapped reads for reads with 1 alignment and max 3 mismatch
    filter_umr_mismatch(ch_star_out_filtering)

    // Qualimap on filtered bam
    filtered_qualimap(filter_umr_mismatch.out.filtered_bam, gtf)

    // Perform featurecount quantification
    feature_counts(filter_umr_mismatch.out.filtered_bam, gtf)

    // Filter the annotated featureCounts bam
    // for only annotated/assigned reads with 1
    // target
    filter_for_annotated(feature_counts.out.feature_counts_out_bam)

    // Produce qualimap output of the annotated bam
    // for metrics
    annotated_qualimap(filter_for_annotated.out.annotated_bam, gtf)

    // Generate input channel containing all the files needed for multiqc per samples. 
    // The final channel structure is [sample_id, [file1, file2, file3, ...]]
    // NB fastqc currently does not integrate multip qualimap
    // outputs successfuly so they are manually carried through
    // to the summary_statistics
    ch_merged_fastp_multiqc
      .mix(ch_post_polyA_fastp_multiqc)
      .mix(ch_io_extract_fastp_multiqc)
      .mix(ch_star_multiqc)
      .mix(feature_counts.out.feature_counts_multiqc)
      .groupTuple(by:0, size: 5)
      .map({it.flatten()}).map({[it[0], it.tail()]})
      .set { ch_multiqc_in }

    // Run multiqc  
    multiqc(ch_multiqc_in)
    ch_multiqc_json = multiqc.out.multiqc_json    
    
    // Sort and index bam file
    sort_index_bam(filter_for_annotated.out.annotated_bam)

    // Perform deduplication
    dedup(sort_index_bam.out.sort_index_bam_out)
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
    // [sample_id, min_nuc_gene_cutoff, h5ad, annotated_qualimap,
    // antisense, dedup.log, filtered_qualimap, multiqc_data, raw_qualimap]
    ch_filtered_count_matrices
    .mix(ch_multiqc_json)
    .mix(sort_index_bam.out.antisense_out)
    .mix(raw_qualimap.out.qualimap_txt)
    .mix(filtered_qualimap.out.qualimap_txt)
    .mix(annotated_qualimap.out.annotated_qualimap)
    .mix(ch_io_dedup_log)
    .groupTuple(by:0, size: 7, sort:{it.name})
    .mix(ch_cell_caller_out)
    .groupTuple(by:0, size:2, sort:{sort:{it.getClass() == nextflow.util.ArrayBag ? 1 : 0}})
    .map({it.flatten()})
    .set({ch_summary_report_in})
    
    // Generate summary statistics
    summary_statistics(ch_summary_report_in)

    ch_summary_metrics_and_plot = summary_statistics.out.metrics_csv.join(ch_cell_caller_plot, by:0)

    // Generate single sample report
    single_summary_report(ch_summary_metrics_and_plot, single_sample_report_template, cs_logo)

    // // Generate multi sample report
    // multi_sample_report(single_summary_report.out.single_sample_metric_out.collect(), multi_sample_report_template)
 
}
