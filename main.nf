#!/usr/bin/env nextflow

/*
* A Nextflow pipeline for processing scRNA-Seq data generated using CS Genetics'
* single-cell kit to produce a genes by barcode count table.
*/

nextflow.enable.dsl=2
include {
  features_file; merge_lanes; merged_fastp; io_extract; io_extract_fastp;
  trim_extra_polya; post_polyA_fastp; star;
  create_valid_empty_bam as create_valid_empty_bam_star;
  run_qualimap as raw_qualimap; run_qualimap as filtered_qualimap; run_qualimap as annotated_qualimap;
  filter_umr_mismatch; feature_counts; filter_for_annotated; multiqc;
  sort_index_bam; dedup; io_count; count_matrix;
  filter_count_matrix; cell_caller; summary_statistics; single_summary_report;
  multi_sample_report
  } from './modules/processes.nf'

def order_integer_first(it){
  try{
        // Will raise exception if not int
        it.isInteger()
        0
      } catch(MissingMethodException e1){
        // path will end here and therefore return 1
        1
      }
}

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
    single_sample_report_template = file(params.single_sample_report_template)
    multi_sample_report_template = file(params.multi_sample_report_template)

    // Create the whitelist object
    whitelist = file(params.whitelist_path)

    // Make a file object of the STAR dir
    star_index = file(params.star_index_dir)

    // Create feature file for count_matrix from GTF
    features_file(gtf)
    feature_file_out = features_file.out.modified_gtf

    // Create empty qualimap output template path object
    empty_qualimap_template = file(params.empty_qualimap_template)

    // This process will merge fastqs split over multiple lanes 
    // and count the number of reads in the merged fastq
    merge_lanes(ch_input)
    ch_merge_lanes_out = merge_lanes.out.merge_lanes_out

    // Set the barcode_pattern 
    barcode_pattern="CCCCCCCCCCCCC"
    merged_fastp(ch_merge_lanes_out, barcode_pattern)
    ch_merged_fastp_multiqc = merged_fastp.out.merged_fastp_multiqc

    // Get the whitelist and extract the IOs from the fastqs using umitools
    io_extract(ch_merge_lanes_out, whitelist, barcode_pattern)
    ch_io_extract_out = io_extract.out.io_extract_out
    ch_io_extract_log = io_extract.out.io_extract_log

    // Trim and remove low quality reads with fastp
    io_extract_fastp(ch_io_extract_out)
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

    // Filter for empty fastq
    // Pipe good to STAR
    // Pipe empty to create_valid_empty_bam_star
    trim_extra_polya.out.trim_extra_polya_out
          .branch { 
            good_fastq: it[1].countFastq() > 0
            empty_fastq: it[1].countFastq() == 0
            }
          .set{polyA_out_ch}

    // Align the good fastqs with STAR
    star(polyA_out_ch.good_fastq, star_index)

    // Create an empty one-line-header bam for the empty fastq samples
    // and send this into each of the qualimap process
    // and dedup (i.e. skip )
    create_valid_empty_bam_star(polyA_out_ch.empty_fastq.map({[it[0], it[1], "_Aligned.sortedByCoord.out"]}))

    // Qualimap on STAR output
    // Annotate the channel objects with dummy counts of either 1 or 0
    // depending on whether the bams are empty are not to either run
    // qualimap or populate an empty qualimap template
    raw_qualimap(star.out.out_bam.map({[it[0], it[1], 1]}).mix(create_valid_empty_bam_star.out.out_bam.map({[it[0], it[1], 0]})), gtf, empty_qualimap_template, "raw")

    // Filter the mapped reads for reads with 1 alignment and max 3 mismatch
    filter_umr_mismatch(star.out.out_bam.mix(create_valid_empty_bam_star.out.out_bam))

    // Qualimap on filtered bam
    filtered_qualimap(filter_umr_mismatch.out.filtered_bam, gtf, empty_qualimap_template, "filtered")

    // Perform featurecount quantification
    feature_counts(filter_umr_mismatch.out.filtered_bam, gtf)

    // Filter the annotated featureCounts bam
    // for only annotated/assigned reads with 1 target
    filter_for_annotated(feature_counts.out.out_bam)

    // Produce qualimap output of the annotated bam for metrics
    annotated_qualimap(filter_for_annotated.out.annotated_bam, gtf, empty_qualimap_template, "annotated")

    // Generate input channel containing all the files needed for multiqc per samples. 
    // The final channel structure is [sample_id, [file1, file2, file3, ...]]
    // NB fastqc currently does not integrate multip qualimap
    // outputs successfuly so they are manually carried through
    // to the summary_statistics
    ch_merged_fastp_multiqc
      .mix(ch_post_polyA_fastp_multiqc)
      .mix(ch_io_extract_fastp_multiqc)
      .groupTuple(by:0, size: 3)
      .map({it.flatten()}).map({[it[0], it.tail()]})
      .set { ch_multiqc_in }

    // Run multiqc  
    multiqc(ch_multiqc_in)
    
    // Sort and index bam file
    sort_index_bam(filter_for_annotated.out.annotated_bam)

    // Perform deduplication
    dedup(sort_index_bam.out.sort_index_bam_out)

    // Generate file for count matrix
    io_count(dedup.out.io_dedup_sam)
    ch_io_count_out = io_count.out.io_count_out

    // Generate raw count matrix
    count_matrix(ch_io_count_out, whitelist, feature_file_out)
    ch_h5ad = count_matrix.out.h5ad

    // Run cell caller
    cell_caller(ch_h5ad)
    ch_cell_caller_out = cell_caller.out.cell_caller_out //[val(sample_id), int(cell_caller_nuc_gene_threshold)]

    // Sort the groupTuple so that the int is always
    // first and then flatten the tuple list to return a 3mer
    // N.B. We were originally sorting by class (sort:{it.getClass() == sun.nio.fs.UnixPath ? 1 : 0})
    // But for some reason this only worked locally and not on NextFlow tower
    ch_filter_count_matrix_in = ch_cell_caller_out.mix(ch_h5ad)
    .groupTuple(by: 0, size:2, sort:{order_integer_first(it)})
    .map{[it[0], it[1][0], it[1][1]]}

    // Output filtered (cells only) count tables
    filter_count_matrix(ch_filter_count_matrix_in)

    // structure of ch_summary_report_in is
    // [sample_id, min_nuc_gene_cutoff, h5ad, annotated_qualimap,
    // antisense, dedup.log, filtered_qualimap, multiqc_data, raw_qualimap]
    filter_count_matrix.out.cell_only_count_matrix
    .mix(multiqc.out.multiqc_json)
    .mix(sort_index_bam.out.antisense_out)
    .mix(raw_qualimap.out.qualimap_txt)
    .mix(filtered_qualimap.out.qualimap_txt)
    .mix(annotated_qualimap.out.qualimap_txt)
    .mix(dedup.out.io_dedup_log)
    .groupTuple(by:0, size: 7, sort:{it.name})
    .mix(ch_cell_caller_out)
    .groupTuple(by:0, size:2, sort:{sort:{order_integer_first(it)}})
    .map({it.flatten()})
    .set({ch_summary_report_in})
    
    // Generate summary statistics
    summary_statistics(ch_summary_report_in)

    ch_summary_metrics_and_plot = summary_statistics.out.metrics_csv.join(cell_caller.out.cell_caller_plot, by:0)

    // Generate single sample report
    single_summary_report(ch_summary_metrics_and_plot, single_sample_report_template)

    // Generate multi sample report
    multi_sample_report(single_summary_report.out.single_sample_metric_out.collect(), multi_sample_report_template)
 
}
