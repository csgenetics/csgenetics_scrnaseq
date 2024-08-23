#!/usr/bin/env nextflow

/*
* A Nextflow pipeline for processing scRNA-Seq data generated using CS Genetics'
* single-cell kit to produce a genes by barcode count table.
*/

nextflow.enable.dsl=2
include {
  download_star_index; download_gtf; download_input_csv; download_barcode_list; download_public_fastq;
  features_file; merge_lanes; merged_fastp; io_extract; io_extract_fastp;
  trim_extra_polya; post_polyA_fastp; star;
  create_valid_empty_bam as create_valid_empty_bam_star;
  run_qualimap as raw_qualimap; run_qualimap as annotated_qualimap;
  gtf2bed; run_rseqc as raw_rseqc; run_rseqc as annotated_rseqc;
  feature_counts; multiqc;
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

  // Some users have issues accessing the S3 resources we have hosted publicly in our
  // s3://csgx.public.readonly bucket due to their local AWS configurations.
  // One workaround is to use anonymous access to the S3 resources
  // using aws.client.anonymous = true configuration.
  // However, this may then prevent the users from accessing their own S3 resources.
  // To solve this issue, to access public S3 resources, we will run a process that 
  // uses the AWS CLI but with no user configuration using the --no-sign-request flag.
  // We will perform this for all resources starting with 's3://csgx.public.readonly'
  // The following paths need to be checked:
  // params.star_index
  // params.gtf
  // params.input_csv
  
  // Check whether params.star_index starts with s3://csgx.public.readonly
  // and if it does, download the file in a process and set the star_index to the downloaded file
  if (params.star_index.startsWith("s3://csgx.public.readonly")){
    star_index = download_star_index()
  } else {
    star_index = file(params.star_index)
  }
  
  // Check whether params.gtf starts with s3://csgx.public.readonly
  // and if it does, download the file in a process and set the gtf to the downloaded file
  if (params.gtf.startsWith("s3://csgx.public.readonly")){
    gtf = download_gtf()
  } else {
    gtf = file(params.gtf)
  }

  // Check whether params.input_csv starts with s3://csgx.public.readonly
  // and if it does, download the file in a process and set the input_csv to the downloaded file
  if (params.input_csv.startsWith("s3://csgx.public.readonly")){
    input_csv = download_input_csv()
  } else {
    input_csv = Channel.fromPath(params.input_csv)
  }
  
  // Finally, we need to check the fastq files to see if they
  // are hosted in the s3://csgx.public.readonly bucket.
  // If so, download via a process, else directly make a file object
  // so that locally configured AWS credentials can be used.
  // We make the assumption that if R1 and R2 for each row are hosted in the same bucket.
  input_csv
    .splitCsv(header:true, sep:',', strip: true)
    .map {row ->
        // To comply with nf-core standards, we will move from working with 'sample_id'
        // as the first column, to 'sample'. To maintain reverse compatibility with
        // e.g. the test resource input_csv, we will make it so that the objects of the row
        // are accessed by index rather than specific key (column name).
        // The .splitCsv method returns a map object, so we will convert it to a list of its values.
        def rowList = row.values().toList()
        def fastq_1 = rowList[1]
        def fastq_2 = rowList[2]
        def download = fastq_1.startsWith('s3://csgx.public.readonly') && fastq_2.startsWith('s3://csgx.public.readonly') ? 'download' : 'no_download'
        return [ rowList[0], fastq_1, fastq_2, download]
    }
    .branch {
        download: it[3] == 'download'
        no_download: it[3] == 'no_download'
    }
    .set { split_ch }

  split_ch.download
    .map { it[0..2] } // remove the 'download' string from the tuple
    .set { to_download_ch }

  split_ch.no_download
    .map { it[0..2] } // remove the 'no_download' string from the tuple
    .map { [it[0], file(it[1]), file(it[2])] } // convert the paths to File objects
    .set { no_download_ch }

  // Download the fastq files from the s3://csgx.public.readonly bucket
  download_public_fastq(to_download_ch)

  download_public_fastq.out.downloaded_fastqs
    .mix(no_download_ch)
    .groupTuple(by:0)
    .set { ch_input }

  // The following paths will always need to be downloaded from the s3://csgx.public.readonly bucket:
  // params.barcode_list_path
  barcode_list = download_barcode_list()

  // Create path objects to HTML report templates
  single_sample_report_template = file("${baseDir}/templates/single_sample_report_template.html.jinja2")
  multi_sample_report_template = file("${baseDir}/templates/multi_sample_report_template.html.jinja2")
  // Create empty qualimap output template path object
  empty_qualimap_template = file("${baseDir}/templates/empty_qualmap.txt")

  // Create feature file for count_matrix from GTF
  features_file(gtf)
  feature_file_out = features_file.out.modified_gtf

  // ch_input is in the form: [sample, [R1_L001, R1_L002], [R2_L001, R2_L002]],  [sample, [R1_L001], [R2_L001]]
  // The number of items in each of the sub arrays depends on how many lanes the sample was sequenced over.
  // We only want to pass those samples that have multiple lanes worth of fastqs into merge_lanes.
  // For those samples that do not have multiple lanes, we will skip merge lanes and merge with the output
  // of merge_lanes.
  
  // Split ch_input base on whether there are >1 R1 fastqs
  ch_input.branch{
    multiple_lanes: it[1].size() > 1
    single_lane: it[1].size() == 1
  }.set{
    ch_input_split
  }

  // After identifying the multiple lanes samples
  // We will split these futher so that they are in the form:
  // [sample, R1, [R1_L001, R1_L002]], [sample, R2, [R2_L001, R2_L002]]
  // This does 2 things. It allows us to futher parallelize the merge_lanes process,
  // and it means we don't have to rely on R1 or R2 being in the file names.
  // E.g. some systems use names like *r_1* and *r_2*.
  merge_lanes_in = ch_input_split.multiple_lanes
    .flatMap{[[it[0], "R1", it[1]], [it[0], "R2", it[2]]]}  

  // This process will merge fastqs split over multiple lanes 
  merge_lanes(merge_lanes_in)
  ch_merge_lanes_out = merge_lanes.out.merge_lanes_out

  // Combine the R1 and R2 reads per sample ensuring R1 comes first
  // --> [sample, *.merged.R1*, *.merged.R2*]
  ch_merge_lanes_out_merged = ch_merge_lanes_out.groupTuple(by:0, size:2).map({it[1][0] == "R1" ? [it[0], it[2][0], it[2][1]] : [it[0], it[2][1], it[2][0]]})

  // Flatten the R1 and R2 in the non-merged fastq pairs
  // [sample, [R1], [R2]] --> [sample, R1, R2]
  ch_input_split_single_lane_flattened = ch_input_split.single_lane
    .map{[it[0], it[1][0], it[2][0]]}

  // Merge the merged and non-merged fastqs
  io_extract_in_ch = ch_merge_lanes_out_merged.mix(ch_input_split_single_lane_flattened)

  // Set the barcode_pattern
  // TODO run separate taks of merged_fastp for each of the R1 and R2 files
  // to increase parallelization.
  barcode_pattern="CCCCCCCCCCCCC"
  merged_fastp(io_extract_in_ch, barcode_pattern)
  ch_merged_fastp_multiqc = merged_fastp.out.merged_fastp_multiqc

  // Get the barcode_list and extract the IOs from the fastqs using umitools
  io_extract_script = "${baseDir}/bin/io_extract.awk"
  io_extract(io_extract_in_ch, barcode_list, io_extract_script)
  ch_io_extract_out = io_extract.out.io_extract_out

  // Trim and remove low quality reads with fastp
  io_extract_fastp(ch_io_extract_out)
  ch_io_extract_fastp_out = io_extract_fastp.out.fastp_out
  ch_io_extract_fastp_multiqc = io_extract_fastp.out.fastp_multiqc

  // Trim extra polyA
  trim_polyA_script = file("${baseDir}/bin/trim_poly_A.awk")
  trim_extra_polya(ch_io_extract_fastp_out, trim_polyA_script)
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

  // Filter star outputs by unique alignment counts
  // If 0 unique alignments need to pass bam into
  // create_valid_empty_bam_star
  star.out.out_bam
        .branch { 
          good_bam: it[2].toInteger() > 0
          bad_bam: it[2].toInteger() == 0
          }
        .set{star_out_ch}

  // Create an empty one-line-header bam that can be read
  // by samtools. We only create this for those samples that
  // had 0 reads after QC (i.e. after trim_extra_polyA) or
  // after mapping i.e. after star.
  create_valid_empty_bam_star(polyA_out_ch.empty_fastq.map({[it[0], "_Aligned.sortedByCoord.out"]}).mix(star_out_ch.bad_bam.map({[it[0], "_Aligned.sortedByCoord.out"]})))

  // Qualimap on STAR output
  // Annotate the channel objects with dummy counts of either 1 or 0
  // depending on whether the bams are empty are not to either run
  // qualimap or populate an empty qualimap template
  raw_qualimap(star_out_ch.good_bam.map({[it[0], it[1], 1]}).mix(create_valid_empty_bam_star.out.out_bam.map({[it[0], it[1], 0]})), gtf, empty_qualimap_template, "raw")

  // Process to convert input GTF to gene model bed for rseqc
  gtf2bed_script = file("${baseDir}/bin/gtf2bed")
  gtf2bed(gtf, gtf2bed_script)
  // RSeQC read distribution on STAR output
  // Set gene_model based on genome (note -genome option isn't required...)
  raw_rseqc(star_out_ch.good_bam.map({[it[0], it[1], 1]}).mix(create_valid_empty_bam_star.out.out_bam.map({[it[0], it[1], 0]})), gtf2bed.out.bed, "raw")
  ch_raw_rseqc_multiqc = raw_rseqc.out.rseqc_log
  // Perform featurecount quantification
  // The 1 and 0 being added in the map represent bams that contain (1)
  // or do not (0) contain alignments.
  // If alignments are present featureCounts is run,
  // else the empty bam is simply copied for collection from the process.
  feature_counts(star_out_ch.good_bam.map({[it[0], it[1], 1]}).mix(create_valid_empty_bam_star.out.out_bam.map({[it[0], it[1], 0]})), gtf, "${baseDir}/bin/assign_multi_mappers.gawk")

  // Produce qualimap output of the annotated bam for metrics
  annotated_rseqc(feature_counts.out.out_bam, gtf2bed.out.bed, "annotated")
  ch_annotated_rseqc_multiqc = annotated_rseqc.out.rseqc_log
  annotated_qualimap(feature_counts.out.out_bam, gtf, empty_qualimap_template, "annotated")

  // Generate input channel containing all the files needed for multiqc per samples. 
  // The final channel structure is [sample, [file1, file2, file3, ...]]
  // NB fastqc currently does not integrate multip qualimap
  // outputs successfuly so they are manually carried through
  // to the summary_statistics
  ch_merged_fastp_multiqc
    .mix(ch_post_polyA_fastp_multiqc)
    .mix(ch_io_extract_fastp_multiqc)
    .mix(ch_raw_rseqc_multiqc)
    .mix(ch_annotated_rseqc_multiqc)
    .groupTuple(by:0, size: 5)
    .map({it.flatten()}).map({[it[0], it.tail()]})
    .set { ch_multiqc_in }

  // ch_multiqc_in.view()
  // Run multiqc  
  multiqc(ch_multiqc_in)
  
  // Sort and index bam file
  sort_index_bam(feature_counts.out.out_bam)

  // Perform deduplication
  dedup(sort_index_bam.out.sort_index_bam_out)

  // Generate file for count matrix
  io_count(dedup.out.io_dedup_sam)
  ch_io_count_out = io_count.out.io_count_out

  // Generate raw count matrix
  count_matrix(ch_io_count_out, barcode_list, feature_file_out)
  ch_h5ad = count_matrix.out.h5ad

  // Run cell caller
  cell_caller(ch_h5ad)
  ch_cell_caller_out = cell_caller.out.cell_caller_out //[val(sample), int(cell_caller_nuc_gene_threshold)]

  // Sort the groupTuple so that the int is always
  // first and then flatten the tuple list to return a 3mer
  // N.B. We were originally sorting by class (sort:{it.getClass() == sun.nio.fs.UnixPath ? 1 : 0})
  // but for some reason this only worked locally and not on Seqera Platform
  ch_filter_count_matrix_in = ch_cell_caller_out.mix(ch_h5ad)
  .groupTuple(by: 0, size:2, sort:{order_integer_first(it)})
  .map{[it[0], it[1][0], it[1][1]]}

  // Output filtered (cells only) count tables
  filter_count_matrix(ch_filter_count_matrix_in)

  // structure of ch_summary_statistics_in is
  // [sample, min_nuc_gene_cutoff, raw_h5ad, annotated_qualimap,
  // antisense, dedup.log, multiqc_data, raw_qualimap]
  ch_cell_caller_out
  .join(filter_count_matrix.out.raw_count_matrix)
  .join(annotated_qualimap.out.qualimap_txt)
  .join(sort_index_bam.out.antisense_out)
  .join(dedup.out.io_dedup_log)
  .join(multiqc.out.multiqc_json)
  .join(raw_qualimap.out.qualimap_txt)
  .join(raw_rseqc.out.rseqc_log)
  .join(annotated_rseqc.out.rseqc_log)
  .set({ch_summary_statistics_in})

  // Generate summary statistics
  summary_statistics(ch_summary_statistics_in)

  ch_summary_metrics_and_plot = summary_statistics.out.metrics_csv.join(cell_caller.out.cell_caller_plot, by:0)

  // Generate single sample report
  single_summary_report(ch_summary_metrics_and_plot, single_sample_report_template)

  // Generate multi sample report
  multi_sample_report(single_summary_report.out.single_sample_metric_out.collect(), multi_sample_report_template)
 
}
