#!/usr/bin/env nextflow

/*
* A Nextflow pipeline for processing scRNA-Seq data generated using CS Genetics'
* single-cell kit to produce a genes by barcode count table.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genomes config e.g. gtf, star_index
// Follows nf-core pattern from igenomes.config
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

def requireSafeSampleId(value) {
    def sample = value == null ? '' : value.toString()
    if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/) || sample in ['.', '..']) {
        throw new IllegalArgumentException(
            "Invalid sample ID '${sample}': use 1-128 letters, numbers, '.', '_' " +
            "or '-', beginning with a letter or number; '.' and '..' are not allowed"
        )
    }
    return sample
}

def canonicalManualThreshold(value, label) {
    def raw = value == null ? '' : value.toString().trim()
    if (raw.isEmpty()) {
        return 'nan'
    }
    def threshold
    try {
        threshold = new BigDecimal(raw)
    } catch (NumberFormatException exc) {
        throw new IllegalArgumentException(
            "${label} manual Cell Caller threshold must be a finite non-negative number"
        )
    }
    if (threshold.signum() < 0) {
        throw new IllegalArgumentException(
            "${label} manual Cell Caller threshold must be a finite non-negative number"
        )
    }
    if (threshold > 308) {
        throw new IllegalArgumentException(
            "${label} manual Cell Caller threshold is too large"
        )
    }
    return threshold.stripTrailingZeros().toPlainString()
}

def canonicalMinimumCountThreshold(value) {
    def raw = value == null ? '' : value.toString().trim()
    if (!(raw ==~ /\d+/)) {
        throw new IllegalArgumentException(
            "minimum_count_threshold must be a non-negative integer"
        )
    }
    def threshold = new BigInteger(raw)
    if (threshold > Integer.MAX_VALUE) {
        throw new IllegalArgumentException(
            "minimum_count_threshold is too large"
        )
    }
    return threshold.intValue()
}

def mergeSampleThresholds(sample, thresholds, mixedSpecies) {
    def parts = thresholds.collect { threshold ->
        mixedSpecies ? threshold.split('_', -1).toList() : [threshold]
    }
    def width = mixedSpecies ? 2 : 1
    if (parts.any { it.size() != width }) {
        throw new IllegalArgumentException(
            "Sample '${sample}' has malformed manual Cell Caller thresholds"
        )
    }
    def merged = (0..<width).collect { index ->
        def supplied = parts.collect { it[index] }.findAll { it != 'nan' }.unique()
        if (supplied.size() > 1) {
            throw new IllegalArgumentException(
                "Sample '${sample}' has conflicting manual Cell Caller thresholds: ${thresholds.unique()}"
            )
        }
        return supplied ? supplied[0] : 'nan'
    }
    return mixedSpecies ? merged.join('_') : merged[0]
}

def requireSamplesheetHeader(row, mixedSpecies) {
    def keys = row.keySet().collect { it.toString() }
    def accepted
    if (mixedSpecies) {
        accepted = [
            ['sample', 'fastq_1', 'fastq_2'],
            ['sample_id', 'fastq_1', 'fastq_2'],
            ['sample', 'fastq_1', 'fastq_2', 'hsap_manual_cell_caller_threshold', 'mmus_manual_cell_caller_threshold'],
            ['sample_id', 'fastq_1', 'fastq_2', 'hsap_manual_cell_caller_threshold', 'mmus_manual_cell_caller_threshold'],
            ['sample', 'fastq_1', 'fastq_2', 'hsap_manual_cellcaller_threshold', 'mmus_manual_cellcaller_threshold'],
            ['sample_id', 'fastq_1', 'fastq_2', 'hsap_manual_cellcaller_threshold', 'mmus_manual_cellcaller_threshold'],
        ]
    } else {
        accepted = [
            ['sample', 'fastq_1', 'fastq_2'],
            ['sample_id', 'fastq_1', 'fastq_2'],
            ['sample', 'fastq_1', 'fastq_2', 'manual_cell_caller_threshold'],
            ['sample_id', 'fastq_1', 'fastq_2', 'manual_cell_caller_threshold'],
            ['sample', 'fastq_1', 'fastq_2', 'manual_cellcaller_threshold'],
            ['sample_id', 'fastq_1', 'fastq_2', 'manual_cellcaller_threshold'],
        ]
    }
    if (!accepted.contains(keys)) {
        throw new IllegalArgumentException(
            "Invalid input CSV header ${keys}; expected one of ${accepted}"
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESOLVE GENOME PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Resolve genome attributes into params
// These can be overridden by command-line parameters
params.star_index                   = getGenomeAttribute('star_index')
params.gtf                          = getGenomeAttribute('gtf')
params.mitochondria_chromosome      = getGenomeAttribute('mitochondria_chromosome')
params.mixed_species                = getGenomeAttribute('mixed_species')
params.hsap_mitochondria_chromosome = getGenomeAttribute('hsap_mitochondria_chromosome')
params.mmus_mitochondria_chromosome = getGenomeAttribute('mmus_mitochondria_chromosome')
params.hsap_gene_prefix             = getGenomeAttribute('hsap_gene_prefix')
params.mmus_gene_prefix             = getGenomeAttribute('mmus_gene_prefix')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { save_resolved_configuration } from './modules/local/save_resolved_configuration/main.nf'
include { download_star_index } from './modules/local/download_star_index/main.nf'
include { download_gtf } from './modules/local/download_gtf/main.nf'
include { download_input_csv } from './modules/local/download_input_csv/main.nf'
include { download_barcode_list } from './modules/local/download_barcode_list/main.nf'
include { download_barcode_correction_list } from './modules/local/download_barcode_correction_list/main.nf'
include { download_public_fastq } from './modules/local/download_public_fastq/main.nf'
include { features_file } from './modules/local/features_file/main.nf'
include { merge_lanes } from './modules/local/merge_lanes/main.nf'
include { qc } from './modules/local/qc/main.nf'
include { star } from './modules/local/star/main.nf'
include { create_valid_empty_bam as create_valid_empty_bam_star } from './modules/local/create_valid_empty_bam/main.nf'
include { gtf2bed } from './modules/local/gtf2bed/main.nf'
include { run_rseqc as raw_rseqc } from './modules/local/run_rseqc/main.nf'
include { run_rseqc as annotated_rseqc } from './modules/local/run_rseqc/main.nf'
include { initial_feature_count } from './modules/local/initial_feature_count/main.nf'
include { filter_for_UMRs_mismatch } from './modules/local/filter_for_UMRs_mismatch/main.nf'
include { umr_transcript_assignment } from './modules/local/umr_transcript_assignment/main.nf'
include { umr_exon_assignment } from './modules/local/umr_exon_assignment/main.nf'
include { multimapper_assignment } from './modules/local/multimapper_assignment/main.nf'
include { merge_transcript_exon_umr_bams } from './modules/local/merge_transcript_exon_umr_bams/main.nf'
include { merge_annotated_UMRs_with_annotated_multimappers } from './modules/local/merge_annotated_UMRs_with_annotated_multimappers/main.nf'
include { count_high_conf_annotated_umr_multimap } from './modules/local/count_high_conf_annotated_umr_multimap/main.nf'
include { single_sample_multiqc } from './modules/local/single_sample_multiqc/main.nf'
include { multi_sample_multiqc } from './modules/local/multi_sample_multiqc/main.nf'
include { sort_index_bam } from './modules/local/sort_index_bam/main.nf'
include { dedup } from './modules/local/dedup/main.nf'
include { io_count } from './modules/local/io_count/main.nf'
include { count_matrix } from './modules/local/count_matrix/main.nf'
include { filter_count_matrix } from './modules/local/filter_count_matrix/main.nf'
include { cell_caller } from './modules/local/cell_caller/main.nf'
include { categorize_reads } from './modules/local/categorize_reads/main.nf'
include { summary_statistics } from './modules/local/summary_statistics/main.nf'
include { qc_cascade_plot_single } from './modules/local/qc_cascade_plot_single/main.nf'
include { qc_cascade_plot_multi } from './modules/local/qc_cascade_plot_multi/main.nf'
include { consolidated_report } from './modules/local/consolidated_report/main.nf'

workflow {

  // Canonicalise the only CLI numeric interpolated into a task shell. Nextflow
  // params are immutable once the workflow starts, so carry the validated
  // integer explicitly to every consumer instead of attempting to reassign the
  // params entry. The JSON schema is Launch UI metadata, not runtime validation.
  def minimum_count_threshold = canonicalMinimumCountThreshold(
    params.minimum_count_threshold
  )

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

  // Nextflow currently doesn't have functionality to output the resolved configuration
  // https://github.com/nextflow-io/nextflow/issues/1515
  // They may develop this functionality in the future, but for now we will use a process
  // to output the resolved configuration to a file.
  save_resolved_configuration(minimum_count_threshold)

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
    input_csv = channel.fromPath(params.input_csv)
  }
  
  // Finally, we need to check the fastq files to see if they
  // are hosted in the s3://csgx.public.readonly bucket.
  // If so, download via a process, else directly make a file object
  // so that locally configured AWS credentials can be used.
  // We make the assumption that if R1 and R2 for each row are hosted in the same bucket.
  input_csv
    .splitCsv(header:true, sep:',', strip: true)
    .map {row ->
        requireSamplesheetHeader(row, params.mixed_species)
        def rowList = row.values().toList()
        def sample = requireSafeSampleId(rowList[0])
        def fastq_1 = rowList[1]
        def fastq_2 = rowList[2]
        if (!fastq_1 || !fastq_2 || fastq_1 == fastq_2) {
          throw new IllegalArgumentException(
            "Sample '${sample}' must provide distinct, non-empty fastq_1 and fastq_2 paths"
          )
        }
        def download = fastq_1.startsWith('s3://csgx.public.readonly') && fastq_2.startsWith('s3://csgx.public.readonly') ? 'download' : 'no_download'
        return [sample, fastq_1, fastq_2, download]
    }
    .branch { row ->
        download: row[3] == 'download'
        no_download: row[3] == 'no_download'
    }
    .set { split_ch }

  split_ch.download
    .map { row -> row[0..2] } // remove the 'download' string from the tuple
    .set { to_download_ch }

  split_ch.no_download
    .map { row -> row[0..2] } // remove the 'no_download' string from the tuple
    .map { row -> [row[0], file(row[1]), file(row[2])] } // convert the paths to File objects
    .set { no_download_ch }

  // Download the fastq files from the s3://csgx.public.readonly bucket
  download_public_fastq(to_download_ch)

  download_public_fastq.out.downloaded_fastqs
    .mix(no_download_ch)
    .groupTuple(by:0)
    .set { ch_input }

  // Check whether params.barcode_list starts with s3://csgx.public.readonly
  // and if it does, download the file in a process and set the barcode_list to the downloaded file
  if (params.barcode_list_path.startsWith("s3://csgx.public.readonly")) {
    barcode_list = download_barcode_list()
  } else {
    barcode_list = file(params.barcode_list_path)
  }

  // Check whether params.barcode_correction_list_path starts with s3://csgx.public.readonly
  // and if it does, download the file in a process and set the barcode_correction_list to the downloaded file
  if (params.barcode_correction_list_path.startsWith("s3://csgx.public.readonly")) {
    barcode_correction_list = download_barcode_correction_list()
  } else {
    barcode_correction_list = file(params.barcode_correction_list_path)
  }

  // Create path objects to HTML report template + vendored (offline-safe) assets
  consolidated_report_template = file("${baseDir}/templates/consolidated_report_template.html.jinja2")
  report_vendor_dir = file("${baseDir}/assets/vendor")
  // Create empty rseqc output template path object
  empty_rseqc_template = file("${baseDir}/templates/rseqc_empty_template.txt")
  // Create feature file for count_matrix from GTF
  features_file(gtf)
  feature_file_out = features_file.out.modified_gtf

  // ch_input is in the form: [sample, [R1_L001, R1_L002], [R2_L001, R2_L002]],  [sample, [R1_L001], [R2_L001]]
  // The number of items in each of the sub arrays depends on how many lanes the sample was sequenced over.
  // We only want to pass those samples that have multiple lanes worth of fastqs into merge_lanes.
  // For those samples that do not have multiple lanes, we will skip merge lanes and merge with the output
  // of merge_lanes.
  
  // Split ch_input base on whether there are >1 R1 fastqs
  ch_input.branch{ sample_data ->
    multiple_lanes: sample_data[1].size() > 1
    single_lane: sample_data[1].size() == 1
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
    .flatMap{ sample_data -> [[sample_data[0], "R1", sample_data[1]], [sample_data[0], "R2", sample_data[2]]]}  

  // This process will merge fastqs split over multiple lanes 
  merge_lanes(merge_lanes_in)
  ch_merge_lanes_out = merge_lanes.out.merge_lanes_out

  // Combine the R1 and R2 reads per sample ensuring R1 comes first
  // --> [sample, *.merged.R1*, *.merged.R2*]
  ch_merge_lanes_out_merged = ch_merge_lanes_out.groupTuple(by:0, size:2).map({ grouped -> grouped[1][0] == "R1" ? [grouped[0], grouped[2][0], grouped[2][1]] : [grouped[0], grouped[2][1], grouped[2][0]]})

  // Flatten the R1 and R2 in the non-merged fastq pairs
  // [sample, [R1], [R2]] --> [sample, R1, R2]
  ch_input_split_single_lane_flattened = ch_input_split.single_lane
    .map{ sample_data -> [sample_data[0], sample_data[1][0], sample_data[2][0]]}

  // Merge the merged and non-merged fastqs
  io_extract_in_ch = ch_merge_lanes_out_merged.mix(ch_input_split_single_lane_flattened)

  // Run unified QC process (replaces io_extract + io_extract_fastp + trim_extra_polya + post_polyA_fastp + merged_fastp)
  // This single process:
  // 1. Extracts and corrects barcodes
  // 2. Performs SSS trimming and polyX trimming
  // 3. Trims internal polyA
  // 4. Calculates Q30 metrics for R2 barcode and R1 output
  // 5. Outputs JSON files for MultiQC
  qc(io_extract_in_ch, barcode_correction_list)
  ch_qc_log = qc.out.qc_log
  ch_qc_multiqc = qc.out.qc_multiqc

  // Filter for empty fastq
  // Pipe good to STAR
  // Pipe empty to create_valid_empty_bam_star
  qc.out.qc_out
        .branch { qc_result ->
          good_fastq: qc_result[1].countFastq() > 0
          empty_fastq: qc_result[1].countFastq() == 0
          }
        .set{qc_out_filtered_ch}

  // Align the good fastqs with STAR
  star(qc_out_filtered_ch.good_fastq, star_index)

  // Route STAR output by all accepted alignment evidence, not only unique
  // alignments. A sample with zero unique reads but accepted multi-locus reads
  // is real data and must continue through UMR/multimapper assignment.
  star.out.out_bam
        .branch { star_result ->
          aligned_bam: star_result[3].toLong() > 0
          unaligned_bam: star_result[3].toLong() == 0
          }
        .set{star_out_ch}

  // Create an empty one-line-header bam that can be read
  // by samtools. We only create this for those samples that
  // had 0 reads after QC or no unique/multi-locus alignments after STAR.
  create_valid_empty_bam_star(qc_out_filtered_ch.empty_fastq.map({ qc_result -> [qc_result[0], "_Aligned.sortedByCoord.out"]}).mix(star_out_ch.unaligned_bam.map({ star_result -> [star_result[0], "_Aligned.sortedByCoord.out"]})))

  // Canonical STAR BAM channel for every sample, exactly once. Aligned samples
  // retain STAR's BAM; QC-empty and truly unaligned samples use a valid
  // header-only BAM. Reuse this tuple shape everywhere that needs raw STAR
  // evidence so an empty sample cannot disappear from later joins.
  star_out_ch.aligned_bam
    .map({ star_result -> [star_result[0], star_result[1], true] })
    .mix(create_valid_empty_bam_star.out.out_bam.map({ bam_result -> [bam_result[0], bam_result[1], false] }))
    .set { ch_all_star_bams }

  ch_all_star_bams
    .branch { sample_id, bam, has_alignments ->
      aligned: has_alignments
      empty: !has_alignments
    }
    .set { ch_all_star_bams_split }

  // Process to convert input GTF to gene model bed for RSeQC
  gtf2bed(gtf)

  // RSeQC read distribution on STAR output
  // If alignments are present RSeQC is run; otherwise the empty RSeQC
  // template is populated in the process.
  raw_rseqc_in_ch = ch_all_star_bams.map { sample_id, bam, has_alignments ->
    [sample_id, bam, has_alignments ? 1L : 0L]
  }
  raw_rseqc(raw_rseqc_in_ch, gtf2bed.out.bed, empty_rseqc_template, "raw")
  ch_raw_rseqc_multiqc = raw_rseqc.out.rseqc_log

  // Split feature counting into multiple processes to take advantage of parallel processing
  // Perform featurecount quantification
  // If alignments are present, featureCounts is run,
  // else the empty bam is simply copied for collection from the process.
  initial_feature_count(ch_all_star_bams, gtf)

  // Make a featureCounts BAM channel containing only samples with alignments.
  // Essentially performs an inner join to limit input to samples present in
  // the canonical aligned branch.
  // Using combine instead of join because strict mode has failOnMismatch:true by default
  initial_feature_count_aligned_bam_out_ch = ch_all_star_bams_split.aligned
    .map({ star_result -> [star_result[0]] })
    .combine(initial_feature_count.out.feature_count_bam, by: [0])

  // If alignments are present run UMR and multimapper processing
  // UMRs
  // Generate a bam with only UMRs, and up to 3 mismatches
  filter_for_UMRs_mismatch(initial_feature_count_aligned_bam_out_ch)
  // Get the first set of gene associations based on transcript feature annotations
  umr_transcript_assignment(filter_for_UMRs_mismatch.out.umr_mismatch_bam)
  // Get the second set of gene associations based on exon feature annotations (i.e. exon-tie breaking)
  umr_exon_assignment(filter_for_UMRs_mismatch.out.umr_mismatch_bam, gtf)

  // Multimappers: fused filter + transcript assignment + exon tie-break + merge (one task per sample
  // to cut serial container-starts/staging on the critical path; output unchanged vs the old 4 processes)
  multimapper_assignment(initial_feature_count_aligned_bam_out_ch, gtf, file("${baseDir}/bin/assign_multi_mappers.gawk"))

  // Merge the transcript- and exon-based gene assignments for the umrs
  merge_transcript_exon_umr_bams(umr_transcript_assignment.out.umr_transcript_assigned_bam.combine(umr_exon_assignment.out.umr_exon_assigned_bam, by: 0))

  // Merge the multimapper and UMR bams
  merge_annotated_UMRs_with_annotated_multimappers(merge_transcript_exon_umr_bams.out.high_conf_annotated_umr_bam.combine(multimapper_assignment.out.high_conf_annotated_multimapped_bam, by: 0))
  
  // Re-merge aligned annotations with the canonical header-only BAMs.
  umr_multimapper_annotated_bam_out_ch = merge_annotated_UMRs_with_annotated_multimappers.out.high_conf_annotated_bam
    .mix(ch_all_star_bams_split.empty.map({ star_result -> [star_result[0], star_result[1]] }))

  // Count number of aligned reads
  count_high_conf_annotated_umr_multimap(umr_multimapper_annotated_bam_out_ch)

  // Produce RSeQC output of the annotated bam for metrics
  // The aligned count is kept numeric here and in raw_rseqc_in_ch. A previous
  // Boolean-only process contract silently treated every positive annotated
  // count (for example "1") as false and emitted the empty metrics template.
  annotated_rseqc_in_ch = umr_multimapper_annotated_bam_out_ch.combine(count_high_conf_annotated_umr_multimap.out.aligned_count, by: 0)
  annotated_rseqc(annotated_rseqc_in_ch, gtf2bed.out.bed, empty_rseqc_template, "annotated")
  ch_annotated_rseqc_multiqc = annotated_rseqc.out.rseqc_log

  // Generate input channel containing all the files needed for multiqc per samples.
  // The final channel structure is [sample, [file1, file2, file3, ...]]
  ch_qc_multiqc
    .mix(ch_raw_rseqc_multiqc)
    .mix(ch_annotated_rseqc_multiqc)
    .groupTuple(by:0, size: 3)
    .map({ grouped -> grouped.flatten()}).map({ flattened -> [flattened[0], flattened.tail()]})
    .set { ch_ss_multiqc_in }

  // Run single-sample multiqc  
  single_sample_multiqc(ch_ss_multiqc_in)

  // Collect all files from single sample multiqc input into a single list containing all samples
  // But removing the sample_id
  ch_ms_multiqc_in = ch_ss_multiqc_in.map({ sample_files -> sample_files[1] }).collect()

  // Run multi-sample multiqc
  multi_sample_multiqc(ch_ms_multiqc_in)
  
  // Sort and index bam file
  sort_index_bam(umr_multimapper_annotated_bam_out_ch)

  // Perform deduplication
  dedup_in_ch = sort_index_bam.out.sort_index_bam_out.combine(count_high_conf_annotated_umr_multimap.out.aligned_count, by: 0)
  dedup(dedup_in_ch)

  // Generate file for count matrix
  io_count(dedup.out.io_dedup_sam)
  ch_io_count_out = io_count.out.io_count_out

  // Generate raw count matrix
  count_matrix(ch_io_count_out, barcode_list, feature_file_out)
  ch_h5ad = count_matrix.out.h5ad

  // Extract any user specified thresholds for cell calling from the input_csv and format them for input to the cell caller process
  input_csv
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row ->
        requireSamplesheetHeader(row, params.mixed_species)
        def rowList = row.values().toList()
        def sample = requireSafeSampleId(rowList[0])

        if (params.mixed_species) {
          // In the mixed species case, any user specified hsap and mmus thresholds are expected in columns 4 and 5 of the input csv respectively
          // When specifying thresholds for mixed species, both columns must be present (even if they are left empty for some samples)
          if (rowList.size() < 5) {
            // If one or both columns are missing, no user thresholds are used so the input for cell calling is "nan_nan"
            return [sample, "nan_nan"]
          } else {
            // If both columns are present, we will set the user thresholds to "nan" if they are empty for a particular sample, otherwise we carry through the input value 
            def hsap_threshold
            def mmus_threshold
            if (rowList[3].isEmpty()) {
              hsap_threshold = "nan"
            } else {
              hsap_threshold = canonicalManualThreshold(rowList[3], 'human')
            }
            if (rowList[4].isEmpty()) {
              mmus_threshold = "nan"
            } else {
              mmus_threshold = canonicalManualThreshold(rowList[4], 'mouse')
            }
            return [sample, "${hsap_threshold}_${mmus_threshold}"]
          }
          
        } else {
          // In the single species case, the user specified threshold is expected in column 4
          if (rowList.size() < 4) {
            // If the column is missing, we will set the user threshold to "nan"
            return [sample, "nan"]
          } else {
            // If the column exists but is empty for a particular sample, we will set the threshold to "nan", otherwise we take the value from the input csv
            if (rowList[3].isEmpty()) {
              return [sample, "nan"]
            } else {
              return [sample, canonicalManualThreshold(rowList[3], 'single-species')]
            }
          }
        }
    }
    // Thresholds are sample-level metadata while rows are lanes. Collapse
    // identical values and reject conflicting values before scheduling tasks.
    .groupTuple(by: 0)
    .map { sample, thresholds ->
        return [sample, mergeSampleThresholds(sample, thresholds, params.mixed_species)]
    }
    .set { user_specified_cell_caller_thresholds_ch }

  
  ch_cell_caller = ch_h5ad.combine(user_specified_cell_caller_thresholds_ch, by: 0)

  // Run cell caller
  cell_caller(ch_cell_caller, minimum_count_threshold)
  ch_cell_caller_out = cell_caller.out.cell_caller_out //[val(sample), int(cell_caller_nuc_gene_threshold)]

  // cell_caller_out is [sample_id, threshold] and ch_h5ad is [sample_id, raw_h5ad]; join by sample_id
  // gives the [sample_id, threshold, raw_h5ad] tuple that filter_count_matrix expects, deterministically.
  // (This previously used mix + groupTuple(size:2) with a sort closure that put the int threshold before
  // the path. The sort relied on the path throwing MissingMethodException on .isInteger(); but on Seqera
  // Platform / Fusion the path object responds to isInteger() without throwing, so both elements tied at
  // sort key 0 and groupTuple fell back to arrival order -- a race between cell_caller and count_matrix
  // that silently put the threshold in the path slot on heavy samples, crashing filter_count_matrix with
  // "Not a valid path value: '<threshold>'". join keys on sample_id and is order-deterministic.)
  ch_filter_count_matrix_in = ch_cell_caller_out.join(ch_h5ad, by: 0)

  // Output filtered (cells only) count tables
  filter_count_matrix(ch_filter_count_matrix_in)

  // Create input channel for categorize_reads process
  // Need to combine STAR BAM, raw count matrix H5AD, and qc JSON (main qc.json file)
  ch_categorize_reads_in = ch_all_star_bams
    .map({ star_result -> [star_result[0], star_result[1]]})  // [sample_id, star_bam]
    .combine(filter_count_matrix.out.raw_count_matrix, by: [0])  // [sample_id, star_bam, raw_h5ad]
    .combine(ch_qc_multiqc.map({ qc_result -> [qc_result[0], qc_result[1]]}), by: [0])  // [sample_id, star_bam, raw_h5ad, qc_json]

  // Run categorize_reads to calculate read and count metrics
  categorize_reads(ch_categorize_reads_in)

  // structure of ch_summary_statistics_in is
  // [sample, min_nuc_gene_cutoff, raw_h5ad,
  // antisense, dedup.log, multiqc_data, raw_seqc, annotated_rseqc, read_categorization_csv, qc_log]
  ch_cell_caller_out
  .combine(filter_count_matrix.out.raw_count_matrix, by: 0)
  .combine(sort_index_bam.out.antisense_out, by: 0)
  .combine(dedup.out.io_dedup_log, by: 0)
  .combine(single_sample_multiqc.out.multiqc_json, by: 0)
  .combine(raw_rseqc.out.rseqc_log, by: 0)
  .combine(annotated_rseqc.out.rseqc_log, by: 0)
  .combine(categorize_reads.out.read_categories.map({ read_cat -> [read_cat[0], read_cat[1]]}), by: 0)
  .combine(ch_qc_log, by: 0)
  .set({ch_summary_statistics_in})

  // Generate summary statistics
  summary_statistics(ch_summary_statistics_in)

  // Generate single-sample QC cascade plots. Each task emits a small internal
  // fragment for consolidated_report and a separately published offline page.
  qc_cascade_plot_single(summary_statistics.out.metrics_csv)

  // Generate the equivalent internal + standalone outputs across all samples.
  qc_cascade_plot_multi(summary_statistics.out.metrics_csv.map { metrics_tuple -> metrics_tuple[1] }.collect())

  // Collect, flat, all per-sample inputs for the single consolidated report.
  // Metrics csvs are also published per-sample by summary_statistics.
  ch_all_metrics_csvs = summary_statistics.out.metrics_csv.map { it[1] }.collect()
  // cell_caller_plots = tuple(sample_id, counts_pdf_html, barnyard_html); keep the file paths.
  ch_all_cell_caller_plots = cell_caller.out.cell_caller_plots
    .flatMap { sample_id, counts_pdf, barnyard -> [counts_pdf, barnyard] }
    .collect()
  ch_all_qc_cascade_fragments = qc_cascade_plot_single.out.qc_cascade_fragment
    .map { qc_cascade_tuple -> qc_cascade_tuple[1] }
    .collect()

  // Run-provenance metadata for the report header/footer (all values already known
  // to the pipeline; factual only -- no quality judgement). Serialised to JSON and
  // passed to the consolidated_report process as a single value.
  def provenance_json = groovy.json.JsonOutput.toJson([
    genome:        params.genome,
    annotation:    params.gtf ? file(params.gtf).name : 'N/A',
    mixed:         params.mixed_species,
    pipeline_ver:  workflow.manifest.version ?: 'N/A',
    commit:        workflow.commitId ?: 'N/A',
    revision:      workflow.revision ?: 'N/A',
    run_name:      workflow.runName,
    session_id:    workflow.sessionId.toString(),
    start:         workflow.start.toString(),
    nf_version:    workflow.nextflow.version.toString(),
    outdir:        params.outdir,
    barcode_kit:   params.barcode_list_path ? file(params.barcode_list_path).name : 'N/A',
    count_threshold: minimum_count_threshold,
    homepage:      workflow.manifest.homePage ?: 'https://github.com/csgenetics/csgenetics_scrnaseq'
  ])
  // Never interpolate raw JSON into a task shell or env export: customer paths
  // and run names can contain quotes/metacharacters. Base64 is shell-inert and
  // is decoded with strict validation by create_consolidated_report.py.
  def provenance_base64 = provenance_json.getBytes('UTF-8').encodeBase64().toString()

  // Generate the single consolidated, self-contained experiment report.
  consolidated_report(
    ch_all_metrics_csvs,
    ch_all_cell_caller_plots,
    ch_all_qc_cascade_fragments,
    qc_cascade_plot_multi.out.qc_cascade_fragment,
    consolidated_report_template,
    report_vendor_dir,
    provenance_base64
  )

}
