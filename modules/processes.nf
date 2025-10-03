#!/usr/bin/env nextflow

nextflow.enable.dsl=2
import groovy.json.JsonOutput

process save_resolved_configuration{
  tag "Save resolved configuration"
  publishDir "${params.outdir}/pipeline_info/", mode: 'copy'

  output:
  path("resolved_configuration.txt")

  script:
  json_str = JsonOutput.toJson(params)
  json_indented = JsonOutput.prettyPrint(json_str)
  // NOTE: single quotes are critical here;
  """
  echo '${json_indented}' > resolved_configuration.txt
  """
}

/*
* Download the star index from the s3://csgx.public.readonly bucket
*/
process download_star_index {
  tag "Download STAR index"

  output:
  path("star"), emit: star_index

  script:
  """
  mkdir star
  aws s3 cp --no-sign-request ${params.star_index} ./star --recursive
  """
}

/* Download the gtf file from the s3://csgx.public.readonly bucket */
process download_gtf {
  tag "Download GTF"

  output:
  path("*.gtf")

  script:
  """
  aws s3 cp --no-sign-request ${params.gtf} .
  """
}

/* Download the input_csv file from the s3://csgx.public.readonly bucket */
/* Can also be used to download the user specified cell caller threshold csv, if specified */
process download_input_csv {
  tag "Download input csvs"

  output:
  path("*.csv")

  script:
  """
  aws s3 cp --no-sign-request ${params.input_csv} .
  """
}

/* Download the csgx hosted barcode csv
* from the s3://csgx.public.readonly bucket.
*/
process download_barcode_list {
  tag "Get barcode list"

  output:
  path("*.csv")

  script:
  """
  aws s3 cp --no-sign-request ${params.barcode_list_path} .
  """
}

/* Download the csgx hosted barcode correction list tsv
* from the s3://csgx.public.readonly bucket.
*/
process download_barcode_correction_list {
  tag "Get barcode correction list"

  output:
  path("*.tsv")

  script:
  """
  aws s3 cp --no-sign-request ${params.barcode_correction_list_path} .
  """
}

/* Download a fastq hosted on the s3://csgx.public.readonly bucket */
process download_public_fastq {
  tag "${sample_id} ${fastq_1.tokenize('/').last()}"

  input:
  tuple val(sample_id), val(fastq_1), val(fastq_2)

  output:
  tuple val(sample_id), path("${fastq_1.tokenize('/').last()}"), path("${fastq_2.tokenize('/').last()}"), emit: downloaded_fastqs

  script:
  def fastq_1_name = fastq_1.tokenize('/').last()
  def fastq_2_name = fastq_2.tokenize('/').last()
  """
  echo $fastq_1_name
  aws s3 cp --no-sign-request $fastq_1 $fastq_1_name
  aws s3 cp --no-sign-request $fastq_2 $fastq_2_name
  """
}

/*
* Convert GTF into features file
*/
process features_file {
  input:
  path(gtf)

  output:
  path("${gtf.baseName}_features_names.tsv"), emit: modified_gtf

  shell:
  '''
  features_names.py !{gtf} !{gtf.baseName}_features_names.tsv
  '''
}

/*
* Merge Lanes
*/
process merge_lanes {
  tag "$sample_id $read_num"

  input:
  tuple val(sample_id), val(read_num), path(fastqs)

  output:
  tuple val(sample_id), val(read_num), path("${sample_id}.merged.${read_num}.fastq.gz"), emit: merge_lanes_out

  script:
  """
  # cat ${fastqs} > ${sample_id}.merged.${read_num}.fastq.gz

  # Globs are ordered so lane merging will happen in same order for R1 and R2
  cat *.f*q.gz > ${sample_id}.merged.${read_num}.fastq.gz
  """
}

process merged_fastp{
  tag "$sample_id"

  publishDir "${params.outdir}/fastp", mode: 'copy'

  input:
  tuple val(sample_id), path("${sample_id}.R1.merged.fastq.gz"), path("${sample_id}.R2.merged.fastq.gz")
  val(barcode_pattern)

  output:
  tuple val(sample_id), path("${sample_id}_R*.merged_fastp.html"), path("${sample_id}_R*.merged_fastp.json"), emit: merged_fastp_multiqc

  shell:
  barcode_length = barcode_pattern.size()
  '''
  # Get q30_rate for barcode present in fastq2 file 
  # maximum length which be equal to barcode length 
  barcode_length=$(echo !{barcode_pattern} | tr -cd 'C' | wc -c)

  # R1 fastp
  fastp -i !{sample_id}.R1.merged.fastq.gz \
    -A -G \
    -j !{sample_id}_R1.merged_fastp.json \
    -h !{sample_id}_R1.merged_fastp.html

  # R2 fastp
  fastp -i !{sample_id}.R2.merged.fastq.gz \
    -A -G \
    -b !{barcode_length} \
    -l 13 \
    -j !{sample_id}_R2.merged_fastp.json \
    -h !{sample_id}_R2.merged_fastp.html
  '''
}

/*
* Unified QC process - replaces io_extract + io_extract_fastp + trim_extra_polya
* This process combines:
* 1. Barcode extraction and correction (io_extract functionality)
* 2. SSS trimming and polyX trimming (io_extract_fastp functionality)
* 3. Internal polyA trimming (trim_extra_polya functionality)
* 4. Q30 calculation for both R2 barcode and R1 output
* All in a single Rust binary for maximum performance
*/
process qc {
  tag "$sample_id"

  container 'quay.io/csgenetics/qc:0.3'

  input:
  tuple val(sample_id), path(r1), path(r2)
  path(corrected_barcodelist)

  output:
  tuple val(sample_id), path("${sample_id}.qc.R1.fastq.gz"), emit: qc_out
  tuple val(sample_id), path("${sample_id}.qc.log"), emit: qc_log
  tuple val(sample_id), path("${sample_id}.R1.preQC.fastp.json"), path("${sample_id}.R2.preQC.fastp.json"), path("${sample_id}.R1.postQC.fastp.json"), emit: qc_multiqc

  script:
  """
  # Run unified qc Rust binary (available in PATH from container)
  qc \\
    ${r1} ${r2} ${corrected_barcodelist} \\
    --sample-id ${sample_id} \\
    --output-dir . \\
    --no-filtered

  # Check if output exists, create empty file if not
  if [ ! -f "${sample_id}.qc.R1.fastq.gz" ]; then
    touch ${sample_id}.qc.R1.fastq && gzip ${sample_id}.qc.R1.fastq
  fi
  """
}

/*
* Align reads with STAR
* We do the sorting separately in samtools to avoid the sorting RAM error
* thrown by STAR
*/
process star {
  tag "$sample_id"

  publishDir "${params.outdir}/STAR", mode: 'copy'

  input:
  tuple val(sample_id), path(r1)
  path(index)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.out.bam"), env(uniquely_mapped_reads), emit: out_bam

  script:
  """
      STAR --runThreadN 8 \
        --genomeDir ${index} \
        --readFilesIn ${r1} \
        --outFileNamePrefix ${sample_id}_ \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted \
        --readFilesCommand zcat \
        --outSAMattributes Standard \
        --outFilterMultimapNmax 1000

      # Get number of uniquely aligned reads
      uniquely_mapped_reads=\$(grep "Uniquely mapped reads number" ${sample_id}_Log.final.out | cut -d "|" -f 2 | xargs)

  """

}

/*
* Create an empty bam file that has a single line header
* so that the downstream samtools filtering doesn't break.
* Otherwise STAR and featureCount produce an empty bam that is an empty file
* that breaks samtools view: [main_samview] fail to read the header from "247intergenic_Aligned.sortedByCoord.out.bam".
*/
process create_valid_empty_bam{
  tag "$sample_id"

  publishDir "${params.outdir}/STAR", mode: 'copy'

  input:
  tuple val(sample_id), val(prefix)

  output:
  tuple val(sample_id), path("${sample_id}${prefix}.bam"), emit: out_bam

  script:
  """
  echo "@HD	VN:1.4	SO:coordinate" | samtools view -h -b > ${sample_id}${prefix}.bam
  """

}

/*
Process to convert the input GTF to a gene model bed file for rseqc read distribution
*/
process gtf2bed {
  input:
  path(gtf)

  output:
  path("gene_model.bed"), emit: bed

  shell:
  '''
  gtf2bed !{gtf} > gene_model.bed
  '''
}

/*
Generic process for running RSeQC read distribution
*/

process run_rseqc {
  tag "$sample_id"
  
  publishDir "${params.outdir}/RSeQC/read_distribution", mode: 'copy', pattern: "*_rseqc_results.txt", saveAs: {"${sample_id}.${prefix}_RSeQC.txt"}

  input:
  tuple val(sample_id), path(bam), val(count)
  path(bed)
  path(empty_rseqc_template)
  val(prefix)

  output:
  tuple val(sample_id), path("*_rseqc_results.txt"), emit: rseqc_log

  shell:
  '''
  if [[ !{count} > 0 ]]
    then
      read_distribution.py  -i !{bam} -r !{bed} > !{sample_id}_!{prefix}_rseqc_results.txt
    else      
      cat !{empty_rseqc_template} | envsubst > !{sample_id}_!{prefix}_rseqc_results.txt
  fi
  '''
}

process initial_feature_count {
  tag "$sample_id"

  publishDir "${params.outdir}/featureCounts", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), val(aligned_count)
  path(gtf)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam"), emit: feature_count_bam

  script:
  """
  if [[ $aligned_count > 0 ]] # If the bam is not empty
  then
    samtools sort $bam -o ${sample_id}_Aligned.sortedByCoord.out.bam
    # Start by running feature counts on the star output
    # including strandedness and annotation of multimappers
    featureCounts -a $gtf -o ${sample_id}.star.featureCounts.gene.txt -R BAM ${sample_id}_Aligned.sortedByCoord.out.bam -T 4 -t transcript -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
  else
    # Simply rename the input bam so that it can be collected
    cp $bam ${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam
  fi
  """
}

// Filter alignments to allow for at most 3 mismatches to the reference.
process filter_for_UMRs_mismatch {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(featurecount_bam)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.bam.featureCounts.bam"), emit: umr_mismatch_bam

  script:
  """
  samtools view -h -b -e '[NH]==1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${featurecount_bam} > ${sample_id}.UMRs.bam.featureCounts.bam
  """
}

// Collect those UMRs that have unambiguous transcript annotations based on transcript feature (i.e. no exon tie-breaking required)
process umr_transcript_assignment {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_mismatch_bam)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.transcript.assigned.bam"), emit: umr_transcript_assigned_bam

  script:
  """
  # Filter for the reads that were 'Assigned' a transcript target
  samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${umr_mismatch_bam} > ${sample_id}.UMRs.transcript.assigned.bam
  """
}

process umr_exon_assignment {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_mismatch_filtered_bam)
  path(gtf)

  output:
  tuple val(sample_id), path("${sample_id}.UMRs.exon.assigned.bam"), emit: umr_exon_assigned_bam

  script:
  """
  # Filter out the reads that were that were classified as Unassigned_Ambiguity and run them through
  # featureCounts using the exon tag to see if the ambiguity can be cleared up based on exon mapping.
  samtools view -h -b -e '[XS]=="Unassigned_Ambiguity"' -b $umr_mismatch_filtered_bam > ${sample_id}.UMRs.transcript.unassigned_ambiguity.bam

  # We have to remove the XS tag from the *.transcript.unassigned_ambiguity.bam because featureCounts adds
  # an additional tag, that prevents proper fitltering with samtools
  samtools view -h ${sample_id}.UMRs.transcript.unassigned_ambiguity.bam | sed 's/\\tXS\\:Z\\:[^\\t]*//' | samtools view -h -b > ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam

  # Do exon tie-breaking and filter to those that are assigned.
  featureCounts -a $gtf -o ${sample_id}.UMRs.exon.assigned.txt -R BAM ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam -T 4 -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1
  samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${sample_id}.UMRs.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam > ${sample_id}.UMRs.exon.assigned.bam
  """
}

process filter_for_multimappers_mismatch {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(feature_count_bam)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.bam.featureCounts.bam"), emit: multimap_mismatch_bam

  script:
  """
  samtools view -h -b -e '[NH]>1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${feature_count_bam} > ${sample_id}.multimapped.bam.featureCounts.bam
  """
}

process multimapper_transcript_assignment{
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_mismatch_filtered_bam)
  path(multi_mapper_script)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.transcript.assigned.bam"), emit: assigned_bam
  tuple val(sample_id), path("${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam"), emit: unassigned_bam

  script:
  """
  # Sort the bam file by query name in preparation for running through the gawk program for the first time.
  # Run the multimapped, annotated alignments through the gawk program that pulls
  # out alignments that can be associated directly as 'Assigned' and alignments that are ambiguous and should be passed onto
  # exon tie breaking.
  # See the script's header for more information on how it works.
  samtools sort -n $multimapper_mismatch_filtered_bam | samtools view | gawk -f $multi_mapper_script

  # multi_mapper_script produces assigned_reads.sam_body and ambiguous_reads.sam_body corresponding to the Assigned and still ambigous reads, respectively.
  # These files will only be produced if there were reads of the respective type identified.
  # If they exist, we need to convert both of these files back into valid bam format by adding a header and converting from sam->bam
  # If they don't exist then we should create a valid empty bam by the same name.

  if [ -f assigned_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    cat <(samtools view -H $multimapper_mismatch_filtered_bam) assigned_reads.sam_body | samtools view -b -h > ${sample_id}.multimapped.transcript.assigned.bam
  else
    samtools view -H -b $multimapper_mismatch_filtered_bam > ${sample_id}.multimapped.transcript.assigned.bam
  fi

  # If the ambigous_reads.sam_body exists then we convert this back to a valid bam.
  if [ -f ambiguous_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    # before rerunning through featureCounts for exon tie-breaking
    cat <(samtools view -H $multimapper_mismatch_filtered_bam) ambiguous_reads.sam_body | samtools view -h -b > ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam
  else
    samtools view -H -b $multimapper_mismatch_filtered_bam > ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam
  fi
  """
}

process multimapper_exon_assignment{
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_unassigned_bam)
  path(gtf)
  path(multi_mapper_script)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.exon.assigned.bam"), emit: assigned_bam

  script:
  """
  # Run featureCounts using the exon feature to do tie-breaking and then run back through the gawk script to pull out those
  # reads that have a single Assigned alignment.
  featureCounts -a $gtf -o ${sample_id}.multimapped.exon.assigned.txt -R BAM $multimapper_unassigned_bam -T 4 -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
  samtools sort -n ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam | samtools view | gawk -f $multi_mapper_script

  # If the assigned_reads.sam_body file exists then we were successfuly able to pull out further assigned reads
  if [ -f assigned_reads.sam_body ]; then
    # Cat with the headers of the featureCounts bam
    cat <(samtools view -H ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam) assigned_reads.sam_body | samtools view -b -h > ${sample_id}.multimapped.exon.assigned.bam;
  else
    # There were no further reads successfuly annotated
    # Create a valid empty bam
    samtools view -H -b ${sample_id}.multimapped.transcript.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam > ${sample_id}.multimapped.exon.assigned.bam
  fi
  """
}

process merge_transcript_exon_umr_bams {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_transcript_bam), path(umr_exon_bam)

  output:
  tuple val(sample_id), path("${sample_id}.umr.annotated.bam"), emit: high_conf_annotated_umr_bam

  script:
  """
  samtools merge -o ${sample_id}.umr.annotated.bam $umr_transcript_bam $umr_exon_bam
  """
}

process merge_transcript_exon_multimapper_bams {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(multimapper_transcript_bam), path(multimapper_exon_bam)

  output:
  tuple val(sample_id), path("${sample_id}.multimapped.annotated.bam"), emit: high_conf_annotated_multimapped_bam

  script:
  """
  samtools merge -o ${sample_id}.multimapped.annotated.bam $multimapper_transcript_bam $multimapper_exon_bam
  """
}

process merge_annotated_UMRs_with_annotated_multimappers {
  tag "$sample_id"

  publishDir "${params.outdir}/featureCounts", mode: 'copy', pattern: "${sample_id}.annotated.bam"

  input:
  tuple val(sample_id), path(umr_annotated_bam), path(multimapper_annotated_bam)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), emit: high_conf_annotated_bam

  script:
  """
  samtools merge -o ${sample_id}.mapped.sorted.filtered.annotated.bam $umr_annotated_bam $multimapper_annotated_bam
  """
}

process count_high_conf_annotated_umr_multimap {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(umr_multimapper_annotated_bam)

  output:
  tuple val(sample_id), env(alignment_count), emit: aligned_count

  script:
  """
  alignment_count=\$(samtools view -c ${umr_multimapper_annotated_bam})
  """
}

/*
* Run multiqc to single-sample report and data
*/
// TODO Pick up the multiqc_data.json directly and rename it in summary_statistics
// NOTE I am going to see if we can pick this up without a double star glob.
process single_sample_multiqc {
  publishDir "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}", mode: 'copy', pattern: "*_data"
  publishDir "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}", mode: 'copy', pattern: "*_multiqc.html"
  publishDir "${params.outdir}/multiqc/single_sample_multiqc/${sample_id}", mode: 'copy', pattern: "**/multiqc_data.json", saveAs: {"${sample_id}.multiqc.data.json"}

  input:
  tuple val(sample_id), path(multiqc_in_files)

  output:
  path "*_multiqc.html"
  path "*_data"
  tuple val(sample_id), path("**/multiqc_data.json"), emit: multiqc_json

  script:
  """
  multiqc . \
    -f \
    --title "${sample_id} multiqc" \
    --filename "${sample_id}_multiqc.html" \
    -m fastp \
    -m rseqc
  """
}

/*
* Run multiqc to generate multisample report and data
*/
process multi_sample_multiqc {
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc_data"
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc.html"
  publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multisample_multiqc_data/multiqc_data.json", saveAs: {"multisample.multiqc.data.json"}

  input:
  path(multiqc_in_files)

  output:
  path "multisample_multiqc.html"
  path "multisample_multiqc_data"
  path("multisample_multiqc_data/multiqc_data.json")

  script:
  """
  multiqc . \
    -f \
    --title "multisample multiqc" \
    --filename "multisample_multiqc.html" \
    -m fastp \
    -m rseqc
  """
}

/*
* Sort and index bam file using samtools
*/
process sort_index_bam {
  tag "$sample_id"

  input:
  tuple val (sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.antisense.txt"), emit: antisense_out
  tuple val (sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: sort_index_bam_out

  script:
  """
  samtools view -c -f 16 ${bam} > ${sample_id}.antisense.txt
  samtools sort ${bam} -O BAM -o ${sample_id}_sorted.bam
  samtools index ${sample_id}_sorted.bam
  """
}

/*
* Deduplicate reads (currently based on the combination of IO and SSS start site)
*/
process dedup{
  tag "$sample_id"

  publishDir "${params.outdir}/deduplication", mode: 'copy'

  input:
  tuple val (sample_id), path(bam), path(bai), val(alignment_count)

  output:
  tuple val(sample_id), path("${sample_id}.dedup.log"), emit: io_dedup_log
  tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: io_dedup_sam
  
  script:
  """
  # The effect of this will be to deduplicate based on the combination of cell
  # barcode and SSS start site.
  if [[ $alignment_count > 0 ]]
    then
      umi_tools dedup \
        --per-cell \
        --in-sam -I ${bam} \
        --log=${sample_id}.dedup.log > ${sample_id}.dedup.bam
    else
    
    # Then there were no alignments and we should output a dummy .dedup.log
    # and the original bam renamed
      echo "INFO Reads: Input Reads: 0\nINFO Number of reads out: 0\n" > ${sample_id}.dedup.log
      cp $bam ${sample_id}.dedup.bam
  fi
  """
}

/*
* Generate various count file precursors
*/
process io_count {
  tag "$sample_id"

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*_bcGeneSummary.txt'), emit: io_count_out

  shell:
  '''
  ## 1. Generate file for count matrix
  # alternative to umi_tools count
  # command: View the input dedup.sam file, grep for reads with gene assignment (which have the 'XT:' tag)
  # then keep the 1st (read_ID) and 18th fields (gene name, should be the 'XT:' tag itself) only
  # then split the string by '_' (-d '_'), keep fields 2,3,4 which correspond to 'io_sequence', 'gene_name' (4th field is so that we don't lose gene_names containing one '_')
  # use sed to replace the first '_' in each line, and any 'XT:Z:' strings with empty string with sed
  
  # output has 2 columns: io_sequence and gene_name for every deduplicated alignments with gene assignment
  samtools view !{f} | awk '/XT:/ {match($1, /_[A-Z]+_$/); printf substr($0,RSTART+1,RLENGTH-2); match($0, /XT:Z:[A-Za-z0-9_]+/); print "\t" substr($0,RSTART+5,RLENGTH-5)}' > !{sample_id}_bcGeneSummary.txt
  '''

}

/*
* Generate a h5ad count matrix
* output raw_count_matrix (unfiltered for single-cells only)
* Will output ${sample_id}.count_matrix.empty.h5ad if the input file is blank
* Otherwise the file will be called ${sample_id}.count_matrix.h5ad
* The tripartite count table files will only be output if the input
* file is not empty.
*/
process count_matrix {
  tag "$sample_id"

  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "features.tsv.gz"

  input:
  tuple val(sample_id), path(input_file)
  path(barcode_list)
  path(features_file)

  output:
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}.raw_feature_bc_matrix.*h5ad"), emit: h5ad
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true

  script:
  def mixed_args = params.mixed_species ? "--mixed_species True --hsap_mito_chr ${params.hsap_mitochondria_chromosome} --mmus_mito_chr ${params.mmus_mitochondria_chromosome} --hsap_gene_prefix ${params.hsap_gene_prefix} --mmus_gene_prefix ${params.mmus_gene_prefix}" : "--mixed_species False --mito_chr ${params.mitochondria_chromosome}"
  """
  count_matrix.py --barcode_list ${barcode_list} --count_table ${input_file} --gene_list ${features_file} --sample ${sample_id} $mixed_args
  """
}

/*
* Run cell caller - this determines a threshold number of counts to call a cell.
* The threshold is output to stdout and output as cell_caller_out
* If the h5ad matrix is empty, an empty .png will be output and checked for
* in the summary statistic script causing the Cell Caller plot to be hidden.
*/
process cell_caller {
  tag "$sample_id"

  // Publish the plots; the glob pattern is used to collect the mixed species and single species plots
  // Single species are named: {self.sample_name}_pdf_with_cutoff.png
  // Mixed species are named: {self.sample_name}_hsap_pdf_with_cutoff.png and {self.sample_name}_mmus_pdf_with_cutoff.png
  publishDir "${params.outdir}/plots", pattern: "${sample_id}*_pdf_with_cutoff.png", mode: 'copy'

  input:
  tuple val(sample_id), path(count_matrix_h5ad), val(manual_threshold_str)

  output:
  tuple val(sample_id), stdout, emit: cell_caller_out
  tuple val(sample_id), path("${sample_id}_counts_pdf_with_threshold.html"), path("${sample_id}_barnyard_plot.html"), emit: cell_caller_plots
  // publiDir statement only works on files output in the ouput directive
  tuple val(sample_id), path("${sample_id}*_pdf_with_cutoff.png"), emit: cell_caller_png

  script:
  """
  cell_caller.py --sample_name ${sample_id} --minimum_count_threshold ${params.minimum_count_threshold} --count_matrix ${count_matrix_h5ad} --single_species ${!params.mixed_species} --manual_threshold_str $manual_threshold_str
  """
}

/*
* Filter count table - produce a count table that has been filtered for barcodes that pass the cell caller threshold
* If an empty file is input as the .h5ad (due to fail at count_matrix), an "${sample_id}.*.cell_only.count_matrix.empty.h5ad"
* file will be created and carried into summary_metrics.py. The tripartite matrix files will not be ouput in this case.
*/
process filter_count_matrix{
  tag "$sample_id"

  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "*.filtered_feature_bc_matrix.h5ad"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "features.tsv.gz"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "*.raw_feature_bc_matrix.h5ad"

  input:
  tuple val(sample_id), val(count_threshold), path(h5ad_raw_count_matrix)

  output:
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}*.filtered_feature_bc_matrix.*h5ad")
  // Output and emit the raw h5ad matrix for use in summary statistics
  tuple val(sample_id), path("${sample_id}*.raw_feature_bc_matrix.*h5ad"), emit: raw_count_matrix

  script:
  def mixed_args = params.mixed_species ? "TRUE" : "FALSE"
  """
  filter_count_matrix.py ${count_threshold} ${h5ad_raw_count_matrix} ${sample_id} $mixed_args
  """
}

/*
* Categorize reads based on whether they are associated with cellular or non-cellular barcodes
* Also calculate counts in cells vs out of cells
*/
process categorize_reads {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(star_bam), path(raw_count_matrix_h5ad), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}.read_categorization.csv"), emit: read_categories

  script:
  def barcode_length = params.barcode_pattern.count('C')
  def mixed_flag = params.mixed_species ? "--mixed_species" : ""
  """
  categorize_reads.py \
    --sample_id ${sample_id} \
    --star_bam ${star_bam} \
    --raw_count_matrix_h5ad ${raw_count_matrix_h5ad} \
    --fastp_json ${fastp_json} \
    --barcode_length ${barcode_length} \
    ${mixed_flag}
  """
}

/*
* Generate the summary statistics required for HTML report
*/
process summary_statistics {
  tag "$sample_id"

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy', pattern: "*.csv"

  input:
  tuple val(sample_id), val(minimum_count_threshold), path(raw_h5ad), path(antisense), path(dedup), path("${sample_id}.multiqc.data.json"), path("${sample_id}_raw_rseqc_results.txt"), path("${sample_id}_annotated_rseqc_results.txt"), path(read_categorization_csv), path(qc_log)
  output:
  tuple val(sample_id), path("${sample_id}.metrics.csv"), emit: metrics_csv

  script:
  def mixed_flag = params.mixed_species ? "--mixed-species" : ""
  """
  summary_statistics.py \
    --sample-id ${sample_id} \
    --h5ad ${raw_h5ad} \
    --multiqc-json ${sample_id}.multiqc.data.json \
    --antisense ${antisense} \
    --dedup-log ${dedup} \
    --raw-rseqc ${sample_id}_raw_rseqc_results.txt \
    --annotated-rseqc ${sample_id}_annotated_rseqc_results.txt \
    --read-categorization ${read_categorization_csv} \
    --qc-log ${qc_log} \
    ${mixed_flag}
  """
}

/*
* Generate a per sample html report 
*/
process single_summary_report {
  tag "$sample_id"

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy'
  
  input:
  tuple val(sample_id), path(metrics_csv), path(pdf_plot_html), path(barnyard_plot_html)
  path(html_template)

  output:
  tuple val(sample_id), path("${sample_id}_report.html")
  path("${sample_id}.metrics.csv"), emit: single_sample_metric_out

  script:
  """
  create_single_sample_report.py $sample_id $pdf_plot_html $barnyard_plot_html $metrics_csv $html_template ${params.mixed_species}
  """
}

/*
* Generate an experiment summary report containing statistics for multiple samples
*/
process multi_sample_report {
  publishDir "${params.outdir}/report/", mode: 'copy'

  input:
  path(csvs)
  path(template)

  output:
  path('multisample_report.html')
  path('multisample_out.csv')

  script:
  """
  create_multi_sample_report.py $template ${params.mixed_species}
  """
}

