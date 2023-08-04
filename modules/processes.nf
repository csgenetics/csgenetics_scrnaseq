#!/usr/bin/env nextflow

/*
* This script contains all of the processes used by the daily R&D pipeline
*/

nextflow.enable.dsl=2

/*
* Convert GTF into features file
*/
process features_file {
  label 'c2m8'

  input:
  path(gtf)

  output:
  path("${gtf.baseName}_features_names.tsv"), emit: modified_gtf

  shell:
  '''
  features_names.py !{gtf} !{gtf.baseName}_features_names_tmp.tsv
  # Modile 10x triple underscore for mm10 as this breaks the pipeline
  sed 's/mm10___/mm10_/g' !{gtf.baseName}_features_names_tmp.tsv > !{gtf.baseName}_features_names.tsv
  '''
}

/*
* Merge Illumina Lanes
*/
process merge_lanes {
  tag "$sample_id"
  label 'c1m4'

  input:
  tuple val(sample_id), path(fastq_1), path(fastq_2)

  output:
  tuple val(sample_id), path("${sample_id}_R1.merged.fastq.gz"), path("${sample_id}_R2.merged.fastq.gz"), emit: merge_lanes_out

  shell:
  '''
  cat !{fastq_1} > !{sample_id}_R1.merged.fastq.gz
  cat !{fastq_2} > !{sample_id}_R2.merged.fastq.gz
  '''
}

process merged_fastp{
  tag "$sample_id"
  label 'c4m2'

  publishDir "${params.outdir}/fastp", mode: 'copy'

  input:
  tuple val(sample_id), path(r1), path(r2)
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
  fastp -i !{r1} \
    -A -G \
    -j !{sample_id}_R1.merged_fastp.json \
    -h !{sample_id}_R1.merged_fastp.html

  # R2 fastp
  fastp -i !{r2} \
    -A -G \
    -b !{barcode_length} \
    -l 13 \
    -j !{sample_id}_R2.merged_fastp.json \
    -h !{sample_id}_R2.merged_fastp.html
  '''
}

/*
* Use umitools to label reads with IOs from the parameter provided whitelist
*/
process io_extract {
  tag "$sample_id"
  label 'c2m4'

  input:
  tuple val(sample_id), path(r1), path(r2)
  path(whitelist)
  val(barcode_pattern)

  output:
  tuple val(sample_id), path("${sample_id}_R1.io_extract.fastq.gz"), path("${sample_id}_R2.io_extract.fastq.gz"), emit: io_extract_out
  tuple val(sample_id), path("extract.log"), emit: io_extract_log

  script:
  """
  cat $whitelist | cut -d ',' -f2 > wl.txt

  umi_tools extract -I $r2 \
    --bc-pattern=${barcode_pattern} \
    --read2-in=$r1 \
    --stdout=${sample_id}_R2.io_extract.fastq.gz \
    --read2-out=${sample_id}_R1.io_extract.fastq.gz \
    --whitelist=wl.txt \
    -L extract.log
  """
}

/*
* Use fastp to trim reads and remove low quality reads
* https://github.com/OpenGene/fastp
* trim 5' sss_nmer
* trim 3' polyA and polyG of length >=15
* remove reads of length <=20
* disable adapter trimming
* disable polyG trimming 
*/
process io_extract_fastp {
  tag "$sample_id"

  label 'c4m2'

  publishDir "${params.outdir}/fastp", pattern: '*.{json,html}', mode: 'copy'

  input: tuple val(sample_id), path(r1), path(r2)

  output:
  tuple val(sample_id), path("${sample_id}_R1.io_extract.fastp.fastq.gz"), emit: fastp_out
  tuple val(sample_id), path("${sample_id}_R1.io_extract.fastp.html"), path("${sample_id}_R1.io_extract.fastp.json"), emit: fastp_multiqc

  shell:
  '''
  # SSS trimming
  # 3' polyX trimming with min length of 15 bases
  # remove reads of length <= 20 bases
  # disable adapter trimming
  # disable polyG trimming 
  
  fastp -i !{r1} \
    -f !{params.sss_nmer} \
    -x --poly_x_min_len 15 \
    -l 20 \
    -A \
    -G \
    -j !{sample_id}_R1.io_extract.fastp.json \
    -h !{sample_id}_R1.io_extract.fastp.html \
    --stdout \
     2> fastp.log | gzip > !{sample_id}_R1.io_extract.fastp.fastq.gz
  '''     
}

/*
* Trim extra polyA tail using a bespoke python scripts
* Fastp only trims polyX tails if it is right at the end of 3' end
* This script finds and trims polyA (of length >=15) regardless of where it is in the read. eg. xxxxxxx[15As]CTG --> xxxxxxx
* Also remove reads of length <=20 post polyA trimming
*/
process trim_extra_polya {
  tag "$sample_id"
  label 'c2m4'

  input:
  tuple val(sample_id), path(fastq)

  output:
  tuple val(sample_id), path("${sample_id}_R1.polyA_trimmed.fastq.gz"), emit: trim_extra_polya_out
  tuple val(sample_id), path("${sample_id}_trim_polyA_metrics.csv"), emit: trim_extra_polya_log1

  script:
  """
  trim_extra_polya.py $fastq $sample_id
  """
}

process post_polyA_fastp{
  tag "$sample_id"

  label 'c4m2'

  publishDir "${params.outdir}/fastp", pattern: '*.{json,html}', mode: 'copy'

  input: tuple val(sample_id), path(r1)

  output:
  tuple val(sample_id), path("${sample_id}_R1.post_polyA_fastp.fastq.gz"), emit: fastp_out
  tuple val(sample_id), path("${sample_id}_R1.post_polyA_fastp.html"), path("${sample_id}_R1.post_polyA_fastp.json"), emit: fastp_multiqc

  script:
  """
  # remove reads of length <= 20 bases
  # disable adapter trimming
  # disable polyG trimming 
  
  fastp -i $r1 \
    -l 20 \
    -A \
    -G \
    -j ${sample_id}_R1.post_polyA_fastp.json \
    -h ${sample_id}_R1.post_polyA_fastp.html \
    --stdout \
     2> fastp.log | gzip > ${sample_id}_R1.post_polyA_fastp.fastq.gz
  """
}

/*
* Align reads with STAR
*/
process star {
  tag "$sample_id"
  label 'c8m40'

  publishDir "${params.outdir}/STAR", mode: 'copy'

  input:
  tuple val(sample_id), path(r1)
  path(index)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), env(uniquely_mapped_reads), emit: out_bam

  script:
  """
      STAR --runThreadN 8 \
        --genomeDir ${index} \
        --readFilesIn ${r1} \
        --outFileNamePrefix ${sample_id}_ \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --outSAMattributes Standard;
      
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

  label 'c1m1'

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
* Generic process for running qualimap.
*/
process run_qualimap {
  tag "$sample_id"
  
  label 'c8m16'

  publishDir "${params.outdir}/qualimap", mode: 'copy'

  input:
  tuple val (sample_id), path(bam), val(count)
  path(gtf)
  path(empty_qualimap_template)
  val(prefix)

  output:
  tuple val (sample_id), path("${sample_id}.${prefix}_qualimap.txt"), emit: qualimap_txt

  shell:
  '''
  if [[ !{count} > 0 ]]
    then
      qualimap rnaseq -outdir !{sample_id}_!{prefix}_qualimap -a proportional -bam !{bam} -p strand-specific-forward -gtf !{gtf} --java-mem-size=16G 
      mv !{sample_id}_!{prefix}_qualimap/rnaseq_qc_results.txt !{sample_id}.!{prefix}_qualimap.txt
    else
      BAMNAME=!{bam}
      GTFNAME=!{gtf}
      export BAMNAME GTFNAME
      cat !{empty_qualimap_template} | envsubst > !{sample_id}.!{prefix}_qualimap.txt
  fi
  '''
}

/*
* Filter for reads with 1 alignment and max of 3 differences from reference seq.
*/
process filter_umr_mismatch{
tag "$sample_id"

label "c2m2"

input:
tuple val(sample_id), path(star_out_bam)

output:
tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.bam"), env(alignment_count), emit: filtered_bam

script:
"""
samtools view -h -b -e '[NH]==1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b $star_out_bam > ${sample_id}.mapped.sorted.filtered.bam
alignment_count=\$(samtools view -c ${sample_id}.mapped.sorted.filtered.bam)
"""
}

/*
* Quantify genes using feature counts
*/
process feature_counts {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/featureCounts", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), val(aligned_count)
  path(gtf)

  output:
  tuple val(sample_id), path('*.bam'), emit: out_bam

  shell:
  """
  if [[ $aligned_count > 0 ]] # If the bam is not empty
    then
      # Run featureCounts as normal
      featureCounts -a $gtf -o ${sample_id}_gene_assigned.txt -R BAM $bam -T 4 -t exon,intron -g gene_id --fracOverlap 0.5 --extraAttributes gene_name
    else
      # Simply rename the input bam so that it can be collected
      cp $bam ${sample_id}.mapped.sorted.filtered.bam.featureCounts.bam
  fi
  """
}

process filter_for_annotated {
  tag "$sample_id"

  label "c2m2"

  input:
  tuple val(sample_id), path(feature_counts_out_bam)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), env(alignment_count), emit: annotated_bam

  script:
  """
  samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b $feature_counts_out_bam > ${sample_id}.mapped.sorted.filtered.annotated.bam 
  alignment_count=\$(samtools view -c ${sample_id}.mapped.sorted.filtered.annotated.bam)
  """
}

/*
* Run multiqc
*/
process multiqc {
  label 'c1m1'

  publishDir "${params.outdir}/multiqc/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(multiqc_in_files)

  output:
  path "*_multiqc.html"
  path "*_data"
  tuple val(sample_id), path("${sample_id}.multiqc.data.json"), emit: multiqc_json

  script:
  """
  multiqc . \
    -f \
    --title "${sample_id} multiqc" \
    --filename "${sample_id}_multiqc.html" \
    -m fastqc \
    -m fastp \
    -m star \
    -m featureCounts
  mv ${sample_id}_multiqc_data/multiqc_data.json ./${sample_id}.multiqc.data.json
  """
}

/*
* Sort and index bam file using samtools
*/
process sort_index_bam {
  tag "$sample_id"
  label 'c1m2'

  input:
  tuple val (sample_id), path(bam), val(count)

  output:
  tuple val(sample_id), path("${sample_id}.antisense.txt"), emit: antisense_out
  tuple val (sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), env(alignment_count), emit: sort_index_bam_out

  script:
  """
  samtools view -c -f 16 ${bam} > ${sample_id}.antisense.txt
  samtools sort ${bam} -O BAM -o ${sample_id}_sorted.bam
  samtools index ${sample_id}_sorted.bam
  alignment_count=\$(samtools view -c ${sample_id}_sorted.bam)
  """
}

/*
* Deduplicate reads (currently based on the combination of IO and SSS start site)
*/
process dedup{
  tag "$sample_id"
  label 'c4m1'

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
  label 'c4m1'
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
  samtools view !{f} | awk '/XT:/ {match($0, /_[A-Z]+_/); printf substr($0,RSTART+1,RLENGTH-2); match($0, /XT:Z:[A-Za-z0-9_]+/); print "\t" substr($0,RSTART+5,RLENGTH-5)}' > !{sample_id}_bcGeneSummary.txt
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
  label 'c4m4'

  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "*.raw_feature_bc_matrix.h5ad"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "features.tsv.gz"

  input:
  tuple val(sample_id), path(input_file)
  path(whitelist)
  path(features_file)

  output:
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}.raw_feature_bc_matrix.*h5ad"), emit: h5ad
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true

  script:
  """
  count_matrix.py --white_list ${whitelist} --count_table ${input_file} --gene_list ${features_file} --sample ${sample_id}
  """
}

/*
* Run cell caller - this determines a threshold number of nuclear genes detected to call a cell.
* The threshold is output to stdout and output as cell_caller_out
* If the h5ad matrix is empty, an empty .png will be output and checked for
* in the summary statistic script causing the Cell Caller plot to be hidden.
*/
process cell_caller {
  tag "$sample_id"
  label 'c4m2'

  publishDir "${params.outdir}/plots", pattern: "${sample_id}_pdf_with_cutoff.png", mode: 'copy'

  input:
  tuple val(sample_id), path(count_matrix_h5ad)

  output:
  tuple val(sample_id), stdout, emit: cell_caller_out
  tuple val(sample_id), path("${sample_id}_pdf_with_cutoff.png"), emit: cell_caller_plot

  script:
  def mixed_args = params.mixed_species ? "--mixed TRUE --mt_prefix ${params.hsap_mitochondria_prefix} --mt_prefix2 ${params.mmus_mitochondria_prefix}" : "--mixed FALSE --mt_prefix ${params.mitochondria_prefix}"
  """
  cell_caller.py --sample ${sample_id} --min_nucGene ${params.min_nuc_gene} --count_matrix $count_matrix_h5ad --mt_chromosome '${params.mt_chromosome}'
  """
}

/*
* Filter count table - produce a count table that has been filtered for barcodes that pass the cell caller threshold
* If an empty file is input as the .h5ad (due to fail at count_matrix), an "${sample_id}.*.cell_only.count_matrix.empty.h5ad"
* file will be created and carried into summary_metrics.py. The tripartite matrix files will not be ouput in this case.
*/
process filter_count_matrix{
  tag "$sample_id"
  label 'c2m2'

  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "*.filtered_feature_bc_matrix.h5ad"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "matrix.mtx.gz"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "barcodes.tsv.gz"
  publishDir "${params.outdir}/count_matrix/filtered_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "features.tsv.gz"

  input:
  tuple val(sample_id), val(nuc_gene_threshold), path(h5ad_raw_count_matrix)

  output:
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("barcodes.tsv.gz"), path("features.tsv.gz"), path("matrix.mtx.gz"), optional: true
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}*.filtered_feature_bc_matrix.*h5ad")
  // Output and emit the raw h5ad matrix for use in summary statistics
  tuple val(sample_id), path("${sample_id}*.raw_feature_bc_matrix.*h5ad"), emit: raw_count_matrix

  script:
  def mixed_args = params.mixed_species ? "TRUE ${params.hsap_mitochondria_prefix} ${params.mmus_mitochondria_prefix} ${params.hsap_gene_prefix} ${params.mmus_gene_prefix} ${params.purity}" : "FALSE ${params.mitochondria_prefix}"
  """
<<<<<<< HEAD
  filter_count_matrix.py ${nuc_gene_threshold} ${h5ad_raw_count_matrix} ${sample_id} $mixed_args
=======
  filter_count_matrix.py ${nuc_gene_threshold} ${h5ad_raw_count_matrix} ${sample_id} '${params.mt_chromosome}'
>>>>>>> add chromosome MT to anndata
  """
}

/*
* Generate the summary statistics required for HTML report
*/
process summary_statistics {
  tag "$sample_id"
  label 'c1m8'

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy', pattern: "*.csv"
  
  input:
  tuple val(sample_id), val(min_nuc_gene_cutoff), path(raw_h5ad), path(annotated_qualimap), path(antisense), path(dedup), path(filtered_qualimap), path(multiqc_data_json), path(raw_qualimap)
  output:
  tuple val(sample_id), path("${sample_id}.metrics.csv"), emit: metrics_csv

  script:
  def mixed_args = params.mixed_species ? "TRUE ${params.hsap_mitochondria_prefix} ${params.mmus_mitochondria_prefix} ${params.hsap_gene_prefix} ${params.mmus_gene_prefix} ${params.purity}" : "FALSE ${params.mitochondria_prefix}"
  """
  summary_statistics.py $sample_id $h5ad $multiqc_data_json $antisense $dedup $raw_qualimap $filtered_qualimap $annotated_qualimap
  """
}

/*
* Generate a per sample html report 
*/
process single_summary_report {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy'
  
  input:
  tuple val(sample_id), path(metrics_csv), path(plot_png)
  path(html_template)

  output:
  tuple val(sample_id), path("${sample_id}_report.html")
  path("${sample_id}.metrics.csv"), emit: single_sample_metric_out

  script:
  """
  create_single_sample_report.py $sample_id $plot_png $metrics_csv $html_template ${params.mixed_species}
  """
}

/*
* Generate an experiment summary report containing statistics for multiple samples
*/
process multi_sample_report {
  label 'c2m4'

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

