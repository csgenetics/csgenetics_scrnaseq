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
  tuple val (sample_id), path(f)
  path(index)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: star_out_bams

  shell:
  '''
  # alignment filtering step to be performed by samtools

  STAR --runThreadN 8 \
    --genomeDir !{index} \
    --readFilesIn !{f} \
    --outFileNamePrefix !{sample_id}_ \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --outSAMattributes Standard
  '''

}

/*
* Run Qualimap, generates a report, also useful for getting related mapping composition statistics
*/
process raw_qualimap {
  tag "$sample_id"
  
  label 'c8m16'

  publishDir "${params.outdir}/qualimap", mode: 'copy'

  input:
  tuple val (sample_id), path(bam)
  path(gtf)

  output:
  tuple val (sample_id), path("${sample_id}.raw_qualimap.txt"), emit: qualimap_txt

  shell:
  '''
  qualimap rnaseq -outdir !{sample_id}_raw_qualimap -a proportional -bam !{bam} -p strand-specific-forward -gtf !{gtf} --java-mem-size=16G 
  mv !{sample_id}_raw_qualimap/rnaseq_qc_results.txt !{sample_id}.raw_qualimap.txt
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
tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.bam"), emit: filtered_bam

script:
"""
samtools view -h -e '[NH]==1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b $star_out_bam > ${sample_id}.mapped.sorted.filtered.bam 
"""
}


/*
* Run Qualimap, generates a report, also useful for getting related mapping composition statistics
*/
process filtered_qualimap {
  tag "$sample_id"
  label 'c8m16'

  publishDir "${params.outdir}/qualimap", mode: 'copy'

  input:
  tuple val(sample_id), path(bam)
  path(gtf)

  output:
  tuple val (sample_id), path("${sample_id}.filtered_qualimap.txt"), emit: qualimap_txt

  shell:
  '''
  qualimap rnaseq -outdir !{sample_id}_filtered_qualimap -a proportional -bam !{bam} -p strand-specific-forward -gtf !{gtf} --java-mem-size=16G 
  mv !{sample_id}_filtered_qualimap/rnaseq_qc_results.txt !{sample_id}.filtered_qualimap.txt
  '''
}

/*
* Quantify genes using feature counts
*/
process feature_counts {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/featureCounts", mode: 'copy'

  input:
  tuple val(sample_id), path(bam)
  path(gtf)

  output:
  tuple val(sample_id), path('*.bam'), emit: feature_counts_out_bam

  shell:
  """
  featureCounts -a $gtf -o ${sample_id}_gene_assigned.txt -R BAM $bam -T 4 -t exon,intron -g gene_id --fracOverlap 0.5 --extraAttributes gene_name
  """
}

process filter_for_annotated {
  tag "$sample_id"

  label "c2m2"

  input:
  tuple val(sample_id), path(feature_counts_out_bam)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), emit: annotated_bam

  script:
  """
  samtools view -h -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b $feature_counts_out_bam > ${sample_id}.mapped.sorted.filtered.annotated.bam 
  """
}

process annotated_qualimap {
  tag "$sample_id"
  label 'c8m16'

  publishDir "${params.outdir}/qualimap", mode: 'copy'

  input:
  tuple val(sample_id), path(bam)
  path(gtf)

  output:
  tuple val (sample_id), path("${sample_id}.annotated_qualimap.txt"), emit: annotated_qualimap

  shell:
  '''
  qualimap rnaseq -outdir !{sample_id}_annotated_qualimap -a proportional -bam !{bam} -p strand-specific-forward -gtf !{gtf} --java-mem-size=16G 
  mv !{sample_id}_annotated_qualimap/rnaseq_qc_results.txt !{sample_id}.annotated_qualimap.txt
  '''
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
  tuple val (sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.antisense.txt"), emit: antisense_out
  tuple val (sample_id), path('*.{bam,bai}'), emit: sort_index_bam_out

  shell:
  '''
  samtools view -c -f 16 !{bam} > !{sample_id}.antisense.txt
  samtools sort !{bam} -o !{sample_id}_sorted.bam
  samtools index !{sample_id}_sorted.bam
  '''
}

/*
* Deduplicate reads (currently based on the combination of IO and SSS start site)
*/
process dedup{
  tag "$sample_id"
  label 'c4m1'

  publishDir "${params.outdir}/io_count", mode: 'copy'

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path("${sample_id}.dedup.log"), emit: io_dedup_log
  tuple val(sample_id), path("${sample_id}.dedup.sam"), emit: io_dedup_sam
  
  shell:
  '''
  # The effect of this will be to deduplicate based on the combination of cell
  # barcode and SSS start site.

  umi_tools dedup \
    --per-cell \
    --in-sam -I !{f} \
    --out-sam -S !{sample_id}.dedup.sam \
    --log=!{sample_id}.dedup.log
  '''
}

/*
* Generate various count file precursors
*/
process io_count {
  tag "$sample_id"
  label 'c4m1'

  publishDir "${params.outdir}/io_count", mode: 'copy'

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*_bcGeneSummary.txt'), emit: io_count_out
  tuple val(sample_id), path('*_goodUMRcount.txt'),  emit: io_goodumr_count

  shell:
  '''
  ## 1. Generate file for count matrix
  # alternative to umi_tools count
  # command: View the input dedup.sam file, grep for reads with gene assignment (which have the 'XT:' tag)
  # then keep the 1st (read_ID) and 18th fields (gene name, should be the 'XT:' tag itself) only
  # then split the string by '_' (-d '_'), keep fields 2,3,4 which correspond to 'io_sequence', 'gene_name' (4th field is so that we don't lose gene_names containing one '_')
  # use sed to replace the first '_' in each line, and any 'XT:Z:' strings with empty string with sed

  # output has 2 columns: io_sequence and gene_name for every deduplicated alignments with gene assignment

  samtools view !{f} | grep 'XT:' | sed 's/mm10___/mm10_/g' |cut -f 1,18 | cut -f 2,3,4 -d '_' | sed 's/_//' | sed 's/XT:Z://g'> !{sample_id}_bcGeneSummary.txt

  ## 2. Generate deduplicated reads per IO file 
  # similar command as above, but include all alignments, instead of just those with gene assignment
  samtools view !{f} | cut -f 1,3| cut -f 2,3 -d '_' | sed 's/_//g' | sort | uniq -c > !{sample_id}_goodUMRcount.txt
  '''

}

/*
* Generate a h5ad count matrix
* output raw_count_matrix (unfiltered for single-cells only)
*/
process count_matrix {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/count_matrix/raw_count_matrix/${sample_id}", mode: 'copy'
  
  input:
  tuple val(sample_id), path(f)
  path(w)
  path(features_file)


  output:
  tuple val(sample_id), path("${sample_id}.count_matrix.h5ad"), emit: h5ad
  tuple val(sample_id), path("${sample_id}.matrix.mtx.gz"), emit: raw_matrix
  tuple val(sample_id), path("${sample_id}.barcodes.tsv.gz"), emit: raw_barcodes
  tuple val(sample_id), path("${sample_id}.features.tsv.gz"), emit: raw_features

  script:
  """
  count_matrix.py --white_list ${w} --count_table ${f} --gene_list ${features_file} --sample ${sample_id}
  """
}

/*
* Run cell caller - this determines a threshold number of nuclear genes detected to call a cell.
* The threshold is output to stdout and output as cell_caller_out
*/
process cell_caller {
  tag "$sample_id"
  label 'c4m2'

  publishDir "${params.outdir}/plots", pattern: '*.png', mode: 'copy'

  input:
  tuple val(sample_id), path(count_matrix_h5ad)

  output:
  // returns the call_caller-calculated nuclear gene threshold
  tuple val(sample_id), stdout, emit: cell_caller_out 
  tuple val(sample_id), path('*.png'), emit: cell_caller_plot

  script:
  """
  cell_caller.py --sample ${sample_id} --min_nucGene ${params.min_nuc_gene} --count_matrix $count_matrix_h5ad
  """
}

/*
* Filter count table - produce a count table that has been filtered for barcodes that pass the cell caller threshold
*/
process filter_count_matrix{
  tag "$sample_id"
  label 'c2m2'

  publishDir "${params.outdir}/count_matrix/cell_only_count_matrix/${sample_id}/", mode: 'copy'

  input:
  tuple val(sample_id), val(nuc_gene_threshold), path(h5ad_raw_count_matrix)

  output:
  // Output the 3 part barcode, features, matrix 
  tuple val(sample_id), path("${sample_id}.*.cell_only.barcodes.tsv.gz"), path("${sample_id}.*.cell_only.features.tsv.gz"), path("${sample_id}.*.cell_only.matrix.mtx.gz")
  // Output the h5ad matrix
  tuple val(sample_id), path("${sample_id}.*.cell_only.count_matrix.h5ad"), emit: cell_only_count_matrix

  script:
  """
  filter_count_matrix.py ${nuc_gene_threshold} ${h5ad_raw_count_matrix} ${sample_id}
  """
}

/*
* Generate the summary statistics required for HTML report
*/
process summary_statistics {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy', pattern: "*.csv"
  
  input:
  tuple val(sample_id), val(min_nuc_gene_cutoff), path(h5ad), path(annotated_qualimap), path(antisense), path(dedup), path(filtered_qualimap), path(multiqc_data_json), path(raw_qualimap)

  output:
  tuple val(sample_id), path("${sample_id}.metrics.csv"), emit: metrics_csv

  script:
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
  path(cs_logo)

  output:
  tuple val(sample_id), path("${sample_id}_report.html")
  path("${sample_id}.metrics.csv"), emit: single_sample_metric_out

  script:
  """
  create_single_sample_report.py $sample_id $plot_png $metrics_csv $html_template
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
  create_multi_sample_report.py $template
  """
}

