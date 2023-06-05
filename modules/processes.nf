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
* Run BaseSpace Download
*/
process basespace {
  label 'c1m4'

  output:
  path('*.fastq.gz')

  when:
  !params.fastq

  shell:
  '''
  bs download project \
    --extension=fastq.gz \
    --api-server=!{params.bs_api_server} \
    --access-token=!{params.bs_access_token} \
    -i !{params.bs_project_id} -o .

  # mv fastq files from subfolders to one level higher
  find . -mindepth 2 -type f -print -exec mv {} . \\;
  '''
}

/*
* Merge Illumina Lanes
*/

process merge_lanes {
  tag "$sample_id"
  label 'c1m4'

  input:
  tuple val(sample_id), file(fastq_1), file(fastq_2)

  output:
  tuple val(sample_id), env(numreads), path('merged'), emit: merge_lanes_out
  tuple val(sample_id), env(numreads), path('numreads.txt'), emit: numreads_log

  shell:
  '''
  mkdir merged
  # merging R1 and R2
  zcat !{sample_id}_L*_R1*.f*q.gz | gzip > merged/!{sample_id}_R1.fastq.gz
  zcat !{sample_id}_L*_R2*.f*q.gz | gzip > merged/!{sample_id}_R2.fastq.gz
  
  # count number of raw reads
  numreads=$(( $(zcat merged/!{sample_id}_R1.f*q.gz | wc -l) / 4))
  echo $numreads > numreads.txt
  '''
}


/*
* Run FASTQC - generates a report, also useful for getting related base composition statistics
*/
process fastqc {
  tag "$sample_id"
  label 'c2m4'

  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  tuple val(sample_id), path(f)

  output:
  path ('*.{zip,html}'), emit: fastqc_multiqc
  tuple val(sample_id), path("${sample_id}_R1_fastqc.zip"), emit: fastqc_out

  shell:
  '''
  in=merged
  fastqc -q -o . $in/*.fastq.gz
  '''
}

/*
* Use fastp to trim R2 reads equal to barcode length to get q30 rate for barcode region only 
* https://github.com/OpenGene/fastp
*/
process barcode {
  tag "$sample_id"
  label 'c4m2'

  input:
  tuple val (sample_id), path(f)
  val (barcode_pattern)

  output:
  path('*.{json,html}'), emit: barcode_multiqc

  shell:
  '''
  mkdir fastp_out
     
  # Get q30_rate for barcode present in fastq2 file 
  # maximum length which be equal to barcode length 
  barcode_length=$(echo !{barcode_pattern} | tr -cd 'C' | wc -c)
  in=merged
  fastp -i $in/!{sample_id}_R2.fastq.gz \
    -l $barcode_length \
    -b $barcode_length \
    -j !{sample_id}_R2_fastp.json \
    -h !{sample_id}_R2_fastp.html 
  '''
}

/*
* Use umitools to label reads with IOs from the parameter provided whitelist
*/
process io_extract {
  tag "$sample_id"
  label 'c2m4'

  input:
  tuple val (sample_id), path(f)
  path(w)
  val (barcode_pattern)

  output:
  tuple val (sample_id), path('io_extract_out'), emit: io_extract_out
  tuple val (sample_id), path('extract.log'), emit: io_extract_log

  shell:
  '''
  # check input directory name  
  in=merged
  
  mkdir io_extract_out

  cat !{w} | cut -d ',' -f2 > wl.txt
  umi_tools extract -I $in/!{sample_id}_R2.fastq.gz \
    --bc-pattern=!{barcode_pattern} \
    --read2-in=$in/!{sample_id}_R1.fastq.gz \
    --stdout=io_extract_out/!{sample_id}_R2.fastq.gz \
    --read2-out=io_extract_out/!{sample_id}_R1.fastq.gz \
    --whitelist=wl.txt \
    -L extract.log
  
  '''
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
process fastp {
  tag "$sample_id"
  label 'c4m2'

  publishDir "${params.outdir}/fastp", pattern: '*.{json,html}', mode: 'copy'
  input: tuple val (sample_id), path(f)

  output:
  tuple val (sample_id), path('fastp_out'), emit: fastp_out
  path('*.{json,html}'), emit: fastp_multiqc
  tuple val (sample_id), path('fastp.log'), emit: fastp_log

  shell:
  '''
  mkdir fastp_out
  # SSS trimming
  # 3' polyX trimming with min length of 15 bases
  # remove reads of length <= 20 bases
  # disable adapter trimming
  # disable polyG trimming 
  
  fastp -i io_extract_out/!{sample_id}_R1.fastq.gz \
    -f !{params.sss_nmer} \
    -x --poly_x_min_len 15 \
    -l 20 \
    -A \
    -G \
    -j !{sample_id}_R1_fastp.json \
    -h !{sample_id}_R1_fastp.html \
    --stdout \
     2> fastp.log | gzip > fastp_out/!{sample_id}_R1.fastq.gz
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
  tuple val (sample_id), path('fastp_out')

  output:
  tuple val (sample_id), path('*.fastq.gz'), emit: trim_extra_polya_out
  tuple val (sample_id), path('trim_polyA_metrics.csv'), emit: trim_extra_polya_log1
  tuple val (sample_id), path('read_lengths_post_trimming.csv'), emit:trim_extra_polya_log2

  shell:
  '''
  trim_extra_polya.py fastp_out/!{sample_id}_R1.fastq.gz !{sample_id}
  '''
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
  tuple val (sample_id), path('*Aligned.sortedByCoord.out.bam'), emit: star_out
  path ("*.{out,tab,mate1}"), emit: star_multiqc

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
process qualimap {
  tag "$sample_id"
  label 'c8m16'

  publishDir "${params.outdir}/qualimap", mode: 'copy'

  input:
  tuple val (sample_id), path(f)
  path(gtf)

  output:
  tuple val (sample_id), path("*_qualimap.txt"), emit: qualimap_txt

  shell:
  '''
  qualimap rnaseq -outdir !{sample_id} -a proportional -bam !{f} -p strand-specific-forward -gtf !{gtf} --java-mem-size=16G 
  mv !{sample_id}/rnaseq_qc_results.txt !{sample_id}_qualimap.txt
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
  tuple val (sample_id), path(f)
  path(gtf)

  output:
  tuple val(sample_id), path('*.bam'), emit: feature_counts_out
  path("*.{txt,summary}"), emit: feature_counts_multiqc

  shell:
  '''
  featureCounts -a !{gtf} -o !{sample_id}_gene_assigned.txt -R BAM !{f} -T 4 -t exon,intron,intergenic -g gene_id --fracOverlap 0.5 --extraAttributes gene_name
  '''
}


/*
* Run multiqc
*/
process multiqc {
  label 'c1m1'

  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
  path(f)

  output:
  file "*_multiqc.html"
  file "*_data"
  path("*/*.json") , emit: multiqc_json

  script:
  def filename = "${params.bs_project_id}_multiqc.html"
  def reporttitle = "${params.bs_project_id} (csgx/rnaseq)"
  """
  multiqc . \
    -f \
    --title "$reporttitle" \
    --filename "$filename" \
    -m fastqc \
    -m fastp \
    -m fastp \
    -m star \
    -m featureCounts
  """
}


/*
* Sort and index bam file using samtools
*/
process sort_index_bam {
  tag "$sample_id"
  label 'c1m2'

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*.txt'), emit: antisense_out
  tuple val (sample_id), path('*.{bam,bai}'), emit: sort_index_bam_out

  shell:
  '''
  samtools view -c -f 16 !{f} > !{sample_id}_antisense.txt
  samtools sort !{f} -o !{sample_id}_sorted.bam
  samtools index !{sample_id}_sorted.bam
  '''
}

/*
* Group reads by cell ids/IOs using umitools
*/
process group {
  tag "$sample_id"

  label 'c4m2'

  publishDir "${params.outdir}/io_count", mode: 'copy'

  input:
  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*_group.sam'), emit: io_group_sam
  tuple val(sample_id), path('*_group_filtered.sam'), emit: io_group_filtered_sam
  tuple val(sample_id), path('*_group.tsv'), emit: io_group_tsv
  tuple val(sample_id), path('*_group.log'), emit: io_group_log

  shell:

  if (params.remove_singletons)
  '''
  # Originally this process was just umi_tools count.
  # We had to change it to remove singleton reads from the bam file.
  # We followed instructions from here:
  # https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html#step-6-counting-molecules

  # 1. run umi_tools group first
  # remove --per-gene  --gene-tag=XT  --assigned-status-tag=XS  to group reads by start position only, not by gene
  # add --read-length to group reads by start and end position

  umi_tools group \
    --read-length \
    --per-cell \
    -I !{sample_id}_sorted.bam \
    --group-out=!{sample_id}_group.tsv \
    --output-bam --out-sam -S !{sample_id}_group.sam \
    --log=!{sample_id}_group.log

    # 2. filter out singleton reads
    # https://github.com/CGATOxford/UMI-tools/issues/274
    grep @ !{sample_id}_group.sam > sam_header
    grep -v @ !{sample_id}_group.sam > sam_reads
    awk '$6 !~ /^1$/' !{sample_id}_group.tsv | awk '{print $1}' | tail -n+2 > selected_reads
    awk 'NR==FNR{a[$0]; next} $1 in a' selected_reads sam_reads > sam_selected_reads
    cat sam_header sam_selected_reads > !{sample_id}_group_filtered.sam

  '''
  else

  '''
  # Originally this process was just umi_tools count.
  # We had to change it to remove singleton reads from the bam file.
  # We followed instructions from here:
  # https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html#step-6-counting-molecules

  # run umi_tools group first
  # remove --per-gene  --gene-tag=XT  --assigned-status-tag=XS  to group reads by start position only, not by gene
  # add --read-length to group reads by start and end position

  umi_tools group \
    --read-length \
    --per-cell \
    -I !{sample_id}_sorted.bam \
    --group-out=!{sample_id}_group.tsv \
    --output-bam --out-sam -S !{sample_id}_group.sam \
    --log=!{sample_id}_group.log

  cp !{sample_id}_group.sam !{sample_id}_group_filtered.sam
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
  tuple val(sample_id), path('*_dedup.log'), emit: io_dedup_log
  tuple val(sample_id), path('*_dedup.sam'), emit: io_dedup_sam
  
  shell:
  '''
  # If we are NOT using unique IOs, remove the --per-gene and --gene-tag flags.
  # The effect of this will be to deduplicate based on the combination of cell
  # barcode and SSS start site. This has the caveat that we expect multiple SSS
  # products per input mRNA, so this will inflate some counts. Additionally, the
  # maximum counts detectable per gene will be equal to the gene length. This
  # is problematic for highly expressed genes, which will have deflated counts.

  umi_tools dedup \
    --read-length \
    --per-cell \
    --in-sam -I !{f} \
    --out-sam -S !{sample_id}_dedup.sam \
    --log=!{sample_id}_dedup.log
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

  publishDir "${params.outdir}/count_matrix/", mode: 'copy'
  
  input:
  tuple val (sample_id), path(f)
  path(w)
  path(features_file)


  output:
  tuple val(sample_id), path('*.h5ad'), emit: h5ad
  tuple val(sample_id), path('raw_count_matrix/*/matrix.mtx.gz'), emit: raw_matrix
  tuple val(sample_id), path('raw_count_matrix/*/barcodes.tsv.gz'), emit: raw_barcodes
  tuple val(sample_id), path('raw_count_matrix/*/features.tsv.gz'), emit: raw_features

  shell:
  '''
  mkdir -p raw_count_matrix/!{sample_id}
  count_matrix.py --white_list !{w} --count_table !{f} --gene_list !{features_file} --sample !{sample_id}
  '''
}



/*
* Run cell caller - this determines a threshold number of Nuclear gene detected per single cell, above which reagent_bead is considered to have successfully indexed a single cell
*/
process cell_caller {
  tag "$sample_id"
  label 'c4m2'

  publishDir "${params.outdir}/plots", pattern: '*.png', mode: 'copy'

  input:

  tuple val (sample_id), path(f)

  output:
  tuple val(sample_id), path('*.txt'), emit: cell_caller_out
  tuple val(sample_id), path('*.png'), emit: cell_caller_plot


  shell:
    '''
    cell_caller.py --sample !{sample_id} --min_nucGene !{params.min_nuc_gene} | cut -f2 -d, > !{sample_id}_cell_caller.txt
    '''
}



/*
* Generate a Summary report
* output filter_count_matrix and html summary report 
*/
process summary_report {
  tag "$sample_id"
  label 'c4m4'

  publishDir "${params.outdir}/report/", mode: 'copy'
  
  input:
  tuple val (sample_id), path(f)
  tuple val (sample_id), path(w)
  tuple val (sample_id), path(g)
  path(m)
  tuple val (sample_id), path(a)
  file(q)
  tuple val (sample_id), path(p)
  tuple val (sample_id), path(r)

  output:
  tuple val(sample_id), path('*.html'), emit: report_html
  tuple val(sample_id), path('filter_count_matrix/*/*.mtx.gz') 
  tuple val(sample_id), path('filter_count_matrix/*/*.tsv.gz')

  shell:
  '''
  mkdir -p filter_count_matrix/!{sample_id}/
  
  web_summary.R  --matrix !{f} --barcodes !{w} --features !{g} --sample !{sample_id} --multiqc_json !{m} --antisense !{a} --qualimap_report !{sample_id}_qualimap.txt --plot !{p} --cell_caller !{r}
  '''
}
