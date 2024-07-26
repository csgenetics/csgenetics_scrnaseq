#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
process download_input_csv {
  tag "Download input csv"

  output:
  path("*.csv")

  script:
  """
  aws s3 cp --no-sign-request ${params.input_csv} .
  """
}

/* Download the csgx hosted templates and barcode list
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
  tag "$sample_id"

  input:
  tuple val(sample_id), path(fastq_1), path(fastq_2)

  output:
  tuple val(sample_id), path("${sample_id}_R1.merged.fastq.gz"), path("${sample_id}_R2.merged.fastq.gz"), emit: merge_lanes_out

  shell:
  '''
  #cat !{fastq_1} > !{sample_id}_R1.merged.fastq.gz
  #cat !{fastq_2} > !{sample_id}_R2.merged.fastq.gz


  # merging R1 and R2
  # globs are ordered so lane merging will happen in same order for R1 and R2
  # If there's only one file for R1 or R2, we just rename it
  R1_files=$(ls !{fastq_1} | wc -l)
  if [ "$R1_files" -gt 1 ]; then
    cat *R1*.f*q.gz > !{sample_id}_R1.merged.fastq.gz
  elif [ "$R1_files" -eq 1 ]; then
    mv *R1*.f*q.gz !{sample_id}_R1.merged.fastq.gz
  fi

  # For R2 files
  R2_files=$(ls !{fastq_2} | wc -l)
  if [ "$R2_files" -gt 1 ]; then
    cat *R2*.f*q.gz > !{sample_id}_R2.merged.fastq.gz
  elif [ "$R2_files" -eq 1 ]; then
    mv *R2*.f*q.gz !{sample_id}_R2.merged.fastq.gz
  fi
  '''
}

process merged_fastp{
  tag "$sample_id"

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
* Extract the 13bp barcode from the R2 read and append it to the header of the R1 read.
* Remove any reads where the barcode does not exactly match a barcode in the barcode list.
*/
process io_extract {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(r1), path(r2)
  path(barcode_list)
  path(io_extract_script)

  output:
  tuple val(sample_id), path("${sample_id}.io_extract.R1.fastq.gz"), emit: io_extract_out

  script:
  """
  cat $barcode_list | cut -d ',' -f2 > barcode_list.txt
  gawk -v r2=${r2} -v sample_id=${sample_id} -v bc_length=13 -f $io_extract_script barcode_list.txt <(zcat ${r1})

  # If it doesn't exist, create an empty file to collect
  if [ ! -f "${sample_id}.io_extract.R1.fastq.gz" ]; then
    touch ${sample_id}.io_extract.R1.fastq && gzip ${sample_id}.io_extract.R1.fastq
  fi
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

  publishDir "${params.outdir}/fastp", pattern: '*.{json,html}', mode: 'copy'

  input: tuple val(sample_id), path(r1)

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
    -l 5 \
    -A \
    -G \
    -j !{sample_id}_R1.io_extract.fastp.json \
    -h !{sample_id}_R1.io_extract.fastp.html \
    --stdout \
     2> fastp.log | gzip > !{sample_id}_R1.io_extract.fastp.fastq.gz
  '''     
}

/*
* Trims reads from positions that match the regexs: A{15,} or A{13,}CG 
* I.e. A homopolymers >= 15 bp or >=13bp + CG are identified and trimmed along with any following bps.
* Only reads >20bp in length (after trimming) are retained .
*/
process trim_extra_polya {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(fastq)
  path(trim_polyA_script)

  output:
  tuple val(sample_id), path("${sample_id}_R1.polyA_trimmed.fastq.gz"), emit: trim_extra_polya_out

  script:
  """
  zcat $fastq | awk -f $trim_polyA_script
  # Check if the output file exists and rename it to the sample_id
  # If it doesn't exist, create an empty file
  if [ -f "good.fastq.gz" ]; then
    mv good.fastq.gz ${sample_id}_R1.polyA_trimmed.fastq.gz
  else
    touch ${sample_id}_R1.polyA_trimmed.fastq && gzip ${sample_id}_R1.polyA_trimmed.fastq
  fi
  """
}

process post_polyA_fastp{
  tag "$sample_id"

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
* Generic process for running qualimap.
*/
process run_qualimap {
  tag "$sample_id"
  
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
      qualimap rnaseq -outdir !{sample_id}_!{prefix}_qualimap -a proportional -bam !{bam} -p strand-specific-forward -gtf !{gtf} --java-mem-size=!{task.memory.toGiga()}G
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
* Annotate the bam alignment with gene feature annotations
* and use these to filter the reads so that only reads with
* unambiguous annotations are carried through.
* We implement a strategy to reassign ambiguous initial annotations.
* We handle uniquely mapped reads (UMRs) and multimapped reads separately.
* UMRs
* 1 - Reads are first annotated using the gene feature. Unambigously annotated reads
* (tagged "Assigned") are extracted and retained to merge with other retained reads.
* 2 - Those reads that were classified as Unassigned_Ambiguity are passed into a second
* round of featureCount, this time annotating against the exon feature. 'Assigned' reads
* are retained to merge with other retained reads.
* Multimapped reads
* 1 - Reads that are asociated unambiguously to only 1 gene
* are retained as UMRs to the associated alignment.
* 2 - Reads that are ambiguously assigned are annotated with the exon feature
* and assigned reads from this annotation are retained.
* See the top of the gawk script for further details on how the multimappers are handled.
*/
process feature_counts {
  tag "$sample_id"

  publishDir "${params.outdir}/featureCounts", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), val(aligned_count)
  path(gtf)
  path(multi_mapper_resolution_script)

  output:
  tuple val(sample_id), path("${sample_id}.mapped.sorted.filtered.annotated.bam"), env(alignment_count), emit: out_bam

  shell:
  """
  if [[ $aligned_count > 0 ]] # If the bam is not empty
    then
      samtools sort $bam -o ${sample_id}_Aligned.sortedByCoord.out.bam
      # Start by running feature counts on the star output
      # including strandedness and annotation of multimappers
      featureCounts -a $gtf -o ${sample_id}.star.featureCounts.gene.txt -R BAM ${sample_id}_Aligned.sortedByCoord.out.bam -T 4 -t transcript -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M

      # Process the multimapped and uniquely mapped reads separately.
      # Generate the UMRs
      samtools view -h -b -e '[NH]==1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam > ${sample_id}.UMRs.bam.featureCounts.bam
      # Generate the multimapped reads
      samtools view -h -b -e '[NH]>1 && ([nM]==0 || [nM]==1 || [nM]==2 || [nM]==3)' -b ${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam > ${sample_id}.multimapped.bam.featureCounts.bam

      ##### UMR processing ####
      # Filter out the reads that were 'Assigned' a gene target
      samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${sample_id}.UMRs.bam.featureCounts.bam > ${sample_id}.UMRs.gene.assigned.bam

      # Filter out the reads that were that were classified as Unassigned_Ambiguity and run them through
      # featureCounts using the exon tag to see if the ambiguity can be cleared up based on exon mapping.
      samtools view -h -b -e '[XS]=="Unassigned_Ambiguity"' -b ${sample_id}.UMRs.bam.featureCounts.bam > ${sample_id}.UMRs.gene.unassigned_ambiguity.bam

      # We have to remove the XS tag from the *.gene.unassigned_ambiguity.bam because featureCounts adds
      # an additional tag, that prevents proper fitltering with samtools
      samtools view -h ${sample_id}.UMRs.gene.unassigned_ambiguity.bam | sed 's/\\tXS\\:Z\\:[^\\t]*//' | samtools view -h -b > ${sample_id}.UMRs.gene.unassigned_ambiguity.no_xs_tag.bam

      # Do exon tie-breaking and filter to those that are assigned.
      featureCounts -a $gtf -o ${sample_id}.UMRs.exon.assigned.txt -R BAM ${sample_id}.UMRs.gene.unassigned_ambiguity.no_xs_tag.bam -T 4 -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1
      samtools view -h -b -e '[XN]==1 && [XT] && [XS]=="Assigned"' -b ${sample_id}.UMRs.gene.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam > ${sample_id}.UMRs.exon.assigned.bam  

      ##### Multimapper processing ####
      # Sort the bam file by query name in preparation for running through the gawk program for the first time.
      # Run the multimapped, annotated alignments through the gawk program that pulls
      # out alignments that can be associated directly as 'Assigned' and alignments that are ambiguous and should be passed onto
      # exon tie breaking.
      # See the script's header for more information on how it works.
      samtools sort -n ${sample_id}.multimapped.bam.featureCounts.bam | samtools view | gawk -f $multi_mapper_resolution_script

      # Produces assigned_reads.sam_body and ambiguous_reads.sam_body corresponding to the Assigned and still ambigous reads, respectively.
      # These files will only be produced if there were reads of the respective type identified.
      # The assigned_reads.sam_body should be converted into a valid bam by adding the headers back in
      if [ -f assigned_reads.sam_body ]; then
        # Cat with the headers of the featureCounts bam
        cat <(samtools view -H ${sample_id}.multimapped.bam.featureCounts.bam) assigned_reads.sam_body | samtools view -b -h > ${sample_id}.multimapped.gene.assigned.bam;
        # Remove the file so that it doesn't get confused with the next round of the gawk script's outputs
        rm assigned_reads.sam_body
      fi

      # If the ambigous_reads.sam_body exists then we convert this back to a genuine bam, by adding the headers and then rerun it through featureCounts using the -t exon
      # assignment to do the exon tie break. NB. the featureCount-derived bam tags e.g. XS, XT and XN have been removed in the gawk script.
      if [ -f ambiguous_reads.sam_body ]; then
        # Cat with the headers of the featureCounts bam
        # before rerunning through featureCounts for exon tie-breaking
        cat <(samtools view -H ${sample_id}.multimapped.bam.featureCounts.bam) ambiguous_reads.sam_body | samtools view -h -b > ${sample_id}.multimapped.gene.unassigned_ambiguity.no_xs_tag.bam;
        # Do exon tie-breaking and then run back through the gawk script to pull out those
        # reads that have a single Assigned alignment.
        featureCounts -a $gtf -o ${sample_id}.multimapped.exon.assigned.txt -R BAM ${sample_id}.multimapped.gene.unassigned_ambiguity.no_xs_tag.bam -T 4 -t exon -g gene_id --fracOverlap 0.5 --extraAttributes gene_name -s 1 -M
        samtools sort -n ${sample_id}.multimapped.gene.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam | samtools view | gawk -f $multi_mapper_resolution_script
        # If the assigned_reads.sam_body file exists then we were successfuly able to pull out further assigned reads
        if [ -f assigned_reads.sam_body ]; then
          # Cat with the headers of the featureCounts bam
          cat <(samtools view -H ${sample_id}.multimapped.gene.unassigned_ambiguity.no_xs_tag.bam.featureCounts.bam) assigned_reads.sam_body | samtools view -b -h > ${sample_id}.multimapped.exon.assigned.bam;
        fi
        rm ambiguous_reads.sam_body
        # If the assigned_reads.sam_body doesn't exist then we weren't able to pull out any further assigned reads. 
      fi
      
      # Finally, merge together any of the four *.assigned.bam files that exist
      if [ -f ${sample_id}.UMRs.exon.assigned.bam ]; then
        samtools merge -o ${sample_id}.UMRs.annotated.bam ${sample_id}.UMRs.gene.assigned.bam ${sample_id}.UMRs.exon.assigned.bam;
      else
        mv ${sample_id}.UMRs.gene.assigned.bam ${sample_id}.UMRs.annotated.bam;
      fi

      if [ -f ${sample_id}.multimapped.gene.assigned.bam ]; then
        samtools merge -o ${sample_id}.multimapped.temp.annotated.bam ${sample_id}.UMRs.annotated.bam ${sample_id}.multimapped.gene.assigned.bam;
      else
        mv ${sample_id}.UMRs.annotated.bam ${sample_id}.multimapped.temp.annotated.bam;
      fi

      if [ -f ${sample_id}.multimapped.exon.assigned.bam ]; then
        samtools merge -o ${sample_id}.mapped.sorted.filtered.annotated.bam ${sample_id}.multimapped.temp.annotated.bam ${sample_id}.multimapped.exon.assigned.bam;
      else
        mv ${sample_id}.multimapped.temp.annotated.bam ${sample_id}.mapped.sorted.filtered.annotated.bam;
      fi

    else
      # Simply rename the input bam so that it can be collected
      cp $bam ${sample_id}.mapped.sorted.filtered.annotated.bam
  fi

  alignment_count=\$(samtools view -c ${sample_id}.mapped.sorted.filtered.annotated.bam)
  """
}

/*
* Run multiqc
*/
process multiqc {
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

  publishDir "${params.outdir}/count_matrix/raw_feature_bc_matrix/${sample_id}/", mode: 'copy', pattern: "*.raw_feature_bc_matrix.h5ad"
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
  def mixed_args = params.mixed_species ? "--mixed_species True --hsap_mito_chr ${params.hsap_mitochondria_chromosome} --mmus_mito_chr ${params.mmus_mitochondria_chromosome} --hsap_gene_prefix ${params.hsap_gene_prefix} --mmus_gene_prefix ${params.mmus_gene_prefix}" : "--mixed_species False --mito_symbol ${params.mitochondria_chromosome}"
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

  publishDir "${params.outdir}/plots", pattern: "${sample_id}_pdf_with_cutoff.png", mode: 'copy'

  input:
  tuple val(sample_id), path(count_matrix_h5ad)

  output:
  tuple val(sample_id), stdout, emit: cell_caller_out
  tuple val(sample_id), path("${sample_id}*_pdf_with_cutoff.png"), emit: cell_caller_plot

  script:
  """
  cell_caller.py --sample_name ${sample_id} --minimum_count_threshold ${params.minimum_count_threshold} --count_matrix ${count_matrix_h5ad} --single_species ${!params.mixed_species}
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
* Generate the summary statistics required for HTML report
*/
process summary_statistics {
  tag "$sample_id"

  publishDir "${params.outdir}/report/${sample_id}", mode: 'copy', pattern: "*.csv"
  
  input:
  tuple val(sample_id), val(minimum_count_threshold), path(raw_h5ad), path(annotated_qualimap), path(antisense), path(dedup), path(multiqc_data_json), path(raw_qualimap)
  output:
  tuple val(sample_id), path("${sample_id}.metrics.csv"), emit: metrics_csv

  script:
  def mixed_args = params.mixed_species ? "TRUE" : "FALSE"
  """
  summary_statistics.py ${sample_id} ${raw_h5ad} ${multiqc_data_json} ${antisense} ${dedup} ${raw_qualimap} ${annotated_qualimap} $mixed_args
  """
}

/*
* Generate a per sample html report 
*/
process single_summary_report {
  tag "$sample_id"

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

