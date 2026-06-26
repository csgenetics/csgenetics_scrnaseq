# External Pipeline Map — csgenetics_scrnaseq

Source: `/nssd2/humebc/external/csgenetics_scrnaseq`. Single-cell RNA-seq -> genes-by-barcode count table. DSL2, single anonymous `workflow {}` in `main.nf` (no named sub-workflows). All process bodies live in `modules/processes.nf`. Manifest `nextflowVersion = '!>=25.10.2'`, strict syntax enabled, `process.resourceLimits` used (no check_max).

## 1. Workflow DAG (`main.nf`)

### Setup / param resolution (main.nf:18-42, 75-185)
- `getGenomeAttribute()` (main.nf:18) reads `params.genomes[params.genome][attr]`, returns null if absent.
- Genome attrs hoisted to params at parse time (main.nf:35-42): `star_index, gtf, mitochondria_chromosome, mixed_species, hsap_mitochondria_chromosome, mmus_mitochondria_chromosome, hsap_gene_prefix, mmus_gene_prefix`. Overridable on CLI.
- `save_resolved_configuration()` invoked unconditionally first (main.nf:94) — dumps resolved `params` JSON.

### S3 public-bucket conditional download pattern (main.nf:98-176)
For each public resource, `if (path.startsWith("s3://csgx.public.readonly"))` -> run a `download_*` process (uses `aws s3 cp --no-sign-request`); else build a `file()`/`channel.fromPath` directly. Applies to:
- `params.star_index` -> `download_star_index` (main.nf:98)
- `params.gtf` -> `download_gtf` (main.nf:106)
- `params.input_csv` -> `download_input_csv` (main.nf:114)
- `params.barcode_list_path` -> `download_barcode_list` (main.nf:164)
- `params.barcode_correction_list_path` -> `download_barcode_correction_list` (main.nf:172)

### FASTQ input handling (main.nf:125-160)
- `input_csv.splitCsv(header,sep=',',strip)` -> `.map` accesses **by index** not key (rowList[0]=sample, [1]=fastq_1, [2]=fastq_2) for nf-core `sample` rename compatibility (main.nf:127-138).
- Per-row `download` flag set if BOTH fastqs start with `s3://csgx.public.readonly`. `.branch{download/no_download}` (main.nf:139-152).
- `download` rows -> `download_public_fastq`; `no_download` rows -> `file()` objects.
- `download_public_fastq.out.downloaded_fastqs.mix(no_download_ch).groupTuple(by:0)` -> `ch_input` of `[sample, [R1...], [R2...]]` (main.nf:157-160).

### Lane merging (main.nf:194-224)
- `ch_input.branch{ multiple_lanes: R1.size()>1; single_lane: R1.size()==1 }` (main.nf:194).
- multiple_lanes `flatMap` into `[sample,"R1",[fastqs]]` and `[sample,"R2",[fastqs]]` -> `merge_lanes` (parallelizes R1/R2; avoids reliance on R1/R2 in filenames).
- `merge_lanes.out` `groupTuple(by:0,size:2)` then `.map` reorders to ensure R1 first -> `ch_merge_lanes_out_merged`.
- single_lane flattened to `[sample,R1,R2]`.
- `io_extract_in_ch = ch_merge_lanes_out_merged.mix(single_lane_flattened)`.

### QC + alignment (main.nf:233-264)
- `qc(io_extract_in_ch, barcode_correction_list)` — unified Rust binary (barcode extract/correct, SSS/polyX trim, internal polyA trim, Q30, MultiQC JSON). Outputs `qc_out`, `qc_log`, `qc_multiqc`.
- `qc.out.qc_out.branch{ good_fastq: countFastq()>0; empty_fastq: countFastq()==0 }` (main.nf:240).
- `star(good_fastq, star_index)` -> emits `[sample, bam, env(uniquely_mapped_reads)]`.
- `star.out.out_bam.branch{ good_bam: uniq>0; bad_bam: uniq==0 }` (main.nf:253).
- `create_valid_empty_bam_star` (aliased) fed by `empty_fastq` + `bad_bam`, makes a 1-line-header BAM so samtools downstream does not break.

### RSeQC raw + feature counting (main.nf:266-289)
- `gtf2bed(gtf)` -> gene model bed.
- `raw_rseqc` (alias of `run_rseqc`) fed by good_bam (count flag 1) mixed with empty bams (count flag 0); also takes bed, empty_rseqc_template, literal `"raw"`.
- `initial_feature_count(good_bam[flag1] + empty[flag0], gtf)` — featureCounts `-t transcript -g gene_id --fracOverlap 0.5 -s 1 -M`.
- `initial_feature_count_good_bam_out_ch` = inner-join of good_bam samples with feature_count_bam via `.combine(by:[0])` (uses combine not join because strict mode `failOnMismatch:true`) (main.nf:289).

### UMR + multimapper annotation fork (main.nf:291-320)
Two parallel branches off `initial_feature_count_good_bam_out_ch`:
- **UMR**: `filter_for_UMRs_mismatch` (NH==1, nM<=3) -> `umr_transcript_assignment` (transcript-feature Assigned) and `umr_exon_assignment` (exon tie-break on Unassigned_Ambiguity; takes gtf).
- **Multimapper**: `filter_for_multimappers_mismatch` (NH>1, nM<=3) -> `multimapper_transcript_assignment` (gawk script `bin/assign_multi_mappers.gawk`; emits assigned + unassigned) -> `multimapper_exon_assignment` (gtf + gawk; exon tie-break).
- Merge per branch: `merge_transcript_exon_umr_bams` (combine by:0), `merge_transcript_exon_multimapper_bams`.
- `merge_annotated_UMRs_with_annotated_multimappers` (combine by:0) -> `${sample}.mapped.sorted.filtered.annotated.bam` (published).
- `umr_multimapper_annotated_bam_out_ch` = that output **mixed with** `create_valid_empty_bam_star.out` (re-introduces 0-alignment samples) (main.nf:317).

### Counts of annotated reads + annotated RSeQC (main.nf:320-325)
- `count_high_conf_annotated_umr_multimap` -> `[sample, env(alignment_count)]`.
- `annotated_rseqc` (alias) fed by annotated bam combined with alignment_count, bed, template, literal `"annotated"`.

### MultiQC (main.nf:329-344)
- `ch_ss_multiqc_in` = `qc_multiqc` mixed with raw + annotated rseqc logs, `groupTuple(by:0,size:3)`, flattened to `[sample,[files]]`.
- `single_sample_multiqc(ch_ss_multiqc_in)`.
- `ch_ms_multiqc_in` = collect of all sample file-lists (sample id stripped) -> `multi_sample_multiqc`.

### Dedup + count matrix (main.nf:347-359)
- `sort_index_bam(umr_multimapper_annotated_bam_out_ch)` -> antisense count + sorted/indexed bam.
- `dedup` = sort_index output combined with alignment_count -> `umi_tools dedup --per-cell` (skips when count==0).
- `io_count(dedup.out.io_dedup_sam)` -> `*_bcGeneSummary.txt` via samtools+awk.
- `count_matrix(io_count_out, barcode_list, feature_file_out)` -> raw h5ad + tripartite mtx. `features_file(gtf)` (main.nf:184) produces the features tsv.

### Cell calling + filtering (main.nf:362-429)
- Re-split `input_csv` to extract user manual thresholds (main.nf:362-405): mixed_species expects hsap/mmus in cols 4/5 -> `"hsap_mmus"` string (else `"nan_nan"`); single-species col 4 -> value or `"nan"`.
- `ch_cell_caller = ch_h5ad.combine(thresholds, by:0)` -> `cell_caller` -> `[sample, stdout threshold]` + plots.
- `ch_filter_count_matrix_in`: `cell_caller_out.mix(ch_h5ad).groupTuple(by:0,size:2, sort by order_integer_first)` so int threshold comes first, then path (main.nf:64-73, 418-420). `order_integer_first` distinguishes int vs path (replaces a class-based sort that failed on Seqera Platform).
- `filter_count_matrix` -> filtered + raw h5ad (`raw_count_matrix` emit).

### Read categorization + reports (main.nf:427-471)
- `ch_categorize_reads_in` = good_bam STAR bam + filtered raw h5ad + qc json (combine by 0 twice).
- `categorize_reads` -> `${sample}.read_categorization.csv`.
- `summary_statistics`: big combine of cell_caller threshold, raw h5ad, antisense, dedup log, single_sample_multiqc json, raw rseqc, annotated rseqc, read categorization, qc log (main.nf:438-450) -> `${sample}.metrics.csv`.
- `qc_cascade_plot_single(metrics_csv)`.
- `single_summary_report(metrics_csv + cell_caller plots + qc cascade, template)` -> per-sample html + emits metrics csv.
- `qc_cascade_plot_multi(collect of all single_sample metric csvs)`.
- `multi_sample_report(collect metric csvs, template, multi qc cascade)` -> multisample html/csv/plots.

## 2. Process table

Container resolved from `conf/images.config`. All script sections are **`script:`** (none use `shell:` or `template:`). Language = bash unless noted. Several `script:` blocks build Groovy vars before the heredoc.

| Process | Key inputs | Key outputs | Container | Script lang / notes |
|---|---|---|---|---|
| save_resolved_configuration | (none) | resolved_configuration.txt | ubuntu:bionic | bash; Groovy `JsonOutput` builds string before echo |
| download_star_index | params.star_index | star/ dir (emit star_index) | quay.io/csgenetics/aws_download:0.1 | bash; `aws s3 cp --no-sign-request` |
| download_gtf | params.gtf | *.gtf | quay.io/csgenetics/aws_download:0.1 | bash |
| download_input_csv | params.input_csv | *.csv | quay.io/csgenetics/aws_download:0.1 | bash |
| download_barcode_list | params.barcode_list_path | *.csv | quay.io/csgenetics/aws_download:0.1 | bash |
| download_barcode_correction_list | params.barcode_correction_list_path | *.tsv | quay.io/csgenetics/aws_download:0.1 | bash |
| download_public_fastq | tuple(sample, fastq_1 val, fastq_2 val) | downloaded_fastqs tuple | quay.io/csgenetics/aws_download:0.1 | bash; Groovy `.tokenize('/').last()` |
| features_file | gtf | `${gtf.baseName}_features_names.tsv` (emit modified_gtf) | quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0 | bash calling `features_names.py` |
| merge_lanes | tuple(sample, read_num, fastqs) | `${sample}.merged.${read_num}.fastq.gz` | ubuntu:bionic | bash; `cat *.f*q.gz` (glob-ordered) |
| qc | tuple(sample,r1,r2); corrected_barcodelist | qc_out R1 fastq; qc_log; qc_multiqc (3 fastp jsons) | quay.io/csgenetics/qc:0.1 | bash calling Rust `qc` binary |
| star | tuple(sample,r1); index | bam + env(uniquely_mapped_reads) | quay.io/csgenetics/samtools_star:0.0.2 | bash; STAR `--runThreadN 8`, `--outFilterMultimapNmax 1000`, unsorted BAM. publishDir STAR |
| create_valid_empty_bam (alias create_valid_empty_bam_star) | tuple(sample,prefix) | `${sample}${prefix}.bam` | quay.io/biocontainers/samtools:1.17--hd87286a_1 | bash; samtools 1-line header. publishDir STAR |
| gtf2bed | gtf | gene_model.bed | quay.io/csgenetics/gtf2bed:2.0 | bash `gtf2bed` |
| run_rseqc (aliases raw_rseqc / annotated_rseqc) | tuple(sample,bam,count); bed; empty_template; val(prefix) | *_rseqc_results.txt (emit rseqc_log) | quay.io/biocontainers/rseqc:5.0.3--py39hf95cd2a_0 | bash; if count>0 read_distribution.py else envsubst template. publishDir RSeQC |
| initial_feature_count | tuple(sample,bam,aligned_count); gtf | featureCounts.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools sort + featureCounts. publishDir featureCounts |
| filter_for_UMRs_mismatch | tuple(sample,featurecount_bam) | UMRs.bam.featureCounts.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools `-e '[NH]==1 && nM<=3'` |
| umr_transcript_assignment | umr_mismatch_bam | UMRs.transcript.assigned.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools filter XS==Assigned |
| umr_exon_assignment | umr_mismatch_bam; gtf | UMRs.exon.assigned.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools + featureCounts exon tie-break |
| filter_for_multimappers_mismatch | feature_count_bam | multimapped.bam.featureCounts.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools `-e '[NH]>1 && nM<=3'` |
| multimapper_transcript_assignment | multimap bam; gawk script | assigned_bam + unassigned_bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; samtools + `gawk -f assign_multi_mappers.gawk` |
| multimapper_exon_assignment | unassigned bam; gtf; gawk script | multimapped.exon.assigned.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash; featureCounts exon + gawk |
| merge_transcript_exon_umr_bams | umr transcript+exon bams | umr.annotated.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash samtools merge |
| merge_transcript_exon_multimapper_bams | mm transcript+exon bams | multimapped.annotated.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash samtools merge |
| merge_annotated_UMRs_with_annotated_multimappers | umr + mm annotated bams | mapped.sorted.filtered.annotated.bam | quay.io/csgenetics/featurecount_samtools:0.2 | bash samtools merge. publishDir featureCounts |
| count_high_conf_annotated_umr_multimap | annotated bam | env(alignment_count) | quay.io/csgenetics/featurecount_samtools:0.2 | bash samtools view -c |
| single_sample_multiqc | tuple(sample,multiqc files) | html, *_data, multiqc_json | quay.io/csgenetics/multiqc:0.1 | bash multiqc `-m unified_qc -m rseqc`. publishDir multiqc |
| multi_sample_multiqc | collected multiqc files | multisample html/data/json | quay.io/csgenetics/multiqc:0.1 | bash multiqc. publishDir multiqc |
| sort_index_bam | tuple(sample,bam) | antisense.txt; sorted bam+bai | quay.io/biocontainers/samtools:1.17--hd87286a_1 | bash samtools |
| dedup | sorted bam+bai+alignment_count | dedup.log; dedup.bam | quay.io/csgenetics/umi-tools-csgx:0.9 | bash `umi_tools dedup --per-cell` (skip if 0). publishDir deduplication |
| io_count | dedup bam | *_bcGeneSummary.txt | quay.io/biocontainers/samtools:1.17--hd87286a_1 | bash samtools+awk |
| count_matrix | bcGeneSummary; barcode_list; features_file | raw h5ad; tripartite mtx | quay.io/csgenetics/scanpy_anndata:0.0.4 | bash `count_matrix.py`; Groovy builds mixed/single `--mixed_species` args. publishDir count_matrix/raw |
| cell_caller | h5ad; manual_threshold_str | stdout threshold; plots html | quay.io/csgenetics/scanpy_anndata:0.0.4 | bash `cell_caller.py` (uses params.minimum_count_threshold, !mixed). publishDir plots |
| filter_count_matrix | tuple(sample,count_threshold,raw h5ad) | filtered h5ad; raw h5ad (emit raw_count_matrix); mtx | quay.io/csgenetics/scanpy_anndata:0.0.4 | bash `filter_count_matrix.py`. publishDir count_matrix filtered+raw |
| categorize_reads | star bam; raw h5ad; fastp json | read_categorization.csv | quay.io/csgenetics/scanpy_anndata:0.0.4 | bash `categorize_reads.py`; barcode_length = count of 'C' in params.barcode_pattern |
| summary_statistics | sample, threshold, h5ad, antisense, dedup, multiqc json, raw+annotated rseqc, read cat csv, qc log | metrics.csv | quay.io/csgenetics/scanpy_anndata:0.0.4 | bash `summary_statistics.py`. publishDir report |
| qc_cascade_plot_single | metrics_csv | qc_cascade.html | quay.io/csgenetics/html_build:0.1.0 | bash `qc_cascade_plot.py --mode single`. publishDir report |
| qc_cascade_plot_multi | collected metric csvs | multisample_qc_cascade.html | quay.io/csgenetics/html_build:0.1.0 | bash `qc_cascade_plot.py --mode multi`. publishDir report |
| single_summary_report | metrics csv, plots, qc cascade, template | report.html; metric csv (emit) | quay.io/csgenetics/html_build:0.1.0 | bash `create_single_sample_report.py`. publishDir report |
| multi_sample_report | collected csvs, template, qc cascade | multisample html/csv/plots | quay.io/csgenetics/html_build:0.1.0 | bash `create_multi_sample_report.py`. publishDir report |

Note: every container is referenced only inside the `docker`/`singularity`/`test*` profiles (images.config is included only there). `local`, `conda`, `aws` profiles do NOT include images.config.

## 3. Genome / reference handling

`conf/genomes.config` defines `params.genomes` as a pure data map (no if statements), selected by `params.genome`. Each entry supplies the 8 attrs hoisted in main.nf. All references hosted on `s3://csgx.public.readonly` (STAR_2_7_11b indices + GTFs), triggering the no-sign-request download path.

| genome key | mixed | mito chr | gene prefixes |
|---|---|---|---|
| GRCh38 | false | MT | - |
| GRCm39 | false | MT | - |
| BDGP6 (fly) | false | mitochondrion_genome | - |
| Sscrofa11 (pig) | false | MT | - |
| mouse_human_mix | **true** | (uses hsap=GRCh38_MT, mmus=GRCm39_MT) | GRCh38_ / GRCm39_ |
| GRCz11 (zebrafish) | false | MT | - |
| GRCh38_test (test.config only) | false | MT | - (reduced index/gtf) |

`mixed_species` drives mixed-arg branches in `count_matrix`, `cell_caller`, `filter_count_matrix`, `categorize_reads`, `summary_statistics`, reports, the input_csv threshold parsing (main.nf:367), and STAR memory (60GB vs 40GB, base.config:53).

### Key default params (nextflow.config:22-50)
- `input_csv = null`, `outdir = './results'`.
- `barcode_list_path`, `barcode_correction_list_path` -> default IDT_IO_kit_v2 on public S3.
- `barcode_pattern = "CCCCCCCCCCCCC"` (13bp).
- `minimum_count_threshold = 100`, `sss_nmer = 8`.
- `awsregion = null`, `awsqueue = null`.
- `max_memory = 256.GB`, `max_cpus = 16`, `max_time = 240.h` -> feed `process.resourceLimits` (nextflow.config:176).
- Global: `nextflow.enable.strict = true`, `process.shell = ['/bin/bash','-euo','pipefail']`, env `PYTHONNOUSERSITE/R_PROFILE_USER/R_ENVIRON_USER`. timeline/report/trace/dag all enabled under pipeline_info.

## 4. Profiles (nextflow.config:53-156)

Each profile re-`includeConfig 'conf/base.config'` inside its own block (comment: required for resource allocation to work). base.config sets `errorStrategy='retry'`, `maxRetries=5`, plus per-process cpus/memory (`memory = { N.GB * task.attempt }`).

| Profile | genomes | test cfg | images | engine |
|---|---|---|---|---|
| local | yes | - | no | executor=local |
| docker | yes | - | yes | docker (runOptions `-u $(id -u):$(id -g)`) |
| singularity | yes | - | yes | singularity |
| conda | yes | - | conda_envs.config | conda (createTimeout 60min) |
| test | - | test.config | yes | docker |
| test_pbmc_4_sample_full | - | test_pbmc_4_sample_full.config | yes | docker |
| test_hsap_mmus_2_sample_full | - | test_hsap_mmus_2_sample_full.config | yes | docker |
| test_singularity | yes | test.config | yes | singularity |
| test_conda | yes | test.config | conda_envs.config | conda |
| aws | yes | - | no (uses aws.config) | aws.config: executor=awsbatch, queue=params.awsqueue, region, batch.cliPath |

`conf/test.config`: 2-sample 1M-read minimal dataset; genome `GRCh38_test` with reduced index/gtf; input_csv + barcode_list on public S3.
`conf/aws.config`: `aws.region`, `aws.batch.cliPath='/home/ec2-user/miniconda/bin/aws'`, `process.executor='awsbatch'`, `process.queue=params.awsqueue`.

## 5. Notable for overhaul (HOW, not WHAT)
- `qc` already a unified Rust binary replacing multiple legacy fastp/io_extract steps (main.nf:226-235 comments).
- `create_valid_empty_bam` and the count-flag (1/0) mixing pattern exist solely to keep 0-read/0-alignment samples flowing without breaking samtools — heavy branch/mix bookkeeping.
- `.combine(by:0)` used instead of `.join` deliberately because strict mode defaults `failOnMismatch:true` (main.nf:288 comment).
- `order_integer_first` sort workaround for Seqera Platform (class-based sort failed) (main.nf:64-73).
- STAR runs `--runThreadN 8` hardcoded while base.config gives `cpus=16` — mismatch worth flagging.
- featureCounts thread flags hardcoded `-T 4` across feature processes vs base.config cpus.
