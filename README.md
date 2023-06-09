# CS Genetics scRNA-Seq pipeline
**A Nextflow pipeline for processing scRNA-Seq data generated using CS Genetics' single-cell kit to produce a genes by barcode count table.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Contents
- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Locally with dockerized processes](#locally-with-dockerized-processes) 
  - [Dockerized](#dockerized)
- [Running the pipeline on MacOS](#running-the-pipeline-on-macos)
- [Specifying input sequencing files](#specifying-input-sequencing-files)
- [Testing the pipeline](#testing-the-pipeline)
- [Launching the pipeline directly from the Github repo](#launching-the-pipeline-directly-from-the-csgeneticcsgenetics_scrnaseq-github-repo)
- [Updating the pipeline](#updating-the-pipeline)
- [Specifying a pipeline version](#specifying-a-pipeline-version)
- [Available profiles](#available-profiles)
  - [test](#test)
  - [docker](#docker)
- [Configurable parameters](#configurable-parameters)
  - [profile](#profile)
  - [outdir](#outdir)
  - [star_index_dir](#star_index_dir)
  - [gtf_path](#gtf_path)
  - [whitelist_path](#whitelist_path)
  - [min_nuc_gene](#min_nuc_gene)
  - [depth_min](#depth_min)
- [Outputs](#outputs)
  - [count_matrix](#count_matrix)
  - [report](#report)
  - [pipeline_info](#pipeline_info)
  - [fastp](#fastp)
  - [fastqc](#fastqc)
  - [featureCounts](#fastqc)
  - [io_count](#io_count)
  - [multiqc](#multiqc)
  - [plots](#plots)
  - [qualimap](#qualimap)
  - [STAR](#star)
- [Log files](#log-files)
- [Examples](#examples)
  - [Example 1](#example-1)
  - [Example 2](#example-2)

## Introduction

**CS Genetics' scRNA-Seq pipeline** is a bioinformatics best-practice analysis pipeline for processing single-cell RNA-Seq data from their single-cell RNA-Seq kits.

It runs on a Unix-like operating system (E.g. Linux).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable manner. The pipeline uses Docker containers to run the main pipeline instance and its constituent processes making installation trivial and results reproducible.

It processes FASTQ files generated by CS Genetics' scRNA-Seq library kit to count matrices that can be loaded directly into [Seurat](https://satijalab.org/seurat/index.html) or [scanpy](https://scanpy.readthedocs.io/en/stable/) for further analyses.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Running the pipeline

The pipeline is run natively in your environment but using docker containers for each of the pipeline's processes.

To run the pipeline locally, you must have Nextflow ([installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)) and git ([installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)) installed. It is recommended that you work with either the Docker of Singularity profiles
(see [Available profiles](#available-profiles) below) to enable the use of preconfigured containers for each of the pipeline processes.

For each profile you must have the respective program installed: [Docker installation](https://docs.docker.com/get-docker/), [Singularity installation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

With Nextflow, Docker and git installed, clone the csgenetics_scrnaseq repository to a specified directory and change into that directory.

```bash
git clone https://github.com/csgenetics/csgenetics_scrnaseq.git $HOME/analysis && cd $HOME/analysis
```

The pipeline is then run using the Nextflow executable and in this case the `docker` profile.

For specification of input sequences see [Specifying input sequencing files](#specifying-input-sequencing-files).

```bash
nextflow run main.nf -profile docker --input_csv $HOME/analysis/input_csv/input_csv.csv
```

For a full list of the configurable parameters that can be can be supplied to the pipeline
and other options for configuration [Configurable parameters](#configurable-parameters).
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Running the pipeline on MacOS
While the pipeline can be launched on MacOS, some of the processes are RAM intensive.

E.g. the STAR mapping process is currently configured to run in a container that is allocated 40GB of RAM.
The qualimap process is configured to use 16GB of RAM.

The actual resources utilized will depend on the charachter of the samples being analysed.

However, resource allocations may exceed those available on a standard MacOS laptop/desktop.

As such it is recommened to run the pipeline on an HPC system.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Specifying input sequencing files
The sequencing files to be analysed are specified using an input csv file.

An example template can be found [here](input_csv/template.csv).

The header row must be present and the full paths to files should be given.

The `fastq_1` should contain the sequencing data that will be mapped to the genome. `fastq_2` should contain the CS Genetics barcode.

The full path to the input csv should be supplied to the pipeline using the `--input_csv` flag.

E.g.
```bash
nextflow run main.nf -profile docker --input_csv <path/to/csv_dir/input_csv.csv>
```

Multiple sets of sequencing files (e.g. from multiple lanes of sequencing)
can be merged by the pipeline and used for a single sample
by supplying the same sample name but with different sequencing file sets
on separate lines.

E.g.
```bash
sample_id,fastq_1,fastq_2
Sample1,/home/example_user/analysis/raw_reads/example_Sample1_L001_R1_001.fastq.gz,/home/example_user/analysis/raw_reads/example_Sample1_L001_R1_001.fastq.gz
Sample1,/home/example_user/analysis/raw_reads/example_Sample1_L002_R1_001.fastq.gz,/home/example_user/analysis/raw_reads/example_Sample1_L002_R1_001.fastq.gz
Sample2,/home/example_user/analysis/raw_reads/example_Sample2_L001_R1_001.fastq.gz,/home/example_user/analysis/raw_reads/example_Sample2_L001_R1_001.fastq.gz
```
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Testing the pipeline
To test that your environment is setup correctly the pipeline can be run using the test profile:

```bash
nextflow run main.nf -profile test
```

The test proile will run using a set of remotely hosted resources. By default the work and results directories will be created in the current working directory at `./work` and `./results`, respectively.

For a full list of the configurable parameters that can be can be supplied to the pipeline
and other options for configuration see [Configurable parameters](#configurable-parameters).
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Launching the pipeline directly from the csgenetic/csgenetics_scrnaseq Github repo
In the above examples, the csgenetics_scrnaseq git repostory was cloned locally
and the pipeline was launched specifying the main.nf Nextflow script.

Alternatively the pipeline can be launched directly from the GitHub repository specifying its qualified name: `csgenetics/csgenetics_scrnaseq`.

See the [Nextflow documentation on Pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html) for further details.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

### Updating the pipeline
When launching the pipeline from a local clone of the GitHub repository, you will need to keep the pipeline up-to-date by refularly pulling down updates:
```
git pull
```
When launching the pipeline by specifying the qualified name of the pipeline, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull csgenetics/csgenetics_scrnaseq
```
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

### Specifying a pipeline version

It's a good idea to specify a pipeline version (alternatively referred to as a tag or a release) when launching the pipeline. This ensures that a specific version of the pipeline code and software are used when you run your pipeline, thus ensuring reproducibilty.

See the [Nextflow documentation on Pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html) for further details.

The available tags can be displayed in a cloned repository using:
`git tag`

Alternatively, the releases can be viewed [online](https://github.com/csgenetics/csgenetics_scrnaseq/releases).

For reference, the version will be logged in reports when you run the pipeline.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Available profiles
Nextflow pipeline configurable parameters can be set in groups by specifying profiles.

See the [Config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) section of the Netflow documentation for further details.

There are four profiles available for the CS Genetics scRNA-Seq pipeline:
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
  - Runs Docker containers for each process using Docker
- `test_singularity`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
  - Runs Docker containers for each process using Singularity
- `docker`
  - A generic configuration profile that enables use of pre-configured Docker containers for each process run using Docker.
- `singularity`
  - A generic configuration profile that enables use of pre-configured Docker containers for each process run using Singularity.
  - See the Nextflow documentation on [Singularity](https://www.nextflow.io/docs/latest/container.html#singularity) for further details.

If no profile is set, then the [local executor](https://www.nextflow.io/docs/latest/executor.html#local) will be used and software required for each process to run should be pre-installed on your local system.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

### test
Launching the pipeline with this profile will set configuration parameters so that remotely hosted small, human fastq files and a remotely hosted GRCh38 set of resources are used. The pipeline will use docker containers for each of the processes.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results_test    # Finished results
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
E.g.
```bash
nextflow run main.nf -profile test
```

When running the test profile, do not supply the `--input-csv` argument. A remotely hosted input csv is used.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

To run the `test` profile using Singualrity to launch the Docker containers use the profile `test_singularity`

E.g.
```bash
nextflow run main.nf -profile singularity_test
```
### docker
Launching the pipeline with this profile configures the pipeline to use pre-specified Docker containers for each of the processes. It is recommended to run the pipeline using this profile.

E.g.
```bash
nextflow run main.nf -profile docker --input_csv $HOME/analysis/input_csv/input_csv.csv
```
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Configurable parameters

Please refer to the [Nextflow documentation on configuration](https://www.nextflow.io/docs/latest/config.html)
for a general introduction to configuring Nextflow pipelines.

### `profile`

Use this parameter to choose a configuration profile. See [Available profiles](#available-profiles).

```bash
-profile docker
```
N.B. note the single hyphen.

### `outdir`

The output directory where the results will be saved.

```bash
--outdir <path/to/output/dir>
```

### `star_index_dir`
Specify the path of the STAR index directory. Required for mapping.

By default the remotely hosted Human STAR index is used see [below](#premade-star-indexes).

```bash
--star_index_dir s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/star/
```

<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

#### Premade STAR indexes
There are premade remotely hosted STAR indexes for the following species (remotely hosted path given):
- Human: s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/star/

If you are working with one of these species, you can provide the remotely hosted directory
to the `--star_index_dir` parameter. The pipeline will automatically download the resource.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

#### Generating a STAR index
If you are working with a different species or wish to create your own indexes for a different genome,
please follow the instructions for creating a STAR index [here](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) (section 'Generating genome indexes').
 
### `gtf_path` 

Path to the gtf annotation file.

There are remotely hosted GTF files for the following species (remotely hosted path given):
- Human: s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf

By default the remotely hosted Human gtf is used.

```bash
--gtf_path s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf
```
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

### `whitelist_path`
Specify the whitelist path to use.
The following whitelists are hosted remotely at: 

- IDT_IO_kit_v2.csv: s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v2.csv

By default the IDT_IO_kit_v2.csv whitelist is used.
```bash
--whitelist_path s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v2.csv
```

### `min_nuc_gene`

The minimum number of nuclear genes that must be detected for a given barcode to be considered a cell.

```bash
--min_nuc_gene 100
```

### `depth_min`

Fastqs with less than this number of raw reads will be removed from the analysis.

```bash
--depth_min 100000
```

<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Outputs
The output directory is specified using the `--outdir` flag. In this directory the pipeline outputs a number of useful results files organised within the following subdirectories:

### `count_matrix`
Contains the count matrices. The matrices are output in two different formats:
- .h5ad (compatible with scanpy in Python)
- the tripartite barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz (compatible with Seurat in R)

The matrices are output as 'raw' (containing all cellular barcodes) and 'cell_only' (containing a subset of the barcodes that were classified as cells through meeting the minimum nuclear genes detected threshold).

### `report`
Contains the per sample .html reports detailing key statistics for each of the samples.
Metrics are also provided in .csv format.

Contains the experiment report that consolidates metrics for all samples.

### `pipeline_info`
Contains trace files related to the execution of the pipeline

### `fastp`
Contains the html and json files output from the fastp QC tasks per sample.

### `fastqc`
Contains the html and .zip files output from the fastqc QC tasks per sample.

### `featureCounts`
Contains the files output from the featureCounts process including .bam files and summaries of feature assignment.

### `io_count`
Contains the files associated with deduplication and grouping of reads.

### `multiqc`
Contains files related to the MultiQC output.

### `plots`
Plots of the Cell Caller profiles used to generate the minimum detected nuclear genes threshold
for cell calling.

### `qualimap`
Contains the qualimap output logs used for assessing mapping metrics.

### `STAR`
Contains the bam files output from STAR and associated mapping log files.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Log files
As standard, Nextflow produces a `.nextflow.log` file in the directory from which the pipeline was run.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

## Examples
Below are some examples of launching the pipeline with explanations of the commands.

### Example 1
```bash
nextflow run main.nf -profile docker --input_csv <path/to/input/csv> --outdir ./results --star_index_dir </local/path/to/STAR_index/dir> --gtf_path </local/path/to/gtf> --whitelist s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v2.csv
```

A pipeline is launched using a locally installed version of
Nextflow from a locally cloned copy of the csgenetics/csgenetics_scrnaseq repository.

A local path to the STAR index directory and a
corresponding gtf are provided. The remotely hosted v2 whitelist is selected explicitly.

The output directory is set to `./results`.

The docker profile is selected, configuring the pipeline to use pre-specified Docker containers for each of the processes.

The 'work' directory will by created by default at `./work`.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>

### Example 2
```bash
nextflow run csgenetics/csgenetics_scrnaseq -r 0.0.1 -profile docker --input_csv <path/to/input/csv>
```

A pipeline is launched by directly specifying the GitHub repository.
Version 0.0.1 of the pipeline is used.

The docker profile is selected, configuring the pipeline to use pre-specified Docker containers for each of the processes.

The remotely hosted Human STAR index and gtf file are used by default. The remotely hosted v2 whitelist is used by default.
<div style="text-align: right"><a href="#cs-genetics-scrna-seq-pipeline">top</a></div>
