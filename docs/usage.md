# Configuration of the CS Genetics scRNA-Seq pipeline

## Table of contents

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Specifying a pipeline version](#specifying-a-pipeline-version)
  - [Available profiles](#available-profiles)
    - [test](#test)
    - [docker](#docker)
- [Configurable parameters](#configurable-parameters)
  - [`-profile`](#-profile)
  - [`--outdir`](#--outdir)
  - [`--star_index_dir`](#star-index-dir)
  - [`--gtf_path`](#gtf_path)
  - [`--whitelist_path`](#whitelist_path)
  - [`--min_nuc_gene`](#min_nuc_gene)
  - [`--depth_min`](#depth_min)
  - [`--remove_singletons`](#remove_singletons)

## Introduction
This document details:
- How to launch the pipeline
- How to update the pipeline
- How to run specific versions of the pipeline
- Available profiles
- The configurable parameters

This document details the configurable parameters of the pipeline.

Please refer to the [Nextflow documentation on configuration](https://www.nextflow.io/docs/latest/config.html)
for a general introduction to configuring Nextflow pipelines.

## Running the pipeline

If launching the pipeline from a local clone of the GitHub repository, the typical command for running the pipeline is:
```bash
nextflow run main.nf -profile <optional profile> --input_csv $HOME/analysis/input_csv/input_csv.csv
```

or directly from the GitHub repository:

```bash
nextflow run csgenetics/csgenetics_scrnaseq -profile <optional profile> --input_csv $HOME/analysis/input_csv/input_csv.csv
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

### Updating the pipeline
When launching the pipeline from a local clone of the GitHub repository, you will need to keep the pipeline up-to-date by refularly pulling down updates:
```
git pull
```
When launching the pipeline by specifying the qualified name of the pipeline, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull basepair/csgenetics
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

### Specifying a pipeline version

It's a good idea to specify a pipeline version (alternatively referred to as a tag or a release) when launching the pipeline. This ensures that a specific version of the pipeline code and software are used when you run your pipeline, thus ensuring reproducibilty.

See the [Nextflow documentation on Pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html) for further details.

The available tags can be displayed in a cloned repository using:
`git tag`

Alternatively, the releases can be viewed [online](https://github.com/csgenetics/csgenetics_scrnaseq/releases).

For reference, the version will be logged in reports when you run the pipeline.
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

### Available profiles
Nextflow pipeline configurable parameters can be set in groups by specifying profiles.

See the [Config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) section of the Netflow documentation for further details.

There are two profiles available for the CS Genetics scRNA-Seq pipeline:
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile that enables use of pre-configured Docker containers for each process

If no profile is set, then the [local executor](https://www.nextflow.io/docs/latest/executor.html#local) will be used and software required for each process to run should be pre-installed on your local system.
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

#### test
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
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

#### docker
Launching the pipeline with this profile includes configures the pipeline to use pre-specified Docker containers for each of the processes. It is recommended to run the pipeline using this profile.

E.g.
```bash
nextflow run main.nf -profile docker --input_csv $HOME/analysis/input_csv/input_csv.csv
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

## Configurable parameters

### `-profile`
N.B. note the single hyphen

Use this parameter to choose a configuration profile. See [Available profiles](#available-profiles).

### `--outdir`

The output directory where the results will be saved.

### `--star_index_dir`
Specify the path of the STAR index directory. Required for mapping.

```bash
--star_index_dir s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/star/
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

#### Premade STAR indexes
There are premade remotely hosted STAR indexes for the following species (remotely hosted path given):
- Human: s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/star/

If you are working with one of these species, you can provide the remotely hosted directory
to the `--star_index_dir` parameter. The pipeline will automatically download the resource.
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

#### Generating a STAR index
If you are working with a different species or wish to create your own indexes for a different genome,
please follow the instructions for creating a STAR index [here](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) (section 'Generating genome indexes').
 
### `--gtf_path` 

Path to the gtf annotation file.

There are remotely hosted GTF files for the following species (remotely hosted path given):
- Human: s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf

```bash
--gtf_path s3://csgx.public.readonly/resources/references/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>

### `--whitelist_path`
Specify the whitelist path to use. The whitelist to use will be specific to the kit version used.
The whitelists are hosted remotely at: 

- IDT_IO_kit_v1.csv: s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v1.csv
- IDT_IO_kit_v2.csv: s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v1.csv

```bash
--whitelist_path s3://csgx.public.readonly/resources/whitelists/IDT_IO_kit_v1.csv
```

### `--min_nuc_gene`

The minimum number of nuclear genes that must be detected for a given barcode to be considered a cell.

```bash
--min_nuc_gene 100
```

### `--depth_min`

Fastqs with less than this number of raw reads will be removed from the analysis.

```bash
--depth_min 100000
```

### `--remove_singletons`

If set to `true` UMRs with only 1 count per sample will be filtered out of the analysis.

```bash
--remove_singletons false
```
<div style="text-align: right"><a href="#configuration-of-the-cs-genetics-scrna-seq-pipeline">top</a></div>
