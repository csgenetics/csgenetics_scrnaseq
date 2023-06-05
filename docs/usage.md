# Configuration of the CS Genetics scRNA-Seq pipeline

## Table of contents

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [`-profile`](#-profile)
  - [`--fastq`](#--fastq)
  - [`--fastq_path`](#--fastq_path)
  - [`--access-token`](#--access-token)
  - [`--bs_project_id`](#--bs_project_id)
  - [`--species_path`](#--species_path)
  - [`--genome_path`](#--genome_path)
  - [`--gtf_path`](#--gtf_path)
  - [`--io_whitelist`](#--io_whitelist)
  - [`--whitelist_path`](#--whitelist_path)
  - [`--min_nuc_gene`](#--min_nuc_gene)
  - [`--purity`](#--purity)
  - [`--depth_min`](#--depth_min)
  - [`--remove_singletons`](#--remove_singletons)
  - [`--dedup`](#--dedup)
  
- [Other command line parameters](#other-command-line-parameters)
  - [`--outdir`](#--outdir)
  - [`-resume`](#-resume)
  - [`-c`](#-c)

      <!-- TOC END -->

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
nextflow run main.nf -profile <optional profile> --input-csv $HOME/analysis/input_csv/input_csv.csv
```

or directly from the GitHub repository:

```bash
nextflow run csgenetics/csgenetics_scrnaseq -profile <optional profile> --input-csv $HOME/analysis/input_csv/input_csv.csv
```

### Updating the pipeline
When launching the pipeline from a local clone of the GitHub repository, you will need to keep the pipeline up-to-date by refularly pulling down updates:
```
git pull
```
When launching the pipeline by specifying the qualified name of the pipeline, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull basepair/csgenetics
```

### Specifying a pipeline version

It's a good idea to specify a pipeline version (alternatively referred to as a tag or a release) when launching the pipeline. This ensures that a specific version of the pipeline code and software are used when you run your pipeline, thus ensuring reproducibilty.

See the [Nextflow documentation on Pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html) for further details.

The available tags can be displayed in a cloned repository using:
`git tag`

Alternatively, the releases can be viewed [online](https://github.com/csgenetics/csgenetics_scrnaseq/releases).

For reference, the version will be logged in reports when you run the pipeline.


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

#### docker
Launching the pipeline with this profile includes configures the pipeline to use pre-specified Docker containers for each of the processes. It is recommended to run the pipeline using this profile.

E.g.
```bash
nextflow run main.nf -profile docker --input-csv $HOME/analysis/input_csv/input_csv.csv
```

## Configurable parameters

### `-profile`
N.B. note the single hyphen

Use this parameter to choose a configuration profile. See [Available profiles](#available-profiles).

### `--outdir`

The output directory where the results will be saved.

### `--species_path`
<!-- TODO remove -->
Use this to specify the location of your species refrences and index file 
Already built index are stored in s3 bucket for human (GRCh38) and mouse (GMCh38) genome . 
For Building custom Genome Indexes please refer to Generating Genome Indexes section 

```bash
--species_path `s3://bp-publc/reflib/csgenetics/GRCh38.ensembl.release_103/`
```

### `--genome_path`
<!-- TODO rename -->
Specify the path of star index file 

```bash
--genome_path `s3://bp-publc/reflib/csgenetics/GRCh38.ensembl.release_103/star_index`
```
 
### `--gtf_path` 

Path to the gtf annotation file

```bash
--gtf_path `s3://bp-publc/reflib/csgenetics/GRCh38.ensembl.release_103/`
```

### `--io_whitelist`
<!-- TODO remove -->
To Run UMI-tools and generate an inferred whitelist based on IOs detected in reads
As we provide our own whitelist, this step is currently not run in the pipeline by default

```bash
--io_whitelist false
```

### `--whitelist_path`
<!-- TODO host whitelist remotely. -->
Specify the whitelist path to use. The whitelist to use will be specific to the kit version used.
The whitelists are hosted remotely at: 

Defaults to bin/IDT_IO_kit_v2.csv.

```bash
--whitelist_path `/bin/IDT_IO_kit_v1.csv`
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

### `--dedup`
<!-- TODO remove -->
Deduplicate reads using UMI-tools

```bash
--dedup true
```
