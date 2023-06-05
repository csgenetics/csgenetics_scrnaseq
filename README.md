# CS Genetics scRNA-Seq pipeline
**A Nextflow pipeline for processing scRNA-Seq data generated using CS Genetics' single-cell kit to produce a genes by barcode count table.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40csgenetics-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/csgenetics)

## Contents
- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
    - [Dockerized](#dockerized)
    - [Locally with dockerized processes](#locally-with-dockerized-processes)
- [Testing the pipeline](#testing-the-pipeline)


## Introduction

**CS Genetics' scRNA-Seq pipeline** is a bioinformatics best-practice analysis pipeline for processing single-cell RNA-Seq data.

It is designed to run on a Unix-like operating system (Linux, macOS, etc)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable manner. The pipeline uses Docker containers to run the main pipeline instance and its constituent processes making installation trivial and results highly reproducible.

It processes FASTQ files generated by CS Genetics' scRNA-Seq library kit to count matrices that can be loaded directly into [Seurat](https://satijalab.org/seurat/index.html) or [scanpy](https://scanpy.readthedocs.io/en/stable/) for further analyses.

## Running the pipeline

The pipeline can be run in one of two ways:
- **Dockerized (recommended)**
- **Locally with dockerized processes**

### Dockerized

The pipeline is run in a Docker container that has Nextflow installed. Within the container, the pipeline uses dedicated Docker containers for each of the executed processes.

To run the pipeline dockerized, you must have Docker ([installation instructions](https://docs.docker.com/get-docker/)) and git ([installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)) installed.

With Docker and git installed, clone the csgenetics_scrnaseq repository to a specified directory and change into that directory.

```bash
git clone https://github.com/csgenetics/csgenetics_scrnaseq.git $HOME/analysis && cd $HOME/analysis
```

The pipeline is then run by starting a docker image and supplying
it with commandline arguments. including the docker profile. E.g.

```bash
docker run -it --rm -e USER=$USER \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $HOME/analysis:$HOME/analysis \
  public.ecr.aws/csgenetics/csgenetics_scrnaseq_base:0.0.1 \ 
  nextflow run main.nf -profile docker -w $HOME/analysis/work/ --outdir $HOME/analysis/results/ --input-csv $HOME/analysis/input_csv/input_csv.csv
```

where:

- `-it --rm ` runs the container interactively and cleans up after exiting.

- `-e USER=$USER` makes your USER environmental variable visible instide the container.

- `-v /var/run/docker.sock:/var/run/docker.sock` mounts `/var/run/docker.sock` at `/var/run/docker.sock` enabling Docker to spawn sister containers.

- `-v $HOME/analysis:$HOME/analysis` mounts `$HOME/analysis` at `$HOME/analysis` setting up a persistent working directory for the pipeline run.

- `nextflow run main.nf ...` is the command for launching the pipeline including commandline-suopplied arguments (...)


For a full list of the configurable parameters that can be can be supplied to the pipeline
and other options for configuration see the [usage docs](docs/usage.md).


### Locally with dockerized processes
The pipeline is run natively in your environment using docker containers for each of the pipeline's processes.

To run the pipeline locally with dockerized processes, you must have Nextflow ([installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)), Docker ([installation instructions](https://docs.docker.com/get-docker/)) and git ([installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)) installed.

With Nextflow, Docker and git installed, clone the csgenetics_scrnaseq repository to a specified directory and change into that directory.

```bash
git clone https://github.com/csgenetics/csgenetics_scrnaseq.git $HOME/analysis && cd $HOME/analysis
```

The pipeline is then run using the Nextflow executable and the 'docker' profile.

```bash
nextflow run main.nf -profile docker --input-csv $HOME/analysis/input_csv/input_csv.csv
```

For a full list of the configurable parameters that can be can be supplied to the pipeline
and other options for configuration see the [usage docs](docs/usage.md).

## Testing the pipeline
To test that your environment is setup correctly the pipeline can be run using the test profiles like so

```bash
nextflow run main.nf -profile test
```

The test proile will run using a set of remotely hosted resources. By default the work and results directories will be created in the current working directory at `./work` and `./results`, respectively.

For a full list of the configurable parameters that can be can be supplied to the pipeline
and other options for configuration see the [usage docs](docs/usage.md).

For more details on running 

## Launching the pipeline directly from the csgenetic/csgenetics_scrnaseq Github repo
In the above examples, the csgenetics_scrnaseq git repostory was cloned locally
and the pipeline was launched specifying the main.nf Nextflow script.

Alternatively the pipeline can be launched directly from the GitHub repository specifying its qualified name: `csgenetics/csgenetics_scrnaseq`.

See the [Nextflow documentation on Pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html) for further details.








### Prebuilt genome refresources
After intial QC the seuqencing files will be mapped against your specified genome using [STAR](https://github.com/alexdobin/STAR).

For mapping to be conducted, a STAR genome index is required.

The path to the directory containing the STAR index must be supplied to the pipeline using the species_path parameter either by modifying the config files or as a command line-supplied argument.

For more details see the [usage docs](docs/usage.md).



or through supplying the 

For mouse and human, pre-built 

### a. Generating Genome Indexes

STAR source code and binaries can be downloaded from GitHub: named releases from https://github.com/alexdobin/STAR/releases

Strongly recommended that users generate their own genome indexes with most up-to-date assemblies and annotations

Examples of acceptable genome sequence files:
- **ENSEMBL**: Download the FASTA file containing all the chromosomes together in the [Human Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index) genome, which has primary assembly in the filename. Navigate to the Gene annotation section of the Ensembl website and click on the Download [Human GTF link](https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/). This takes you to an FTP site with a list of GTF files available.

- **GENCODE**: Files marked with PRI (primary) should be used as a refrences files. Strongly recommended for mouse and human:http://www.gencodegenes.org/

The basic options to generate genome indices using STAR are as follows:

```bash 
STAR \
    --runMode genomeGenerate \
    --genomeDir /path/to/genomeDir \
    --genomeFastaFiles /path/to/genome/fasta.fa \
    --runThreadN 4 \
    --sjdbGTFfile /path/to/annotations.gtf
```

## Pipeline-Environments

Different options are available to run the pipeline using various configurations in this repo:

### Nextflow with docker images 

The following are required to run the CSgenetics scRNAseq pipeline;
1. Prerequisite software
    - [`Nextflow`](https://nf-co.re/usage/installation) (preferably v21.x or higher)
    - [`Docker`](https://docs.docker.com)

2. Reference genome
    - Already build index files are stored for Human (GRCh38) and Mouse (GRCm39) genome in [S3 bucket](https://s3.console.aws.amazon.com/s3/buckets/bptest0?prefix=csgenetics/&region=us-east-1).
    - Please refer to Generating Genome Indexes section for custom genome downloading and building STAR index.

      **NOTE**:- Please make sure to change the species_path option in the nextflow.config based on genome used.

3. Clone the csgenetics repository by running the following command:

```bash
git clone https://github.com/basepair/csgenetics.git && cd csgenetics
```

4. Test it on a minimal dataset with a single command

```bash
nextflow run main.nf -profile test
```

5. Start running your own analysis!
    - See [usage docs](docs/usage.md) for all of the available options when running the pipeline.
   
```bash
nextflow run main.nf -profile <docker>
```


### Nextflow Local/Manual Installation

If the software is not already available, you will need to install it

#### Dependencies

##### Core tools

* Unix-like operating system (Linux, macOS, etc)
* [Nextflow](https://www.nextflow.io/)  >=21.10.3
* [JAVA](https://www.oracle.com/java/technologies/downloads/) >=11
* [Basespace-cli](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview) 0.1 

##### Bio tools

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.9
* [Fastp](https://github.com/OpenGene/fastp/releases/tag/v0.23.1) 0.23.1
* [Umi-tools-csgx](https://pypi.org/project/umi-tools-csgx/) 1.1.1
* [Samtools](http://www.htslib.org/download/) 1.14
* [Star](https://github.com/alexdobin/STAR/releases/tag/2.7.9a) 2.7.9a
* [Subread](https://subread.sourceforge.net/) 2.0.1
* [Multiqc](https://multiqc.info) 1.11
* [Qualimap](https://bioconda.github.io/recipes/qualimap/README.html) 2.2.2d

##### Python packages

* [Python](https://www.python.org/downloads/) 3.9.7
* [Pysam](https://pypi.org/project/pysam/0.19.1/) 0.19.1
* [Scanpy](https://pypi.org/project/scanpy/1.7.2/) 1.7.2

##### R packages

* [R](https://www.r-project.org/) 4.2.1
* optparse 1.7.3              
* R2HTML 2.3.3               
* rjson 0.2.21                
* scales 1.2.1               
* ggplot2 3.4.1               
* DropletUtils 1.18.1        
* SingleCellExperiment 1.20.0 
* SummarizedExperiment 1.28.0
* Biobase 2.58.0             
* GenomicRanges 1.50.2       
* GenomeInfoDb 1.34.9         
* IRanges 2.32.0             
* S4Vectors 0.36.1            
* BiocGenerics 0.44.0        
* MatrixGenerics 1.10.0       
* matrixStats 0.63.0         
* Matrix 1.5-3     
  


Start running your own analysis!

- See [usage docs](docs/usage.md) for all of the available options when running the pipeline.
   
```bash
nextflow run main.nf -profile <standard>
```

