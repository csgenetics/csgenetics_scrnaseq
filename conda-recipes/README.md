# CS Genetics Conda Recipes

This directory contains the Conda recipes whose source is maintained in this public
repository:

- **csgenetics-gtf2bed** - a Rust GTF-to-BED converter
- **multiqc-csgenetics** - a platform-independent MultiQC plugin for pipeline reports

The pipeline also consumes the published `csgenetics-qc` package. Its source and package
recipe are maintained outside this repository, so there is deliberately no local
`conda-recipes/qc` recipe. The customer-facing version pins remain in
[`conda_envs/qc.yml`](../conda_envs/qc.yml) and
[`conda_envs/multiqc.yml`](../conda_envs/multiqc.yml); those environment files are the
authoritative dependencies used by the pipeline's Conda profile.

The pipeline uses the standard `umi_tools` package from Bioconda, not a custom build.

## Prerequisites

Install [Pixi](https://pixi.sh/) and create an anaconda.org account with upload access to
the `cs_genetics` organization. The checked-in `pixi.toml` supplies `conda-build` and the
Anaconda client.

## Build the GTF-to-BED converter

The recipe packages a Linux MUSL binary built from the tracked Rust crate at
`images/gtf2bed`. With the MUSL target installed, run these commands from the repository
root:

```bash
cd images/gtf2bed
cargo build --locked --release --target x86_64-unknown-linux-musl
cd ../../conda-recipes
pixi run conda-build gtf2bed
```

## Build the MultiQC plugin

From the repository root:

```bash
cd conda-recipes
pixi run conda-build multiqc-csgenetics
```

The recipe builds `multiqc-csgenetics` from the tracked source at `images/multiqc`.
Conda-build prints the exact output package path when either recipe completes.

## Upload a package

Authenticate using the Anaconda client without writing credentials into this public
repository, then upload the package path reported by conda-build:

```bash
cd conda-recipes
pixi run anaconda upload <path-reported-by-conda-build> --user cs_genetics
```

Verify the published packages at
[anaconda.org/cs_genetics/csgenetics-gtf2bed](https://anaconda.org/cs_genetics/csgenetics-gtf2bed)
and
[anaconda.org/cs_genetics/multiqc-csgenetics](https://anaconda.org/cs_genetics/multiqc-csgenetics).

## Update a recipe

1. Update the version in the source package and its matching `meta.yaml` together.
2. Reset the recipe build number to `0` for a new version, or increment it when rebuilding
   unchanged source.
3. Build and test the package with the applicable command above.
4. Upload the exact artifact produced by conda-build.
