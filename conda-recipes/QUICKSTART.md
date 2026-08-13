# Quick Start: Building the CS Genetics MultiQC Conda Package

## Step 1: Install Pixi

Install [Pixi](https://pixi.sh/). The `pixi.toml` in this directory supplies
`conda-build` and the Anaconda client.

## Step 2: Build the MultiQC plugin

From the repository root:

```bash
cd conda-recipes
pixi run conda-build multiqc-csgenetics
```

Conda-build tests the package and prints the exact output artifact path.

## Step 3: Upload to Anaconda.org

```bash
cd conda-recipes
pixi run anaconda upload <path-reported-by-conda-build> --user cs_genetics
```

Verify the published package at
[anaconda.org/cs_genetics/multiqc-csgenetics](https://anaconda.org/cs_genetics/multiqc-csgenetics).

The pipeline also installs the separately maintained, published `csgenetics-qc` package.
Do not try to build it from this directory: its authoritative customer pin is
[`conda_envs/qc.yml`](../conda_envs/qc.yml), and this repository intentionally contains
no QC source recipe.

---
