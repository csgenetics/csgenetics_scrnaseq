# Quick Start: Building CS Genetics Conda Packages

## Step 1: Install conda-build

```bash
# Install conda-build and anaconda-client
conda install -n base conda-build anaconda-client

# Verify installation
conda build --version
anaconda --version
```


## Step 3: Build Your First Package - MultiQC Plugin

Start with the easiest package (Python, no compilation):

```bash
cd csgenetics_scrnaseq/conda-recipes

# Build it
pixi run conda-build multiqc-csgenetics

# Should see:
# BUILD START: multiqc-csgenetics-0.1.0-py_0
# ...
# TEST START: multiqc-csgenetics-0.1.0-py_0
# ...
# BUILD END: multiqc-csgenetics-0.1.0-py_0
```

**Output location:**
```
.pixi/envs/default/conda-bld/noarch/multiqc-csgenetics-0.1.0-py_0.conda
```

## Step 4: Upload to Anaconda.org

```bash
cd csgenetics_scrnaseq/conda-recipes

# Upload
pixi run anaconda upload .pixi/envs/default/conda-bld/noarch/multiqc-csgenetics-0.1.0-py_0.conda --user cs_genetics

# Verify on web
# Visit: https://anaconda.org/cs_genetics/multiqc-csgenetics
```


## Step 5: Build QC Binary

**Important:** QC binary is platform-specific! The recipe is already configured to build from local source.

```bash
cd csgenetics_scrnaseq/conda-recipes

# Build (Rust toolchain installed automatically by conda-build)
pixi run conda-build qc

# Upload
pixi run anaconda upload .pixi/envs/default/conda-bld/linux-64/csgenetics-qc-0.3.4-0.conda --user cs_genetics
```

**Notes:**
- The recipe uses `cargo install` which handles compilation and installation
- Build takes about 2 minutes (15 seconds for Rust compilation)
- Output: `.pixi/envs/default/conda-bld/linux-64/csgenetics-qc-0.3.4-0.conda`

### Building for Multiple Platforms

Build on each target platform for best compatibility:

```bash
# On Linux server
cd conda-recipes
pixi run conda-build qc
pixi run anaconda upload .pixi/envs/default/conda-bld/linux-64/csgenetics-qc-0.3.4-0.conda --user cs_genetics

# On Mac Intel
cd conda-recipes
pixi run conda-build qc
pixi run anaconda upload .pixi/envs/default/conda-bld/osx-64/csgenetics-qc-0.3.4-0.conda --user cs_genetics

# On Mac ARM (M1/M2/M3)
cd conda-recipes
pixi run conda-build qc
pixi run anaconda upload .pixi/envs/default/conda-bld/osx-arm64/csgenetics-qc-0.3.4-0.conda --user cs_genetics
```

---

## You're Done!

Your packages are now available to anyone:

```bash
# Users can install with:
conda install -c cs_genetics csgenetics-qc multiqc-csgenetics
conda install -c bioconda umi_tools

# Or in environment.yml:
channels:
  - cs_genetics
  - conda-forge
  - bioconda
dependencies:
  - csgenetics-qc=0.3.4
  - multiqc-csgenetics=0.1.0
  - umi_tools  # from bioconda
```

**Verify your packages:**
- https://anaconda.org/cs_genetics/csgenetics-qc
- https://anaconda.org/cs_genetics/multiqc-csgenetics
```

---
