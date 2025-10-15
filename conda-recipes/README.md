# CS Genetics Conda Recipes

This directory contains conda recipes for building CS Genetics custom packages for the scRNA-seq pipeline.

## Overview

Two custom packages are provided:

1. **csgenetics-qc** - Rust binary for barcode extraction and quality control
2. **multiqc-csgenetics** - Custom MultiQC plugins for pipeline reporting

Note: The pipeline uses standard `umi_tools` from bioconda, not a custom version.

## Prerequisites

```bash
# Install conda-build
conda install conda-build anaconda-client

# Create anaconda.org account and CS Genetics organization
# Visit: https://anaconda.org/
```

## Building Packages

### 1. QC Binary (Platform-Specific)

The QC binary is a compiled Rust executable and must be built separately for each platform.

**Prerequisites:**
- Rust toolchain (automatically installed by conda-build via the `rust` compiler package)
- Pixi environment with conda-build (see Prerequisites section above)

**Build for Linux (most common for HPC):**

```bash
cd conda-recipes
pixi run conda-build qc

# Output package location:
# .pixi/envs/default/conda-bld/linux-64/csgenetics-qc-0.3.4-0.conda
```

**Notes:**
- The recipe uses `cargo install` which handles binary compilation and installation
- Currently builds from local source at `../../../rnaseq/images/qc`
- The meta.yaml can be updated to build from GitHub once the repo has proper tags and LICENSE

**Build for Mac:**

Building on the target platform is recommended for best compatibility:

```bash
# On Mac Intel
cd conda-recipes
pixi run conda-build qc
# Output: .pixi/envs/default/conda-bld/osx-64/csgenetics-qc-0.3.4-0.conda

# On Mac ARM (M1/M2/M3)
cd conda-recipes
pixi run conda-build qc
# Output: .pixi/envs/default/conda-bld/osx-arm64/csgenetics-qc-0.3.4-0.conda
```

**Troubleshooting:**

If you encounter "cannot find qc binary" errors, the build.sh script uses `cargo install --path . --root $PREFIX` which properly handles binary installation even with `strip = true` in Cargo.toml.

The test phase verifies:
1. Binary exists and is executable
2. Running `qc` without arguments shows usage message

### 2. MultiQC Plugin (Platform-Independent)

This is a Python package (noarch) that provides a custom MultiQC module for CS Genetics QC output.

**Prerequisites:**
- Pixi environment with conda-build (see Prerequisites section above)

**Build:**

```bash
cd conda-recipes
pixi run conda-build multiqc-csgenetics

# Output package location:
# .pixi/envs/default/conda-bld/noarch/multiqc-csgenetics-0.1.0-py_0.conda
```

**Notes:**
- The recipe builds from local source at `../../images/multiqc`
- Contains the `unified_qc` module that parses CS Genetics QC JSON output
- Registers as a MultiQC plugin via Python entry points

**What this plugin does:**
- Provides custom MultiQC module for CS Genetics unified QC binary output
- Automatically discovered by MultiQC when installed
- Displays barcode extraction stats, trimming metrics, and quality distributions


## Uploading to Anaconda.org

```

### Upload Packages

**Using Pixi environment:**

```bash
cd conda-recipes

# Upload QC binary (do for each platform you built)
pixi run anaconda upload .pixi/envs/default/conda-bld/linux-64/csgenetics-qc-0.3.4-0.conda --user cs_genetics

# Upload MultiQC plugin
pixi run anaconda upload .pixi/envs/default/conda-bld/noarch/multiqc-csgenetics-0.1.0-0.conda --user cs_genetics
```

**Authentication:**
The upload command uses the `ANACONDA_API_TOKEN` environment variable for authentication. Make sure it's set:
```bash
export ANACONDA_API_TOKEN=your_token_here
# Or add to ~/.bashrc for persistence
```

**Verify uploads:**
- QC binary: https://anaconda.org/cs_genetics/csgenetics-qc
- MultiQC plugin: https://anaconda.org/cs_genetics/multiqc-csgenetics


## Using the Packages

### Add CS Genetics Channel

Users add the channel to their conda configuration:

```bash
# Add channel globally
conda config --add channels cs_genetics

# Or specify in environment file
```

### In Conda Environment Files

```yaml
name: my_environment
channels:
  - cs_genetics
  - conda-forge
  - bioconda
dependencies:
  - csgenetics-qc=0.3.4
  - multiqc-csgenetics=0.1.0
  - umi_tools  # from bioconda
```

### Direct Installation

```bash
# Install from cs_genetics channel
conda install -c cs_genetics csgenetics-qc
conda install -c cs_genetics multiqc-csgenetics

# Install umi_tools from bioconda
conda install -c bioconda umi_tools
```


## Updating Packages

### Version Bumping

1. Update version in source code
2. Update version in `meta.yaml`
3. Commit and tag (for git sources)
4. Rebuild package
5. Upload new version

### Example: Updating QC to v0.3.5

```bash
# 1. Update qc source and tag
cd /path/to/qc
# ... make changes ...
git tag v0.3.5
git push --tags

# 2. Update conda recipe
cd conda-recipes/qc
# Edit meta.yaml: change version to 0.3.5
# Edit meta.yaml: increment build number OR reset to 0 for new version

# 3. Rebuild
conda build .

# 4. Upload
anaconda upload ~/miniconda3/conda-bld/linux-64/csgenetics-qc-0.3.5-0.tar.bz2 --user csgenetics
```
