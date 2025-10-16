# AWS EC2 Testing Guide for CS Genetics scRNA-Seq Pipeline

This guide provides step-by-step instructions for testing the pipeline on a fresh AWS EC2 instance, including testing both Docker and Conda execution profiles.

## Table of Contents
- [EC2 Instance Specifications](#ec2-instance-specifications)
- [Step 1: Launch EC2 Instance](#step-1-launch-ec2-instance)
- [Step 2: Connect to Your Instance](#step-2-connect-to-your-instance)
- [Step 3: Install System Dependencies](#step-3-install-system-dependencies)
- [Step 4: Install Docker](#step-4-install-docker)
- [Step 5: Install Conda](#step-5-install-conda)
- [Step 6: Install Nextflow](#step-6-install-nextflow)
- [Step 7: Clone the Pipeline](#step-7-clone-the-pipeline)
- [Step 8: Configure AWS Credentials](#step-8-configure-aws-credentials)
- [Step 9: Test with Docker Profile](#step-9-test-with-docker-profile)
- [Step 10: Test with Conda Profile](#step-10-test-with-conda-profile)
- [Troubleshooting](#troubleshooting)

## EC2 Instance Specifications

### Recommended Instance Type
- **Instance Type**: `m5.4xlarge` (recommended)
- **vCPUs**: 16 (minimum for STAR alignment)
- **Memory**: 64 GB (provides comfortable headroom for STAR)
- **Storage**: 100 GB (EBS gp3 volume recommended)

### Why These Specs?
- The STAR alignment process requires 16 CPUs and 40 GB RAM (see [base.config](conf/base.config:55-58))
- The `m5.4xlarge` provides 64 GB RAM, giving ample headroom for STAR alignment
- The test profile processes 2 samples with 1M raw reads each
- Docker images and conda environments require significant disk space (~30-40 GB combined)
- Additional space needed for STAR index, GTF files, and results

### Alternative Instance Types
- `c5.4xlarge`: 16 vCPUs, 32 GB RAM (may be insufficient for STAR - not recommended)
- `c5.9xlarge`: 36 vCPUs, 72 GB RAM (more expensive, good for very large datasets)
- `m5.2xlarge`: 8 vCPUs, 32 GB RAM (cheaper but may be slow for STAR alignment)

## Step 1: Launch EC2 Instance

1. Log into your [AWS Console](https://console.aws.amazon.com/)
2. Navigate to **EC2** → **Instances** → **Launch Instance**

### Configure Instance Settings:

**Name and tags:**
- Name: `csgenetics-pipeline-test`

**Application and OS Images (Amazon Machine Image):**
- **AMI**: Ubuntu Server 24.04 LTS (HVM), EBS General Purpose (SSD) Volume Type
- **Architecture**: 64-bit (x86)

**Instance type:**
- Select: `m5.4xlarge` (recommended - see specifications above)

**Key pair (login):**
- Create new key pair or select existing
- **Type**: RSA
- **Format**: `.pem` (for SSH from Mac/Linux) or `.ppk` (for PuTTY on Windows)
- Download and save the key file securely

**Network settings:**
- **VPC**: Select your VPC (or use default)
- **Auto-assign public IP**: Enable
- **Firewall (security groups)**:
  - Create new security group or select existing
  - **Required rule**: Allow SSH (port 22) from your IP address
  - Click "Add security group rule" if needed:
    - Type: SSH
    - Source: My IP

**Configure storage:**
- **Size**: 100 GB
- **Volume type**: gp3 (recommended) or gp2
- **Delete on termination**: Checked (optional, for cleanup)

**Advanced details:**
- Leave defaults unless you have specific requirements
- **IAM instance profile**: If accessing private S3 buckets, attach an IAM role with S3 read permissions

3. Click **Launch Instance**
4. Wait for instance to reach "Running" state

## Step 2: Connect to Your Instance

### Option A: Using AWS Console (Session Manager)
1. Select your instance in the EC2 console
2. Click **Connect** → **Session Manager** → **Connect**

### Option B: Using SSH (Recommended)
1. Note your instance's **Public IPv4 address** from the EC2 console
2. Set permissions on your key file:
   ```bash
   chmod 400 /path/to/your-key.pem
   ```
3. Connect via SSH:
   ```bash
   ssh -i /path/to/your-key.pem ubuntu@<PUBLIC_IP_ADDRESS>
   ```

## Step 3: Create a Non-Root User with Sudo Access

**Note**: If you're using the Ubuntu AMI, you're already logged in as the `ubuntu` user which has sudo access. You can skip this step and proceed to Step 4.

If you're logged in as root or want to create a dedicated user for pipeline work:

```bash
# Create a new user (replace 'pipelineuser' with your desired username)
sudo adduser pipelineuser

# Add the user to the sudo group
sudo usermod -aG sudo pipelineuser

# Switch to the new user
su - pipelineuser

# Verify sudo access
sudo whoami
# Should output: root
```

**For SSH access with the new user**: You'll need to copy the SSH authorized keys:

```bash
# As root or ubuntu user, copy SSH keys to new user
sudo mkdir -p /home/pipelineuser/.ssh
sudo cp ~/.ssh/authorized_keys /home/pipelineuser/.ssh/
sudo chown -R pipelineuser:pipelineuser /home/pipelineuser/.ssh
sudo chmod 700 /home/pipelineuser/.ssh
sudo chmod 600 /home/pipelineuser/.ssh/authorized_keys

# Now you can SSH directly as the new user:
# ssh -i /path/to/your-key.pem pipelineuser@<PUBLIC_IP_ADDRESS>
```

All subsequent installation steps should be performed as a non-root user with sudo privileges.

## Step 4: Install Git

Git is required to clone the pipeline repository.

```bash
# Update package lists
sudo apt-get update

# Install git
sudo apt-get install -y git

# Verify installation
git --version
```

## Step 5: Install Docker

### Install Docker Engine

Follow the official Docker installation instructions for Ubuntu:

**[Docker Engine Installation for Ubuntu](https://docs.docker.com/engine/install/ubuntu/)**

After installation, verify Docker is working:

```bash
# Verify Docker installation
sudo docker run hello-world
```

### Configure Docker for Non-Root User

Follow the official Docker post-installation steps to run Docker as a non-root user:

**[Docker Post-Installation Steps for Linux](https://docs.docker.com/engine/install/linux-postinstall/)**

At minimum, you need to add your user to the docker group and reload the group membership:

```bash
# Add your user to the docker group
sudo usermod -aG docker $USER

# Apply group changes (or logout/login, or use newgrp)
newgrp docker

# Verify Docker works without sudo
docker run hello-world
```

## Step 6: Install Java (Required for Nextflow)

Nextflow requires Java 11 or newer. **We strongly recommend using SDKMAN to install Java.**

**Important**: During testing, we found that installing Java 25 directly from Oracle (following the Java installation docs) caused compatibility issues with Nextflow, even though the Nextflow documentation states it supports Java 25. SDKMAN installation of Java 17 LTS avoids these issues and provides the most reliable Nextflow compatibility.

```bash
# Install prerequisites
sudo apt-get install -y unzip zip curl

# Install SDKMAN
curl -s https://get.sdkman.io | bash

# Initialize SDKMAN in your current shell
source "$HOME/.sdkman/bin/sdkman-init.sh"

# Install Java 17 (LTS version - tested and confirmed working with Nextflow)
sdk install java 17.0.10-tem

# Verify Java installation
java --version
```

Expected output should show OpenJDK 17.0.10 or newer.

**Why SDKMAN?**
- Avoids version compatibility issues (Java 25 from Oracle caused problems during testing)
- Easy to manage multiple Java versions if needed
- Well-tested with Nextflow
- Recommended over `apt-get install default-jre` or direct Oracle JDK downloads

## Step 7: Install Nextflow

Follow the official Nextflow installation instructions:

**[Nextflow Installation Documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation)**

Quick installation steps:

```bash
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Make it executable
chmod +x nextflow

# Create a local bin directory if it doesn't exist
mkdir -p $HOME/.local/bin/

# Move nextflow to the local bin directory
mv nextflow $HOME/.local/bin/

# Add local bin to PATH in .bashrc
echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.bashrc

# Reload shell configuration
source ~/.bashrc

# Verify installation
which nextflow
nextflow -version
```

Expected output should show Nextflow version 24.04.3 or newer.

## Step 8: Install Conda

Follow the official Miniconda installation instructions for Linux:

**[Miniconda Installation for Linux](https://docs.anaconda.com/miniconda/install/#linux-2)**

After installation, you must accept the Anaconda Terms of Service for the default channels:

```bash
# Accept Terms of Service for default channels (required for pipeline)
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Verify conda installation
conda --version
```

**Important**: The pipeline requires Conda 25.9.1 or newer. If you have an older version, update it:

```bash
# Update conda to latest version
conda update -n base -c defaults conda

# Verify conda version is 25.9.1 or newer
conda --version
```

## Step 9: Clone the Pipeline

```bash
# Clone the pipeline repository
git clone https://github.com/csgenetics/csgenetics_scrnaseq.git
cd csgenetics_scrnaseq

# By default, you'll be on the main branch
# For development/testing, you can checkout other branches if needed:
# git checkout devel

# Verify current branch
git branch
```

**Note**: Most users should work with the `main` branch, which contains stable releases.

## Step 10: Test with Conda Profile

The conda profile is tested first because it will reveal any conda-specific setup issues.

### Run the Test Profile with Conda

```bash
# Ensure you're in the pipeline directory
cd ~/csgenetics_scrnaseq

# Run the test_conda profile
nextflow run main.nf -profile test_conda
```

### Expected Behavior on First Run

1. **Conda Environment Creation**:
   - Nextflow creates conda environments for each process (this takes time on first run)
   - Environments are cached in `work/conda/` directory
   - Subsequent runs will be much faster due to caching
   - You'll see messages like: "Creating Conda env: ..."

2. **Processing Phase**:
   - Downloads STAR index, GTF, test data from S3 (uses `--no-sign-request`, no AWS credentials needed)
   - Runs QC, STAR alignment, feature counting, cell calling
   - Generates reports

3. **Completion**:
   - Results in `./results_test` directory
   - Success message: "Completed at: [timestamp]"

### Verify Test Success

```bash
# Check for successful completion
tail -n 20 .nextflow.log

# Check results directory
ls -lh results_test/

# Expected directories:
# - report/          (HTML QC report)
# - count_matrix/    (gene expression matrices)
# - multiqc/         (MultiQC report)
```

**Expected Results**: The test profile processes 2 samples with 1M raw reads each, resulting in approximately 300-500 cells per sample with ~60 genes detected per cell.

## Step 11: Test with Docker Profile

### Clean Previous Run

```bash
# Remove previous work directory to ensure clean test
rm -rf work/
rm -rf results_test/
rm -f .nextflow.log*
```

### Run the Test Profile with Docker

```bash
# Ensure you're in the pipeline directory
cd ~/csgenetics_scrnaseq

# Run the test profile (uses Docker containers)
nextflow run main.nf -profile test
```

### Expected Behavior

1. **Download Phase**:
   - Downloads STAR index, GTF, test data from S3
   - Downloads Docker images automatically

2. **Processing Phase**:
   - QC processes run
   - STAR alignment (most time-consuming step)
   - Feature counting, cell calling, report generation

3. **Completion**:
   - Results in `./results_test` directory
   - Success message: "Completed at: [timestamp]"

### Verify Test Success

```bash
# Check for successful completion
tail -n 20 .nextflow.log

# Check results directory
ls -lh results_test/

# Expected directories:
# - report/          (HTML QC report)
# - count_matrix/    (gene expression matrices)
# - multiqc/         (MultiQC report)
```

## Step 12: Test with Full Dataset (Optional)

To test the pipeline with a more realistic, production-scale dataset, use the `test_pbmc_4_sample_full` profile. This processes 4 PBMC samples with 54-66M reads per sample.

### Clean Previous Run

```bash
# Remove previous test results
rm -rf work/
rm -rf results_test/
rm -f .nextflow.log*
```

### Run Full Dataset Test

```bash
# Ensure you're in the pipeline directory
cd ~/csgenetics_scrnaseq

# Run the full PBMC test profile with Docker
nextflow run main.nf -profile test_pbmc_4_sample_full
```

### Expected Behavior

1. **Download Phase**:
   - Downloads full STAR index
   - Downloads 4 large FASTQ files (54-66M reads each)
   - Total download size: Several GB

2. **Processing Phase**:
   - Processes 4 samples in parallel
   - Each sample results in 16-24k raw reads per cell
   - STAR alignment is the most time-consuming step
   - This will take significantly longer than the basic test profile

3. **Completion**:
   - Results in `./results_test` directory
   - Much larger output files than basic test

### Verify Full Test Success

```bash
# Check for successful completion
tail -n 20 .nextflow.log

# Check results directory
ls -lh results_test/

# Expected: Similar structure but with 4 samples worth of data
```

**Note**: This test provides a realistic assessment of pipeline performance with production-scale data.

