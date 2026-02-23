# Installation Guide

## Prerequisites

### System Requirements
- Linux/Unix system (Ubuntu 20.04+, CentOS 8+, or macOS)
- 16GB+ RAM recommended
- 100GB+ disk space
- 8+ CPU cores recommended

### Required Software

#### 1. Conda (Miniconda or Anaconda)

**Miniconda (recommended):**
```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts - recommend "yes" to initialize conda
```

**macOS:**
```bash
# Using Homebrew
brew install miniconda
conda init zsh  # or bash
```

#### 2. Nextflow (optional, for workflow execution)

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Move to PATH
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -v
```

## Installation Steps

### Method 1: Using Conda Environment (Recommended)

```bash
# Clone or download this pipeline
cd /path/to/pipeline

# Create environment from file
conda env create -f environment.yml

# Activate environment
conda activate variant-tracing

# Verify installations
fastqc --version
bwa
samtools --version
gatk --version
vep --help
iqtree2 --version
```

### Method 2: Manual Installation

If you prefer to install tools individually:

```bash
# Create new environment
conda create -n variant-tracing python=3.11
conda activate variant-tracing

# Quality Control
conda install -c bioconda fastqc trimmomatic

# Alignment
conda install -c bioconda bwa samtools

# Variant Calling
conda install -c bioconda gatk4 bcftools

# Annotation
conda install -c bioconda ensembl-vep

# Population & Ancestry
conda install -c bioconda plink admixture haplogrep

# Phylogenetics
conda install -c bioconda iqtree
```

### Download Reference Genome

```bash
# Create data directory
mkdir -p data/reference

# Download human GRCh38 from NCBI
cd data/reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p12/GCF_000001405.39_GRCh38.p12_genomic.fna.gz

# Rename
gunzip GCF_000001405.39_GRCh38.p12_genomic.fna.gz
mv GCF_000001405.39_GRCh38.p12_genomic.fna GRCh38.fa

# Index with BWA (will be done automatically by pipeline)
```

### Download Test Data (Optional)

```bash
# Create reads directory
mkdir -p data/reads

# Download NA12878 sample from NIST (optional)
# See: https://www.nist.gov/programs-projects/genomic-reference-materials
```

## Verification

After installation, verify all tools are available:

```bash
# Activate environment
conda activate variant-tracing

# Test each tool
fastqc --version
# Expected: FastQC v0.12.1

bwa 2>&1 | head -3
# Expected: Program: bwa (alignment via Burrows-Wheeler transform)

samtools --version
# Expected: samtools 1.21

gatk --version
# Expected: The Genome Analysis Toolkit (GATK) v4.5.0.0

vep --help | head -5
# Expected: Usage: vep [options]

iqtree2 --version
# Expected: IQ-TREE 2.3.3

# Test python scripts
python3 bin/vcf_to_phy.py --help
# Expected: usage: vcf_to_phy.py [-h] vcf output
```

## Troubleshooting

### Java Memory Issues
If you encounter Java memory errors with GATK:
```bash
# Set Java heap size
export JAVA_TOOL_OPTIONS="-Xmx8g"
```

### VEP Cache
VEP requires cache files. If not present, they will be downloaded:
```bash
# Pre-download VEP cache (optional)
vep_install -a cf -s homo_sapiens -y GRCh38
```

### Permission Issues
```bash
# Make scripts executable
chmod +x scripts/*.sh
chmod +x run_pipeline.sh
```

## Next Steps

After installation, see [USAGE.md](USAGE.md) for running the pipeline.
