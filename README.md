# Genomic Variant Origin Tracing Pipeline

A bioinformatics pipeline for tracing the origin of genomic variants from raw sequencing data (FASTQ) to ancestry, population, and phylogenetic analysis.

## Overview

This pipeline processes raw FASTQ files through a complete analysis workflow to identify variants, annotate their functional significance, determine population origins, and place them in an evolutionary context.

## Pipeline Stages

### Stage 1: Quality Control (QC)
**Tools:** FastQC, Trimmomatic

- **FastQC**: Generates comprehensive quality reports on raw sequencing reads, including per-base quality scores, GC content, sequence length distribution, and adapter content
- **Trimmomatic**: Removes adapter sequences and low-quality bases; filters reads below minimum length threshold

**Why it matters:** Poor-quality reads can cause false positive variant calls. QC ensures downstream analysis uses high-quality data.

### Stage 2: Alignment
**Tools:** BWA-MEM, SAMtools

- **BWA-MEM**: Aligns trimmed reads to the human GRCh38 reference genome using the Burrows-Wheeler Transform algorithm
- **SAMtools**: Converts SAM to sorted BAM format and creates indices for efficient access

**Why it matters:** Accurate alignment is critical for variant calling. BWA-MEM is optimal for Illumina reads >70bp.

### Stage 3: Variant Calling
**Tools:** GATK HaplotypeCaller, VariantFiltration, SelectVariants

- **HaplotypeCaller**: Performs local de novo assembly to identify germline variants (SNPs and indels)
- **VariantFiltration**: Applies hard filters to remove low-quality calls based on statistical thresholds
- **SelectVariants**: Separates SNPs and indels into distinct files

**Why it matters:** GATK HaplotypeCaller is the gold standard for germline variant discovery, with high sensitivity and specificity.

### Stage 4: Functional Annotation
**Tools:** Ensembl VEP (Variant Effect Predictor)

- **VEP**: Annotates variants with functional consequences including:
  - Gene and transcript information
  - Protein changes (missense, nonsense, frameshift)
  - Conservation scores (PhyloP, CADD)
  - Pathogenicity predictions
  - Population frequencies (gnomAD, TOPMed)

**Why it matters:** Understanding functional impact helps prioritize variants for downstream analysis and clinical interpretation.

### Stage 5: Population Frequency Analysis
**Tools:** BCFtools, PLINK

- **BCFtools stats**: Calculates variant statistics including Ti/Tv ratio, allele frequencies, and transition/transversion counts
- **PLINK**: Converts VCF to PLINK format (BED/BIM/FAM) for population genetics analyses

**Why it matters:** Population frequency helps distinguish common benign variants from rare pathogenic mutations.

### Stage 6: Ancestry Inference
**Tools:** ADMIXTURE, HaploGrep

- **ADMIXTURE**: Estimates ancestral population proportions using maximum likelihood
- **HaploGrep**: Determines mitochondrial DNA haplogroup from variants

**Why it matters:** Ancestry context is essential for interpreting variant pathogenicity and understanding evolutionary origins.

### Stage 7: Phylogenetic Analysis
**Tools:** IQ-TREE

- **IQ-TREE**: Builds maximum-likelihood phylogenetic trees from SNP alignments
- Uses ModelFinder to select optimal substitution model
- Performs ultrafast bootstrap for branch support

**Why it matters:** Phylogenetic placement reveals evolutionary relationships and can identify shared ancestry between samples.

## Directory Structure

```
genomic-variant-tracing/
├── main.nf                 # Nextflow workflow definition
├── nextflow.config         # Nextflow configuration
├── environment.yml         # Conda environment specification
├── run_pipeline.sh         # Master pipeline runner
├── bin/
│   └── vcf_to_phy.py       # VCF to Phylip converter
├── scripts/
│   ├── 01_qc.sh            # Quality control
│   ├── 02_alignment.sh     # Read alignment
│   ├── 03_variant_calling.sh  # Variant discovery
│   ├── 04_annotation.sh    # Functional annotation
│   ├── 05_population.sh    # Population frequency
│   ├── 06_ancestry.sh      # Ancestry inference
│   └── 07_phylogeny.sh     # Phylogenetic analysis
├── docs/
│   └── INSTALL.md          # Installation guide
└── data/
    ├── reads/              # Input FASTQ files
    └── reference/          # Reference genome
```

## Quick Start

### 1. Install Dependencies
```bash
# Create conda environment
conda env create -f environment.yml
conda activate variant-tracing

# Install Nextflow (optional)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Prepare Input Data
```bash
# Place FASTQ files in data/reads/
# Naming convention: sample_1.fastq.gz, sample_2.fastq.gz

# Ensure reference genome exists in data/reference/GRCh38.fa
```

### 3. Run Pipeline
```bash
# Using Nextflow (recommended)
nextflow run main.nf \
    --sample_name sample1 \
    --reference data/reference/GRCh38.fa \
    --reads "data/reads/sample1_{1,2}.fastq.gz"

# Or using shell scripts
bash run_pipeline.sh sample1 8
```

## Output Files

| Stage | Directory | Files |
|-------|----------|-------|
| QC | `results/01_qc/` | FastQC HTML reports, trimmed FASTQs |
| Alignment | `results/02_alignment/` | Sorted BAM, BAI index |
| Variants | `results/03_variants/` | Raw/filtered VCF, SNPs, indels |
| Annotation | `results/04_annotation/` | Annotated VCE, VEP summary |
| Population | `results/05_population/` | Stats, PLINK files |
| Ancestry | `results/06_ancestry/` | ADMIXTURE Q-matrix, haplogroup |
| Phylogeny | `results/07_phylogeny/` | Phylip alignment, tree file |

## Pipeline Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Input FASTQ glob pattern | `data/reads/*_{1,2}.fastq.gz` |
| `--reference` | Reference genome FASTA | `data/reference/GRCh38.fa` |
| `--sample_name` | Sample identifier | `sample` |
| `--output_dir` | Output directory | `results` |
| `--threads` | Number of threads | `8` |
| `--K` | ADMIXTURE K value | `3` |

## Testing with Sample Data

To verify the pipeline works correctly:

```bash
# Create test directory structure
mkdir -p data/reads data/reference

# Download a small test dataset (NA12878 subset)
# Place test files in data/reads/

# Run pipeline on test data
bash run_pipeline.sh test_sample 4
```

## Performance Notes

- **Runtime**: ~4-8 hours for 30x WGS sample (8 threads)
- **Storage**: ~50GB for intermediate files, ~1GB final output
- **Memory**: 16GB RAM recommended, 32GB for GATK

## Troubleshooting

See [INSTALL.md](docs/INSTALL.md) for detailed troubleshooting guide.

Common issues:
- Java memory errors: Set `export JAVA_TOOL_OPTIONS="-Xmx8g"`
- Missing VEP cache: Auto-downloads on first run
- Tool not found: Ensure conda environment is activated

## Citation

If you use this pipeline, please cite the individual tools:

- Li H., Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics.
- McKenna A. et al. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res.
- Danecek P. et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience.
- McLaren W. et al. (2016) The Ensembl Variant Effect Predictor. Genome Biology.
- Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics.
- Min B et al. (2024) IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference. Molecular Biology and Evolution.

## License

MIT License - See LICENSE file for details.
