# Usage Guide

## Input Requirements

### FASTQ Files
- Format: FASTQ (`.fastq.gz` or `.fq.gz`)
- Naming: `{sample_name}_1.fastq.gz` and `{sample_name}_2.fastq.gz`
- Single-end: `{sample_name}.fastq.gz`
- Quality: Illumina phred+33 or phred+64

### Reference Genome
- Format: FASTA (`.fa` or `.fasta`)
- Recommended: GRCh38.primary assembly
- Must be indexed with BWA before use

## Running the Pipeline

### Option 1: Using Nextflow (Recommended)

```bash
# Basic usage
nextflow run main.nf \
    --sample_name NA12878 \
    --reference data/reference/GRCh38.fa \
    --reads "data/reads/NA12878_{1,2}.fastq.gz"

# With custom parameters
nextflow run main.nf \
    --sample_name sample1 \
    --reference /path/to/GRCh38.fa \
    --reads "data/reads/*_{1,2}.fastq.gz" \
    --output_dir my_results \
    --threads 16 \
    -resume
```

The `-resume` flag continues from the last completed step if the pipeline was interrupted.

### Option 2: Using Shell Scripts

```bash
# Make scripts executable
chmod +x scripts/*.sh run_pipeline.sh

# Run full pipeline
bash run_pipeline.sh sample1 8

# Or run individual stages
bash scripts/01_qc.sh sample1
bash scripts/02_alignment.sh sample1 data/reference/GRCh38.fa
bash scripts/03_variant_calling.sh sample1 data/reference/GRCh38.fa
# ... etc
```

## Pipeline Parameters

### Nextflow Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Input FASTQ glob pattern | `data/reads/*_{1,2}.fastq.gz` |
| `--reference` | Reference genome FASTA | `data/reference/GRCh38.fa` |
| `--sample_name` | Sample identifier | `sample` |
| `--output_dir` | Output directory | `results` |
| `--threads` | Threads per process | `8` |
| `--platform` | Sequencing platform | `illumina` |
| `-resume` | Resume from last checkpoint | false |

### Shell Script Parameters

```bash
# Stage 1: QC
bash scripts/01_qc.sh <sample_name> [reads_dir] [output_dir] [threads]

# Stage 2: Alignment  
bash scripts/02_alignment.sh <sample_name> [reference] [trimmed_dir] [output_dir] [threads]

# Stage 3: Variant Calling
bash scripts/03_variant_calling.sh <sample_name> [reference] [bam_dir] [output_dir] [threads]

# Stage 4: Annotation
bash scripts/04_annotation.sh <sample_name> [vcf_dir] [output_dir]

# Stage 5: Population
bash scripts/05_population.sh <sample_name> [vcf_dir] [output_dir] [threads]

# Stage 6: Ancestry
bash scripts/06_ancestry.sh <sample_name> [plink_dir] [vcf_dir] [output_dir] [K] [threads]

# Stage 7: Phylogeny
bash scripts/07_phylogeny.sh <sample_name> [vcf_dir] [output_dir] [threads]
```

## Sample Commands

### Running on Multiple Samples

```bash
# Process multiple samples with Nextflow
nextflow run main.nf \
    --reads "data/reads/{sample1,sample2,sample3}_{1,2}.fastq.gz" \
    --reference data/reference/GRCh38.fa
```

### Custom ADMIXTURE K Value

```bash
# Run ancestry analysis with K=5 populations
nextflow run main.nf \
    --K 5 \
    --sample_name sample1 \
    -profile conda
```

### Using Singularity Container

```bash
# If Singularity is available
nextflow run main.nf -profile singularity
```

## Output Interpretation

### Variant Calling Results

The filtered VCF contains high-confidence variant calls:

```bash
# View variant summary
bcftools stats results/03_variants/sample1_snps.vcf

# Count variants
grep -v "^#" results/03_variants/sample1_snps.vcf | wc -l
```

### Annotation Results

VEP annotation adds CSQ field with predictions:

```bash
# Extract coding consequences
grep "missense" results/04_annotation/sample1_annotated.vcf | head

# Extract high-impact variants
grep "HIGH" results/04_annotationsample1_annotated.vcf | head
```

### Ancestry Results

ADMIXTURE Q-matrix shows ancestry proportions:

```bash
# View ancestry proportions
cat results/06_ancestry/sample1.3.Q

# View haplogroup
cat results/06_ancestry/sample1_haplogroup.txt
```

### Phylogenetic Results

IQ-TREE output includes tree file:

```bash
# View tree (requires figtree or similar)
# results/07_phylogeny/sample1_iqtree.treefile
```

## Batch Processing

For processing multiple samples:

```bash
# Create list of samples
ls data/reads/*_1.fastq.gz | sed 's/_1.fastq.gz//' | sed 's/data\/reads\///' > samples.txt

# Run in parallel with GNU parallel
cat samples.txt | parallel -j 4 'nextflow run main.nf --sample_name {}'
```

## Performance Optimization

### Memory Settings

```bash
# For GATK (in nextflow.config)
withName: GATK_HAPLOTYPECALLER {
    memory = '16 GB'
}

# For shell scripts
export JAVA_TOOL_OPTIONS="-Xmx8g"
```

### Thread Settings

```bash
# Adjust based on available cores
# In nextflow.config
process {
    withName: BWA_MEM {
        cpus = 16
    }
}
```

## Common Issues

### Error: "No reads found"

Check FASTQ naming:
```bash
# Should match pattern: sample_1.fastq.gz, sample_2.fastq.gz
ls data/reads/
```

### Error: "Reference not indexed"

BWA index is created automatically. If issues:
```bash
bwa index data/reference/GRCh38.fa
samtools faidx data/reference/GRCh38.fa
```

### Error: VCF empty after variant calling

Check BAM quality:
```bash
samtools flagstat results/02_alignment/sample1_sorted.bam
```

## Validation

To validate pipeline works correctly:

```bash
# Test with known sample
# Download GIAB sample and run pipeline
# Compare results to known truth set
```

## Support

For issues and questions:
1. Check tool documentation
2. Review log files in `results/.nextflow/`
3. Verify input data quality with FastQC
