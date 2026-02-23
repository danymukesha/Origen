#!/bin/bash
# =============================================================================
# Stage 3: Variant Calling - GATK HaplotypeCaller
# =============================================================================
# This script calls variants from aligned reads:
# 1. GATK HaplotypeCaller - Calls germline variants
# 2. GATK VariantFiltration - Filters low-quality variants
# 3. GATK SelectVariants - Separates SNPs and indels
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
REFERENCE="${2:-data/reference/GRCh38.fa}"
BAM_DIR="${3:-results/02_alignment}"
OUTPUT_DIR="${4:-results/03_variants}"
THREADS="${5:-8}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Input file
BAM="${BAM_DIR}/${SAMPLE_NAME}_sorted.bam"

if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM file not found: $BAM"
    echo "Please run 02_alignment.sh first"
    exit 1
fi

echo "=========================================="
echo "Stage 3: Variant Calling"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Reference: $REFERENCE"
echo "BAM: $BAM"
echo "Threads: $THREADS"
echo "=========================================="

# Step 1: GATK HaplotypeCaller - Call germline variants
echo "[1/3] Running GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$BAM" \
    -O "${OUTPUT_DIR}/${SAMPLE_NAME}_raw.vcf" \
    --native-pair-hmm-threads "$THREADS" \
    -ERC GVCF

echo "Raw VCF: ${OUTPUT_DIR}/${SAMPLE_NAME}_raw.vcf"

# Step 2: VariantFiltration - Apply hard filters
echo "[2/3] Running GATK VariantFiltration..."
gatk VariantFiltration \
    -V "${OUTPUT_DIR}/${SAMPLE_NAME}_raw.vcf" \
    -O "${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf" \
    --filter-name "QD2" --filter "QD < 2.0" \
    --filter-name "QUAL30" --filter "QUAL < 30.0" \
    --filter-name "SOR3" --filter "SOR > 3.0" \
    --filter-name "MQ40" --filter "MQ < 40.0" \
    --filter-name "MQRankSum-12.5" --filter "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum-8" --filter "ReadPosRankSum < -8.0"

echo "Filtered VCF: ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf"

# Step 3: SelectVariants - Separate SNPs and indels
echo "[3/3] Running GATK SelectVariants..."

# Select SNPs
gatk SelectVariants \
    -V "${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf" \
    -select-type SNP \
    -O "${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf"

# Select Indels
gatk SelectVariants \
    -V "${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf" \
    -select-type INDEL \
    -O "${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf"

echo "SNPs: ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf"
echo "Indels: ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf"

# Summary statistics
echo ""
echo "Variant Statistics:"
bcftools stats -s - "${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf" | head -20

echo ""
echo "=========================================="
echo "Variant Calling Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_raw.vcf"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf"
echo ""
echo "Next step: Run 04_annotation.sh"
