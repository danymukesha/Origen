#!/bin/bash
# =============================================================================
# Stage 5: Population Frequency Analysis - BCFtools and PLINK
# =============================================================================
# This script analyzes population allele frequencies:
# 1. BCFtools stats - Calculate variant statistics
# 2. PLINK - Convert VCF to PLINK format for downstream analysis
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
VCF_DIR="${2:-results/03_variants}"
OUTPUT_DIR="${5:-results/05_population}"
THREADS="${6:-8}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Input file
VCF="${VCF_DIR}/${SAMPLE_NAME}_snps.vcf"

if [[ ! -f "$VCF" ]]; then
    echo "Error: VCF file not found: $VCF"
    echo "Please run 03_variant_calling.sh first"
    exit 1
fi

echo "=========================================="
echo "Stage 5: Population Frequency Analysis"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Input VCF: $VCF"
echo "=========================================="

# Step 1: BCFtools stats
echo "[1/2] Running bcftools stats..."
bcftools stats -s - "$VCF" > "${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt"

echo "Stats: ${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt"

# Step 2: Convert VCF to PLINK format
echo "[2/2] Converting VCF to PLINK format..."
plink \
    --vcf "$VCF" \
    --make-bed \
    --out "${OUTPUT_DIR}/${SAMPLE_NAME}" \
    --threads "$THREADS" \
    --double-id \
    --allow-extra-chr

echo "PLINK files: ${OUTPUT_DIR}/${SAMPLE_NAME}.bed/bim/fam"

# Summary
echo ""
echo "Population Analysis Summary:"
echo "Total SNPs: $(wc -l < "${OUTPUT_DIR}/${SAMPLE_NAME}.bim")"
echo ""
echo "Allele frequency distribution:"
awk 'NR>1 {print $5}' "${OUTPUT_DIR}/${SAMPLE_NAME}.bim" | sort | uniq -c | head -10

echo ""
echo "=========================================="
echo "Population Analysis Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.bed"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.bim"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.fam"
echo ""
echo "Next step: Run 06_ancestry.sh"
