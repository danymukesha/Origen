#!/bin/bash
# =============================================================================
# Stage 6: Ancestry Analysis - ADMIXTURE and HaploGrep
# =============================================================================
# This script performs ancestry inference:
# 1. ADMIXTURE - Estimates ancestral population proportions
# 2. HaploGrep - Determines mitochondrial haplogroup
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
PLINK_DIR="${2:-results/05_population}"
VCF_DIR="${3:-results/03_variants}"
OUTPUT_DIR="${4:-results/06_ancestry}"
K="${5:-3}"  # Number of ancestral populations
THREADS="${6:-8}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Input files
BED="${PLINK_DIR}/${SAMPLE_NAME}.bed"
BIM="${PLINK_DIR}/${SAMPLE_NAME}.bim"
FAM="${PLINK_DIR}/${SAMPLE_NAME}.fam"
VCF="${VCF_DIR}/${SAMPLE_NAME}_snps.vcf"

echo "=========================================="
echo "Stage 6: Ancestry Analysis"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "K (populations): $K"
echo "=========================================="

# Step 1: ADMIXTURE
if [[ -f "$BED" && -f "$BIM" && -f "$FAM" ]]; then
    echo "[1/2] Running ADMIXTURE..."
    admixture --cv="$THREADS" "$BED" "$K" > "${OUTPUT_DIR}/${SAMPLE_NAME}_admixture.log"
    
    mv "${SAMPLE_NAME}.${K}.Q" "${OUTPUT_DIR}/"
    mv "${SAMPLE_NAME}.${K}.P" "${OUTPUT_DIR}/"
    
    echo "ADMIXTURE Q-matrix: ${OUTPUT_DIR}/${SAMPLE_NAME}.${K}.Q"
    echo "ADMIXTURE P-matrix: ${OUTPUT_DIR}/${SAMPLE_NAME}.${K}.P"
else
    echo "Warning: PLINK files not found, skipping ADMIXTURE"
fi

# Step 2: HaploGrep for mtDNA haplogroup
if [[ -f "$VCF" ]]; then
    echo "[2/2] Running HaploGrep..."
    haplogrep2 \
        --vcf "$VCF" \
        --out "${OUTPUT_DIR}/${SAMPLE_NAME}_haplogroup.txt" \
        --format vcf
    
    echo "Haplogroup: $(tail -5 "${OUTPUT_DIR}/${SAMPLE_NAME}_haplogroup.txt" | head -1)"
else
    echo "Warning: VCF file not found, skipping HaploGrep"
fi

echo ""
echo "=========================================="
echo "Ancestry Analysis Complete!"
echo "=========================================="
echo "Output files:"
[[ -f "$BED" ]] && echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.${K}.Q"
[[ -f "$BED" ]] && echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.${K}.P"
[[ -f "$VCF" ]] && echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_haplogroup.txt"
echo ""
echo "Next step: Run 07_phylogeny.sh"
