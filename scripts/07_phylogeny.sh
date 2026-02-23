#!/bin/bash
# =============================================================================
# Stage 7: Phylogenetic Analysis - IQ-TREE
# =============================================================================
# This script performs phylogenetic analysis:
# 1. VCF to Phylip conversion
# 2. IQ-TREE - Maximum likelihood phylogenetic tree
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
VCF_DIR="${2:-results/03_variants}"
OUTPUT_DIR="${3:-results/07_phylogeny}"
THREADS="${4:-8}"

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
echo "Stage 7: Phylogenetic Analysis"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Input VCF: $VCF"
echo "=========================================="

# Step 1: Convert VCF to Phylip
echo "[1/2] Converting VCF to Phylip format..."
python3 "$(dirname "$0")/../bin/vcf_to_phy.py" "$VCF" "${OUTPUT_DIR}/${SAMPLE_NAME}.phy"

echo "Phylip: ${OUTPUT_DIR}/${SAMPLE_NAME}.phy"

# Count sites
NUM_SITES=$(tail -n +2 "${OUTPUT_DIR}/${SAMPLE_NAME}.phy" | head -1 | awk '{print length($2)}')
echo "Alignment length: $NUM_SITES sites"

if [[ "$NUM_SITES" -lt 100 ]]; then
    echo "Warning: Low number of sites for phylogenetic analysis"
    echo "Results may not be reliable with <100 sites"
fi

# Step 2: IQ-TREE phylogenetic reconstruction
echo "[2/2] Running IQ-TREE..."
iqtree2 -s "${OUTPUT_DIR}/${SAMPLE_NAME}.phy" \
    -m MFP \
    -bb 1000 \
    -nt "$THREADS" \
    -pre "${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree"

echo "IQ-TREE complete."

# Summary
echo ""
echo "Phylogenetic Analysis Summary:"
echo "Best model: $(grep "Best model" "${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree.iqtree" | awk '{print $3}')"
echo "Tree file: ${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree.treefile"

echo ""
echo "=========================================="
echo "Phylogenetic Analysis Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}.phy"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree.treefile"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree.iqtree"
echo ""
echo "To visualize the tree:"
echo "  iqtree2 -s ${OUTPUT_DIR}/${SAMPLE_NAME}.phy -t ${OUTPUT_DIR}/${SAMPLE_NAME}_iqtree.treefile --trees"
