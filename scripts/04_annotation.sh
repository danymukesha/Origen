#!/bin/bash
# =============================================================================
# Stage 4: Functional Annotation - Ensembl VEP
# =============================================================================
# This script annotates variants with functional information:
# - Gene and transcript information
# - Protein changes (amino acid substitutions)
# - Conservation scores
# - Pathogenicity predictions
# - Population frequency databases
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
VCF_DIR="${2:-results/03_variants}"
OUTPUT_DIR="${4:-results/04_annotation}"

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
echo "Stage 4: Functional Annotation"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Input VCF: $VCF"
echo "=========================================="

# Run VEP annotation
echo "[1/1] Running Ensembl VEP..."
vep \
    -i "$VCF" \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}_annotated.vcf" \
    --format vcf \
    --species homo_sapiens \
    --assembly GRCh38 \
    --html \
    --stats_file "${OUTPUT_DIR}/${SAMPLE_NAME}_vep_summary.html" \
    --cache --offline \
    --numbers \
    --total_length \
    --show_ref_allele \
    --canonical \
    --transcript_version \
    --uniprot \
    --domains \
    --regulatory \
    --celldata \
    --gene_phenotype \
    --variant_class

echo "VEP annotation complete."

# Summary of annotation
echo ""
echo "Annotation Summary:"
echo "Total variants annotated: $(grep -v "^#" "${OUTPUT_DIR}/${SAMPLE_NAME}_annotated.vcf" | wc -l)"

echo ""
echo "=========================================="
echo "Annotation Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_annotated.vcf"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_vep_summary.html"
echo ""
echo "Next steps:"
echo "  - Run 05_population.sh for population frequency"
echo "  - Run 06_ancestry.sh for ancestry analysis"
echo "  - Run 07_phylogeny.sh for phylogenetic analysis"
