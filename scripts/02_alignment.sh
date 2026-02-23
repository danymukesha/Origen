#!/bin/bash
# =============================================================================
# Stage 2: Alignment - BWA-MEM and SAMtools
# =============================================================================
# This script aligns trimmed reads to the reference genome:
# 1. BWA-MEM - Aligns reads using BWA-MEM algorithm
# 2. SAMtools - Converts SAM to sorted BAM and indexes
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
REFERENCE="${2:-data/reference/GRCh38.fa}"
TRIMMED_DIR="${3:-results/01_qc/trimmomatic}"
OUTPUT_DIR="${4:-results/02_alignment}"
THREADS="${5:-8}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Input files
R1_PAIRED="${TRIMMED_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
R2_PAIRED="${TRIMMED_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"

if [[ ! -f "$R1_PAIRED" ]]; then
    echo "Error: Trimmed read 1 not found: $R1_PAIRED"
    echo "Please run 01_qc.sh first"
    exit 1
fi

if [[ ! -f "$R2_PAIRED" ]]; then
    echo "Error: Trimmed read 2 not found: $R2_PAIRED"
    echo "Please run 01_qc.sh first"
    exit 1
fi

echo "=========================================="
echo "Stage 2: Alignment"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Reference: $REFERENCE"
echo "Threads: $THREADS"
echo "=========================================="

# Check if BWA index exists, create if not
BASENAME=$(basename "$REFERENCE")
INDEX_PREFIX="${REFERENCE%.*}"

if [[ ! -f "${INDEX_PREFIX}.bwt" ]]; then
    echo "[0/3] Creating BWA index..."
    bwa index -a bwtsw "$REFERENCE"
else
    echo "[0/3] BWA index found, skipping..."
fi

# Step 1: BWA-MEM alignment
echo "[1/3] Running BWA-MEM..."
bwa mem -t "$THREADS" \
    -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA" \
    "$REFERENCE" \
    "$R1_PAIRED" "$R2_PAIRED" \
    > "${OUTPUT_DIR}/${SAMPLE_NAME}.sam"

echo "BWA-MEM complete: ${OUTPUT_DIR}/${SAMPLE_NAME}.sam"

# Step 2: Convert SAM to BAM and sort
echo "[2/3] Converting SAM to sorted BAM..."
samtools sort -@ "$THREADS" \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}.sam"

echo "Sorted BAM: ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"

# Step 3: Index BAM
echo "[3/3] Indexing BAM..."
samtools index -@ "$THREADS" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"

# Remove SAM file to save space
rm "${OUTPUT_DIR}/${SAMPLE_NAME}.sam"

# Summary statistics
echo ""
echo "Alignment Statistics:"
samtools flagstat "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"

echo ""
echo "=========================================="
echo "Alignment Stage Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"
echo "  - ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam.bai"
echo ""
echo "Next step: Run 03_variant_calling.sh"
