#!/bin/bash
# =============================================================================
# Stage 1: Quality Control - FastQC and Trimmomatic
# =============================================================================
# This script performs quality control on raw FASTQ files:
# 1. FastQC - Generates quality reports on raw reads
# 2. Trimmomatic - Trims adapters and low-quality bases
# =============================================================================

set -e

# Configuration
SAMPLE_NAME="${1:-sample}"
READS_DIR="${2:-data/reads}"
OUTPUT_DIR="${3:-results/01_qc}"
THREADS="${4:-8}"

# Create output directories
mkdir -p "${OUTPUT_DIR}/fastqc"
mkdir -p "${OUTPUT_DIR}/trimmomatic"

# Find input files
R1="${READS_DIR}/${SAMPLE_NAME}_1.fastq.gz"
R2="${READS_DIR}/${SAMPLE_NAME}_2.fastq.gz"

if [[ ! -f "$R1" ]]; then
    echo "Error: Read 1 not found: $R1"
    exit 1
fi

if [[ ! -f "$R2" ]]; then
    echo "Error: Read 2 not found: $R2"
    exit 1
fi

echo "=========================================="
echo "Stage 1: Quality Control"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Reads: $R1, $R2"
echo "Threads: $THREADS"
echo "=========================================="

# Step 1: FastQC on raw reads
echo "[1/2] Running FastQC..."
fastqc -t "$THREADS" \
    -o "${OUTPUT_DIR}/fastqc" \
    "$R1" "$R2"

echo "FastQC complete. Reports saved to ${OUTPUT_DIR}/fastqc"

# Step 2: Trimmomatic - Trim adapters and low-quality bases
echo "[2/2] Running Trimmomatic..."
trimmomatic PE -threads "$THREADS" \
    "$R1" "$R2" \
    "${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_1_paired.fastq.gz" \
    "${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_1_unpaired.fastq.gz" \
    "${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_2_paired.fastq.gz" \
    "${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Trimmomatic complete."
echo "Trimmed reads saved to ${OUTPUT_DIR}/trimmomatic"

# Summary
echo "=========================================="
echo "QC Stage Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/fastqc/${SAMPLE_NAME}_1_fastqc.html"
echo "  - ${OUTPUT_DIR}/fastqc/${SAMPLE_NAME}_2_fastqc.html"
echo "  - ${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_1_paired.fastq.gz"
echo "  - ${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_2_paired.fastq.gz"
echo ""
echo "Next step: Run alignment.sh"
