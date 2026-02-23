#!/bin/bash
# =============================================================================
# Master Pipeline Runner
# =============================================================================
# This script runs all pipeline stages sequentially.
# Usage: ./run_pipeline.sh <sample_name> [threads]
# =============================================================================

set -e

SAMPLE_NAME="${1:-sample}"
THREADS="${2:-8}"
REFERENCE="${3:-data/reference/GRCh38.fa}"

echo "=========================================="
echo "GENOMIC VARIANT ORIGIN TRACING PIPELINE"
echo "=========================================="
echo "Sample: $SAMPLE_NAME"
echo "Threads: $THREADS"
echo "Reference: $REFERENCE"
echo "=========================================="

# Check dependencies
command -v nextflow >/dev/null 2>&1 && USE_NEXTFLOW=1 || USE_NEXTFLOW=0

if [[ "$USE_NEXTFLOW" -eq 1 ]]; then
    echo "Using Nextflow workflow"
    echo ""
    
    nextflow run main.nf \
        --sample_name "$SAMPLE_NAME" \
        --reference "$REFERENCE" \
        --reads "data/reads/${SAMPLE_NAME}_*{1,2}.fastq.gz" \
        --threads "$THREADS" \
        -profile conda
    
else
    echo "Nextflow not found, running shell scripts sequentially"
    echo ""
    
    # Stage 1: QC
    echo "=========================================="
    echo "Starting Stage 1: Quality Control"
    echo "=========================================="
    bash scripts/01_qc.sh "$SAMPLE_NAME" "data/reads" "results" "$THREADS"
    
    # Stage 2: Alignment
    echo "=========================================="
    echo "Starting Stage 2: Alignment"
    echo "=========================================="
    bash scripts/02_alignment.sh "$SAMPLE_NAME" "$REFERENCE" "results/01_qc/trimmomatic" "results" "$THREADS"
    
    # Stage 3: Variant Calling
    echo "=========================================="
    echo "Starting Stage 3: Variant Calling"
    echo "=========================================="
    bash scripts/03_variant_calling.sh "$SAMPLE_NAME" "$REFERENCE" "results/02_alignment" "results" "$THREADS"
    
    # Stage 4: Annotation
    echo "=========================================="
    echo "Starting Stage 4: Annotation"
    echo "=========================================="
    bash scripts/04_annotation.sh "$SAMPLE_NAME" "results/03_variants" "results" "results/04_annotation"
    
    # Stage 5: Population Analysis
    echo "=========================================="
    echo "Starting Stage 5: Population Analysis"
    echo "=========================================="
    bash scripts/05_population.sh "$SAMPLE_NAME" "results/03_variants" "results" "results/05_population" "$THREADS"
    
    # Stage 6: Ancestry Analysis
    echo "=========================================="
    echo "Starting Stage 6: Ancestry Analysis"
    echo "=========================================="
    bash scripts/06_ancestry.sh "$SAMPLE_NAME" "results/05_population" "results/03_variants" "results/06_ancestry" "3" "$THREADS"
    
    # Stage 7: Phylogenetic Analysis
    echo "=========================================="
    echo "Starting Stage 7: Phylogenetic Analysis"
    echo "=========================================="
    bash scripts/07_phylogeny.sh "$SAMPLE_NAME" "results/03_variants" "results/07_phylogeny" "$THREADS"
fi

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE!"
echo "=========================================="
echo "Results available in: results/"
echo "  01_qc/         - Quality control reports"
echo "  02_alignment/  - Aligned BAM files"
echo "  03_variants/   - Variant calls (VCF)"
echo "  04_annotation/ - Functional annotations"
echo "  05_population/ - Population frequency"
echo "  06_ancestry/   - Ancestry inference"
echo "  07_phylogeny/  - Phylogenetic trees"
echo "=========================================="
