#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.help = false
params.reads = "data/reads/*_{1,2}.fastq.gz"
params.reference = "data/reference/GRCh38.fa"
params.output_dir = "results"
params.sample_name = "sample"
params.threads = 8
params.platform = "illumina"

if (params.help) {
    log.info """
    ============================================
    Genomic Variant Origin Tracing Pipeline
    ============================================
    
    Usage:
        nextflow run main.nf --reads "data/reads/*_{1,2}.fastq.gz" \\
                             --reference data/reference/GRCh38.fa \\
                             --sample_name sample1
    
    Options:
        --reads              Input FASTQ files (glob pattern)
        --reference          Reference genome FASTA
        --sample_name        Sample name for output files
        --output_dir         Output directory
        --threads            Number of threads per process
        
    ============================================
    """
    exit 0
}

log.info "========================================"
log.info "Genomic Variant Origin Tracing Pipeline"
log.info "========================================"
log.info "Sample: ${params.sample_name}"
log.info "Reads: ${params.reads}"
log.info "Reference: ${params.reference}"
log.info "Output: ${params.output_dir}"
log.info "Threads: ${params.threads}"
log.info "========================================"

// Validate input files
reference_file = file(params.reference)
if (!reference_file.exists()) {
    error "Reference genome not found: ${params.reference}"
}

// ============================================
// STAGE 1: Quality Control
// ============================================

process FASTQC {
    tag "FastQC on ${sample_id}"
    publishDir "${params.output_dir}/01_qc/fastqc", mode: 'copy'
    conda 'bioconda::fastqc=0.12.1'
    
    input:
        tuple val(sample_id), path(reads)
    
    output:
        path("*_fastqc.zip"), emit: fastqc_zip
        path("*_fastqc.html"), emit: fastqc_html
    
    script:
    """
    fastqc -t ${params.threads} -o . $reads
    """
}

process TRIMMOMATIC {
    tag "Trimmomatic on ${sample_id}"
    publishDir "${params.output_dir}/01_qc/trimmomatic", mode: 'copy'
    conda 'bioconda::trimmomatic=0.39'
    
    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("*_paired.fastq.gz"), emit: trimmed_reads
        path("*_unpaired.fastq.gz"), emit: unpaired_reads
        path("*_trimming_report.txt"), emit: trimming_report
    
    script:
    def (r1, r2) = reads
    """
    trimmomatic PE -threads ${params.threads} \\
        $r1 $r2 \\
        ${sample_id}_1_paired.fastq.gz ${sample_id}_1_unpaired.fastq.gz \\
        ${sample_id}_2_paired.fastq.gz ${sample_id}_2_unpaired.fastq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// ============================================
// STAGE 2: Alignment
// ============================================

process BWA_INDEX {
    label 'bwa_index'
    publishDir "${params.output_dir}/00_reference", mode: 'copy'
    conda 'bioconda::bwa=0.7.18'
    
    input:
        path(reference)
    
    output:
        path("*.amb"), emit: amb
        path("*.ann"), emit: ann
        path("*.bwt"), emit: bwt
        path("*.pac"), emit: pac
        path("*.sa"), emit: sa
    
    script:
    """
    bwa index -a bwtsw $reference
    """
}

process BWA_MEM {
    tag "BWA-MEM on ${sample_id}"
    publishDir "${params.output_dir}/02_alignment", mode: 'copy'
    conda 'bioconda::bwa=0.7.18'
    
    input:
        tuple val(sample_id), path(reads)
        path(reference)
    
    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    
    script:
    def (r1, r2) = reads
    """
    bwa mem -t ${params.threads} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${params.platform}" \\
        $reference $r1 $r2 > ${sample_id}.sam
    """
}

process SAMTOOLS_SORT {
    tag "SAMtools sort on ${sample_id}"
    publishDir "${params.output_dir}/02_alignment", mode: 'copy'
    conda 'bioconda::samtools=1.21'
    
    input:
        tuple val(sample_id), path(sam)
    
    output:
        tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
        tuple val(sample_id), path("${sample_id}_sorted.bam.bai"), emit: bai
    
    script:
    """
    samtools sort -@ ${params.threads} -o ${sample_id}_sorted.bam $sam
    samtools index -@ ${params.threads} ${sample_id}_sorted.bam
    """
}

// ============================================
// STAGE 3: Variant Calling
// ============================================

process GATK_HAPLOTYPECALLER {
    tag "GATK HaplotypeCaller on ${sample_id}"
    publishDir "${params.output_dir}/03_variants", mode: 'copy'
    conda 'bioconda::gatk4=4.5.0.0'
    
    input:
        tuple val(sample_id), path(bam), path(bai)
        path(reference)
    
    output:
        tuple val(sample_id), path("${sample_id}_raw.vcf"), emit: vcf
    
    script:
    """
    gatk HaplotypeCaller \\
        -R $reference \\
        -I $bam \\
        -O ${sample_id}_raw.vcf \\
        --native-pair-hmm-threads ${params.threads} \\
        -ERC GVCF
    """
}

process GATK_VARIANTFILTRATION {
    tag "GATK VariantFiltration on ${sample_id}"
    publishDir "${params.output_dir}/03_variants", mode: 'copy'
    conda 'bioconda::gatk4=4.5.0.0'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        tuple val(sample_id), path("${sample_id}_filtered.vcf"), emit: filtered_vcf
    
    script:
    """
    gatk VariantFiltration \\
        -V $vcf \\
        -O ${sample_id}_filtered.vcf \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    """
}

process SELECT_VARIANTS {
    tag "SelectVariants on ${sample_id}"
    publishDir "${params.output_dir}/03_variants", mode: 'copy'
    conda 'bioconda::gatk4=4.5.0.0'
    
    input:
        tuple val(sample_id), path(filtered_vcf)
    
    output:
        tuple val(sample_id), path("${sample_id}_snps.vcf"), emit: snps_vcf
        tuple val(sample_id), path("${sample_id}_indels.vcf"), emit: indels_vcf
    
    script:
    """
    gatk SelectVariants \\
        -V $filtered_vcf \\
        -select-type SNP \\
        -O ${sample_id}_snps.vcf
    
    gatk SelectVariants \\
        -V $filtered_vcf \\
        -select-type INDEL \\
        -O ${sample_id}_indels.vcf
    """
}

// ============================================
// STAGE 4: Functional Annotation
// ============================================

process VEP_ANNOTATION {
    tag "VEP annotation on ${sample_id}"
    publishDir "${params.output_dir}/04_annotation", mode: 'copy'
    conda 'bioconda::ensembl-vep=2024.3'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        tuple val(sample_id), path("${sample_id}_annotated.vcf"), emit: annotated_vcf
        tuple val(sample_id), path("${sample_id}_vep_summary.html"), emit: vep_html
    
    script:
    """
    vep -i $vcf \\
        -o ${sample_id}_annotated.vcf \\
        --format vcf \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        --html \\
        --stats_file ${sample_id}_vep_summary.html \\
        --cache --offline \\
        --numbers \\
        --total_length \\
        --show_ref_allele
    """
}

// ============================================
// STAGE 5: Population Frequency Analysis
// ============================================

process BCFTOOLS_STATS {
    tag "bcftools stats on ${sample_id}"
    publishDir "${params.output_dir}/05_population", mode: 'copy'
    conda 'bioconda::bcftools=1.21'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        path("${sample_id}_stats.txt"), emit: stats
    
    script:
    """
    bcftools stats -s - $vcf > ${sample_id}_stats.txt
    """
}

process PLINK_CONVERSION {
    tag "PLINK conversion on ${sample_id}"
    publishDir "${params.output_dir}/05_population", mode: 'copy'
    conda 'bioconda::plink=1.9b'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        tuple val(sample_id), path("${sample_id}.bed"), emit: bed
        tuple val(sample_id), path("${sample_id}.bim"), emit: bim
        tuple val(sample_id), path("${sample_id}.fam"), emit: fam
    
    script:
    """
    plink --vcf $vcf \\
        --make-bed \\
        --out ${sample_id} \\
        --threads ${params.threads}
    """
}

// ============================================
// STAGE 6: Ancestry Analysis
// ============================================

process ADMIXTURE {
    tag "ADMIXTURE on ${sample_id}"
    publishDir "${params.output_dir}/06_ancestry", mode: 'copy'
    conda 'bioconda::admixture=1.3.0'
    
    input:
        tuple val(sample_id), path(bed), path(bim), path(fam)
    
    output:
        path("${sample_id}.${params.K}.Q"), emit: q_matrix
        path("${sample_id}.${params.K}.P"), emit: p_matrix
    
    script:
    """
    admixture --cv=${params.threads} ${bed}.bed ${params.K} > ${sample_id}_admixture.log
    mv ${sample_id}.${params.K}.Q ${params.output_dir}/06_ancestry/
    mv ${sample_id}.${params.K}.P ${params.output_dir}/06_ancestry/
    """
}

process HAPLOGREP {
    tag "HaploGrep on ${sample_id}"
    publishDir "${params.output_dir}/06_ancestry", mode: 'copy'
    conda 'bioconda::haplogrep=2.1.1'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        path("${sample_id}_haplogroup.txt"), emit: haplogroup
    
    script:
    """
    haplogrep2 --vcf $vcf --out ${sample_id}_haplogroup.txt --format vcf
    """
}

// ============================================
// STAGE 7: Phylogenetic Analysis
// ============================================

process VCF_TO_PHYLIP {
    tag "VCF to Phylip on ${sample_id}"
    publishDir "${params.output_dir}/07_phylogeny", mode: 'copy'
    
    input:
        tuple val(sample_id), path(vcf)
    
    output:
        path("${sample_id}.phy"), emit: phylip
    
    script:
    """
    python3 ${projectDir}/bin/vcf_to_phy.py $vcf ${sample_id}.phy
    """
}

process IQTREE {
    tag "IQ-TREE on ${sample_id}"
    publishDir "${params.output_dir}/07_phylogeny", mode: 'copy'
    conda 'bioconda::iqtree=2.3.3'
    
    input:
        tuple val(sample_id), path(phylip)
    
    output:
        path("${sample_id}_tree.nwk"), emit: tree
        path("${sample_id}_iqtree"), emit: iqtree_out
    
    script:
    """
    iqtree2 -s ${phylip} \\
        -m MFP \\
        -bb 1000 \\
        -nt ${params.threads} \\
        -pre ${sample_id}_iqtree
    """
}

// ============================================
// WORKFLOW DEFINITION
// ============================================

workflow {
    // Parse input reads
    Channel.fromFilePairs(params.reads, flat: true)
        .ifEmpty { error "No reads found matching: ${params.reads}" }
        .set { reads_ch }
    
    sample_id = params.sample_name
    
    // Stage 1: QC
    FASTQC(reads_ch)
    TRIMMOMATIC(reads_ch)
    trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    
    // Stage 2: Alignment
    bwa_index_ch = BWA_INDEX(params.reference)
    BWA_MEM(trimmed_reads, params.reference)
    SAMTOOLS_SORT(BWA_MEM.out.sam)
    
    // Stage 3: Variant Calling
    GATK_HAPLOTYPECALLER(SAMTOOLS_SORT.out.bam, params.reference)
    GATK_VARIANTFILTRATION(GATK_HAPLOTYPECALLER.out.vcf)
    SELECT_VARIANTS(GATK_VARIANTFILTRATION.out.filtered_vcf)
    
    // Stage 4: Annotation
    VEP_ANNOTATION(SELECT_VARIANTS.out.snps_vcf)
    
    // Stage 5: Population Frequency
    BCFTOOLS_STATS(SELECT_VARIANTS.out.snps_vcf)
    PLINK_CONVERSION(SELECT_VARIANTS.out.snps_vcf)
    
    // Stage 6: Ancestry
    ADMIXTURE(PLINK_CONVERSION.out.bed.join(PLINK_CONVERSION.out.bim).join(PLINK_CONVERSION.out.fam))
    HAPLOGREP(SELECT_VARIANTS.out.snps_vcf)
    
    // Stage 7: Phylogeny
    VCF_TO_PHYLIP(SELECT_VARIANTS.out.snps_vcf)
    IQTREE(VCF_TO_PHYLIP.out.phylip)
}
