#!/usr/bin/env python3
"""
Convert VCF to Phylip format for phylogenetic analysis.
This script extracts SNP positions from a VCF file and creates a phylip alignment.
"""

import sys
import argparse
from collections import defaultdict


def parse_vcf(vcf_file):
    """Parse VCF file and extract SNP information."""
    samples = []
    snps = defaultdict(dict)
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    samples = parts[9:]  # Sample columns start at index 9
            else:
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                
                # Only process biallelic SNPs
                if len(ref) == 1 and len(alt) == 1 and ',' not in alt:
                    for i, sample in enumerate(samples):
                        geno = parts[9 + i].split(':')[0]
                        if geno in ['0/0', '0/1', '1/0', '1/1', './.']:
                            if geno in ['0/0', '0/1', '1/0']:
                                snps[pos][sample] = ref if geno in ['0/0'] else alt
                            elif geno == '1/1':
                                snps[pos] = alt
                            else:
                                snps[pos][sample] = 'N'
    
    return samples, snps


def write_phylip(samples, snps, output_file):
    """Write phylip format alignment."""
    num_sites = len(snps)
    num_seqs = len(samples)
    
    with open(output_file, 'w') as f:
        f.write(f"{num_seqs} {num_sites}\n")
        
        for sample in samples:
            seq = ''.join([snps[pos].get(sample, 'N') for pos in sorted(snps.keys())])
            # Phylip format: max 10 chars for sequence name, space-padded
            f.write(f"{sample[:10]:<10}{seq}\n")


def main():
    parser = argparse.ArgumentParser(description='Convert VCF to Phylip format')
    parser.add_argument('vcf', help='Input VCF file')
    parser.add_argument('output', help='Output Phylip file')
    
    args = parser.parse_args()
    
    samples, snps = parse_vcf(args.vcf)
    write_phylip(samples, snps, args.output)
    
    print(f"Converted {len(snps)} SNPs for {len(samples)} samples to {args.output}")


if __name__ == '__main__':
    main()
