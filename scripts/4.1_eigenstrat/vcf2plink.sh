#!/usr/bin/env bash
set -euo pipefail

# Purpose: convert a high-confidence VCF to binary PLINK format for PCA/EIGENSTRAT processing.
# Input: high-confidence biallelic VCF file.
# Output: PLINK .bed, .bim, and .fam files.
# Software: PLINK

INPUT_VCF="<path/to/highconfidence_biallelic_filtered.vcf.gz>"
OUTPUT_DIR="<path/to/plink_output_dir>"
BASENAME="<output_prefix>"
THREADS=16

mkdir -p "$OUTPUT_DIR"

plink \
  --vcf "$INPUT_VCF" \
  --allow-extra-chr \
  --threads "$THREADS" \
  --make-bed \
  --out "${OUTPUT_DIR}/${BASENAME}"
