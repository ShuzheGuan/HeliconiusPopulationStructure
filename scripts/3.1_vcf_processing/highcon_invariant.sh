#!/usr/bin/env bash
set -euo pipefail

# Purpose: create a high-confidence invariant-site callset.
# Input: invariant-site VCF file.
# Output: high-confidence invariant-site VCF after low-depth genotype masking and missingness filtering.
# Software: GATK VariantFiltration, bcftools, tabix

INPUT_VCF="<path/to/invariant_input.vcf.gz>"
OUTPUT_DIR="<path/to/highcon_invariant_output_dir>"
BASENAME="<output_prefix>"
THREADS=32

DP_THRESHOLD=8
MISSINGNESS_THRESHOLD=50

mkdir -p "$OUTPUT_DIR"

VCF_STEP1="${OUTPUT_DIR}/${BASENAME}_DP${DP_THRESHOLD}_GATKfiltered.vcf.gz"
VCF_STEP2="${OUTPUT_DIR}/${BASENAME}_DP${DP_THRESHOLD}_max${MISSINGNESS_THRESHOLD}missing.vcf.gz"

gatk VariantFiltration \
  -V "$INPUT_VCF" \
  --genotype-filter-expression "DP < ${DP_THRESHOLD}" \
  --genotype-filter-name "lowDP" \
  --set-filtered-genotype-to-no-call \
  -O "$VCF_STEP1"

tabix -p vcf "$VCF_STEP1"
bcftools stats --threads "$THREADS" "$VCF_STEP1" > "${VCF_STEP1%.vcf.gz}.stats.txt"

bcftools filter \
  -i "F_MISSING <= ${MISSINGNESS_THRESHOLD}/100" \
  --threads "$THREADS" \
  -Oz \
  -o "$VCF_STEP2" \
  "$VCF_STEP1"

tabix -p vcf "$VCF_STEP2"
bcftools stats --threads "$THREADS" "$VCF_STEP2" > "${VCF_STEP2%.vcf.gz}.stats.txt"
