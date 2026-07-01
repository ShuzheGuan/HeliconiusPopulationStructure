#!/usr/bin/env bash
set -euo pipefail

# Purpose: subset samples from a VCF to ingroup only and split sites by variant class.
# Input: input VCF file and sample ID list.
# Output: sample-subset VCF plus indel, biallelic SNP, multiallelic SNP, overlap, and invariant VCF files.
# Software: GATK SelectVariants, bcftools, tabix

INPUT_VCF="<path/to/input.vcf.gz>"
SAMPLE_LIST="<path/to/sample_id_list.txt>"
OUTPUT_DIR="<path/to/split_vcf_output_dir>"
BASENAME="<output_prefix>"
THREADS=32

mkdir -p "$OUTPUT_DIR"

SELECTED_VCF="${OUTPUT_DIR}/${BASENAME}.vcf.gz"
OUT_INDEL="${OUTPUT_DIR}/${BASENAME}_indel.vcf.gz"
OUT_BIALLELIC="${OUTPUT_DIR}/${BASENAME}_biallelic.vcf.gz"
OUT_MULTIALLELIC="${OUTPUT_DIR}/${BASENAME}_multiallelic.vcf.gz"
OUT_OVERLAP="${OUTPUT_DIR}/${BASENAME}_overlap.vcf.gz"
OUT_INVARIANT="${OUTPUT_DIR}/${BASENAME}_invariant.vcf.gz"

SAMPLE_ARGS=()
while read -r SAMPLE; do
  [[ -n "$SAMPLE" ]] || continue
  SAMPLE_ARGS+=(--sample-name "$SAMPLE")
done < "$SAMPLE_LIST"

gatk SelectVariants \
  -V "$INPUT_VCF" \
  -O "$SELECTED_VCF" \
  --remove-unused-alternates \
  "${SAMPLE_ARGS[@]}"

tabix -p vcf "$SELECTED_VCF"

bcftools view --threads "$THREADS" -v indels "$SELECTED_VCF" -Oz -o "$OUT_INDEL"
tabix -p vcf "$OUT_INDEL"
bcftools stats --threads "$THREADS" "$OUT_INDEL" > "${OUT_INDEL%.vcf.gz}.stats.txt"

bcftools view --threads "$THREADS" -v snps -m2 -M2 "$SELECTED_VCF" -Oz -o "$OUT_BIALLELIC"
tabix -p vcf "$OUT_BIALLELIC"
bcftools stats --threads "$THREADS" "$OUT_BIALLELIC" > "${OUT_BIALLELIC%.vcf.gz}.stats.txt"

bcftools view --threads "$THREADS" -v snps -m3 -e 'TYPE~"indel" || TYPE~"overlap"' "$SELECTED_VCF" -Oz -o "$OUT_MULTIALLELIC"
tabix -p vcf "$OUT_MULTIALLELIC"
bcftools stats --threads "$THREADS" "$OUT_MULTIALLELIC" > "${OUT_MULTIALLELIC%.vcf.gz}.stats.txt"

bcftools view --threads "$THREADS" -i 'TYPE~"overlap" && TYPE!~"indel"' "$SELECTED_VCF" -Oz -o "$OUT_OVERLAP"
tabix -p vcf "$OUT_OVERLAP"
bcftools stats --threads "$THREADS" "$OUT_OVERLAP" > "${OUT_OVERLAP%.vcf.gz}.stats.txt"

bcftools view --threads "$THREADS" -i 'ALT="."' "$SELECTED_VCF" -Oz -o "$OUT_INVARIANT"
tabix -p vcf "$OUT_INVARIANT"
bcftools stats --threads "$THREADS" "$OUT_INVARIANT" > "${OUT_INVARIANT%.vcf.gz}.stats.txt"
