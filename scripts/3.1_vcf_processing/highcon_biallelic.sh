#!/usr/bin/env bash
set -euo pipefail

# Purpose: create a high-confidence biallelic SNP callset.
# Input: biallelic SNP VCF.
# Output: PASS, max-missingness, and MAC-filtered VCF.
# Software: GATK VariantFiltration, bcftools, tabix

INPUT_VCF="<path/to/biallelic_snp_input.vcf.gz>"
OUTPUT_DIR="<path/to/highcon_biallelic_output_dir>"
PREFIX="<output_prefix>"
THREADS=32

GQ=30
AD=3
MAX_MISSING=0.50
MAC=2

mkdir -p "$OUTPUT_DIR"

GATK_VCF="${OUTPUT_DIR}/${PREFIX}_GATKfiltering_GQ${GQ}_AD${AD}.vcf.gz"
PASS_VCF="${OUTPUT_DIR}/${PREFIX}_filtered_GQ${GQ}_AD${AD}.vcf.gz"
MISS_VCF="${OUTPUT_DIR}/${PREFIX}_filtered_GQ${GQ}_AD${AD}_max50missing.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/${PREFIX}_filtered_GQ${GQ}_AD${AD}_max50missing_MAC${MAC}.vcf.gz"

gatk VariantFiltration \
  -V "$INPUT_VCF" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  --genotype-filter-expression "GQ < ${GQ}.0" --genotype-filter-name "lowGQ" \
  --genotype-filter-expression "isHet == 1 && (AD[0] < ${AD} || AD[1] < ${AD})" --genotype-filter-name "lowSupport_Het" \
  --genotype-filter-expression "isHomRef == 1 && AD[0] < ${AD}" --genotype-filter-name "lowSupport_HomRef" \
  --genotype-filter-expression "isHomVar == 1 && AD[1] < ${AD}" --genotype-filter-name "lowSupport_HomAlt" \
  --set-filtered-genotype-to-no-call \
  -O "$GATK_VCF"
tabix -p vcf "$GATK_VCF"

bcftools view -f PASS --threads "$THREADS" "$GATK_VCF" -Oz -o "$PASS_VCF"
tabix -p vcf "$PASS_VCF"

bcftools filter -i "F_MISSING <= ${MAX_MISSING}" --threads "$THREADS" "$PASS_VCF" -Oz -o "$MISS_VCF"
tabix -p vcf "$MISS_VCF"

bcftools view -i "MAC>=${MAC}" --threads "$THREADS" "$MISS_VCF" -Oz -o "$FINAL_VCF"
tabix -p vcf "$FINAL_VCF"
