#!/usr/bin/env bash
set -euo pipefail

# Purpose: create a high-confidence indel callset.
# Input: indel/complex VCF and reference genome.
# Output: normalized, PASS, max-missingness, and MAC-filtered indel VCF.
# Software: bcftools, GATK VariantFiltration, tabix

INPUT_VCF="<path/to/indel_or_complex_input.vcf.gz>"
REFERENCE="<path/to/reference_genome.fasta>"
OUTPUT_DIR="<path/to/highcon_indel_output_dir>"
PREFIX="<output_prefix>"
THREADS=32

GQ=30
AD=3
MAX_MISSING=0.50
MAC=5

mkdir -p "$OUTPUT_DIR"

NORM_VCF="${OUTPUT_DIR}/${PREFIX}_norm.vcf.gz"
INDEL_VCF="${OUTPUT_DIR}/${PREFIX}_norm_indel.vcf.gz"
GATK_VCF="${OUTPUT_DIR}/${PREFIX}_norm_indel_GATKfiltering_GQ${GQ}_AD${AD}.vcf.gz"
PASS_VCF="${OUTPUT_DIR}/${PREFIX}_norm_indel_filtered_GQ${GQ}_AD${AD}.vcf.gz"
MISS_VCF="${OUTPUT_DIR}/${PREFIX}_norm_indel_filtered_GQ${GQ}_AD${AD}_max50missing.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/${PREFIX}_norm_indel_filtered_GQ${GQ}_AD${AD}_max50missing_MAC${MAC}.vcf.gz"

bcftools norm -m -both -f "$REFERENCE" "$INPUT_VCF" -Oz -o "$NORM_VCF"
tabix -p vcf "$NORM_VCF"

bcftools view -i 'TYPE~"indel"' --threads "$THREADS" "$NORM_VCF" -Oz -o "$INDEL_VCF"
tabix -p vcf "$INDEL_VCF"

gatk VariantFiltration \
  -V "$INDEL_VCF" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
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
