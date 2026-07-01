#!/usr/bin/env bash
set -euo pipefail

# Purpose: convert a PLINK dataset to population-stratified allele frequencies for TreeMix.
# Input: PLINK .bed/.bim/.fam files and a two-column sample-population TSV.
# Output: filtered PLINK files and a .frq.strat file.
# Software: PLINK, awk

PLINK_PREFIX="<path/to/plink_prefix_without_extension>"
SAMPLE_POP_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_DIR="<path/to/plink_frequency_output_dir>"
OUTPUT_PREFIX="<dataset_label>"
CHR_SET="<species_chr_set>"

MAC=2
GENO=0.05

mkdir -p "$OUTPUT_DIR"

KEEP_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.keep"
WITHIN_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.within"
FILTERED_PREFIX="${OUTPUT_DIR}/${OUTPUT_PREFIX}_filtered"
FREQ_PREFIX="${OUTPUT_DIR}/${OUTPUT_PREFIX}_freq"

awk 'BEGIN {OFS="\t"} NR > 1 && $1 != "" && $2 != "" {print $1, $1}' \
  "$SAMPLE_POP_TSV" > "$KEEP_FILE"

awk 'BEGIN {OFS="\t"} NR > 1 && $1 != "" && $2 != "" {print $1, $1, $2}' \
  "$SAMPLE_POP_TSV" > "$WITHIN_FILE"

plink \
  --bfile "$PLINK_PREFIX" \
  --keep "$KEEP_FILE" \
  --mac "$MAC" \
  --geno "$GENO" \
  --chr-set "$CHR_SET" \
  --make-bed \
  --out "$FILTERED_PREFIX"

plink \
  --bfile "$FILTERED_PREFIX" \
  --freq \
  --within "$WITHIN_FILE" \
  --chr-set "$CHR_SET" \
  --out "$FREQ_PREFIX"
