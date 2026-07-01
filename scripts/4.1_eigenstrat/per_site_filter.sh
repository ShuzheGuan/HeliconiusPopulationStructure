#!/usr/bin/env bash
set -euo pipefail

# Purpose: filter a PLINK dataset by per-site missingness, with and without an additional MAF filter.
# Input: PLINK .bed, .bim, and .fam files.
# Output: per-site filtered PLINK files with geno0.10 and geno0.10_maf05 suffixes.
# Software: PLINK

INPUT_PREFIX="<path/to/plink_input_prefix_without_extension>"
OUTPUT_DIR="<path/to/per_site_filtered_plink_dir>"
BASENAME="<output_prefix>"
# Set this to the chromosome count used by the dataset.
CHR_SET="<species_chr_set>"
THREADS=16

GENO_THRESHOLD=0.10
MAF_THRESHOLD=0.05

mkdir -p "$OUTPUT_DIR"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --geno "$GENO_THRESHOLD" \
  --maf "$MAF_THRESHOLD" \
  --make-bed \
  --out "${OUTPUT_DIR}/${BASENAME}_geno0.10_maf05" \
  --threads "$THREADS"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --geno "$GENO_THRESHOLD" \
  --make-bed \
  --out "${OUTPUT_DIR}/${BASENAME}_geno0.10" \
  --threads "$THREADS"
