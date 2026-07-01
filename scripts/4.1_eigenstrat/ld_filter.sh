#!/usr/bin/env bash
set -euo pipefail

# Purpose: create LD-controlled PLINK datasets by physical thinning and LD pruning.
# Input: per-site filtered PLINK .bed, .bim, and .fam files.
# Output: PLINK files thinned by 2 kb spacing and by LD pruning.
# Software: PLINK

INPUT_PREFIX="<path/to/per_site_filtered_plink_prefix_without_extension>"
OUTPUT_DIR="<path/to/ld_filtered_plink_dir>"
BASENAME="<output_prefix>"
# Set this to the chromosome count used by the dataset.
CHR_SET="<species_chr_set>"
THREADS=16

BP_SPACE=2000
LD_WINDOW="5kb"
LD_STEP=5
LD_R2=0.2

mkdir -p "$OUTPUT_DIR"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --bp-space "$BP_SPACE" \
  --make-bed \
  --out "${OUTPUT_DIR}/${BASENAME}_2kb_ld" \
  --threads "$THREADS"

TMP_PREFIX="${OUTPUT_DIR}/${BASENAME}_rsquare_ld_tmp"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
  --out "$TMP_PREFIX" \
  --threads "$THREADS"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --extract "${TMP_PREFIX}.prune.in" \
  --make-bed \
  --out "${OUTPUT_DIR}/${BASENAME}_rsquare_ld" \
  --threads "$THREADS"
