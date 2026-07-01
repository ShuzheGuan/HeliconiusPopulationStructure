#!/usr/bin/env bash
set -euo pipefail

# Purpose: run ADMIXTURE with cross-validation across K values for one dataset.
# Input: LD-filtered binary PLINK .bed, .bim, and .fam files.
# Output: ADMIXTURE .Q, .P, and per-K log files with CV errors.
# Software: ADMIXTURE

OUTPUT_DIR="<path/to/admixture_output_dir>"
DATASET_LABEL="<dataset_label>"
# Use a full path here because ADMIXTURE writes .Q and .P files to OUTPUT_DIR.
PLINK_PREFIX="<full/path/to/ld_filtered_plink_prefix_without_extension>"
THREADS=16
K_MIN=2
K_MAX=15

mkdir -p "$OUTPUT_DIR"

INPUT_BED="${PLINK_PREFIX}.bed"
INPUT_BIM="${PLINK_PREFIX}.bim"
INPUT_FAM="${PLINK_PREFIX}.fam"

if [[ ! -s "$INPUT_BED" || ! -s "$INPUT_BIM" || ! -s "$INPUT_FAM" ]]; then
  exit 1
fi

for K in $(seq "$K_MIN" "$K_MAX"); do
  (
    cd "$OUTPUT_DIR"

    admixture \
      --cv \
      "$INPUT_BED" \
      "$K" \
      -j"$THREADS" \
      > "${DATASET_LABEL}_admixture_k${K}.log"
  )
done
