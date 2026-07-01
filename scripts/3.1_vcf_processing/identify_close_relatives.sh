#!/usr/bin/env bash
set -euo pipefail

# Purpose: identify close relatives from a QC-thinned PLINK dataset.
# Input: PLINK BED/BIM/FAM files.
# Output: KING IBS output files for close-relative assessment.
# Software: KING

PLINK_PREFIX="<path/to/plink_dataset_prefix>"
OUTPUT_PREFIX="<path/to/king_output_prefix>"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

king \
  -b "${PLINK_PREFIX}.bed" \
  --ibs \
  --prefix "$OUTPUT_PREFIX"
