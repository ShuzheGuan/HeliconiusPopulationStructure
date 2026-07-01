#!/usr/bin/env bash
set -euo pipefail

# Purpose: run TreeMix across a range of migration-edge values.
# Input: gzipped TreeMix frequency file.
# Output: TreeMix result files for each migration-edge value.
# Software: TreeMix

INPUT_FRQ_GZ="<path/to/input_treemix.frq.gz>"
OUTPUT_DIR="<path/to/treemix_output_dir>"
OUTPUT_PREFIX="<dataset_label>"
ROOT_POP="<root_population_label>"

M_START=0
M_END=10
K=50

mkdir -p "$OUTPUT_DIR"

for M in $(seq "$M_START" "$M_END"); do
  treemix \
    -i "$INPUT_FRQ_GZ" \
    -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_m${M}" \
    -root "$ROOT_POP" \
    -m "$M" \
    -k "$K" \
    -global \
    -se \
    -bootstrap
done
