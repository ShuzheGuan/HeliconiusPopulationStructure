#!/usr/bin/env bash
set -euo pipefail

# Purpose: calculate windowed pi, FST, and dxy with pixy.
# Input: VCF containing invariant and biallelic sites, plus a two-column sample-population TSV file.
# Output: pixy windowed pi, FST, and dxy result tables.
# Software: pixy

VCF="<path/to/highconfidence_invariant_and_biallelic_sites.vcf.gz>"
POPULATIONS_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_DIR="<path/to/pixy_output_dir>"
OUTPUT_PREFIX="<output_prefix>"
WINDOW_SIZE=50000
THREADS=32

mkdir -p "$OUTPUT_DIR"

[[ -s "$VCF" ]] || exit 1
[[ -s "${VCF}.tbi" || -s "${VCF}.csi" ]] || exit 1
[[ -s "$POPULATIONS_TSV" ]] || exit 1

pixy \
  --vcf "$VCF" \
  --populations "$POPULATIONS_TSV" \
  --window_size "$WINDOW_SIZE" \
  --n_cores "$THREADS" \
  --stats pi fst dxy \
  --output_folder "$OUTPUT_DIR" \
  --output_prefix "$OUTPUT_PREFIX"
