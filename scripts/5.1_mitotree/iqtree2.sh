#!/usr/bin/env bash
set -euo pipefail

# Purpose: infer a mitochondrial maximum-likelihood tree from the trimmed alignment.
# Input: trimAl-trimmed mitochondrial alignment FASTA.
# Output: IQ-TREE2 output files, including treefile, model, log, and support values.
# Software: IQ-TREE2

INPUT_ALIGNMENT="<path/to/mito_concat_aligned_trimAI_trimmed.fasta>"
OUTPUT_PREFIX="<path/to/iqtree_output_prefix>"
THREADS="AUTO"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

[[ -s "$INPUT_ALIGNMENT" ]] || exit 1

iqtree2 \
  -s "$INPUT_ALIGNMENT" \
  -m MFP \
  -bb 1000 \
  -bnni \
  -alrt 1000 \
  -abayes \
  -nt "$THREADS" \
  -pre "$OUTPUT_PREFIX"
