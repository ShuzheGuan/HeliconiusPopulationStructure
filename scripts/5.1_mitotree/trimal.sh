#!/usr/bin/env bash
set -euo pipefail

# Purpose: trim mitochondrial alignment columns with trimAl automated trimming.
# Input: MAFFT-aligned mitochondrial multi-FASTA file.
# Output: trimmed alignment FASTA and HTML trimming summary.
# Software: trimAl

INPUT_ALIGNMENT="<path/to/mito_concat_aligned.fasta>"
OUTPUT_TRIMMED_ALIGNMENT="<path/to/mito_concat_aligned_trimAI_trimmed.fasta>"
OUTPUT_HTML_SUMMARY="<path/to/mito_concat_aligned_trimAI_summary.html>"

mkdir -p "$(dirname "$OUTPUT_TRIMMED_ALIGNMENT")" "$(dirname "$OUTPUT_HTML_SUMMARY")"

[[ -s "$INPUT_ALIGNMENT" ]] || exit 1

trimal \
  -in "$INPUT_ALIGNMENT" \
  -out "$OUTPUT_TRIMMED_ALIGNMENT" \
  -automated1 \
  -htmlout "$OUTPUT_HTML_SUMMARY"
