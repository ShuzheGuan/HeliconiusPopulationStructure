#!/usr/bin/env bash
set -euo pipefail

# Purpose: concatenate mitochondrial FASTA files and align them with MAFFT.
# Input: directory of mitochondrial FASTA files.
# Output: aligned mitochondrial multi-FASTA file.
# Software: MAFFT

INPUT_FASTA_DIR="<path/to/mitochondrial_fasta_dir>"
OUTPUT_ALIGNMENT="<path/to/mito_concat_aligned.fasta>"
THREADS=32

OUTPUT_DIR="$(dirname "$OUTPUT_ALIGNMENT")"
mkdir -p "$OUTPUT_DIR"

[[ -d "$INPUT_FASTA_DIR" ]] || exit 1

shopt -s nullglob
FASTA_FILES=( "$INPUT_FASTA_DIR"/*.fasta )
[[ "${#FASTA_FILES[@]}" -gt 0 ]] || exit 1

TMP_CONCAT="$(mktemp "${OUTPUT_DIR}/mito_concat_XXXXXX.fasta")"
trap 'rm -f "$TMP_CONCAT"' EXIT

cat "${FASTA_FILES[@]}" > "$TMP_CONCAT"

mafft \
  --maxiterate 1000 \
  --localpair \
  --adjustdirectionaccurately \
  --thread "$THREADS" \
  "$TMP_CONCAT" \
> "$OUTPUT_ALIGNMENT"
