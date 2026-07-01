#!/usr/bin/env bash
set -euo pipefail

# Purpose: build a custom BLAST database from reference mitochondrial genomes.
# Input: reference mitochondrial FASTA files.
# Output: concatenated reference FASTA and nucleotide BLAST database files.
# Software: BLAST+

REF_DIR="<path/to/reference_mitochondrial_fastas>"
OUT_DIR="<path/to/blast_database_dir>"
REFS_FASTA="${OUT_DIR}/mito_refs.fasta"
DB_PREFIX="${OUT_DIR}/mito_refs_db"

mkdir -p "$OUT_DIR"

find "$OUT_DIR" -maxdepth 1 -type f \( -name "mito_refs_db.*" -o -name "mito_refs.fasta" \) -delete

shopt -s nullglob
REF_FILES=( "$REF_DIR"/*.fa "$REF_DIR"/*.fasta )

cat "${REF_FILES[@]}" > "$REFS_FASTA"

makeblastdb \
  -in "$REFS_FASTA" \
  -dbtype nucl \
  -out "$DB_PREFIX"
