#!/usr/bin/env bash
set -euo pipefail

# Purpose: query assembled mitochondrial FASTA files from one dataset against the custom BLAST database.
# Input: assembled mitochondrial FASTA files and a nucleotide BLAST database.
# Output: per-sample BLAST result tables and a dataset-level BLAST summary table.
# Software: BLAST+

FASTA_DIR="<path/to/mitochondrial_fastas>"
DATASET_LABEL="<dataset_label>"
DB_PREFIX="<path/to/blast_database_dir>/mito_refs_db"
OUT_DIR="<path/to/blast_output_dir>"

mkdir -p "$OUT_DIR"

SUMMARY="${OUT_DIR}/${DATASET_LABEL}_blast_summary.tsv"
printf 'Sample\tRank\tTopHit\tPercentIdentity\tAlignLength\n' > "$SUMMARY"

find "$FASTA_DIR" -maxdepth 1 -type f -name "*.fasta" | sort | while read -r SAMPLE; do
  PREFIX="$(basename "$SAMPLE" .fasta)"
  OUT_TSV="${OUT_DIR}/${PREFIX}_vs_mito_refs.tsv"

  blastn \
    -query "$SAMPLE" \
    -db "$DB_PREFIX" \
    -outfmt 6 \
    -max_target_seqs 10 \
    -evalue 1e-3 \
    > "$OUT_TSV"

  if [[ -s "$OUT_TSV" ]]; then
    awk -v prefix="$PREFIX" 'NR <= 2 {printf "%s\t%d\t%s\t%s%%\t%s\n", prefix, NR, $2, $3, $4}' "$OUT_TSV" >> "$SUMMARY"
  else
    printf '%s\t1\tNoHit\tNA\tNA\n' "$PREFIX" >> "$SUMMARY"
    printf '%s\t2\tNoHit\tNA\tNA\n' "$PREFIX" >> "$SUMMARY"
  fi
done
