#!/usr/bin/env bash
set -euo pipefail

# Purpose: add population labels to one EIGENSTRAT .ind file using a sample-label table.
# Input: EIGENSTRAT .ind file and a two-column sample-label TSV file.
# Output: labeled .ind file with the third column replaced by population labels.
# Software: awk

INPUT_IND="<path/to/input.ind>"
LABEL_TSV="<path/to/sample_label.tsv>"
OUTPUT_IND="<path/to/output_with_poplabels.ind>"

mkdir -p "$(dirname "$OUTPUT_IND")"

awk -F'\t' -v OFS=' ' '
  FNR == NR {
    gsub(/\r/, "", $1)
    gsub(/\r/, "", $2)
    LABEL_BY_ID[$1] = $2
    next
  }
  {
    gsub(/\r/, "", $0)
    split($0, FIELD, /[[:space:]]+/)
    SAMPLE_ID = FIELD[1]
    SEX = FIELD[2]
    LABEL = (SAMPLE_ID in LABEL_BY_ID) ? LABEL_BY_ID[SAMPLE_ID] : "unknown"
    if (LABEL == "") LABEL = "unknown"
    print SAMPLE_ID, SEX, LABEL
  }
' "$LABEL_TSV" "$INPUT_IND" > "$OUTPUT_IND"
