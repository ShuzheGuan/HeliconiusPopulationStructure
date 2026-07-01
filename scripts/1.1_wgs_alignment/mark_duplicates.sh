#!/usr/bin/env bash
set -euo pipefail

# Purpose: mark duplicate reads in sorted BAM files.
# Input: sorted BAM files and sample ID list.
# Output: duplicate-marked BAM files and duplicate metrics files.
# Software: GATK MarkDuplicatesSpark

BAM_DIR="<path/to/aligned_sorted_bam_dir>"
OUTPUT_DIR="<path/to/deduplicated_bam_dir>"
ID_LIST="<path/to/sample_id_list.txt>"

mkdir -p "$OUTPUT_DIR"

while read -r SAMPLE_ID; do
  [[ -n "$SAMPLE_ID" ]] || continue

  INPUT_BAM="${BAM_DIR}/${SAMPLE_ID}_sorted_RG.bam"
  OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}_deduplicated.bam"
  METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_ID}_metrics.txt"

  gatk MarkDuplicatesSpark \
    -I "$INPUT_BAM" \
    -O "$OUTPUT_BAM" \
    -M "$METRICS_FILE"
done < "$ID_LIST"
