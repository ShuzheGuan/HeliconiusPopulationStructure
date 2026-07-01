#!/usr/bin/env bash
set -euo pipefail

# Purpose: trim paired-end raw reads with fastp.
# Input: paired-end raw FASTQ files.
# Output: trimmed paired-end FASTQ files and fastp quality-control reports.
# Software: fastp

RAW_DIR="<path/to/raw_reads_dir>"
TRIMMED_DIR="<path/to/trimmed_reads_dir>"
LOG_DIR="<path/to/fastp_log_dir>"

mkdir -p "$TRIMMED_DIR" "$LOG_DIR"

find "$RAW_DIR" -maxdepth 1 -type f -name "*1.fq.gz" | sort | while read -r READ1; do
  base="$(basename "$READ1")"
  prefix="${base%1.fq.gz}"
  RUN_ID="${prefix%[_.]}"

  READ2="${RAW_DIR}/${prefix}2.fq.gz"
  [[ -f "$READ2" ]] || continue

  fastp \
    -i "$READ1" \
    -I "$READ2" \
    -o "${TRIMMED_DIR}/${RUN_ID}_1_trimmed.fastq.gz" \
    -O "${TRIMMED_DIR}/${RUN_ID}_2_trimmed.fastq.gz" \
    -j "${LOG_DIR}/${RUN_ID}_fastp.json" \
    -h "${LOG_DIR}/${RUN_ID}_fastp.html"
done
