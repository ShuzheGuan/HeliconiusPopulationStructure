#!/usr/bin/env bash
set -euo pipefail

# Purpose: discover and mask high-confidence repeats, then merge repeat masks with inaccessible genome masks.
# Input: reference genome FASTA and existing inaccessible BED file from depth/soft-clipping filters.
# Output: RepeatModeler library, RepeatMasker annotation, high-confidence repeat BED, and final inaccessible BED file.
# Software: RepeatModeler, RepeatMasker, bedtools

GENOME="<path/to/reference_genome.fasta>"
GENOME_NAME="<reference_genome_label>"
OUT_DIR="<path/to/repeat_filtering_output_dir>"
INACCESSIBLE_BED="<path/to/DP_softclip_inaccessible_windows.bed>"
THREADS=32

REPEATMODELER_DIR="${OUT_DIR}/repeatmodeler_${GENOME_NAME}"
REPEATMASKER_DIR="${REPEATMODELER_DIR}/repeatmasker_output"
KNOWN_REPEAT_LIB="${REPEATMODELER_DIR}/consensi.fa.classified_known"
UNKNOWN_REPEAT_LIB="${REPEATMODELER_DIR}/consensi.fa.classified_unknown"
HIGH_CONFIDENCE_GFF="${OUT_DIR}/${GENOME_NAME}_repeatmasker_high_confidence.gff"
HIGH_CONFIDENCE_BED="${OUT_DIR}/${GENOME_NAME}_repeatmasker_high_confidence.bed"
FINAL_INACCESSIBLE_BED="${OUT_DIR}/${GENOME_NAME}_inaccessible_final.bed"

mkdir -p "$REPEATMODELER_DIR" "$REPEATMASKER_DIR"

ln -sf "$GENOME" "${REPEATMODELER_DIR}/${GENOME_NAME}.fa"

BuildDatabase \
  -name "${REPEATMODELER_DIR}/${GENOME_NAME}_db" \
  "${REPEATMODELER_DIR}/${GENOME_NAME}.fa"

RepeatModeler \
  -database "${REPEATMODELER_DIR}/${GENOME_NAME}_db" \
  -pa "$THREADS" \
  -LTRStruct

CONSENSI="$(find "$REPEATMODELER_DIR" -name "consensi.fa.classified" | head -n 1)"

awk -v known="$KNOWN_REPEAT_LIB" -v unknown="$UNKNOWN_REPEAT_LIB" '
  BEGIN { RS=">"; ORS="" }
  NR > 1 {
    if ($0 ~ /#Unknown/) {
      print ">" $0 >> unknown
    } else {
      print ">" $0 >> known
    }
  }
' "$CONSENSI"

RepeatMasker \
  -pa "$THREADS" \
  -gff \
  -xsmall \
  -dir "$REPEATMASKER_DIR" \
  -lib "$KNOWN_REPEAT_LIB" \
  "$GENOME"

GFF_FILE="$(find "$REPEATMASKER_DIR" -name "*.gff" | head -n 1)"

awk '$3 == "similarity" && $6 <= 10 && ($5 - $4 + 1) >= 100 && $9 !~ /Unknown/' \
  "$GFF_FILE" \
  > "$HIGH_CONFIDENCE_GFF"

awk '{
  match($9, /Motif:([^" ]+)/, arr)
  print $1"\t"($4 - 1)"\t"$5"\t"arr[1]"\t"$6"\t"$7
}' "$HIGH_CONFIDENCE_GFF" > "$HIGH_CONFIDENCE_BED"

cat "$INACCESSIBLE_BED" "$HIGH_CONFIDENCE_BED" | \
  bedtools sort -i - | \
  bedtools merge -i - \
  > "$FINAL_INACCESSIBLE_BED"
