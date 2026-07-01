#!/usr/bin/env bash
set -euo pipefail

# Purpose: define depth-based accessible genome windows across samples.
# Input: reference genome index, duplicate-marked BAM files, and sample ID list.
# Output: per-sample depth summaries, per-window accessibility flags, and consensus accessible/inaccessible BED files.
# Software: bedtools, samtools, mosdepth, Python/pandas

REFERENCE_FAI="<path/to/reference_genome.fasta.fai>"
BAM_DIR="<path/to/deduplicated_bam_dir>"
SAMPLE_LIST="<path/to/sample_id_list.txt>"
OUT_DIR="<path/to/depth_filtering_output_dir>"

WINDOW_SIZE=200
CONSENSUS_PROPORTION=0.50

WINDOWS_BED="${OUT_DIR}/genome_windows_${WINDOW_SIZE}bp.bed"
AVERAGE_DEPTH_TSV="${OUT_DIR}/sample_average_depth.tsv"
PER_SAMPLE_DEPTH_DIR="${OUT_DIR}/per_sample_depth"
PER_SAMPLE_FLAG_DIR="${OUT_DIR}/per_sample_window_accessibility"
CONSENSUS_BED="${OUT_DIR}/consensus_window_summary_DP.bed"
ACCESSIBLE_BED="${OUT_DIR}/consensus_accessible_DP_windows.bed"
INACCESSIBLE_BED="${OUT_DIR}/consensus_inaccessible_DP_windows.bed"

mkdir -p "$OUT_DIR" "$PER_SAMPLE_DEPTH_DIR" "$PER_SAMPLE_FLAG_DIR"

TOTAL_SAMPLES="$(awk 'NF {n++} END {print n + 0}' "$SAMPLE_LIST")"
CONSENSUS_THRESHOLD="$(awk -v n="$TOTAL_SAMPLES" -v p="$CONSENSUS_PROPORTION" 'BEGIN {print int(n * p + 0.999999)}')"

bedtools makewindows \
  -g "$REFERENCE_FAI" \
  -w "$WINDOW_SIZE" \
  > "$WINDOWS_BED"

printf 'sample_id\taverage_depth\n' > "$AVERAGE_DEPTH_TSV"

while read -r SAMPLE; do
  [[ -n "$SAMPLE" ]] || continue

  BAM="${BAM_DIR}/${SAMPLE}_deduplicated.bam"
  AVG_DEPTH="$(samtools depth -a "$BAM" | awk '{sum += $3} END {if (NR > 0) print sum / NR; else print 0}')"

  printf '%s\t%s\n' "$SAMPLE" "$AVG_DEPTH" >> "$AVERAGE_DEPTH_TSV"

  mosdepth \
    --by "$WINDOWS_BED" \
    -t 4 \
    -n \
    -x "${PER_SAMPLE_DEPTH_DIR}/${SAMPLE}" \
    "$BAM"
done < "$SAMPLE_LIST"

while read -r SAMPLE; do
  [[ -n "$SAMPLE" ]] || continue

  REGIONS_FILE="${PER_SAMPLE_DEPTH_DIR}/${SAMPLE}.regions.bed.gz"
  OUTPUT_FILE="${PER_SAMPLE_FLAG_DIR}/${SAMPLE}_window_accessibility.tsv"

  AVG_DEPTH="$(awk -v s="$SAMPLE" '$1 == s {print $2}' "$AVERAGE_DEPTH_TSV")"
  MIN_DEPTH="$(awk -v d="$AVG_DEPTH" 'BEGIN {printf "%.5f", d * 0.5}')"
  MAX_DEPTH="$(awk -v d="$AVG_DEPTH" 'BEGIN {printf "%.5f", d * 2.0}')"

  zcat "$REGIONS_FILE" | \
  awk -v min="$MIN_DEPTH" -v max="$MAX_DEPTH" -v OFS="\t" '{
    flag = ($4 >= min && $4 <= max) ? 1 : 0;
    print $1, $2, $3, flag;
  }' > "$OUTPUT_FILE"
done < "$SAMPLE_LIST"

export WINDOWS_BED PER_SAMPLE_FLAG_DIR CONSENSUS_BED CONSENSUS_THRESHOLD
python3 <<'PY'
import glob
import os
import pandas as pd

windows = pd.read_csv(
    os.environ["WINDOWS_BED"],
    sep="\t",
    header=None,
    names=["chr", "start", "end"],
)

flag_files = sorted(glob.glob(os.path.join(os.environ["PER_SAMPLE_FLAG_DIR"], "*_window_accessibility.tsv")))
flag_columns = [pd.read_csv(f, sep="\t", header=None, usecols=[3]) for f in flag_files]

flags = pd.concat(flag_columns, axis=1)
threshold = int(os.environ["CONSENSUS_THRESHOLD"])

windows["access_count"] = flags.sum(axis=1).astype(int)
windows["flag"] = windows["access_count"].apply(lambda x: "in" if x >= threshold else "out")
windows.to_csv(os.environ["CONSENSUS_BED"], sep="\t", header=False, index=False)
PY

awk '$5 == "in"  {print $1, $2, $3}' OFS="\t" "$CONSENSUS_BED" > "$ACCESSIBLE_BED"
awk '$5 == "out" {print $1, $2, $3}' OFS="\t" "$CONSENSUS_BED" > "$INACCESSIBLE_BED"
