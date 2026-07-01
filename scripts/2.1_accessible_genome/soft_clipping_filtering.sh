#!/usr/bin/env bash
set -euo pipefail

# Purpose: define soft-clipping-based accessible genome windows and intersect them with depth-accessible windows.
# Input: genome windows, duplicate-marked BAM files, sample ID list, and depth-accessible BED file.
# Output: per-sample soft-clipping summaries, consensus soft-clipping BED files, and depth/soft-clipping intersection BED file.
# Software: samtools, bedtools, Python/pandas

WINDOWS_BED="<path/to/genome_windows_1kb.bed>"
BAM_DIR="<path/to/deduplicated_bam_dir>"
SAMPLE_LIST="<path/to/sample_id_list.txt>"
DP_ACCESSIBLE_BED="<path/to/consensus_accessible_DP_windows.bed>"
OUT_DIR="<path/to/softclip_filtering_output_dir>"

MIN_SOFTCLIP_BP=20
MAX_SOFTCLIP_PROPORTION=0.2
CONSENSUS_PROPORTION=0.80

PER_SAMPLE_RAW_DIR="${OUT_DIR}/per_sample_softclip_counts"
PER_SAMPLE_SUMMARY_DIR="${OUT_DIR}/per_sample_softclip_summary"
CONSENSUS_BED="${OUT_DIR}/consensus_window_summary_softclip.bed"
ACCESSIBLE_BED="${OUT_DIR}/consensus_accessible_softclip_windows.bed"
INACCESSIBLE_BED="${OUT_DIR}/consensus_inaccessible_softclip_windows.bed"
INTERSECTION_BED="${OUT_DIR}/consensus_accessible_DP_softclip_intersection.bed"

mkdir -p "$PER_SAMPLE_RAW_DIR" "$PER_SAMPLE_SUMMARY_DIR"

TOTAL_SAMPLES="$(awk 'NF {n++} END {print n + 0}' "$SAMPLE_LIST")"
CONSENSUS_THRESHOLD="$(awk -v n="$TOTAL_SAMPLES" -v p="$CONSENSUS_PROPORTION" 'BEGIN {print int(n * p + 0.999999)}')"

while read -r SAMPLE; do
  [[ -n "$SAMPLE" ]] || continue

  BAM="${BAM_DIR}/${SAMPLE}_deduplicated.bam"
  TOTAL_READS="${PER_SAMPLE_RAW_DIR}/${SAMPLE}_total_reads.bed"
  SOFTCLIP_BAM="${PER_SAMPLE_RAW_DIR}/${SAMPLE}_softclip${MIN_SOFTCLIP_BP}.bam"
  SOFTCLIP_READS="${PER_SAMPLE_RAW_DIR}/${SAMPLE}_softclip_reads.bed"
  SUMMARY="${PER_SAMPLE_SUMMARY_DIR}/${SAMPLE}_softclip_summary.csv"

  bedtools coverage \
    -a "$WINDOWS_BED" \
    -b "$BAM" \
    -counts \
    > "$TOTAL_READS"

  samtools view -h "$BAM" | \
  awk -v min_bp="$MIN_SOFTCLIP_BP" '{
    if ($0 ~ /^@/) { print; next }
    cigar = $6
    while (match(cigar, /([0-9]+)S/, m)) {
      if (m[1] >= min_bp) {
        print $0
        break
      }
      cigar = substr(cigar, RSTART + RLENGTH)
    }
  }' | samtools view -bS - > "$SOFTCLIP_BAM"

  bedtools coverage \
    -a "$WINDOWS_BED" \
    -b "$SOFTCLIP_BAM" \
    -counts \
    > "$SOFTCLIP_READS"

  paste "$TOTAL_READS" "$SOFTCLIP_READS" | \
  awk -v max_prop="$MAX_SOFTCLIP_PROPORTION" -v OFS="," '
    BEGIN {print "chr,start,end,total_reads,softclip_reads,prop,flag"}
    {
      total = $4
      soft = $8
      prop = (total > 0) ? soft / total : 0
      flag = (prop >= max_prop) ? 1 : 0
      print $1, $2, $3, total, soft, prop, flag
    }
  ' > "$SUMMARY"
done < "$SAMPLE_LIST"

export WINDOWS_BED PER_SAMPLE_SUMMARY_DIR CONSENSUS_BED CONSENSUS_THRESHOLD
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

summary_files = sorted(glob.glob(os.path.join(os.environ["PER_SAMPLE_SUMMARY_DIR"], "*_softclip_summary.csv")))
flag_columns = []
for path in summary_files:
    col = pd.read_csv(path, usecols=[6])
    col = col.applymap(lambda x: 1 if x == 0 else 0)
    flag_columns.append(col)

flags = pd.concat(flag_columns, axis=1)
threshold = int(os.environ["CONSENSUS_THRESHOLD"])

windows["access_count"] = flags.sum(axis=1).astype(int)
windows["flag"] = windows["access_count"].apply(lambda x: "in" if x >= threshold else "out")
windows[["chr", "start", "end", "access_count", "flag"]].to_csv(
    os.environ["CONSENSUS_BED"],
    sep="\t",
    header=False,
    index=False,
)
PY

awk '$5 == "in"  {print $1, $2, $3}' OFS="\t" "$CONSENSUS_BED" > "$ACCESSIBLE_BED"
awk '$5 == "out" {print $1, $2, $3}' OFS="\t" "$CONSENSUS_BED" > "$INACCESSIBLE_BED"

bedtools intersect \
  -a "$ACCESSIBLE_BED" \
  -b "$DP_ACCESSIBLE_BED" \
  > "$INTERSECTION_BED"
