#!/usr/bin/env bash
set -euo pipefail

# Purpose: import per-sample GVCFs into GenomicsDB workspaces by genomic interval.
# Input: sample-name map linking sample IDs to GVCF files, and interval list.
# Output: GenomicsDB workspaces for each interval.
# Software: GATK GenomicsDBImport

MAP_LIST="<path/to/sample_name_map.txt>"
INTERVALS_LIST="<path/to/genome_intervals.list>"
WORKSPACE_DIR="<path/to/genomicsdb_workspace_dir>"
READER_THREADS=4

mkdir -p "$WORKSPACE_DIR"

while read -r INTERVAL; do
  [[ -n "$INTERVAL" ]] || continue

  INTERVAL_WORKSPACE="${WORKSPACE_DIR}/workspace_${INTERVAL}"

  gatk GenomicsDBImport \
    --genomicsdb-workspace-path "$INTERVAL_WORKSPACE" \
    --sample-name-map "$MAP_LIST" \
    --intervals "$INTERVAL" \
    --reader-threads "$READER_THREADS"
done < "$INTERVALS_LIST"
