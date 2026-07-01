#!/usr/bin/env bash
set -euo pipefail

# Purpose: test SMC++ estimate settings across window-size (-w) and regularization-penalty values.
# Input: one directory of .smc.gz files and one population list.
# Output: SMC++ estimate output directories for each population, w, and rp combination.
# Software: SMC++, Singularity

SMCPP_CONTAINER="<path/to/smcpp.sif>"
SMCPP_DIR="<path/to/smc_files_dir>"
POPULATIONS_FILE="<path/to/populations.txt>"
OUTPUT_ROOT="<path/to/smcpp_parameter_test_output_dir>"
DATASET_LABEL="<dataset_label>"

MUTATION_RATE="<mutation_rate_per_site_per_generation>"
W_LIST=(20 50 100)
RP_LIST=(6 5 4)
EM_ITERATIONS=50
KNOTS=12
TIMEPOINT_LOWER=10000
TIMEPOINT_UPPER=1000000
CORES=16

mkdir -p "$OUTPUT_ROOT"

while read -r POP; do
  [[ -n "$POP" ]] || continue

  mapfile -t INPUT_FILES < <(
    for FILE in "$SMCPP_DIR"/*.smc.gz; do
      PID="$(zcat -f "$FILE" 2>/dev/null | head -n 1 | sed -n 's/.*"pids":[[:space:]]*\["\([^"]*\)".*/\1/p')"
      [[ "$PID" == "$POP" ]] && printf '%s\n' "$FILE"
    done
  )

  for W in "${W_LIST[@]}"; do
    for RP in "${RP_LIST[@]}"; do
      OUTDIR="${OUTPUT_ROOT}/${DATASET_LABEL}_${POP}_rp${RP}_w${W}"
      mkdir -p "$OUTDIR"

      singularity exec --pwd "$OUTDIR" "$SMCPP_CONTAINER" \
        smc++ estimate \
          --em-iterations "$EM_ITERATIONS" \
          --spline cubic \
          --regularization-penalty "$RP" \
          --timepoints "$TIMEPOINT_LOWER" "$TIMEPOINT_UPPER" \
          --knots "$KNOTS" \
          -w "$W" \
          --cores "$CORES" \
          -o "$OUTDIR" \
          "$MUTATION_RATE" \
          "${INPUT_FILES[@]}"
    done
  done
done < "$POPULATIONS_FILE"
