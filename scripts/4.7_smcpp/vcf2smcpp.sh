#!/usr/bin/env bash
set -euo pipefail

# Purpose: convert VCF data to SMC++ .smc.gz files for each population, contig, and distinguished individual.
# Input: one VCF, one sample-population TSV, one inaccessible-genome mask BED, and a population-contig job list.
# Output: population/individual/contig-level .smc.gz files.
# Software: SMC++, Singularity

JOB_LIST="<path/to/vcf2smcpp_job.list>"   # population<TAB>contig
VCF="<path/to/invariant_biallelic.vcf.gz>"
POPULATION_TSV="<path/to/sample_population_labels.tsv>"
MASK_BED="<path/to/inaccessible.bed.gz>"
OUTPUT_DIR="<path/to/smcpp_raw_data_output_dir>"
SMCPP_CONTAINER="<path/to/smcpp.sif>"
CORES=8

mkdir -p "$OUTPUT_DIR"

while read -r POP CONTIG; do
  [[ -n "${POP:-}" && -n "${CONTIG:-}" ]] || continue
  [[ "$POP" == "population" ]] && continue

  mapfile -t INDIVIDUALS < <(awk -F '\t' -v pop="$POP" '$2 == pop {print $1}' "$POPULATION_TSV")
  INDIVIDUAL_CSV="$(IFS=,; echo "${INDIVIDUALS[*]}")"

  for INDIVIDUAL in "${INDIVIDUALS[@]}"; do
    singularity exec "$SMCPP_CONTAINER" \
      smc++ vcf2smc \
        --cores "$CORES" \
        -m "$MASK_BED" \
        "$VCF" \
        "${OUTPUT_DIR}/${POP}_${INDIVIDUAL}_${CONTIG}.smc.gz" \
        "$CONTIG" \
        "${POP}:${INDIVIDUAL_CSV}" \
        -d "$INDIVIDUAL" "$INDIVIDUAL"
  done
done < "$JOB_LIST"
