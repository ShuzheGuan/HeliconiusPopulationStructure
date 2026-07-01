#!/usr/bin/env bash
set -euo pipefail

# Purpose: jointly genotype GenomicsDB workspaces into VCF files.
# Input: GenomicsDB workspaces, reference genome FASTA, and interval list.
# Output: joint-genotyped VCF files including invariant sites.
# Software: GATK GenotypeGVCFs

WORKSPACE_DIR="<path/to/genomicsdb_workspace_dir>"
OUTPUT_DIR="<path/to/joint_vcf_dir>"
INTERVALS_LIST="<path/to/genome_intervals.list>"
REFERENCE="<path/to/reference_genome.fasta>"

mkdir -p "$OUTPUT_DIR"

while read -r INTERVAL; do
  [[ -n "$INTERVAL" ]] || continue

  INTERVAL_WORKSPACE="${WORKSPACE_DIR}/workspace_${INTERVAL}"
  OUTPUT_VCF="${OUTPUT_DIR}/${INTERVAL}.vcf.gz"

  gatk GenotypeGVCFs \
    -R "$REFERENCE" \
    -V "gendb://${INTERVAL_WORKSPACE}" \
    -O "$OUTPUT_VCF" \
    --include-non-variant-sites
done < "$INTERVALS_LIST"
