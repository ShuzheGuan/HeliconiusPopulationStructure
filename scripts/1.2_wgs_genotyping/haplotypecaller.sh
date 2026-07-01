#!/usr/bin/env bash
set -euo pipefail

# Purpose: call per-sample GVCFs from duplicate-marked BAM files.
# Input: duplicate-marked BAM files, reference genome FASTA, and sample ID list.
# Output: per-sample GVCF files.
# Software: GATK HaplotypeCaller

BAM_DIR="<path/to/deduplicated_bam_dir>"
GVCF_DIR="<path/to/gvcf_dir>"
ID_LIST="<path/to/sample_id_list.txt>"
REFERENCE="<path/to/reference_genome.fasta>"
THREADS=8

mkdir -p "$GVCF_DIR"

samtools faidx "$REFERENCE"
gatk CreateSequenceDictionary \
  -R "$REFERENCE" \
  -O "${REFERENCE%.*}.dict"

while read -r SAMPLE_ID; do
  [[ -n "$SAMPLE_ID" ]] || continue

  INPUT_BAM="${BAM_DIR}/${SAMPLE_ID}_deduplicated.bam"
  OUTPUT_GVCF="${GVCF_DIR}/${SAMPLE_ID}.g.vcf.gz"

  gatk HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$INPUT_BAM" \
    -O "$OUTPUT_GVCF" \
    --emit-ref-confidence GVCF \
    --native-pair-hmm-threads "$THREADS" \
    --output-mode EMIT_ALL_CONFIDENT_SITES
done < "$ID_LIST"
