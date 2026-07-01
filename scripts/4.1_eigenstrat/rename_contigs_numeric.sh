#!/usr/bin/env bash
set -euo pipefail

# Purpose: replace contig names in a PLINK BIM file with numeric chromosome IDs and rename SNP IDs.
# Input: PLINK .bed, .bim, and .fam files with original contig names.
# Output: PLINK .bed, .bim, and .fam files with numeric chromosome IDs in the BIM file.
# Software: awk

INPUT_PREFIX="<path/to/plink_input_prefix_without_extension>"
OUTPUT_PREFIX="<path/to/plink_output_prefix_without_extension>"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

cp "${INPUT_PREFIX}.bed" "${OUTPUT_PREFIX}.bed"
cp "${INPUT_PREFIX}.fam" "${OUTPUT_PREFIX}.fam"

awk -v OFS="\t" '
BEGIN {
  CONTIG_COUNTER = 0
}
{
  if (!($1 in CONTIG_MAP)) {
    CONTIG_COUNTER++
    CONTIG_MAP[$1] = CONTIG_COUNTER
  }

  ORIGINAL_CONTIG = $1
  $1 = CONTIG_MAP[$1]
  $2 = "snp_" ORIGINAL_CONTIG "_" $4
  print
}
' "${INPUT_PREFIX}.bim" > "${OUTPUT_PREFIX}.bim"
