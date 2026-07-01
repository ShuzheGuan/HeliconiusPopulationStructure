#!/usr/bin/env bash
set -euo pipefail

# Purpose: convert an LD-filtered binary PLINK dataset to EIGENSTRAT format.
# Input: LD-filtered PLINK .bed, .bim, and .fam files.
# Output: EIGENSTRAT .geno, .snp, and .ind files.
# Software: PLINK, EIGENSOFT convertf

INPUT_PREFIX="<path/to/ld_filtered_plink_prefix_without_extension>"
OUTPUT_DIR="<path/to/eigenstrat_output_dir>"
SCRATCH_DIR="<path/to/scratch_dir_for_ped_map_and_convertf_files>"
BASENAME="<output_prefix>"
# Set this to the chromosome count used by the dataset.
CHR_SET="<species_chr_set>"
THREADS=8

mkdir -p "$OUTPUT_DIR" "$SCRATCH_DIR"

PED_PREFIX="${SCRATCH_DIR}/${BASENAME}"
EIGEN_PREFIX="${SCRATCH_DIR}/${BASENAME}"
PARFILE="${SCRATCH_DIR}/${BASENAME}_plink_to_eigenstrat.par"

plink \
  --bfile "$INPUT_PREFIX" \
  --chr-set "$CHR_SET" \
  --recode \
  --out "$PED_PREFIX" \
  --threads "$THREADS"

cat > "$PARFILE" <<EOF
genotypename:    ${PED_PREFIX}.ped
snpname:         ${PED_PREFIX}.map
indivname:       ${PED_PREFIX}.ped
outputformat:    EIGENSTRAT
genotypeoutname: ${EIGEN_PREFIX}.geno
snpoutname:      ${EIGEN_PREFIX}.snp
indivoutname:    ${EIGEN_PREFIX}.ind
familynames:     NO
EOF

convertf -p "$PARFILE"

mv -f "${EIGEN_PREFIX}.geno" "${OUTPUT_DIR}/${BASENAME}.geno"
mv -f "${EIGEN_PREFIX}.snp" "${OUTPUT_DIR}/${BASENAME}.snp"
mv -f "${EIGEN_PREFIX}.ind" "${OUTPUT_DIR}/${BASENAME}.ind"
