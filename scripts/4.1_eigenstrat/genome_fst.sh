#!/usr/bin/env bash
set -euo pipefail

# Purpose: estimate genome-wide pairwise FST among population labels using SMARTPCA.
# Input: EIGENSTRAT .geno, .snp, and .ind files with population labels in the third .ind column.
# Output: SMARTPCA FST parameter file and log file.
# Software: EIGENSOFT smartpca

EIGENSTRAT_PREFIX="<path/to/eigenstrat_prefix_without_extension>"
OUTPUT_DIR="<path/to/genome_fst_output_dir>"
BASENAME="<output_prefix>"
# Set this to the chromosome count used by the dataset.
NUMCHROM="<species_numchrom>"

mkdir -p "$OUTPUT_DIR"

GENO_FILE="${EIGENSTRAT_PREFIX}.geno"
SNP_FILE="${EIGENSTRAT_PREFIX}.snp"
IND_FILE="${EIGENSTRAT_PREFIX}.ind"
PARFILE="${OUTPUT_DIR}/${BASENAME}.fst.par"
LOG_FILE="${OUTPUT_DIR}/${BASENAME}.fst.log"

cat > "$PARFILE" <<EOF
genotypename:    ${GENO_FILE}
snpname:         ${SNP_FILE}
indivname:       ${IND_FILE}
fstonly:         YES
numchrom:        ${NUMCHROM}
blgsize:         0.005
EOF

smartpca -p "$PARFILE" > "$LOG_FILE"
