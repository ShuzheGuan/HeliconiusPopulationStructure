#!/usr/bin/env bash
set -euo pipefail

# Purpose: run standard PCA with SMARTPCA on an EIGENSTRAT dataset.
# Input: EIGENSTRAT .geno, .snp, and .ind files.
# Output: SMARTPCA .evec, .eval, .par, and .log files.
# Software: EIGENSOFT smartpca

EIGENSTRAT_PREFIX="<path/to/eigenstrat_prefix_without_extension>"
OUTPUT_DIR="<path/to/smartpca_output_dir>"
BASENAME="<output_prefix>"
# Set this to the chromosome count used by the dataset.
NUMCHROM="<species_numchrom>"

mkdir -p "$OUTPUT_DIR"

GENO_FILE="${EIGENSTRAT_PREFIX}.geno"
SNP_FILE="${EIGENSTRAT_PREFIX}.snp"
IND_FILE="${EIGENSTRAT_PREFIX}.ind"
EVEC_FILE="${OUTPUT_DIR}/${BASENAME}.pca.evec"
EVAL_FILE="${OUTPUT_DIR}/${BASENAME}.pca.eval"
PARFILE="${OUTPUT_DIR}/${BASENAME}.smartpca.par"
LOG_FILE="${OUTPUT_DIR}/${BASENAME}.smartpca.log"

cat > "$PARFILE" <<EOF
genotypename:    ${GENO_FILE}
snpname:         ${SNP_FILE}
indivname:       ${IND_FILE}
evecoutname:     ${EVEC_FILE}
evaloutname:     ${EVAL_FILE}
altnormstyle:    NO
numoutevec:      10
numoutlieriter:  0
numoutlierevec:  2
outliersigmathresh: 6.0
numchrom:        ${NUMCHROM}
blgsize:         0.005
EOF

smartpca -p "$PARFILE" > "$LOG_FILE"
