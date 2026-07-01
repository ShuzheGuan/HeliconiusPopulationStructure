#!/usr/bin/env bash
set -euo pipefail

# Purpose: run SMARTPCA with projection using a population list.
# Input: EIGENSTRAT .geno, .snp, and .ind files plus a poplist file defining the non-projected populations.
# Output: projected SMARTPCA .evec, .eval, .par, and .log files.
# Software: EIGENSOFT smartpca

EIGENSTRAT_PREFIX="<path/to/eigenstrat_prefix_without_extension>"
POPLIST_FILE="<path/to/poplist.txt>"
OUTPUT_DIR="<path/to/projected_smartpca_output_dir>"
BASENAME="<output_prefix>"
PROJECT_TAG="<projection_label>"
# Set this to the chromosome count used by the dataset.
NUMCHROM="<species_numchrom>"

mkdir -p "$OUTPUT_DIR"

GENO_FILE="${EIGENSTRAT_PREFIX}.geno"
SNP_FILE="${EIGENSTRAT_PREFIX}.snp"
IND_FILE="${EIGENSTRAT_PREFIX}.ind"
EVEC_FILE="${OUTPUT_DIR}/${BASENAME}.${PROJECT_TAG}.pca.evec"
EVAL_FILE="${OUTPUT_DIR}/${BASENAME}.${PROJECT_TAG}.pca.eval"
PARFILE="${OUTPUT_DIR}/${BASENAME}.${PROJECT_TAG}.smartpca.par"
LOG_FILE="${OUTPUT_DIR}/${BASENAME}.${PROJECT_TAG}.smartpca.log"

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
qtmode:          0
poplistname:     ${POPLIST_FILE}
lsqproject:      YES
numchrom:        ${NUMCHROM}
blgsize:         0.005
EOF

smartpca -p "$PARFILE" > "$LOG_FILE"
