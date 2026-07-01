#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/pairwise_msci_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: erato setup for pairwise BPP MSC-I runs.
# Edit POP1 and POP2 for each pairwise comparison.

RAW_BPP_ROOT="<path/to/erato_ncds_bpp_root>"
SAMPLE_POP_TSV="<path/to/erato_sample_population_labels.tsv>"
OUTPUT_ROOT="<path/to/erato_pairwise_msci_output_root>"

POP1="<erato_population_1>"
POP2="<erato_population_2>"
EXCLUDE_DIR="Herato2101"

TARGET_LOCI=2000
MAX_MISSING=0.10
N_REPS=5
RANDOM_SEED_BASE=20260504
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.015 e"
TAUPRIOR="2.1 0.06"
PHIPRIOR="1 1"
