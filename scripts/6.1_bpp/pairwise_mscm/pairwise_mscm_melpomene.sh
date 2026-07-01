#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/pairwise_mscm_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: melpomene setup for pairwise BPP MSC-M runs.
# Edit POP1 and POP2 for each pairwise comparison.

RAW_BPP_ROOT="<path/to/melpomene_ncds_bpp_root>"
SAMPLE_POP_TSV="<path/to/melpomene_sample_population_labels.tsv>"
OUTPUT_ROOT="<path/to/melpomene_pairwise_mscm_output_root>"

POP1="<melpomene_population_1>"
POP2="<melpomene_population_2>"
EXCLUDE_DIR="Hmel221001o"

TARGET_LOCI=2000
N_REPS=5
RANDOM_SEED_BASE=20260504
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.03 e"
TAUPRIOR="2.1 0.03"
WPRIOR="2 0.003"
