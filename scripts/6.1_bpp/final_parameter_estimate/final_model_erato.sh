#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/final_model_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: erato setup for final BPP MSC-I parameter-estimate models.

RAW_BPP_ROOT="<path/to/erato_ncds_bpp_root>"
FOCAL_POP_LIST="<path/to/erato_focal_population_list.txt>"
TRIM_MAP_AUTO="<path/to/erato_two_males_per_population.tsv>"
TRIM_MAP_SEX="<path/to/erato_two_males_per_population.tsv>"
IMAP_MAP_AUTO="<path/to/erato_two_males_per_population.tsv>"
IMAP_MAP_SEX="<path/to/erato_two_males_per_population.tsv>"
OUTPUT_ROOT="<path/to/erato_final_model_output_root>"

SEX_CONTIG="Herato2101"
CHR1_CONTIG="Herato0101"
TARGET_RANDOM_LOCI=4000
RANDOM_SEED=20260504

MODEL_POPS=(Yungas Choco_S him ches)
MODEL_PREFIX="Yungas_vs_Choco_S_vs_him_vs_ches"
MODEL_TREE="(((((Yungas,F[&phi=0.200000])E,(Choco_S,C[&phi=0.500000])D)ANC_erato,B[&phi=0.200000])A,(him,E[&phi=0.200000])F)ANC_era_him, ((ches,D[&phi=0.100000])C,A[&phi=0.200000])B)ROOT;"

DATASETS=(sex_chr chr1 random4k)
NCHAINS=3
THREADS=16
RUN_BPP=false

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.03 e"
TAUPRIOR="2.1 0.03"
PHIPRIOR="1 1"
