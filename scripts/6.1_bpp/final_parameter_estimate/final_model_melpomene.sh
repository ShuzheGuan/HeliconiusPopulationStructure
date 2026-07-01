#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/final_model_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: melpomene setup for final BPP MSC-I parameter-estimate models.

RAW_BPP_ROOT="<path/to/melpomene_ncds_bpp_root>"
FOCAL_POP_LIST="<path/to/melpomene_focal_population_list.txt>"
TRIM_MAP_AUTO="<path/to/melpomene_two_males_per_population.tsv>"
TRIM_MAP_SEX="<path/to/melpomene_two_males_per_population.tsv>"
IMAP_MAP_AUTO="<path/to/melpomene_autosome_imap_sample_population_labels.tsv>"
IMAP_MAP_SEX="<path/to/melpomene_two_males_per_population.tsv>"
OUTPUT_ROOT="<path/to/melpomene_final_model_output_root>"

SEX_CONTIG="Hmel221001o"
CHR1_CONTIG="Hmel201001o"
TARGET_RANDOM_LOCI=4000
RANDOM_SEED=20260504

MODEL_POPS=(cyd_Western_Ecuador tim_Yungas Western_Ecuador Yungas)
MODEL_PREFIX="cyd_Western_Ecuador_vs_tim_Yungas_vs_Western_Ecuador_vs_Yungas"
MODEL_TREE="((((cyd_Western_Ecuador,C[&phi=0.200000])D,(tim_Yungas,F[&phi=0.200000])E)ANC_CTH,A[&phi=0.200000])B, (((Western_Ecuador,D[&phi=0.200000])C,B[&phi=0.200000])A,(Yungas,E[&phi=0.200000])F)ANC_MELP)ROOT;"

DATASETS=(sex_chr chr1 random4k)
NCHAINS=10
THREADS=16
RUN_BPP=false

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.03 e"
TAUPRIOR="2.1 0.03"
PHIPRIOR="1 1"
