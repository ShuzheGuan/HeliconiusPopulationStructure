#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/4tip_topology_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: melpomene setup for four-tip BPP A01 topology runs.

NCDS_BPP_ROOT="<path/to/melpomene_ncds_bpp_root>"
CDS_BPP_ROOT="<path/to/melpomene_cds_bpp_root>"
AUTO_SAMPLE_POP_TSV="<path/to/melpomene_sample_population_labels.tsv>"
SEX_SAMPLE_POP_TSV="<path/to/melpomene_two_males_per_population.tsv>"
POPULATION_LIST="<path/to/melpomene_four_tip_population_list.txt>"
OUTPUT_DIR="<path/to/melpomene_four_tip_topology_output_dir>"

DATASET_LABEL="melpomene"
SEX_CONTIG="Hmel221001o"

BLOCK_SIZE=100
N_CHAINS=2
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=100000

THETAPRIOR="2.1 0.03"
TAUPRIOR="2.1 0.03"

SPECIES_LIST="numata Yungas Western_Ecuador tim_Yungas cyd_Western_Ecuador"
STARTING_TREE="(numata,((Yungas,Western_Ecuador),(tim_Yungas,cyd_Western_Ecuador)));"
