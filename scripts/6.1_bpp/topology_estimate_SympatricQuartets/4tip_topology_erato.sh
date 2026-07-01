#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/4tip_topology_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: erato setup for four-tip BPP A01 topology runs.

NCDS_BPP_ROOT="<path/to/erato_ncds_bpp_root>"
CDS_BPP_ROOT="<path/to/erato_cds_bpp_root>"
AUTO_SAMPLE_POP_TSV="<path/to/erato_sample_population_labels.tsv>"
SEX_SAMPLE_POP_TSV="<path/to/erato_two_males_per_population.tsv>"
POPULATION_LIST="<path/to/erato_four_tip_population_list.txt>"
OUTPUT_DIR="<path/to/erato_four_tip_topology_output_dir>"

DATASET_LABEL="erato"
SEX_CONTIG="Herato2101"

BLOCK_SIZE=100
N_CHAINS=2
THREADS=16

BURNIN=20000
SAMPFREQ=10
NSAMPLE=20000

THETAPRIOR="2.1 0.03"
TAUPRIOR="2.1 0.03"

SPECIES_LIST="herm him ches Choco_S Yungas"
STARTING_TREE="(herm,((ches,him),(Choco_S,Yungas)));"
