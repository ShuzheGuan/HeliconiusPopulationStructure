#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/bigtree_topology_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: erato setup for BPP A01 big-tree topology runs.
# Run once for nCDS and once for CDS by changing RAW_BPP_ROOT and OUTPUT_DIR.

RAW_BPP_ROOT="<path/to/erato_cds_or_ncds_bpp_root>"
SAMPLE_POP_TSV="<path/to/erato_sample_population_labels.tsv>"
SAMPLE_SEX_TSV="<path/to/erato_sample_sex.tsv>"
POPULATION_LIST="<path/to/erato_bigtree_population_list.txt>"
OUTPUT_DIR="<path/to/erato_bigtree_topology_output_dir>"

EXCLUDE_DIR=""
MALE_LABEL="m"
MALES_PER_POP=2
BLOCK_SIZE=100
MAX_MISSING=0.20
N_CHAINS=5
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=100000

THETAPRIOR="2.1 0.03"
TAUPRIOR="2.1 0.03"

SPECIES_LIST="sara clysonymus telesiphe herm him ches Magdalena Pacific_Panama Atlantic Yungas Mosquito Western_Ecuador"
STARTING_TREE="(sara,((clysonymus,telesiphe),(herm,((him,ches),(((Magdalena,Pacific_Panama),(Atlantic,Yungas)),(Mosquito,Western_Ecuador))))));"
