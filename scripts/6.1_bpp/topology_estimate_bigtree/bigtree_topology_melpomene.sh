#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  exec bash "$SCRIPT_DIR/bigtree_topology_run.sh" "${BASH_SOURCE[0]}"
fi

# Purpose: melpomene setup for BPP A01 big-tree topology runs.
# Run once for nCDS and once for CDS by changing RAW_BPP_ROOT and OUTPUT_DIR.

RAW_BPP_ROOT="<path/to/melpomene_cds_or_ncds_bpp_root>"
SAMPLE_POP_TSV="<path/to/melpomene_sample_population_labels.tsv>"
SAMPLE_SEX_TSV="<path/to/melpomene_sample_sex.tsv>"
POPULATION_LIST="<path/to/melpomene_bigtree_population_list.txt>"
OUTPUT_DIR="<path/to/melpomene_bigtree_topology_output_dir>"

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

SPECIES_LIST="doris ismenius numata cyd_Magdalena tim_Napo heu Yungas Guiana Atlantic Magdalena Western_Ecuador Pacific_Panama"
STARTING_TREE="(doris,((ismenius,numata),((cyd_Magdalena,(tim_Napo,heu)),((Yungas,(Guiana,Atlantic)),(Magdalena,(Western_Ecuador,Pacific_Panama))))));"
