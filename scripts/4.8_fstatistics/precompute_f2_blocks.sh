#!/usr/bin/env bash
set -euo pipefail

# Purpose: precompute F2 blocks from an EIGENSTRAT dataset.
# Input: EIGENSTRAT .geno/.snp/.ind files and a one-column population list.
# Output: directory of precomputed F2 blocks.
# Software: R, admixtools

EIGENSTRAT_PREFIX="<path/to/eigenstrat_prefix_without_extension>"
POPULATION_LIST="<path/to/population_list.txt>"
OUTPUT_DIR="<path/to/f2_blocks_output_dir>"
BLGSIZE=5e5

mkdir -p "$OUTPUT_DIR"

export EIGENSTRAT_PREFIX POPULATION_LIST OUTPUT_DIR BLGSIZE

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages(library(admixtools))

data_prefix <- Sys.getenv("EIGENSTRAT_PREFIX")
pop_file <- Sys.getenv("POPULATION_LIST")
out_dir <- Sys.getenv("OUTPUT_DIR")
blgsize <- as.numeric(Sys.getenv("BLGSIZE"))

pops <- readLines(pop_file, warn = FALSE)
pops <- pops[nzchar(pops)]

extract_f2(
  data_prefix,
  out_dir,
  pops = pops,
  auto_only = FALSE,
  blgsize = blgsize
)
RSCRIPT
