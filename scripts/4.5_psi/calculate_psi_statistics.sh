#!/usr/bin/env bash
set -euo pipefail

# Purpose: calculate pairwise PSI statistics from polarized PLINK data and population coordinates.
# Input: filtered polarized PLINK .bed file and rangeExpansion coordinate TSV.
# Output: PSI matrix CSV ordered by population order in the coordinate TSV.
# Software: R, rangeExpansion, readr, tibble

PLINK_BED="<path/to/filtered_polarized_plink.bed>"
COORDS_TSV="<path/to/coords.tsv>"
OUTPUT_PSI_CSV="<path/to/psi_matrix.csv>"
PLOIDY=2

[[ -s "$PLINK_BED" ]] || exit 1
[[ -s "$COORDS_TSV" ]] || exit 1

OUTPUT_DIR="$(dirname "$OUTPUT_PSI_CSV")"
mkdir -p "$OUTPUT_DIR"

export PLINK_BED COORDS_TSV OUTPUT_PSI_CSV PLOIDY

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(rangeExpansion)
  library(readr)
  library(tibble)
})

plink_bed <- Sys.getenv("PLINK_BED")
coords_tsv <- Sys.getenv("COORDS_TSV")
output_psi_csv <- Sys.getenv("OUTPUT_PSI_CSV")
ploidy <- as.integer(Sys.getenv("PLOIDY"))

raw <- rangeExpansion::load.data(
  plink.file = plink_bed,
  coords.file = coords_tsv,
  ploidy = ploidy,
  sep = "\t"
)

pop_obj <- rangeExpansion::make.pop(raw, ploidy = ploidy)
psi <- rangeExpansion::get.all.psi(pop_obj)
pop_order <- unique(as.character(pop_obj$coords$pop))

if (is.null(rownames(psi)) || is.null(colnames(psi))) {
  rownames(psi) <- pop_order
  colnames(psi) <- pop_order
}

psi <- psi[pop_order, pop_order, drop = FALSE]
diag(psi) <- 0

readr::write_csv(
  tibble::rownames_to_column(as.data.frame(psi), var = "pop"),
  output_psi_csv
)
RSCRIPT
