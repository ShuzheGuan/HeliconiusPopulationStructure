#!/usr/bin/env bash
set -euo pipefail

# Purpose: build a full NJ tree and LOWO support from PLINK IBS distances.
# Input: PLINK .bed/.bim/.fam files.
# Output: full/LOWO distance matrices and NJ tree with support labels.
# Software: PLINK, R, ape

PLINK_PREFIX="<path/to/plink_prefix_without_extension>"
OUTPUT_DIR="<path/to/nj_tree_lowo_output_dir>"
OUTPUT_PREFIX="<dataset_label>"
CHR_SET="<species_chr_set>"
THREADS=16
LOWO_BLOCKS=20

DIST_DIR="${OUTPUT_DIR}/distance"
LOWO_DIR="${OUTPUT_DIR}/lowo"
WINDOW_DIR="${OUTPUT_DIR}/windows"
TREE_OUT="${OUTPUT_DIR}/${OUTPUT_PREFIX}_lowo_support.nj.tre"

mkdir -p "$DIST_DIR" "$LOWO_DIR" "$WINDOW_DIR"

FULL_DIST="${DIST_DIR}/${OUTPUT_PREFIX}_full"

plink --bfile "$PLINK_PREFIX" --distance square 1-ibs \
  --out "$FULL_DIST" --allow-extra-chr --chr-set "$CHR_SET" --threads "$THREADS"

TOTAL_SNPS="$(wc -l < "${PLINK_PREFIX}.bim")"
WINDOW_SIZE="$(( (TOTAL_SNPS + LOWO_BLOCKS - 1) / LOWO_BLOCKS ))"

awk -v n="$WINDOW_SIZE" -v b="$LOWO_BLOCKS" -v out="$WINDOW_DIR" '
  {w = int((NR - 1) / n) + 1; if (w > b) w = b; print $2 >> out "/window_" sprintf("%02d", w) ".snp"}
' "${PLINK_PREFIX}.bim"

cut -f2 "${PLINK_PREFIX}.bim" > "${WINDOW_DIR}/all_snps.list"

for BLOCK in $(seq -w 1 "$LOWO_BLOCKS"); do
  EXCLUDE="${WINDOW_DIR}/window_${BLOCK}.snp"
  KEEP="${WINDOW_DIR}/window_${BLOCK}.keep"
  LOWO_DIST="${LOWO_DIR}/${OUTPUT_PREFIX}_LOWO_${BLOCK}"

  grep -vxF -f "$EXCLUDE" "${WINDOW_DIR}/all_snps.list" > "$KEEP"

  plink --bfile "$PLINK_PREFIX" --extract "$KEEP" --distance square 1-ibs \
    --out "$LOWO_DIST" --allow-extra-chr --chr-set "$CHR_SET" --threads "$THREADS"
done

export FULL_DIST LOWO_DIR TREE_OUT OUTPUT_PREFIX LOWO_BLOCKS

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages(library(ape))

read_mdist <- function(prefix) {
  ids <- read.table(paste0(prefix, ".mdist.id"), stringsAsFactors = FALSE)$V2
  mat <- as.matrix(read.table(paste0(prefix, ".mdist"), check.names = FALSE))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- colnames(mat) <- ids
  mat
}

full_tree <- ape::nj(as.dist(read_mdist(Sys.getenv("FULL_DIST"))))
lowo_trees <- list()

for (i in seq_len(as.integer(Sys.getenv("LOWO_BLOCKS")))) {
  prefix <- file.path(Sys.getenv("LOWO_DIR"), sprintf("%s_LOWO_%02d", Sys.getenv("OUTPUT_PREFIX"), i))
  d <- read_mdist(prefix)
  d <- d[full_tree$tip.label, full_tree$tip.label]
  lowo_trees[[i]] <- ape::nj(as.dist(d))
}

full_tree$node.label <- round(100 * ape::prop.clades(full_tree, lowo_trees) / length(lowo_trees), 1)
ape::write.tree(full_tree, file = Sys.getenv("TREE_OUT"))
RSCRIPT
