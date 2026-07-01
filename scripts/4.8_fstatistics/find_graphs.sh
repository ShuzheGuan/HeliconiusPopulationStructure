#!/usr/bin/env bash
set -euo pipefail

# Purpose: search candidate qpGraph models from fixed starting topologies.
# Input: precomputed F2 block directories for erato and melp datasets.
# Output: per-run best-model TSVs, full result RDS files, and top-score summaries.
# Software: R, admixtools, igraph, readr, dplyr, tibble

ERATO_F2_DIR="<path/to/erato_f2_blocks_dir>"
MELP_F2_DIR="<path/to/melp_f2_blocks_dir>"
OUTPUT_DIR="<path/to/find_graphs_output_dir>"

N_SEEDS=20
NUMGRAPHS=100
STOP_GEN2=50
STOP_GEN=300
MAX_ADMIX=8
PLUSMINUS=3

mkdir -p "$OUTPUT_DIR"
TOP_SCORE_DIR="${OUTPUT_DIR}/top_scores"
mkdir -p "$TOP_SCORE_DIR"

RUN=0

for PASS in 1 2; do
  if [[ "$PASS" == 1 ]]; then
    REJECT_F4Z="none"
    Z_LABEL="nozreject"
  else
    REJECT_F4Z="4"
    Z_LABEL="zreject"
  fi

  for SPECIES in erato melp; do
    if [[ "$SPECIES" == "erato" ]]; then
      F2_DIR="$ERATO_F2_DIR"
      OUTPOP="melpomene"
    else
      F2_DIR="$MELP_F2_DIR"
      OUTPOP="erato"
    fi

    for OPT_WORST in TRUE FALSE; do
      if [[ "$OPT_WORST" == "TRUE" ]]; then
        TAG_LABEL="zscore"
      else
        TAG_LABEL="totalscore"
      fi

      for SEED in $(seq 1 "$N_SEEDS"); do
        RUN=$((RUN + 1))
        RUN_ID="$(printf "%03d" "$RUN")"
        SEED_ID="$(printf "%02d" "$SEED")"
        BASE="${SPECIES}_${TAG_LABEL}_${Z_LABEL}_${RUN_ID}_${SEED_ID}"

        export SPECIES F2_DIR OUTPOP OUTPUT_DIR SEED NUMGRAPHS STOP_GEN2 STOP_GEN
        export MAX_ADMIX PLUSMINUS OPT_WORST REJECT_F4Z PASS RUN TAG_LABEL BASE

        Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages(library(admixtools))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

species <- Sys.getenv("SPECIES")
f2_dir <- Sys.getenv("F2_DIR")
outpop <- Sys.getenv("OUTPOP")
out_dir <- Sys.getenv("OUTPUT_DIR")
seed <- as.integer(Sys.getenv("SEED"))
numgraphs <- as.integer(Sys.getenv("NUMGRAPHS"))
stop_gen2 <- as.integer(Sys.getenv("STOP_GEN2"))
stop_gen <- as.integer(Sys.getenv("STOP_GEN"))
max_admix <- as.integer(Sys.getenv("MAX_ADMIX"))
plusminus <- as.integer(Sys.getenv("PLUSMINUS"))
opt_worst <- as.logical(Sys.getenv("OPT_WORST"))
reject_f4z <- Sys.getenv("REJECT_F4Z")
pass <- as.integer(Sys.getenv("PASS"))
run <- as.integer(Sys.getenv("RUN"))
tag_label <- Sys.getenv("TAG_LABEL")
base <- Sys.getenv("BASE")

melp_edges <- rbind(
  c("root", "erato"), c("root", "r1"),
  c("r1", "doris"), c("r1", "r1b"),
  c("r1b", "isnu"), c("r1b", "e2"),
  c("isnu", "ismenius"), c("isnu", "numata"),
  c("e2", "cy_th"),
  c("cy_th", "cydno_Magdalena"), c("cy_th", "th"),
  c("th", "tim_Napo"), c("th", "heu"),
  c("e2", "aroot"),
  c("aroot", "awnene"), c("aroot", "Western_Ecuador"),
  c("awnene", "Yungas"), c("awnene", "ane"),
  c("ane", "Guiana"), c("ane", "Atlantic")
)

erato_edges <- rbind(
  c("root", "melpomene"), c("root", "s1"),
  c("s1", "sara"), c("s1", "s2"),
  c("s2", "ct"), c("ct", "clysonymus"), c("ct", "telesiphe"),
  c("s2", "o"), c("o", "herm"), c("o", "o2"),
  c("o2", "hc"), c("hc", "him"), c("hc", "ches"),
  c("o2", "y"), c("y", "amazon"), c("y", "z"),
  c("amazon", "Atlantic"), c("amazon", "Yungas"),
  c("z", "Mosquito"), c("z", "Western_Ecuador")
)

edges <- if (species == "erato") erato_edges else melp_edges
g0 <- igraph::graph_from_edgelist(as.matrix(edges), directed = TRUE)
ed <- f2_from_precomp(f2_dir)

set.seed(seed)

args <- list(
  data = ed,
  initgraph = g0,
  outpop = outpop,
  numgraphs = numgraphs,
  stop_gen2 = stop_gen2,
  stop_gen = stop_gen,
  max_admix = max_admix,
  opt_worst_residual = opt_worst,
  plusminus_generations = plusminus,
  verbose = TRUE
)

if (reject_f4z != "none") {
  res <- do.call(find_graphs, c(args, list(reject_f4z = as.numeric(reject_f4z))))
} else {
  res <- do.call(find_graphs, args)
}

tbl <- tibble::as_tibble(res)
best <- dplyr::slice_min(tbl, order_by = score, n = 1, with_ties = FALSE)

numadmix <- if ("numadmix" %in% names(best)) best$numadmix[1] else NA_real_
worstZ <- if ("worstZ" %in% names(best)) best$worstZ[1] else NA_real_

summary <- tibble::tibble(
  species = species,
  pass = pass,
  run = run,
  seed = seed,
  tag = tag_label,
  reject_f4z = ifelse(reject_f4z == "none", NA_real_, as.numeric(reject_f4z)),
  models_tested = nrow(tbl),
  score = best$score[1],
  numadmix = numadmix,
  worstZ = worstZ
)

readr::write_tsv(summary, file.path(out_dir, paste0(base, ".tsv")))
saveRDS(tbl, file.path(out_dir, paste0(base, ".rds")))
RSCRIPT
      done
    done
  done
done

shopt -s nullglob

for BATCH in erato_zscore erato_totalscore melp_zscore melp_totalscore; do
  FILES=( "$OUTPUT_DIR"/${BATCH}*.tsv )
  [[ "${#FILES[@]}" -gt 0 ]] || continue

  {
    head -n 1 "${FILES[0]}"
    for FILE in "${FILES[@]}"; do
      tail -n +2 "$FILE"
    done
  } | sort -t$'\t' -k8,8g | head -n 6 > "${TOP_SCORE_DIR}/${BATCH}_top_scores.tsv"
done
