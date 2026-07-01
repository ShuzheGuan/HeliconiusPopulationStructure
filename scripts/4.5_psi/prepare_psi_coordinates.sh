#!/usr/bin/env bash
set -euo pipefail

# Purpose: build the coordinate/population TSV required by rangeExpansion PSI.
# Input: PLINK .fam file, coordinate TSV, and sample-population TSV.
# Output: coords.tsv with id, longitude, latitude, outgroup, region, pop.
# Software: R, readr, dplyr

PLINK_PREFIX="<path/to/focal_polarized_plink_prefix_without_extension>"
COORDINATES_TSV="<path/to/sample_coordinates.tsv>"
POPULATION_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_COORDS="<path/to/coords.tsv>"

mkdir -p "$(dirname "$OUTPUT_COORDS")"

export PLINK_PREFIX COORDINATES_TSV POPULATION_TSV OUTPUT_COORDS

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

fam <- read_table(
  paste0(Sys.getenv("PLINK_PREFIX"), ".fam"),
  col_names = c("fid", "id", "father", "mother", "sex", "phenotype"),
  show_col_types = FALSE
)

coords <- read_tsv(Sys.getenv("COORDINATES_TSV"), show_col_types = FALSE) |>
  mutate(id = as.character(id))

pops <- read_tsv(
  Sys.getenv("POPULATION_TSV"),
  col_names = c("id", "pop"),
  show_col_types = FALSE
) |>
  mutate(id = as.character(id), pop = as.character(pop))

coords |>
  inner_join(fam |> select(id), by = "id") |>
  inner_join(pops, by = "id") |>
  transmute(
    id,
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    outgroup = 0L,
    region = pop,
    pop
  ) |>
  filter(!is.na(longitude), !is.na(latitude)) |>
  write_tsv(Sys.getenv("OUTPUT_COORDS"))
RSCRIPT
