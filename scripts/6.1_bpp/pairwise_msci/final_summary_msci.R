#!/usr/bin/env Rscript

# Purpose: summarize MSC-I estimates across chains and replicates for plotting.
# Input: msci_estimates.tsv and a pair metadata TSV with parapatry/sympatry labels.
# Output: one wide summary table with point, lower, and upper values per parameter.
# Software: R, readr, dplyr, tidyr, stringr

ESTIMATES_TSV <- "<path/to/msci_estimates.tsv>"
PAIR_METADATA_TSV <- "<path/to/pair_metadata.tsv>"
OUTPUT_TSV <- "<path/to/msci_summary.tsv>"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

dat <- read_tsv(ESTIMATES_TSV, show_col_types = FALSE)
pair_metadata <- read_tsv(PAIR_METADATA_TSV, show_col_types = FALSE) %>%
  select(species, pop1, pop2, parapatry) %>%
  distinct()

mean_cols <- names(dat)[str_detect(names(dat), "_mean$")]

dat %>%
  select(species, pop1, pop2, rep, chain, all_of(mean_cols)) %>%
  pivot_longer(
    cols = all_of(mean_cols),
    names_to = "parameter",
    values_to = "chain_mean"
  ) %>%
  mutate(
    parameter = str_remove(parameter, "_mean$"),
    chain_mean = as.numeric(chain_mean)
  ) %>%
  group_by(species, pop1, pop2, parameter, rep) %>%
  summarise(rep_estimate = mean(chain_mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(species, pop1, pop2, parameter) %>%
  summarise(
    point = median(rep_estimate, na.rm = TRUE),
    lower = min(rep_estimate, na.rm = TRUE),
    upper = max(rep_estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(point, lower, upper),
    names_to = "stat",
    values_to = "value"
  ) %>%
  mutate(colname = paste0(parameter, "_", stat)) %>%
  select(species, pop1, pop2, colname, value) %>%
  pivot_wider(names_from = colname, values_from = value) %>%
  left_join(pair_metadata, by = c("species", "pop1", "pop2")) %>%
  relocate(parapatry, .after = pop2) %>%
  arrange(species, pop1, pop2) %>%
  write_tsv(OUTPUT_TSV, na = "NA")
