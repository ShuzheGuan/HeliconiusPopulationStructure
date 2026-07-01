#!/usr/bin/env Rscript

# Purpose: summarize final BPP parameter estimates across chains for plotting.
# Input: final_model_estimates.tsv from summarize_final_model_estimates.py.
# Output: one table with median [HPD low, HPD high] per parameter and dataset.
# Software: R, readr, dplyr, tidyr, stringr

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

ESTIMATES_TSV <- "<path/to/final_model_estimates.tsv>"
OUTPUT_TSV <- "<path/to/final_model_summary.tsv>"

dataset_map <- c(
  "sex_chr" = "sex",
  "chr1" = "chr1",
  "random4k" = "random4k"
)
dataset_order <- c("sex", "chr1", "random4k")

dat <- read_tsv(ESTIMATES_TSV, show_col_types = FALSE)

stat_cols <- names(dat)[
  str_detect(names(dat), "_median$|_2\\.5%HPD$|_97\\.5%HPD$")
]

dat %>%
  mutate(dataset = recode(data_set, !!!dataset_map)) %>%
  select(dataset, chain, all_of(stat_cols)) %>%
  pivot_longer(
    cols = all_of(stat_cols),
    names_to = "param_stat",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    stat = case_when(
      str_detect(param_stat, "_median$") ~ "median",
      str_detect(param_stat, "_2\\.5%HPD$") ~ "hpd_low",
      str_detect(param_stat, "_97\\.5%HPD$") ~ "hpd_high"
    ),
    param = param_stat %>%
      str_remove("_median$") %>%
      str_remove("_2\\.5%HPD$") %>%
      str_remove("_97\\.5%HPD$")
  ) %>%
  filter(!is.na(stat)) %>%
  group_by(dataset, param, stat) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    cell = sprintf("%.6g [%.6g, %.6g]", median, hpd_low, hpd_high),
    dataset = factor(dataset, levels = dataset_order)
  ) %>%
  select(param, dataset, cell) %>%
  pivot_wider(names_from = dataset, values_from = cell) %>%
  select(param, all_of(dataset_order)) %>%
  arrange(param) %>%
  write_tsv(OUTPUT_TSV, na = "NA")
