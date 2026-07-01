#!/usr/bin/env python3

# Purpose: sample random null windows and calculate the same allele-sharing statistic.
# Input: PLINK files, per-chromosome .npz genotype cache, excluded BED intervals, and sample lists.
# Output: one BED of sampled null windows and one TSV with T statistics.
# Software: Python, pandas, numpy

from pathlib import Path
import random

import numpy as np
import pandas as pd

from allele_sharing_cache import (
    MIN_NULL_SNPS,
    allele_sharing_from_cache_window,
    cache_path,
    count_snps_in_window,
    load_chrom_cache,
    load_chrom_lengths,
    load_excluded_intervals,
    sample_null_window,
)
from allele_sharing_lib import load_bim_table, load_ids, plink_file, summary_columns


PLINK_PREFIX = Path("<path/to/plink_prefix_without_extension>")
CACHE_DIR = Path("<path/to/per_chromosome_npz_cache_dir>")
EXCLUDE_BEDS = [Path("<path/to/candidate_segments_to_exclude.bed>")]
TEST_SAMPLE_LIST = Path("<path/to/test_population_samples.list>")
CTH_SAMPLE_LIST = Path("<path/to/CTH_samples.list>")
AMAZONE_SAMPLE_LIST = Path("<path/to/AmazonE_samples.list>")

OUTPUT_BED = Path("<path/to/null_windows.bed>")
OUTPUT_TSV = Path("<path/to/null_allele_sharing_summary.tsv>")
RECIPIENT_LABEL = "<recipient_label>"
WINDOW_BP = 50_000
N_NULL_WINDOWS = 1000
MAX_ATTEMPTS = N_NULL_WINDOWS * 1000
SEED = 42


def main():
    rng = random.Random(SEED)
    bim = load_bim_table(plink_file(PLINK_PREFIX, ".bim"))
    pos_by_chrom = {c: grp["pos"].to_numpy(dtype=np.int64) for c, grp in bim.groupby("chrom")}
    chrom_lengths = load_chrom_lengths(bim)
    excluded = load_excluded_intervals(EXCLUDE_BEDS)

    test = load_ids(TEST_SAMPLE_LIST)
    cth = load_ids(CTH_SAMPLE_LIST)
    amazone = load_ids(AMAZONE_SAMPLE_LIST)
    cols = summary_columns(RECIPIENT_LABEL)
    cache = {}

    OUTPUT_BED.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_BED.open("w") as bed_out, OUTPUT_TSV.open("w") as tsv_out:
        tsv_out.write("\t".join(cols) + "\n")
        placed = 0
        attempts = 0

        while placed < N_NULL_WINDOWS:
            attempts += 1
            if attempts > MAX_ATTEMPTS:
                raise SystemExit(f"only placed {placed} null windows after {attempts} attempts")

            hit = sample_null_window(chrom_lengths, excluded, WINDOW_BP, rng)
            if hit is None:
                continue

            chrom, start, end = hit
            if count_snps_in_window(pos_by_chrom, chrom, start, end) < MIN_NULL_SNPS:
                continue

            if chrom not in cache:
                cache[chrom] = load_chrom_cache(cache_path(CACHE_DIR, chrom))

            pos, geno, iids = cache[chrom]
            n_snps, n_compared, n_only, t_stat = allele_sharing_from_cache_window(
                pos, geno, iids, start, end, test, cth, amazone
            )

            segment = f"null_{placed + 1:06d}"
            bed_out.write(f"{chrom}\t{start}\t{end}\t{segment}\t0\n")
            row = [
                segment,
                chrom,
                str(start),
                str(end),
                str(WINDOW_BP),
                str(n_snps),
                str(n_only),
                str(n_compared),
                "" if pd.isna(t_stat) else f"{t_stat:.6f}",
            ]
            tsv_out.write("\t".join(row) + "\n")
            placed += 1


if __name__ == "__main__":
    main()
