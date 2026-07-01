#!/usr/bin/env python3

# Purpose: build a sample-by-SNP matrix of CTH-like allele dosages.
# Input: CTH-like site table, per-chromosome genotype caches, and sample lists.
# Output: introgression genotype matrix.
# Software: Python, pandas, numpy

from pathlib import Path
import csv

import numpy as np
import pandas as pd

from allele_sharing_cache import cache_path, load_chrom_cache
from allele_sharing_lib import load_ids


CTH_LIKE_SITE_TSV = Path("<path/to/highcon_cth_like_sites.tsv>")
CACHE_DIR = Path("<path/to/per_chromosome_npz_cache_dir>")
SAMPLE_LISTS = [
    Path("<path/to/test_population_samples.list>"),
    Path("<path/to/source_population_samples.list>"),
    Path("<path/to/reference_population_samples.list>"),
]
OUTPUT_TSV = Path("<path/to/introgression_genotype_matrix.tsv>")


def cth_like_dosage(alt_copy, cth_like_allele):
    if np.isnan(alt_copy):
        return "NA"
    alt = int(alt_copy)
    if cth_like_allele == 0:
        return str(2 - alt)
    if cth_like_allele == 1:
        return str(alt)
    return "NA"


def column_label(site):
    if "segment_id" in site._fields:
        return f"{site.segment_id}|{site.chrom}:{int(site.pos)}"
    return f"{site.chrom}:{int(site.pos)}"


def main():
    sites = pd.read_csv(CTH_LIKE_SITE_TSV, sep="\t")
    sites["chrom"] = sites["chrom"].astype(str)
    sites["pos"] = sites["pos"].astype(np.int64)
    sites["cth_like_allele"] = sites["cth_like_allele"].astype(int)
    sites = sites.sort_values(["chrom", "pos"]).reset_index(drop=True)

    samples = sorted(set().union(*(load_ids(path) for path in SAMPLE_LISTS)))
    rows = [[sample] for sample in samples]
    cache = {}

    for site in sites.itertuples(index=False):
        chrom = str(site.chrom)
        pos = int(site.pos)
        allele = int(site.cth_like_allele)

        if chrom not in cache:
            pos_arr, geno, cache_iids = load_chrom_cache(cache_path(CACHE_DIR, chrom))
            cache[chrom] = (pos_arr, geno, {sample: i for i, sample in enumerate(cache_iids)})

        pos_arr, geno, iid_index = cache[chrom]
        idx = np.searchsorted(pos_arr, pos)
        found = idx < pos_arr.size and pos_arr[idx] == pos

        for row, sample in zip(rows, samples):
            if not found or sample not in iid_index:
                row.append("NA")
            else:
                row.append(cth_like_dosage(float(geno[iid_index[sample], idx]), allele))

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="\n") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample"] + [column_label(site) for site in sites.itertuples(index=False)])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
