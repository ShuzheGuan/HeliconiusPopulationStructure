#!/usr/bin/env python3

# Purpose: build per-chromosome genotype caches for introgressed-haplotype analyses.
# Input: PLINK files and sample lists.
# Output: one compressed .npz genotype cache per chromosome.
# Software: Python, pandas, PLINK

from pathlib import Path

import pandas as pd

from allele_sharing_cache import build_chrom_cache, cache_path
from allele_sharing_lib import load_ids, plink_file, write_keep


PLINK_PREFIX = Path("<path/to/plink_prefix_without_extension>")
SAMPLE_LISTS = [
    Path("<path/to/test_population_samples.list>"),
    Path("<path/to/source_population_samples.list>"),
    Path("<path/to/reference_population_samples.list>"),
]
CACHE_DIR = Path("<path/to/per_chromosome_npz_cache_dir>")
CHROMOSOMES = None


def main():
    bim_path = plink_file(PLINK_PREFIX, ".bim")
    fam = pd.read_csv(
        plink_file(PLINK_PREFIX, ".fam"),
        sep=r"\s+",
        header=None,
        usecols=[1],
        names=["iid"],
        dtype={"iid": str},
    )
    fam_iids = set(fam["iid"])
    samples = set().union(*(load_ids(path) for path in SAMPLE_LISTS)) & fam_iids

    chroms = CHROMOSOMES or (
        pd.read_csv(bim_path, sep=r"\s+", header=None, usecols=[0])[0]
        .astype(str)
        .unique()
        .tolist()
    )

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    keep = CACHE_DIR / "_samples.keep"
    write_keep(samples, keep)

    for chrom in chroms:
        out = cache_path(CACHE_DIR, str(chrom))
        build_chrom_cache(PLINK_PREFIX, keep, str(chrom), out)


if __name__ == "__main__":
    main()
