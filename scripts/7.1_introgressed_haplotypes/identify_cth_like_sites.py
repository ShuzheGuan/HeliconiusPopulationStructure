#!/usr/bin/env python3

# Purpose: identify SNPs in high-confidence introgressed intervals where CTH and reference panels differ strongly.
# Input: high-confidence introgressed intervals, PLINK BIM, per-chromosome genotype caches, and CTH/reference sample lists.
# Output: site table with the CTH-like allele labelled as REF or ALT.
# Software: Python, pandas, numpy

from pathlib import Path

import numpy as np
import pandas as pd

from allele_sharing_cache import cache_path, load_chrom_cache
from allele_sharing_lib import load_bim_table, load_ids, plink_file


PLINK_PREFIX = Path("<path/to/plink_prefix_without_extension>")
CACHE_DIR = Path("<path/to/per_chromosome_npz_cache_dir>")
HIGHCONF_INTERVALS = Path("<path/to/highcon_introgressed_haplotype.tsv>")
CTH_SAMPLE_LIST = Path("<path/to/CTH_samples.list>")
REFERENCE_SAMPLE_LIST = Path("<path/to/reference_population_samples.list>")
OUTPUT_TSV = Path("<path/to/highcon_cth_like_sites.tsv>")

AF_DIFF_MIN = 0.80
MIN_CALLED_PER_PANEL = 5


def load_intervals(path):
    table = pd.read_csv(path, sep="\t")
    rename = {"contig": "chrom", "seg_id": "segment_id", "name": "segment_id"}
    table = table.rename(columns=rename)
    if {"chrom", "start", "end", "segment_id"}.issubset(table.columns):
        out = table[["chrom", "start", "end", "segment_id"]].copy()
    else:
        out = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "segment_id"],
        )
    out["chrom"] = out["chrom"].astype(str)
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out["segment_id"] = out["segment_id"].astype(str)
    return out


def intervals_by_chrom(intervals):
    out = {}
    for row in intervals.itertuples(index=False):
        out.setdefault(row.chrom, []).append((row.start, row.end, row.segment_id))
    return out


def positions_in_intervals(pos, intervals):
    keep = np.zeros(pos.shape[0], dtype=bool)
    for start, end, _ in intervals:
        keep |= (pos > start) & (pos <= end)
    return keep


def segment_for_position(pos, intervals):
    for start, end, segment_id in intervals:
        if start < pos <= end:
            return segment_id
    return None


def panel_af(genotypes, rows):
    called = genotypes[rows]
    called = called[~np.isnan(called)]
    if called.size == 0:
        return float("nan"), 0
    return float(np.mean(called / 2.0)), int(called.size)


def main():
    intervals = load_intervals(HIGHCONF_INTERVALS)
    by_chrom = intervals_by_chrom(intervals)
    bim = load_bim_table(plink_file(PLINK_PREFIX, ".bim"))
    cth_ids = load_ids(CTH_SAMPLE_LIST)
    reference_ids = load_ids(REFERENCE_SAMPLE_LIST)

    rows = []
    for chrom in sorted(by_chrom):
        npz = cache_path(CACHE_DIR, chrom)
        pos_arr, geno, cache_iids = load_chrom_cache(npz)
        iid_to_row = {sample: i for i, sample in enumerate(cache_iids)}
        cth_rows = np.array([iid_to_row[s] for s in cth_ids if s in iid_to_row], dtype=np.intp)
        ref_rows = np.array([iid_to_row[s] for s in reference_ids if s in iid_to_row], dtype=np.intp)

        bim_chr = bim[bim["chrom"] == chrom].set_index("pos")
        keep = positions_in_intervals(pos_arr, by_chrom[chrom])

        for idx in np.flatnonzero(keep):
            pos = int(pos_arr[idx])
            if pos not in bim_chr.index:
                continue

            af_cth, n_cth = panel_af(geno[:, idx], cth_rows)
            af_ref, n_ref = panel_af(geno[:, idx], ref_rows)
            if n_cth < MIN_CALLED_PER_PANEL or n_ref < MIN_CALLED_PER_PANEL:
                continue
            if np.isnan(af_cth) or np.isnan(af_ref):
                continue
            if abs(af_cth - af_ref) < AF_DIFF_MIN:
                continue

            brow = bim_chr.loc[pos]
            if isinstance(brow, pd.DataFrame):
                brow = brow.iloc[0]

            rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": str(brow["ref"]),
                    "alt": str(brow["alt"]),
                    "segment_id": segment_for_position(pos, by_chrom[chrom]),
                    "cth_like_allele": 1 if af_cth > af_ref else 0,
                    "af_cth": round(af_cth, 6),
                    "af_reference": round(af_ref, 6),
                    "af_diff": round(abs(af_cth - af_ref), 6),
                }
            )

    out = pd.DataFrame(rows).sort_values(["chrom", "pos"])
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_TSV, sep="\t", index=False)


if __name__ == "__main__":
    main()
