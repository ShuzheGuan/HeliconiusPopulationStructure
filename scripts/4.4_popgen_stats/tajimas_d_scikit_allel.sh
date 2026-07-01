#!/usr/bin/env bash
set -euo pipefail

# Purpose: calculate windowed Tajima's D for each population with scikit-allel.
# Input: VCF and two-column sample-population TSV.
# Output: one Tajima's D table per population and one combined table.
# Software: Python, scikit-allel, pandas

VCF="<path/to/invariant_and_biallelic_sites.vcf.gz>"
POPULATIONS_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_DIR="<path/to/tajimas_d_output_dir>"
OUTPUT_PREFIX="<output_prefix>"
WINDOW_SIZE=50000
MIN_SITES=3

mkdir -p "$OUTPUT_DIR"
export VCF POPULATIONS_TSV OUTPUT_DIR OUTPUT_PREFIX WINDOW_SIZE MIN_SITES

python3 - <<'PY'
import gzip
import os
import re

import allel
import numpy as np
import pandas as pd

vcf = os.environ["VCF"]
pop_tsv = os.environ["POPULATIONS_TSV"]
out_dir = os.environ["OUTPUT_DIR"]
prefix = os.environ["OUTPUT_PREFIX"]
window_size = int(os.environ["WINDOW_SIZE"])
min_sites = int(os.environ["MIN_SITES"])

def clean(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x)).strip("_")

def vcf_samples(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]
    raise RuntimeError("VCF header not found")

pops = pd.read_csv(pop_tsv, sep="\t", header=None, names=["sample", "population"], dtype=str).dropna()
samples_in_vcf = set(vcf_samples(vcf))
rows_all = []

for pop, table in pops.groupby("population", sort=True):
    samples = [s for s in table["sample"] if s in samples_in_vcf]
    rows = []

    callset = allel.read_vcf(vcf, samples=samples, fields=["variants/CHROM", "variants/POS", "calldata/GT"])
    chroms = callset["variants/CHROM"]
    positions = callset["variants/POS"].astype(int)
    genotypes = allel.GenotypeArray(callset["calldata/GT"])

    for chrom in sorted(set(chroms)):
        idx = chroms == chrom
        pos = positions[idx]
        ac = genotypes.compress(idx, axis=0).count_alleles(max_allele=1)
        stop = int(pos.max())
        tajd, windows, counts = allel.windowed_tajima_d(pos, ac, size=window_size, start=1, stop=stop)
        tajd = np.asarray(tajd, dtype=float)
        tajd[np.asarray(counts) < min_sites] = np.nan
        rows.append(pd.DataFrame({
            "population": pop,
            "chrom": chrom,
            "window_start": windows[:, 0],
            "window_end": windows[:, 1],
            "variant_count": counts,
            "tajimasD": tajd,
        }))

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(os.path.join(out_dir, f"{prefix}_{clean(pop)}_{window_size}bp_tajimas_d.tsv"), sep="\t", index=False)
    rows_all.append(out)

pd.concat(rows_all, ignore_index=True).to_csv(
    os.path.join(out_dir, f"{prefix}_{window_size}bp_tajimas_d_all_populations.tsv"),
    sep="\t",
    index=False,
)
PY
