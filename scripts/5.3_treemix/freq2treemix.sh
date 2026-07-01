#!/usr/bin/env bash
set -euo pipefail

# Purpose: convert PLINK .frq.strat output to TreeMix allele-count format.
# Input: PLINK .frq.strat file from plink --freq --within.
# Output: gzipped TreeMix frequency file.
# Software: Python, pandas

INPUT_FRQ_STRAT="<path/to/input.frq.strat>"
OUTPUT_FRQ_GZ="<path/to/output_treemix.frq.gz>"

mkdir -p "$(dirname "$OUTPUT_FRQ_GZ")"

export INPUT_FRQ_STRAT OUTPUT_FRQ_GZ

python3 - <<'PY'
import gzip
import os
import pandas as pd

df = pd.read_csv(os.environ["INPUT_FRQ_STRAT"], sep=r"\s+")
cols = {c.upper(): c for c in df.columns}

snp_col = cols["SNP"]
pop_col = cols["CLST"]
nchr = pd.to_numeric(df[cols["NCHROBS"]], errors="coerce").fillna(0).astype(int)

if "MAC" in cols:
    a1_count = pd.to_numeric(df[cols["MAC"]], errors="coerce").fillna(0).astype(int)
else:
    maf = pd.to_numeric(df[cols["MAF"]], errors="coerce").fillna(0)
    a1_count = (maf * nchr).round().astype(int)

df["counts"] = a1_count.astype(str) + "," + (nchr - a1_count).clip(lower=0).astype(str)

snp_order = df[snp_col].drop_duplicates().tolist()
pop_order = df[pop_col].drop_duplicates().tolist()

matrix = df.pivot(index=snp_col, columns=pop_col, values="counts")
matrix = matrix.reindex(index=snp_order, columns=pop_order).fillna("0,0")

with gzip.open(os.environ["OUTPUT_FRQ_GZ"], "wt") as out:
    out.write("\t".join(map(str, matrix.columns)) + "\n")
    for _, row in matrix.iterrows():
        out.write("\t".join(row.astype(str)) + "\n")
PY
