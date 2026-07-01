#!/usr/bin/env python3

# Purpose: compare candidate allele-sharing scores with length-bin null distributions.
# Input: candidate score TSVs, candidate classifications, and null summary TSVs.
# Output: empirical p-value and null percentile for each candidate segment.
# Software: Python, pandas, numpy

from pathlib import Path

import numpy as np
import pandas as pd


CANDIDATE_SUMMARIES = {
    "<recipient_label>": Path("<path/to/candidate_allele_sharing_summary.tsv>"),
}
CLASSIFIED_TSV = Path("<path/to/candidate_classified.tsv>")
NULL_SUMMARY_DIR = Path("<path/to/null_summary_dir>")
OUTPUT_DIR = Path("<path/to/empirical_p_output_dir>")

MATCH_BY = "density"
SNP_FRAC = 0.50


def null_summary_path(classification):
    pop, rest = classification.split("_", 1)
    return NULL_SUMMARY_DIR / f"{pop}_null_{rest}_null_allele_sharing_summary.tsv"


def match_nulls(null_df, n_snps, length_bp):
    null_snps = pd.to_numeric(null_df["n_snps"], errors="coerce")
    null_len = pd.to_numeric(null_df["length_bp"], errors="coerce")

    if MATCH_BY == "density":
        candidate_value = n_snps / length_bp
        null_value = null_snps / null_len
    elif MATCH_BY == "count":
        candidate_value = n_snps
        null_value = null_snps
    else:
        raise ValueError("MATCH_BY must be 'density' or 'count'")

    lo = candidate_value * (1 - SNP_FRAC)
    hi = candidate_value * (1 + SNP_FRAC)
    return null_df[(null_value >= lo) & (null_value <= hi)]


def empirical_p_upper(candidate_t, null_t):
    n = len(null_t)
    if n == 0:
        return np.nan
    return (1 + np.sum(null_t >= candidate_t)) / (n + 1)


def null_percentile(candidate_t, null_t):
    n = len(null_t)
    if n == 0:
        return np.nan
    return 100 * (1 + np.sum(null_t < candidate_t)) / (n + 1)


def test_population(label, candidate_path, classified):
    candidates = pd.read_csv(candidate_path, sep="\t").rename(columns={"segment": "segment_id"})
    candidates = candidates.merge(
        classified[classified["test_population"] == label][["segment_id", "classification"]],
        on="segment_id",
        how="inner",
    )

    null_cache = {}
    rows = []

    for _, row in candidates.iterrows():
        classification = row["classification"]
        if classification not in null_cache:
            null_cache[classification] = pd.read_csv(null_summary_path(classification), sep="\t")

        candidate_t = float(row["T"])
        matched = match_nulls(null_cache[classification], float(row["n_snps"]), float(row["length_bp"]))
        null_t = pd.to_numeric(matched["T"], errors="coerce").dropna().to_numpy()

        rows.append(
            {
                "segment_id": row["segment_id"],
                "T": candidate_t,
                "n_null": len(null_t),
                "empirical_p": empirical_p_upper(candidate_t, null_t),
                "null_percentile": null_percentile(candidate_t, null_t),
            }
        )

    return pd.DataFrame(rows).sort_values("empirical_p", na_position="last")


def main():
    classified = pd.read_csv(CLASSIFIED_TSV, sep="\t")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for label, candidate_path in CANDIDATE_SUMMARIES.items():
        out = test_population(label, candidate_path, classified)
        out.to_csv(OUTPUT_DIR / f"{label}_allele_sharing_empirical_p.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
