#!/usr/bin/env python3

# Purpose: calculate allele sharing for observed SPrime candidate segments.
# Input: PLINK files, candidate BED, and sample lists for test, CTH, and AmazonE groups.
# Output: one TSV with T statistics for candidate segments.
# Software: Python, pandas, PLINK

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from allele_sharing_lib import (
    allele_sharing_for_window,
    load_bim_table,
    load_ids,
    plink_file,
    summary_columns,
    write_keep,
)


PLINK_PREFIX = Path("<path/to/plink_prefix_without_extension>")
CANDIDATE_BED = Path("<path/to/sprime_candidate_segments.bed>")
TEST_SAMPLE_LIST = Path("<path/to/test_population_samples.list>")
CTH_SAMPLE_LIST = Path("<path/to/CTH_samples.list>")
AMAZONE_SAMPLE_LIST = Path("<path/to/AmazonE_samples.list>")
OUTPUT_TSV = Path("<path/to/candidate_allele_sharing_summary.tsv>")
RECIPIENT_LABEL = "<recipient_label>"


def main():
    bim = load_bim_table(plink_file(PLINK_PREFIX, ".bim"))
    bed = pd.read_csv(
        CANDIDATE_BED,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "segment", "score"],
    )

    test = load_ids(TEST_SAMPLE_LIST)
    cth = load_ids(CTH_SAMPLE_LIST)
    amazone = load_ids(AMAZONE_SAMPLE_LIST)
    cols = summary_columns(RECIPIENT_LABEL)
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(prefix="allele_sharing_candidate_") as tmpdir, OUTPUT_TSV.open("w") as out:
        tmp = Path(tmpdir)
        keep_path = tmp / "samples.keep"
        write_keep(test | cth | amazone, keep_path)
        out.write("\t".join(cols) + "\n")

        for i, seg in bed.iterrows():
            n_snps, n_compared, n_only, t_stat = allele_sharing_for_window(
                PLINK_PREFIX,
                keep_path,
                str(seg["chrom"]),
                int(seg["start"]),
                int(seg["end"]),
                bim,
                test,
                cth,
                amazone,
                tmp / f"segment_{i + 1}",
            )
            row = [
                str(seg["segment"]),
                str(seg["chrom"]),
                str(int(seg["start"])),
                str(int(seg["end"])),
                str(int(seg["end"]) - int(seg["start"])),
                str(n_snps),
                str(n_only),
                str(n_compared),
                "" if pd.isna(t_stat) else f"{t_stat:.6f}",
            ]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
