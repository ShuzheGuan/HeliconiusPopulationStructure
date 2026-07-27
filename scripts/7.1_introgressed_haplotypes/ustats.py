#!/usr/bin/env python3

"""Identify U-statistic candidate alleles from a SNP-major PLINK dataset.

Each PLINK allele is tested separately. An allele is retained when its
frequency is at most W in the reference population and at least Y in the donor
population. The output reports its frequency X in the recipient population.
Frequencies are calculated from non-missing diploid genotypes in each panel.
"""

from pathlib import Path

import numpy as np
import pandas as pd


# Replace these placeholders before running the script.
PLINK_PREFIX = Path("<path/to/plink_prefix_without_extension>")
REFERENCE_SAMPLE_LIST = Path("<path/to/reference_samples.list>")
DONOR_SAMPLE_LIST = Path("<path/to/donor_samples.list>")
RECIPIENT_SAMPLE_LIST = Path("<path/to/recipient_samples.list>")
OUTPUT_TSV = Path("<path/to/Ustats_output.tsv>")

# U-statistic frequency thresholds used in the analysis.
W = 0.01
Y = 0.95

# Leave empty to scan every contig, or provide one PLINK contig name.
CONTIG = ""
CHUNK_SNPS = 50_000


# PLINK BED SNP-major codes:
# 0 = homozygous A1, 1 = missing, 2 = heterozygous, 3 = homozygous A2.
_A1_DOSAGE = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)


def load_ids(path):
    ids = [
        line.strip().split()[0]
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    if not ids:
        raise SystemExit(f"Empty sample list: {path}")
    return set(ids)


def population_a1_frequency(dosage_a1, sample_mask):
    panel = dosage_a1[:, sample_mask]
    called = np.isfinite(panel)
    allele_count = np.nansum(panel, axis=1)
    allele_number = 2.0 * called.sum(axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(
            allele_number > 0,
            allele_count / allele_number,
            np.nan,
        )


def bytes_per_snp(n_samples):
    return (n_samples + 3) // 4


def unpack_bed_chunk(raw, n_samples):
    n_snps, n_bytes = raw.shape
    codes = np.empty((n_snps, n_bytes * 4), dtype=np.uint8)

    for offset, shift in enumerate((0, 2, 4, 6)):
        codes[:, offset::4] = (raw >> shift) & 3

    return _A1_DOSAGE[codes[:, :n_samples]]


def scan_u_sites(
    bed_path,
    bim,
    fam_iids,
    reference_ids,
    donor_ids,
    recipient_ids,
    snp_offset,
):
    n_samples = len(fam_iids)
    n_snps = len(bim)
    output_columns = ["snp_id", "contig", "pos", "X"]

    if n_snps == 0:
        return pd.DataFrame(columns=output_columns)

    reference_mask = np.isin(fam_iids, list(reference_ids))
    donor_mask = np.isin(fam_iids, list(donor_ids))
    recipient_mask = np.isin(fam_iids, list(recipient_ids))

    if (
        not reference_mask.any()
        or not donor_mask.any()
        or not recipient_mask.any()
    ):
        raise SystemExit(
            "A population has no samples in the PLINK FAM file: "
            f"reference={reference_mask.sum()}, "
            f"donor={donor_mask.sum()}, "
            f"recipient={recipient_mask.sum()}"
        )

    bytes_snp = bytes_per_snp(n_samples)
    output = []

    with bed_path.open("rb") as bed_file:
        magic = bed_file.read(3)
        if magic != b"\x6c\x1b\x01":
            raise SystemExit(
                f"Not a SNP-major PLINK BED file: {bed_path}"
            )

        bed_file.seek(3 + snp_offset * bytes_snp)

        for start in range(0, n_snps, CHUNK_SNPS):
            n_chunk = min(CHUNK_SNPS, n_snps - start)
            raw = np.fromfile(
                bed_file,
                dtype=np.uint8,
                count=n_chunk * bytes_snp,
            )

            if raw.size != n_chunk * bytes_snp:
                raise SystemExit(
                    f"Short BED read at SNP offset {snp_offset + start}"
                )

            dosage_a1 = unpack_bed_chunk(
                raw.reshape(n_chunk, bytes_snp),
                n_samples,
            )

            reference_a1 = population_a1_frequency(
                dosage_a1,
                reference_mask,
            )
            donor_a1 = population_a1_frequency(
                dosage_a1,
                donor_mask,
            )
            recipient_a1 = population_a1_frequency(
                dosage_a1,
                recipient_mask,
            )

            bim_chunk = bim.iloc[start : start + n_chunk]

            # Test A1 and A2 separately.
            for reference_frequency, donor_frequency, recipient_frequency in (
                (reference_a1, donor_a1, recipient_a1),
                (1.0 - reference_a1, 1.0 - donor_a1, 1.0 - recipient_a1),
            ):
                keep = (
                    (reference_frequency <= W)
                    & (donor_frequency >= Y)
                    & np.isfinite(recipient_frequency)
                )

                if not np.any(keep):
                    continue

                indices = np.flatnonzero(keep)
                hits = bim_chunk.iloc[indices][
                    ["snp_id", "contig", "pos"]
                ].copy()
                hits["X"] = recipient_frequency[indices]
                output.append(hits)

    if not output:
        return pd.DataFrame(columns=output_columns)

    return pd.concat(output, ignore_index=True)[output_columns]


def main():
    bed_path = Path(f"{PLINK_PREFIX}.bed")
    bim_path = Path(f"{PLINK_PREFIX}.bim")
    fam_path = Path(f"{PLINK_PREFIX}.fam")

    for path in (
        bed_path,
        bim_path,
        fam_path,
        REFERENCE_SAMPLE_LIST,
        DONOR_SAMPLE_LIST,
        RECIPIENT_SAMPLE_LIST,
    ):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}")

    fam = pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        usecols=[1],
        names=["iid"],
        dtype=str,
    )
    fam_iids = fam["iid"].to_numpy()

    bim = pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        names=["contig", "snp_id", "cm", "pos", "a1", "a2"],
        dtype={
            "contig": str,
            "snp_id": str,
            "pos": np.int64,
        },
    )[["contig", "snp_id", "pos", "a1", "a2"]]

    snp_offset = 0

    if CONTIG:
        indices = np.flatnonzero(
            bim["contig"].to_numpy() == CONTIG
        )

        if indices.size == 0:
            OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                columns=["snp_id", "contig", "pos", "X"]
            ).to_csv(OUTPUT_TSV, sep="\t", index=False)
            return

        if int(indices[-1] - indices[0] + 1) != int(indices.size):
            raise SystemExit(
                f"Contig {CONTIG} is not contiguous in the BIM file"
            )

        snp_offset = int(indices[0])
        bim = bim.iloc[
            int(indices[0]) : int(indices[-1]) + 1
        ].reset_index(drop=True)

    reference_ids = load_ids(REFERENCE_SAMPLE_LIST)
    donor_ids = load_ids(DONOR_SAMPLE_LIST)
    recipient_ids = load_ids(RECIPIENT_SAMPLE_LIST)
    fam_set = set(fam_iids.tolist())

    for label, ids in (
        ("reference", reference_ids),
        ("donor", donor_ids),
        ("recipient", recipient_ids),
    ):
        missing = ids - fam_set
        if missing:
            example = next(iter(missing))
            raise SystemExit(
                f"{label}: {len(missing)} sample IDs are absent "
                f"from the FAM file, for example {example}"
            )

    sites = scan_u_sites(
        bed_path,
        bim,
        fam_iids,
        reference_ids,
        donor_ids,
        recipient_ids,
        snp_offset,
    )

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    sites.to_csv(OUTPUT_TSV, sep="\t", index=False)

    print(
        f"Wrote {len(sites):,} U-statistic candidate alleles to "
        f"{OUTPUT_TSV}"
    )


if __name__ == "__main__":
    main()
