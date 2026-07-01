"""Allele sharing: (test∩CTH \\ AmazonE) / alleles in test∪CTH∪AmazonE."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

META_COLS = ["segment", "chrom", "start", "end", "length_bp", "n_snps"]


def summary_columns(recipient_label: str) -> list[str]:
    return META_COLS + [
        f"n_shared_{recipient_label}_CTH_only",
        "n_compared",
        "T",
    ]


def plink_file(prefix: Path, ext: str) -> Path:
    return Path(f"{prefix}{ext}")


def load_ids(path: Path) -> set[str]:
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip()}


def write_keep(samples: set[str], path: Path) -> None:
    with path.open("w") as fh:
        for s in sorted(samples):
            fh.write(f"{s}\t{s}\n")


def load_bim_table(bim_path: Path) -> pd.DataFrame:
    return pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 3, 4, 5],
        names=["chrom", "snp", "pos", "ref", "alt"],
        dtype={"chrom": str, "snp": str, "pos": np.int64, "ref": str, "alt": str},
    )


def snps_in_window(bim: pd.DataFrame, chrom: str, bed_start: int, bed_end: int) -> int:
    from_bp = bed_start + 1
    to_bp = bed_end
    sub = bim[(bim["chrom"] == chrom) & (bim["pos"] >= from_bp) & (bim["pos"] <= to_bp)]
    return len(sub)


def run_plink_recode_a(
    bfile: Path,
    keep: Path,
    chrom: str,
    bed_start: int,
    bed_end: int,
    out_prefix: Path,
) -> Path:
    from_bp = bed_start + 1
    to_bp = bed_end
    range_file = out_prefix.parent / f"{out_prefix.name}.range"
    range_file.write_text(f"{chrom}\t{from_bp}\t{to_bp}\t{out_prefix.name}\n")

    cmd = [
        "plink",
        "--bfile",
        str(bfile),
        "--allow-extra-chr",
        "--keep",
        str(keep),
        "--extract",
        "range",
        str(range_file),
        "--recode",
        "A",
        "--out",
        str(out_prefix),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log_path = Path(f"{out_prefix}.log")
        err = log_path.read_text() if log_path.is_file() else (proc.stderr or proc.stdout or "")
        raise SystemExit(
            f"Error: plink --recode A failed for {chrom}:{bed_start}-{bed_end}\n{err}"
        )

    raw_path = Path(f"{out_prefix}.raw")
    if not raw_path.is_file():
        raise SystemExit(f"Error: missing {raw_path}")
    return raw_path


def _alleles_present(genotypes: np.ndarray) -> set[str]:
    g = genotypes[~np.isnan(genotypes)]
    out: set[str] = set()
    if g.size == 0:
        return out
    if np.any((g == 0) | (g == 1)):
        out.add("ref")
    if np.any((g == 1) | (g == 2)):
        out.add("alt")
    return out


def _raw_snp_columns(raw: pd.DataFrame) -> list[str]:
    meta = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}
    return [c for c in raw.columns if c not in meta]


def allele_sharing_counts(
    raw: pd.DataFrame,
    test: set[str],
    cth: set[str],
    amazone: set[str],
) -> tuple[int, int, float]:
    """Return n_compared (union), n_cth_only (test∩CTH\\E), T."""
    iids = raw["IID"].astype(str)
    test_mask = iids.isin(test).to_numpy()
    cth_mask = iids.isin(cth).to_numpy()
    e_mask = iids.isin(amazone).to_numpy()

    n_compared = 0
    n_cth_only = 0
    n_with_e = 0
    for col in _raw_snp_columns(raw):
        g = pd.to_numeric(raw[col], errors="coerce").to_numpy(dtype=float)
        alleles_test = _alleles_present(g[test_mask])
        alleles_cth = _alleles_present(g[cth_mask])
        alleles_e = _alleles_present(g[e_mask])
        union = alleles_test | alleles_cth | alleles_e
        for allele in ("ref", "alt"):
            if allele not in union:
                continue
            n_compared += 1
            in_test = allele in alleles_test
            in_cth = allele in alleles_cth
            in_e = allele in alleles_e
            if in_test and in_cth and not in_e:
                n_cth_only += 1
            elif in_test and in_cth and in_e:
                n_with_e += 1

    t_stat = float(n_cth_only / n_compared) if n_compared > 0 else float("nan")
    return n_compared, n_cth_only, t_stat


def allele_sharing_for_window(
    bfile: Path,
    keep_path: Path,
    chrom: str,
    bed_start: int,
    bed_end: int,
    bim: pd.DataFrame,
    test: set[str],
    cth: set[str],
    amazone: set[str],
    tmp_prefix: Path,
) -> tuple[int, int, int, float]:
    n_snps = snps_in_window(bim, chrom, bed_start, bed_end)
    if n_snps == 0:
        return 0, 0, 0, float("nan")

    raw_path = run_plink_recode_a(bfile, keep_path, chrom, bed_start, bed_end, tmp_prefix)
    raw = pd.read_csv(raw_path, sep=r"\s+")
    n_compared, n_only, t_stat = allele_sharing_counts(raw, test, cth, amazone)
    return n_snps, n_compared, n_only, t_stat
