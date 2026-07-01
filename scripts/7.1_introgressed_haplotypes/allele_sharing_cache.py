"""Per-chrom genotype cache + null placement for allele-sharing null pools."""

from __future__ import annotations

import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from allele_sharing_lib import (
    allele_sharing_counts,
    load_bim_table,
    load_ids,
    plink_file,
    summary_columns,
    write_keep,
)

MIN_NULL_SNPS = 100

PACIFIC_POOLS = [(10, "Pacific"), (20, "Pacific"), (50, "Pacific"), (100, "Pacific")]
AMAZONW_POOLS = [
    (10, "AmazonW"),
    (20, "AmazonW"),
    (50, "AmazonW"),
    (100, "AmazonW"),
    (200, "AmazonW"),
    (500, "AmazonW"),
    (1000, "AmazonW"),
]
NULL_POOLS = PACIFIC_POOLS + AMAZONW_POOLS
NULL_JOBS_PER_POOL = 100
NULL_TOTAL_PER_POOL = 100_000
NULL_PER_JOB = NULL_TOTAL_PER_POOL // NULL_JOBS_PER_POOL

EXCLUDE_BEDS = [
    "Pacific_AmazonE_sprime_segments.bed",
    "AmazonW_AmazonE_sprime_segments.bed",
]


def pool_stem(bin_kb: int, label: str) -> str:
    return f"{label}_null_{bin_kb}kb"


def pool_from_task_id(task_id: int) -> tuple[int, str, str]:
    """Legacy 1..11 mapping (one job per pool)."""
    if task_id < 1 or task_id > len(NULL_POOLS):
        raise ValueError(f"task_id must be 1..{len(NULL_POOLS)}, got {task_id}")
    bin_kb, label = NULL_POOLS[task_id - 1]
    return bin_kb, label, pool_stem(bin_kb, label)


def pool_batch_from_task_id(task_id: int) -> tuple[int, str, str, int]:
    """SLURM array 1..1100: 100 jobs × 1k nulls per pool (100k total per pool)."""
    n_tasks = len(NULL_POOLS) * NULL_JOBS_PER_POOL
    if task_id < 1 or task_id > n_tasks:
        raise ValueError(f"task_id must be 1..{n_tasks}, got {task_id}")
    pool_idx = (task_id - 1) // NULL_JOBS_PER_POOL
    batch = (task_id - 1) % NULL_JOBS_PER_POOL + 1
    bin_kb, label = NULL_POOLS[pool_idx]
    return bin_kb, label, pool_stem(bin_kb, label), batch


def null_pool_dir(null_root: Path, pool_stem: str) -> Path:
    return null_root / pool_stem


def cache_path(raw_dir: Path, chrom: str) -> Path:
    return raw_dir / f"{chrom}.npz"


def load_chrom_lengths(bim: pd.DataFrame) -> dict[str, int]:
    return bim.groupby("chrom")["pos"].max().astype(int).to_dict()


def load_excluded_intervals(bed_paths: list[Path]) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in bed_paths:
        bed = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "name", "score"],
        )
        for row in bed.itertuples(index=False):
            by_chrom[str(row.chrom)].append((int(row.start), int(row.end)))
    return dict(by_chrom)


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def window_overlaps_excluded(
    chrom: str,
    bed_start: int,
    bed_end: int,
    excluded: dict[str, list[tuple[int, int]]],
) -> bool:
    for ex_start, ex_end in excluded.get(chrom, []):
        if intervals_overlap(bed_start, bed_end, ex_start, ex_end):
            return True
    return False


def count_snps_in_window(
    pos_by_chrom: dict[str, np.ndarray],
    chrom: str,
    bed_start: int,
    bed_end: int,
) -> int:
    from_bp = bed_start + 1
    to_bp = bed_end
    pos = pos_by_chrom.get(chrom)
    if pos is None or pos.size == 0:
        return 0
    return int(np.sum((pos >= from_bp) & (pos <= to_bp)))


def sample_null_window(
    chrom_lengths: dict[str, int],
    excluded: dict[str, list[tuple[int, int]]],
    length_bp: int,
    rng: random.Random,
) -> tuple[str, int, int] | None:
    chroms = [c for c, ln in chrom_lengths.items() if ln >= length_bp]
    if not chroms:
        return None
    for _ in range(5000):
        chrom = rng.choice(chroms)
        max_start = chrom_lengths[chrom] - length_bp
        if max_start <= 0:
            continue
        bed_start = rng.randint(0, max_start)
        bed_end = bed_start + length_bp
        if window_overlaps_excluded(chrom, bed_start, bed_end, excluded):
            continue
        return chrom, bed_start, bed_end
    return None


def build_chrom_cache(
    bfile: Path,
    keep_path: Path,
    chrom: str,
    out_npz: Path,
) -> None:
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cache_chr_") as tmpdir:
        tmp = Path(tmpdir)
        prefix = tmp / chrom
        cmd = [
            "plink",
            "--bfile",
            str(bfile),
            "--allow-extra-chr",
            "--keep",
            str(keep_path),
            "--chr",
            chrom,
            "--recode",
            "A",
            "--out",
            str(prefix),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log = Path(f"{prefix}.log")
            err = log.read_text() if log.is_file() else (proc.stderr or proc.stdout or "")
            raise SystemExit(f"plink failed for {chrom}: {err}")

        raw = pd.read_csv(f"{prefix}.raw", sep=r"\s+")
        meta = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}
        snp_cols = [c for c in raw.columns if c not in meta]
        iids = raw["IID"].astype(str).tolist()
        geno = raw[snp_cols].to_numpy(dtype=np.float32)
        bim_chr = load_bim_table(plink_file(bfile, ".bim"))
        bim_chr = bim_chr[bim_chr["chrom"] == chrom].sort_values("pos")
        if len(snp_cols) != len(bim_chr):
            raise SystemExit(
                f"chr {chrom}: .raw cols {len(snp_cols)} != bim SNPs {len(bim_chr)}"
            )
        pos = bim_chr["pos"].to_numpy(dtype=np.int64)
        np.savez_compressed(out_npz, pos=pos, geno=geno, iids=np.array(iids, dtype=object))


def load_chrom_cache(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    data = np.load(path, allow_pickle=True)
    iids = [str(x) for x in data["iids"].tolist()]
    return data["pos"], data["geno"], iids


def allele_sharing_from_cache_window(
    pos: np.ndarray,
    geno: np.ndarray,
    iids: list[str],
    bed_start: int,
    bed_end: int,
    test: set[str],
    cth: set[str],
    amazone: set[str],
) -> tuple[int, int, int, float]:
    from_bp = bed_start + 1
    to_bp = bed_end
    mask = (pos >= from_bp) & (pos <= to_bp)
    n_snps = int(mask.sum())
    if n_snps == 0:
        return 0, 0, 0, float("nan")

    sub = geno[:, mask]
    raw = pd.DataFrame(sub, columns=[f"snp{i}" for i in range(n_snps)])
    raw.insert(0, "IID", iids)
    n_compared, n_only, t_stat = allele_sharing_counts(raw, test, cth, amazone)
    return n_snps, n_compared, n_only, t_stat
