#!/usr/bin/env python3

"""Summarize BPP locus lengths and informative sites."""

import csv
from collections import Counter
from pathlib import Path


BPP_ROOT = Path("<path/to/bpp_locus_dir>")
OUTPUT_TSV = Path("<path/to/locus_summary.tsv>")


def read_bpp(path):
    with path.open() as handle:
        header = handle.readline().split()
        if len(header) < 2:
            return [], 0
        length = int(header[1])
        seqs = [line.split()[-1].upper() for line in handle if line.strip()]
    return seqs, length


def informative_sites(seqs):
    if not seqs:
        return 0

    total = 0
    for column in zip(*seqs):
        counts = Counter(base for base in column if base in {"A", "C", "G", "T"})
        if len(counts) >= 2 and all(n >= 2 for n in counts.values()):
            total += 1
    return total


OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT_TSV.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["locus_id", "n_samples", "locus_length", "informative_sites"])

    for path in sorted(BPP_ROOT.rglob("*.bpp.txt")):
        seqs, length = read_bpp(path)
        writer.writerow([
            path.stem.replace(".bpp", ""),
            len(seqs),
            length,
            informative_sites(seqs),
        ])
