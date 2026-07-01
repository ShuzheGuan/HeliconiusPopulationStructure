#!/usr/bin/env python3

# Purpose: sample posterior topologies from BPP A01 big-tree MCMC files.
# Input: block-list TSVs and BPP *.mcmc.txt files containing one tree per line.
# Output: one TSV of sampled posterior trees for visualization.
# Software: Python

from pathlib import Path
import csv
import random

DATASET_LABEL = "<dataset_label>"
OUTPUT_TSV = Path("<path/to/sampled_mcmc_topologies.tsv>")

MCMC_ROOTS = {
    "ncds": Path("<path/to/ncds_bigtree_output_root>"),
    "cds": Path("<path/to/cds_bigtree_output_root>"),
}

BLOCK_LISTS = [
    {"path": Path("<path/to/monophyly_blocks.tsv>"), "tree_type": "monophyly"},
    {"path": Path("<path/to/introgression_blocks.tsv>"), "tree_type": "introgression"},
]

N_CHAINS = 5
N_SAMPLES = 100000
RANDOM_SEED = 20260523

TREE_RENAME = {}


def apply_renames(tree):
    out = tree
    for old, new in TREE_RENAME.items():
        out = out.replace(old, new)
    return out


def read_exact_line(path, line_no):
    if not path.exists():
        return "NA"

    with path.open() as handle:
        for i, line in enumerate(handle, start=1):
            if i == line_no:
                return apply_renames(line.rstrip("\n"))
    return "NA"


def read_block_list(path, tree_type, rng):
    rows = []

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            coding_type = row["coding_type"]
            contig = row["contig"]
            block = row["block"]
            chain = rng.randint(1, N_CHAINS)
            sample_index = rng.randint(1, N_SAMPLES)

            mcmc_path = (
                MCMC_ROOTS[coding_type]
                / contig
                / f"{contig}_block{block}_chain{chain}.mcmc.txt"
            )

            rows.append(
                {
                    "dataset": DATASET_LABEL,
                    "contig": contig,
                    "block": block,
                    "coding_type": coding_type,
                    "tree_type": tree_type,
                    "chain": chain,
                    "sample_index": sample_index,
                    "tree": read_exact_line(mcmc_path, sample_index),
                }
            )

    return rows


def main():
    rng = random.Random(RANDOM_SEED)
    rows = []

    for block_list in BLOCK_LISTS:
        rows.extend(
            read_block_list(
                path=block_list["path"],
                tree_type=block_list["tree_type"],
                rng=rng,
            )
        )

    fieldnames = [
        "dataset",
        "contig",
        "block",
        "coding_type",
        "tree_type",
        "chain",
        "sample_index",
        "tree",
    ]

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
