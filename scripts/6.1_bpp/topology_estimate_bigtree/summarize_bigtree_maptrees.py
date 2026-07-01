#!/usr/bin/env python3

# Purpose: extract MAP topologies from BPP A01 big-tree output files.
# Input: BPP output text files containing a '(D)' MAP-tree section.
# Output: one TSV with dataset, coding type, contig, block, chain, MAP tree, and P.
# Software: Python

from pathlib import Path
import csv
import re

RUNS = [
    {
        "root": Path("<path/to/bigtree_output_root>"),
        "species": "<dataset_label>",
        "coding_type": "<cds_or_ncds>",
    },
]

OUTPUT_TSV = Path("<path/to/bigtree_map_tree_summary.tsv>")
POP_RENAME = {}

BLOCK_RE = re.compile(r"_block(?P<block>[0-9]+)_chain(?P<chain>[0-9]+)\.txt$")
P_RE = re.compile(r"\[P\s*=\s*([0-9.eE+-]+)\]")


def apply_renames(tree):
    out = tree
    for old, new in POP_RENAME.items():
        out = out.replace(old, new)
    return out


def parse_map_tree(path):
    found = False

    with path.open(errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("(D)"):
                found = True
                continue

            if not found or not stripped:
                continue

            match = P_RE.search(stripped)
            probability = match.group(1) if match else "NA"
            tree = P_RE.sub("", stripped).strip()
            return apply_renames(tree), probability

    return None


def infer_block_chain(path):
    match = BLOCK_RE.search(path.name)
    if match is None:
        return "NA", "NA"
    return match.group("block"), match.group("chain")


def main():
    fieldnames = ["species", "coding_type", "contig", "block", "chain", "map_tree", "P"]
    rows = []

    for run in RUNS:
        for path in sorted(run["root"].rglob("*_chain*.txt")):
            if path.name.endswith(".mcmc.txt"):
                continue

            result = parse_map_tree(path)
            if result is None:
                continue

            block, chain = infer_block_chain(path)
            map_tree, probability = result
            rows.append(
                {
                    "species": run["species"],
                    "coding_type": run["coding_type"],
                    "contig": path.parent.name,
                    "block": block,
                    "chain": chain,
                    "map_tree": map_tree,
                    "P": probability,
                }
            )

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
