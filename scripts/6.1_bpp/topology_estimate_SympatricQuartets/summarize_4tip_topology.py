#!/usr/bin/env python3

# Purpose: summarize four-tip topology frequencies from BPP MCMC tree files.
# Input: directories containing *_chain*.mcmc.txt files from 4-tip A01 runs.
# Output: topology-frequency TSV and high-confidence subset TSV.
# Software: Python

from pathlib import Path
from collections import Counter, defaultdict
import csv
import re

RUN_ROOT = Path("<path/to/four_tip_topology_output_dir>")
OUTPUT_TSV = Path("<path/to/4tip_topology_estimate.tsv>")
HIGH_CONFIDENCE_TSV = Path("<path/to/4tip_topology_highcon.tsv>")

MIN_TREES = 10000
HIGH_CONFIDENCE_THRESHOLD = 0.5

OUTGROUP = "<outgroup_label>"
TAXON_A_REGION_1 = "<taxon_A_region_1>"
TAXON_B_REGION_1 = "<taxon_B_region_1>"
TAXON_A_REGION_2 = "<taxon_A_region_2>"
TAXON_B_REGION_2 = "<taxon_B_region_2>"

TREE_RENAME = {}
TOPOLOGY_COLUMNS = ["tree1", "tree2a", "tree2b", "tree3", "tree4a", "tree4b", "other"]


def rename_taxon(label):
    out = label
    for old, new in TREE_RENAME.items():
        out = out.replace(old, new)
    return out


def strip_newick_metadata(newick):
    newick = re.sub(r"#[0-9.eE+-]+", "", newick)
    newick = re.sub(r":\s*[0-9.eE+-]+", "", newick)
    return newick.replace(";", "").strip()


def get_clades(newick):
    tokens = re.split(r"([(),])", strip_newick_metadata(newick))
    tokens = [token.strip() for token in tokens if token.strip()]

    stack = []
    clades = []

    for token in tokens:
        if token == "(":
            stack.append(set())
        elif token == ")":
            if not stack:
                continue
            group = stack.pop()
            if stack:
                stack[-1].update(group)

            ingroup = {rename_taxon(x) for x in group if rename_taxon(x) != OUTGROUP}
            if len(ingroup) > 1:
                clades.append(frozenset(ingroup))
        elif token == ",":
            continue
        elif stack:
            stack[-1].add(rename_taxon(token))

    return set(clades)


def classify_tree(newick):
    clades = get_clades(newick)

    a1 = TAXON_A_REGION_1
    b1 = TAXON_B_REGION_1
    a2 = TAXON_A_REGION_2
    b2 = TAXON_B_REGION_2

    geography_1 = frozenset({a1, b1})
    geography_2 = frozenset({a2, b2})
    lineage_a = frozenset({a1, a2})
    lineage_b = frozenset({b1, b2})
    a1_b2 = frozenset({a1, b2})

    trio_a1_a2_b2 = frozenset({a1, a2, b2})
    trio_a1_b1_b2 = frozenset({a1, b1, b2})

    if geography_1 in clades and geography_2 in clades:
        return "tree1"
    if lineage_a in clades and lineage_b in clades:
        return "tree3"
    if geography_2 in clades and trio_a1_a2_b2 in clades:
        return "tree2a"
    if geography_1 in clades and trio_a1_b1_b2 in clades:
        return "tree2b"
    if a1_b2 in clades and trio_a1_a2_b2 in clades:
        return "tree4a"
    if lineage_a in clades and trio_a1_a2_b2 in clades:
        return "tree4b"
    return "other"


def parse_run_directory_name(name):
    parts = name.split("_")
    coding_type = next((x for x in parts if x in {"cds", "ncds"}), "NA")
    if coding_type == "NA":
        return "NA", name

    coding_index = parts.index(coding_type)
    contig = "_".join(parts[coding_index + 1 :])
    return coding_type, contig


def summarize_directory(run_dir):
    coding_type, contig = parse_run_directory_name(run_dir.name)
    block_files = defaultdict(list)

    for path in sorted(run_dir.glob("*_chain*.mcmc.txt")):
        block_id = path.name.split("_chain", 1)[0]
        block_files[block_id].append(path)

    rows = []
    for block_id, paths in sorted(block_files.items()):
        counts = Counter()
        total = 0

        for path in paths:
            with path.open(errors="ignore") as handle:
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    counts[classify_tree(line)] += 1
                    total += 1

        if total < MIN_TREES:
            continue

        row = {
            "chromosome": contig,
            "coding_type": coding_type,
            "block_id": block_id,
            "n_trees": total,
        }
        for topology in TOPOLOGY_COLUMNS:
            row[topology] = counts.get(topology, 0) / total
        rows.append(row)

    return rows


def write_tsv(path, rows):
    fieldnames = ["chromosome", "coding_type", "block_id"] + TOPOLOGY_COLUMNS + ["n_trees"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = []
    for run_dir in sorted(path for path in RUN_ROOT.iterdir() if path.is_dir()):
        rows.extend(summarize_directory(run_dir))

    high_confidence = [
        row
        for row in rows
        if max(row[topology] for topology in TOPOLOGY_COLUMNS) > HIGH_CONFIDENCE_THRESHOLD
    ]

    write_tsv(OUTPUT_TSV, rows)
    write_tsv(HIGH_CONFIDENCE_TSV, high_confidence)


if __name__ == "__main__":
    main()
