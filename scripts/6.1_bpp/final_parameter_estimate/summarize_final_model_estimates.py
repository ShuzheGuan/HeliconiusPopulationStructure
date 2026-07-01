#!/usr/bin/env python3

# Purpose: parse final BPP MSC-I output files into one chain-level estimate table.
# Input: BPP output text files named like <model>_<dataset>_msci_chain<chain>.txt.
# Output: one wide TSV with one row per dataset and chain.
# Software: Python

from pathlib import Path
import csv
import re

RUN_ROOT = Path("<path/to/final_model_output_root>")
OUTPUT_TSV = Path("<path/to/final_model_estimates.tsv>")

DATASETS = ["sex_chr", "chr1", "random4k"]
NODE_LABELS = {}

FILE_RE = re.compile(r"_(?P<dataset>sex_chr|chr1|random4k)_msci_chain(?P<chain>[0-9]+)\.txt$")
HEADER_RE = re.compile(r"^\s*param\s+mean\s+median\s+", re.I)


def clean_param(raw):
    match = re.match(r"^(theta|tau):([0-9]+)$", raw)
    if match:
        prefix, idx = match.group(1), int(match.group(2))
        return f"{prefix}:{idx}:{NODE_LABELS.get(idx, f'Node{idx}')}"

    match = re.match(r"^phi:([0-9]+)<-([0-9]+)$", raw)
    if match:
        target, source = int(match.group(1)), int(match.group(2))
        target_label = NODE_LABELS.get(target, f"Node{target}")
        source_label = NODE_LABELS.get(source, f"Node{source}")
        return f"phi:{target}<-{source}:{target_label}<-{source_label}"

    return raw


def parse_output(path):
    lines = path.read_text(errors="ignore").splitlines()[-500:]
    header_i = None

    for i in range(len(lines) - 1, -1, -1):
        if HEADER_RE.match(lines[i]):
            header_i = i
            break

    if header_i is None:
        return {}

    stats = lines[header_i].split()[1:]
    parsed = {}

    for line in lines[header_i + 1 :]:
        if not line.strip() and parsed:
            break
        parts = line.split()
        if len(parts) < 2:
            continue

        param = clean_param(parts[0])
        values = parts[1 : 1 + len(stats)]
        for stat, value in zip(stats, values):
            parsed[f"{param}_{stat}"] = value

    return parsed


def infer_metadata(path):
    match = FILE_RE.search(path.name)
    if match is None:
        return None
    return {
        "data_set": match.group("dataset"),
        "chain": match.group("chain"),
    }


rows = []
fieldnames = ["data_set", "chain"]

for path in sorted(RUN_ROOT.rglob("*.txt")):
    if path.name.endswith(".mcmc.txt"):
        continue

    metadata = infer_metadata(path)
    if metadata is None:
        continue

    row = dict(metadata)
    row.update(parse_output(path))
    rows.append(row)

    for key in row:
        if key not in fieldnames:
            fieldnames.append(key)

OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT_TSV.open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
