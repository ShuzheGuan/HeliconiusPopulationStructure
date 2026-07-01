#!/usr/bin/env python3

# Purpose: parse pairwise BPP MSC-M run outputs into one estimates table.
# Input: BPP output text files named like <pop1>_vs_<pop2>_chain1.txt.
# Output: one TSV with posterior summaries for theta, tau, and migration parameters.
# Software: Python

from pathlib import Path
import csv
import re

RUN_ROOT = Path("<path/to/pairwise_mscm_run_outputs>")
OUTPUT_TSV = Path("<path/to/mscm_estimates.tsv>")
DATASET_LABEL = "<dataset_label>"

PAIR_SEPARATOR = "_vs_"
POP_RENAME = {}

PARAMS = ["theta:1", "theta:2", "theta:3", "tau:3", "M:1->2", "M:2->1"]
PARAM_NAMES = {
    "theta:1": "theta1",
    "theta:2": "theta2",
    "theta:3": "theta3",
    "tau:3": "tau3",
    "M:1->2": "M1to2",
    "M:2->1": "M2to1",
}

HEADER_RE = re.compile(r"^\s*param\s+mean\s+median\s+", re.I)
PARAM_RE = re.compile(r"^(theta:[1-3]|tau:3|M:[12]->[12])\s+(.+)$")
CHAIN_RE = re.compile(r"(?P<pair>.+)_chain(?P<chain>[0-9]+)\.txt$")
REP_RE = re.compile(r"rep(?:licate)?[_-]?([0-9]+)", re.I)


def rename_pop(label):
    return POP_RENAME.get(label, label)


def clean_column(label):
    return (
        label.lower()
        .replace("(", "")
        .replace(")", "")
        .replace("%", "")
        .replace(".", "")
        .replace("_", "")
    )


def find_column(header, names, fallback):
    cleaned = [clean_column(x) for x in header]
    for name in names:
        target = clean_column(name)
        if target in cleaned:
            return cleaned.index(target)
    return fallback


def parse_run_output(path):
    lines = path.read_text(errors="ignore").splitlines()[-500:]
    header_i = None

    for i in range(len(lines) - 1, -1, -1):
        if HEADER_RE.match(lines[i]):
            header_i = i
            break

    if header_i is None:
        return {}

    header = lines[header_i].split()[1:]
    mean_i = find_column(header, ["mean"], 0)
    median_i = find_column(header, ["median"], 1)
    hpd25_i = find_column(header, ["2.5%HPD", "HPD2.5", "2.5"], 7)
    hpd975_i = find_column(header, ["97.5%HPD", "HPD97.5", "97.5"], 8)
    ess_i = find_column(header, ["ESS"], 9)

    parsed = {}
    for line in lines[header_i + 1 : header_i + 121]:
        if not line.strip() and parsed:
            break
        if line.lstrip().startswith("lnL"):
            break

        match = PARAM_RE.match(line.strip())
        if not match:
            continue

        param = match.group(1)
        values = match.group(2).split()
        if param not in PARAMS:
            continue

        name = PARAM_NAMES[param]
        parsed[f"{name}_mean"] = values[mean_i]
        parsed[f"{name}_median"] = values[median_i]
        parsed[f"{name}_hpd2.5"] = values[hpd25_i]
        parsed[f"{name}_hpd97.5"] = values[hpd975_i]
        parsed[f"{name}_ESS"] = values[ess_i]

    return parsed


def infer_metadata(path):
    chain_match = CHAIN_RE.match(path.name)
    if chain_match is None:
        return None

    pair = chain_match.group("pair").removesuffix("_mscm")
    chain = chain_match.group("chain")
    if PAIR_SEPARATOR not in pair:
        return None

    pop1, pop2 = pair.split(PAIR_SEPARATOR, 1)
    rep = "1"
    for part in path.parts:
        rep_match = REP_RE.search(part)
        if rep_match:
            rep = rep_match.group(1)

    return {
        "species": DATASET_LABEL,
        "pop1": rename_pop(pop1),
        "pop2": rename_pop(pop2),
        "rep": rep,
        "chain": chain,
    }


def main():
    fields = ["species", "pop1", "pop2", "rep", "chain"]
    for param in PARAMS:
        name = PARAM_NAMES[param]
        fields += [
            f"{name}_mean",
            f"{name}_median",
            f"{name}_hpd2.5",
            f"{name}_hpd97.5",
            f"{name}_ESS",
        ]

    rows = []
    for path in sorted(RUN_ROOT.rglob("*_chain*.txt")):
        metadata = infer_metadata(path)
        if metadata is None:
            continue

        row = {field: "" for field in fields}
        row.update(metadata)
        row.update(parse_run_output(path))
        rows.append(row)

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
