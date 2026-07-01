#!/usr/bin/env python3

# Purpose: parse pairwise BPP MSC-I run outputs into one estimates table.
# Input: BPP output text files named like <pop1>_vs_<pop2>_msci_rep1_chain1.txt.
# Output: one TSV with posterior summaries for theta, tau, and phi parameters.
# Software: Python

from pathlib import Path
import csv
import re

RUN_ROOT = Path("<path/to/pairwise_msci_run_outputs>")
OUTPUT_TSV = Path("<path/to/msci_estimates.tsv>")
DATASET_LABEL = "<dataset_label>"

PAIR_SEPARATOR = "_vs_"
POP_RENAME = {}

PARAMS = [
    "theta:1",
    "theta:2",
    "theta:3",
    "theta:4",
    "theta:5",
    "tau:3",
    "tau:4",
    "tau:5",
    "phi:6<-5",
    "phi:7<-4",
]

PARAM_NAMES = {
    "theta:1": "theta1",
    "theta:2": "theta2",
    "theta:3": "theta3",
    "theta:4": "theta4",
    "theta:5": "theta5",
    "tau:3": "tau3",
    "tau:4": "tau4",
    "tau:5": "tau5",
    "phi:6<-5": "phi6to5",
    "phi:7<-4": "phi7to4",
}

HEADER_RE = re.compile(r"^\s*param\s+mean\s+median\s+", re.I)
PARAM_RE = re.compile(r"^(theta:[1-5]|tau:[3-5]|phi:[67]<-[45])\s+(.+)$")
FILE_RE = re.compile(r"(?P<pair>.+)_msci_rep(?P<rep>[0-9]+)_chain(?P<chain>[0-9]+)\.txt$")


def rename_pop(label):
    return POP_RENAME.get(label, label)


def clean_column(label):
    return re.sub(r"[^a-z0-9.%]+", "", label.lower())


def find_column(header, names, fallback):
    cleaned = [clean_column(x) for x in header]
    for name in names:
        target = clean_column(name)
        if target in cleaned:
            return cleaned.index(target)
    return fallback


def value(tokens, index):
    return tokens[index] if index < len(tokens) else "NA"


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
    hpd25_i = find_column(header, ["2.5%HPD", "2.5HPD", "HPD2.5"], 7)
    hpd975_i = find_column(header, ["97.5%HPD", "97.5HPD", "HPD97.5"], 8)
    ess_i = find_column(header, ["ESS*", "ESS", "Neff", "Eff*"], 9)

    parsed = {}
    for line in lines[header_i + 1 : header_i + 121]:
        match = PARAM_RE.match(line.strip())
        if not match:
            if parsed:
                break
            continue

        param = match.group(1)
        tokens = match.group(2).split()
        name = PARAM_NAMES[param]
        parsed[f"{name}_mean"] = value(tokens, mean_i)
        parsed[f"{name}_median"] = value(tokens, median_i)
        parsed[f"{name}_hpd2.5"] = value(tokens, hpd25_i)
        parsed[f"{name}_hpd97.5"] = value(tokens, hpd975_i)
        parsed[f"{name}_ESS"] = value(tokens, ess_i)

    return parsed


def infer_metadata(path):
    match = FILE_RE.match(path.name)
    if match is None:
        return None

    pair = match.group("pair")
    if PAIR_SEPARATOR not in pair:
        return None

    pop1, pop2 = pair.split(PAIR_SEPARATOR, 1)
    return {
        "species": DATASET_LABEL,
        "pop1": rename_pop(pop1),
        "pop2": rename_pop(pop2),
        "rep": f"rep{match.group('rep')}",
        "chain": match.group("chain"),
    }


def main():
    fieldnames = ["species", "pop1", "pop2", "rep", "chain"]
    for param in PARAMS:
        name = PARAM_NAMES[param]
        fieldnames += [
            f"{name}_mean",
            f"{name}_median",
            f"{name}_hpd2.5",
            f"{name}_hpd97.5",
            f"{name}_ESS",
        ]

    rows = []
    for path in sorted(RUN_ROOT.rglob("*_msci_rep*_chain*.txt")):
        metadata = infer_metadata(path)
        if metadata is None:
            continue

        row = {field: "NA" for field in fieldnames}
        row.update(metadata)
        row.update(parse_run_output(path))
        rows.append(row)

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
