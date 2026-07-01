#!/usr/bin/env python3

# Purpose: apply the Savage-Dickey test to pairwise BPP MSC-I MCMC traces.
# Input: BPP *.mcmc.txt files and a pair metadata TSV with parapatry/sympatry labels.
# Output: one TSV with posterior/prior probabilities and BF01 per chain.
# Software: Python

from pathlib import Path
import csv
import re

MCMC_ROOT = Path("<path/to/pairwise_msci_run_outputs>")
PAIR_METADATA_TSV = Path("<path/to/pair_metadata.tsv>")
OUTPUT_TSV = Path("<path/to/savage_dickey_tests.tsv>")
DATASET_LABEL = "<dataset_label>"

THRESHOLD = 0.002886
PRIOR_PROBABILITY = (THRESHOLD * THRESHOLD * 2.0) + THRESHOLD
PAIR_SEPARATOR = "_vs_"
POP_RENAME = {}

FILE_RE = re.compile(r"(?P<pair>.+)_msci_rep(?P<rep>[0-9]+)_chain(?P<chain>[0-9]+)\.mcmc\.txt$")


def rename_pop(label):
    return POP_RENAME.get(label, label)


def load_pair_metadata(path):
    labels = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species = row["species"]
            pop1 = row["pop1"]
            pop2 = row["pop2"]
            label = row["parapatry"]
            labels[(species, pop1, pop2)] = label
            labels[(species, pop2, pop1)] = label
    return labels


def find_column(header, token):
    for i, value in enumerate(header):
        if token in value:
            return i
    return None


def process_mcmc_file(path):
    n1 = n2 = n3 = n_union = n_total = 0
    t4_i = t3_i = p6_i = p7_i = None

    with path.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if not parts:
                continue

            if parts[0] == "Gen":
                t4_i = find_column(parts, "tau:4")
                t3_i = find_column(parts, "tau:3")
                p6_i = find_column(parts, "phi:6")
                p7_i = find_column(parts, "phi:7")
                continue

            if None in (t4_i, t3_i, p6_i, p7_i):
                continue

            try:
                tau4 = float(parts[t4_i])
                tau3 = float(parts[t3_i])
                phi6 = float(parts[p6_i])
                phi7 = float(parts[p7_i])
            except (ValueError, IndexError):
                continue

            criterion1 = phi6 < THRESHOLD and phi7 < THRESHOLD
            criterion2 = phi6 > (1 - THRESHOLD) and phi7 > (1 - THRESHOLD)
            criterion3 = tau3 > 0.0 and (tau4 / tau3) > (1 - THRESHOLD)
            in_union = criterion1 or criterion2 or criterion3

            n1 += int(criterion1)
            n2 += int(criterion2)
            n3 += int(criterion3)
            n_union += int(in_union)
            n_total += 1

    if n_total == 0:
        return None

    posterior = n_union / n_total
    bf01 = posterior / PRIOR_PROBABILITY
    return n1, n2, n3, n_union, n_total, posterior, bf01


def infer_metadata(path):
    match = FILE_RE.match(path.name)
    if match is None:
        return None

    pair = match.group("pair")
    pop1, pop2 = pair.split(PAIR_SEPARATOR, 1)
    return {
        "species": DATASET_LABEL,
        "pop1": rename_pop(pop1),
        "pop2": rename_pop(pop2),
        "rep": f"rep{match.group('rep')}",
        "chain": match.group("chain"),
    }


def main():
    pair_metadata = load_pair_metadata(PAIR_METADATA_TSV)

    fieldnames = [
        "species",
        "pop1",
        "pop2",
        "rep",
        "chain",
        "parapatry",
        "n_criterion1",
        "n_criterion2",
        "n_criterion3",
        "n_union",
        "n_total",
        "posterior",
        "prior",
        "BF01",
    ]

    rows = []
    for path in sorted(MCMC_ROOT.rglob("*_msci_rep*_chain*.mcmc.txt")):
        metadata = infer_metadata(path)
        if metadata is None:
            continue

        result = process_mcmc_file(path)
        if result is None:
            continue

        n1, n2, n3, n_union, n_total, posterior, bf01 = result
        row = dict(metadata)
        row["parapatry"] = pair_metadata.get(
            (metadata["species"], metadata["pop1"], metadata["pop2"]), "NA"
        )
        row.update(
            {
                "n_criterion1": n1,
                "n_criterion2": n2,
                "n_criterion3": n3,
                "n_union": n_union,
                "n_total": n_total,
                "posterior": f"{posterior:.6f}",
                "prior": f"{PRIOR_PROBABILITY:.9f}",
                "BF01": f"{bf01:.6f}",
            }
        )
        rows.append(row)

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
