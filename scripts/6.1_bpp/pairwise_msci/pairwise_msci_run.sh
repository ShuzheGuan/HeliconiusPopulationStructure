#!/usr/bin/env bash
set -euo pipefail

# Purpose: prepare pairwise BPP MSC-I inputs from raw BPP loci.
# Input: raw BPP locus files and a two-column sample-population TSV.
# Output: imap, replicate concatenated sequence files, and two BPP control files per replicate.
# Software: bash, Python, BPP

RAW_BPP_ROOT="<path/to/raw_bpp_locus_root>"
SAMPLE_POP_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_ROOT="<path/to/pairwise_msci_output_root>"

POP1="<population_1_label>"
POP2="<population_2_label>"
EXCLUDE_DIR=""

TARGET_LOCI=2000
MAX_MISSING=0.10
N_REPS=5
RANDOM_SEED_BASE=20260504
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.015 e"
TAUPRIOR="2.1 0.06"
PHIPRIOR="1 1"

if [[ $# -gt 0 ]]; then
  source "$1"
fi

PAIR_NAME="${POP1}_vs_${POP2}"
PAIR_DIR="${OUTPUT_ROOT}/${PAIR_NAME}"

mkdir -p "$PAIR_DIR"

export RAW_BPP_ROOT SAMPLE_POP_TSV PAIR_DIR POP1 POP2 EXCLUDE_DIR
export TARGET_LOCI MAX_MISSING N_REPS RANDOM_SEED_BASE THREADS
export BURNIN SAMPFREQ NSAMPLE THETAPRIOR TAUPRIOR PHIPRIOR

python3 - <<'PY'
from pathlib import Path
import os
import random
import re

raw_root = Path(os.environ["RAW_BPP_ROOT"])
sample_pop_tsv = Path(os.environ["SAMPLE_POP_TSV"])
pair_dir = Path(os.environ["PAIR_DIR"])
pop1 = os.environ["POP1"]
pop2 = os.environ["POP2"]
exclude_dir = os.environ["EXCLUDE_DIR"]

target_loci = int(os.environ["TARGET_LOCI"])
max_missing = float(os.environ["MAX_MISSING"])
n_reps = int(os.environ["N_REPS"])
seed_base = int(os.environ["RANDOM_SEED_BASE"])
threads = int(os.environ["THREADS"])

burnin = os.environ["BURNIN"]
sampfreq = os.environ["SAMPFREQ"]
nsample = os.environ["NSAMPLE"]
thetaprior = os.environ["THETAPRIOR"]
tauprior = os.environ["TAUPRIOR"]
phiprior = os.environ["PHIPRIOR"]
pair_name = f"{pop1}_vs_{pop2}"


def load_samples(imap_path):
    pop1_ids = []
    pop2_ids = []

    with sample_pop_tsv.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2 or parts[0] in {"sample", "sample_id", "id"}:
                continue
            sample_id, pop = parts[0], parts[1]
            if pop == pop1:
                pop1_ids.append(sample_id)
            elif pop == pop2:
                pop2_ids.append(sample_id)

    with imap_path.open("w") as out:
        for sample_id in pop1_ids:
            out.write(f"{sample_id}\t{pop1}\n")
        for sample_id in pop2_ids:
            out.write(f"{sample_id}\t{pop2}\n")

    return pop1_ids, pop2_ids


def list_loci():
    loci = []
    for path in sorted(raw_root.rglob("*.bpp.txt")):
        if exclude_dir and exclude_dir in path.parts:
            continue
        if path.stat().st_size > 0:
            loci.append(path)
    return loci


def subset_locus(path, keep_ids, pop1_ids, pop2_ids):
    with path.open() as handle:
        header = handle.readline().strip()
        header_parts = header.split()
        seq_len = int(header_parts[1])

        kept = []
        missing_count = 0
        seen_pop1 = 0
        seen_pop2 = 0

        for line in handle:
            if not line.startswith("^"):
                continue
            match = re.match(r"^\^(\S+)\s+(.*)$", line.rstrip("\n"))
            if match is None:
                continue
            sample_id = match.group(1)
            if sample_id not in keep_ids:
                continue
            seq = match.group(2).replace(" ", "")
            kept.append((sample_id, seq))
            missing_count += seq.upper().count("N")
            seen_pop1 += sample_id in pop1_ids
            seen_pop2 += sample_id in pop2_ids

    if not kept or seen_pop1 == 0 or seen_pop2 == 0:
        return None

    if missing_count / (seq_len * len(kept)) > max_missing:
        return None

    lines = [f"{len(kept)} {seq_len}"]
    lines.extend(f"^{sample_id}  {seq}" for sample_id, seq in kept)
    return "\n".join(lines) + "\n\n"


def write_ctl(seq_name, imap_name, n_loci, n_pop1, n_pop2, rep, chain):
    job_name = f"{pair_name}_msci_rep{rep}_chain{chain}"
    ctl_path = pair_dir / f"{job_name}.ctl"
    model_d_tree = f"(({pop1},Y[&phi=0.200000])X,({pop2},X[&phi=0.100000])Y)R;"

    ctl_path.write_text(
        f"""seed = -1

seqfile = {seq_name}
Imapfile = {imap_name}

speciesdelimitation = 0
speciestree = 0
usedata = 1
cleandata = 0
clock = 1
nloci = {n_loci}

species&tree = 2 {pop1} {pop2}
 {n_pop1} {n_pop2}
{model_d_tree}

phase = 1 1

thetaprior = {thetaprior}
tauprior = {tauprior}
phiprior = {phiprior}

burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}
finetune = 1

jobname = {job_name}
Threads = {threads} 1 1
print = 1 0 0 0
"""
    )


imap_path = pair_dir / f"{pair_name}.imap.txt"
pop1_ids, pop2_ids = load_samples(imap_path)
keep_ids = set(pop1_ids) | set(pop2_ids)
locus_files = list_loci()

for rep in range(1, n_reps + 1):
    seq_path = pair_dir / f"{pair_name}_2000loci_concat_rep{rep}.txt"
    accepted = 0
    shuffled = locus_files[:]
    random.Random(seed_base + rep).shuffle(shuffled)

    with seq_path.open("w") as out:
        for locus_path in shuffled:
            block = subset_locus(locus_path, keep_ids, set(pop1_ids), set(pop2_ids))
            if block is None:
                continue
            out.write(block)
            accepted += 1
            if accepted == target_loci:
                break

    for chain in (1, 2):
        write_ctl(
            seq_name=seq_path.name,
            imap_name=imap_path.name,
            n_loci=accepted,
            n_pop1=len(pop1_ids),
            n_pop2=len(pop2_ids),
            rep=rep,
            chain=chain,
        )
PY
