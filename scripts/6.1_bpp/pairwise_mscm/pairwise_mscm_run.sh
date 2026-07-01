#!/usr/bin/env bash
set -euo pipefail

# Purpose: prepare pairwise BPP MSC-M inputs from raw BPP loci.
# Input: raw BPP locus files and a two-column sample-population TSV.
# Output: per-replicate seq.txt, Imap.tsv, and two BPP control files.
# Software: bash, Python, BPP

RAW_BPP_ROOT="<path/to/raw_bpp_locus_root>"
SAMPLE_POP_TSV="<path/to/sample_population_labels.tsv>"
OUTPUT_ROOT="<path/to/pairwise_mscm_output_root>"

POP1="<population_1_label>"
POP2="<population_2_label>"
EXCLUDE_DIR=""

TARGET_LOCI=2000
N_REPS=5
RANDOM_SEED_BASE=20260504
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.03 e"
TAUPRIOR="2.1 0.03"
WPRIOR="2 0.003"

if [[ $# -gt 0 ]]; then
  source "$1"
fi

PAIR_NAME="${POP1}_vs_${POP2}"
PAIR_DIR="${OUTPUT_ROOT}/${PAIR_NAME}"

mkdir -p "$PAIR_DIR"

export RAW_BPP_ROOT SAMPLE_POP_TSV PAIR_DIR POP1 POP2 EXCLUDE_DIR
export TARGET_LOCI N_REPS RANDOM_SEED_BASE THREADS
export BURNIN SAMPFREQ NSAMPLE THETAPRIOR TAUPRIOR WPRIOR

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
n_reps = int(os.environ["N_REPS"])
seed_base = int(os.environ["RANDOM_SEED_BASE"])
threads = int(os.environ["THREADS"])

burnin = os.environ["BURNIN"]
sampfreq = os.environ["SAMPFREQ"]
nsample = os.environ["NSAMPLE"]
thetaprior = os.environ["THETAPRIOR"]
tauprior = os.environ["TAUPRIOR"]
wprior = os.environ["WPRIOR"]
pair_name = f"{pop1}_vs_{pop2}"


def load_samples():
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

    return pop1_ids, pop2_ids


def list_loci():
    loci = []
    for path in sorted(raw_root.rglob("*.bpp.txt")):
        if exclude_dir and exclude_dir in path.parts:
            continue
        if path.stat().st_size > 0:
            loci.append(path)
    return loci


def subset_locus(path, keep_ids):
    with path.open() as handle:
        header = handle.readline().strip()
        header_parts = header.split()
        seq_len = int(header_parts[1])

        kept = []
        for line in handle:
            if not line.startswith("^"):
                continue
            match = re.match(r"^\^(\S+)\s+(.*)$", line.rstrip("\n"))
            if match is None:
                continue
            sample_id = match.group(1)
            if sample_id in keep_ids:
                kept.append((sample_id, match.group(2).replace(" ", "")))

    if not kept:
        return None

    lines = [f"{len(kept)} {seq_len}"]
    lines.extend(f"^{sample_id}  {seq}" for sample_id, seq in kept)
    return "\n".join(lines) + "\n\n"


def write_ctl(rep_dir, seq_name, imap_name, n_loci, n_pop1, n_pop2, chain):
    job_name = f"{pair_name}_chain{chain}"
    ctl_path = rep_dir / f"{job_name}.ctl"

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
({pop1},{pop2});

phase = 1 1

thetaprior = {thetaprior}
tauprior = {tauprior}
wprior = {wprior}

migration = 2
{pop1} {pop2}
{pop2} {pop1}

burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}
finetune = 1

jobname = {job_name}
Threads = {threads} 1 1
print = 1 0 0 0
"""
    )


pop1_ids, pop2_ids = load_samples()
keep_ids = set(pop1_ids) | set(pop2_ids)
locus_files = list_loci()

for rep in range(1, n_reps + 1):
    rep_dir = pair_dir / f"rep{rep}"
    rep_dir.mkdir(parents=True, exist_ok=True)

    imap_path = rep_dir / "Imap.tsv"
    with imap_path.open("w") as out:
        for sample_id in pop1_ids:
            out.write(f"{sample_id}\t{pop1}\n")
        for sample_id in pop2_ids:
            out.write(f"{sample_id}\t{pop2}\n")

    seq_path = rep_dir / "seq.txt"
    accepted = 0
    shuffled = locus_files[:]
    random.Random(seed_base + rep).shuffle(shuffled)

    with seq_path.open("w") as out:
        for locus_path in shuffled:
            block = subset_locus(locus_path, keep_ids)
            if block is None:
                continue
            out.write(block)
            accepted += 1
            if accepted == target_loci:
                break

    for chain in (1, 2):
        write_ctl(
            rep_dir=rep_dir,
            seq_name=seq_path.name,
            imap_name=imap_path.name,
            n_loci=accepted,
            n_pop1=len(pop1_ids),
            n_pop2=len(pop2_ids),
            chain=chain,
        )
PY
