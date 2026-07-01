#!/usr/bin/env bash
set -euo pipefail

# Purpose: prepare BPP A01 big-tree topology runs from raw BPP loci.
# Input: raw BPP loci, sample-population TSV, sample-sex TSV, and focal population list.
# Output: blocked BPP sequence files, Imap.tsv, and A01 control files.
# Software: bash, Python, BPP

RAW_BPP_ROOT="<path/to/raw_bpp_locus_root>"
SAMPLE_POP_TSV="<path/to/sample_population_labels.tsv>"
SAMPLE_SEX_TSV="<path/to/sample_sex.tsv>"
POPULATION_LIST="<path/to/focal_population_list.txt>"
OUTPUT_DIR="<path/to/bigtree_topology_output_dir>"

EXCLUDE_DIR=""

MALE_LABEL="m"
MALES_PER_POP=2
BLOCK_SIZE=100
MAX_MISSING=0.20
N_CHAINS=5
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=100000

THETAPRIOR="2.1 0.03"
TAUPRIOR="2.1 0.03"

SPECIES_LIST="<space_separated_population_labels_in_tree_order>"
STARTING_TREE="(<newick_tree_with_population_labels>);"

if [[ $# -gt 0 ]]; then
  source "$1"
fi

mkdir -p "$OUTPUT_DIR"

export RAW_BPP_ROOT SAMPLE_POP_TSV SAMPLE_SEX_TSV POPULATION_LIST OUTPUT_DIR
export EXCLUDE_DIR MALE_LABEL
export MALES_PER_POP BLOCK_SIZE MAX_MISSING N_CHAINS THREADS
export BURNIN SAMPFREQ NSAMPLE THETAPRIOR TAUPRIOR
export SPECIES_LIST STARTING_TREE

python3 - <<'PY'
from pathlib import Path
import os
import re

raw_root = Path(os.environ["RAW_BPP_ROOT"])
sample_pop_tsv = Path(os.environ["SAMPLE_POP_TSV"])
sample_sex_tsv = Path(os.environ["SAMPLE_SEX_TSV"])
population_list = Path(os.environ["POPULATION_LIST"])
output_dir = Path(os.environ["OUTPUT_DIR"])
exclude_dir = os.environ["EXCLUDE_DIR"]
male_label = os.environ["MALE_LABEL"]

males_per_pop = int(os.environ["MALES_PER_POP"])
block_size = int(os.environ["BLOCK_SIZE"])
max_missing = float(os.environ["MAX_MISSING"])
n_chains = int(os.environ["N_CHAINS"])
threads = int(os.environ["THREADS"])

burnin = os.environ["BURNIN"]
sampfreq = os.environ["SAMPFREQ"]
nsample = os.environ["NSAMPLE"]
thetaprior = os.environ["THETAPRIOR"]
tauprior = os.environ["TAUPRIOR"]
species_list = os.environ["SPECIES_LIST"].split()
starting_tree = os.environ["STARTING_TREE"]


def read_two_column_tsv(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            if parts[0] in {"sample", "sample_id", "id"}:
                continue
            rows.append((parts[0], parts[1]))
    return rows


def select_keep_samples():
    focal_pops = [
        line.strip()
        for line in population_list.open()
        if line.strip() and not line.startswith("#")
    ]
    focal_set = set(focal_pops)
    sex = dict(read_two_column_tsv(sample_sex_tsv))

    keep_rows = []
    seen = {pop: 0 for pop in focal_pops}
    for sample_id, pop in read_two_column_tsv(sample_pop_tsv):
        if pop not in focal_set:
            continue
        if sex.get(sample_id) != male_label:
            continue
        if seen[pop] >= males_per_pop:
            continue
        seen[pop] += 1
        keep_rows.append((sample_id, pop))

    return keep_rows


def contig_dirs():
    dirs = [
        path
        for path in sorted(raw_root.iterdir())
        if path.is_dir() and not (exclude_dir and path.name == exclude_dir)
    ]
    return dirs if dirs else [raw_root]


def parse_locus(path, keep_ids):
    with path.open() as handle:
        header = handle.readline().strip()
        parts = header.split()
        seq_len = int(parts[1])

        seqs = []
        missing = 0
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
            missing += seq.upper().count("N")
            seqs.append((sample_id, seq))

    if not seqs:
        return None
    if missing / (seq_len * len(seqs)) > max_missing:
        return None

    lines = [f"{len(seqs)} {seq_len}"]
    lines.extend(f"^{sample_id}  {seq}" for sample_id, seq in seqs)
    return "\n".join(lines) + "\n\n"


def write_block(contig_out, contig_name, block_index, blocks):
    block_path = contig_out / f"{contig_name}_block{block_index}.bpp.txt"
    block_path.write_text("".join(blocks))
    return block_path


def write_controls(contig_out, block_files, keep_rows):
    counts = {pop: 0 for pop in species_list}
    for _, pop in keep_rows:
        counts[pop] = counts.get(pop, 0) + 1

    imap_path = contig_out / "Imap.tsv"
    with imap_path.open("w") as out:
        for sample_id, pop in keep_rows:
            out.write(f"{sample_id}\t{pop}\n")

    nspecies = len(species_list)
    count_line = " ".join(str(counts.get(pop, 0)) for pop in species_list)
    phase_flags = " ".join(["1"] * nspecies)

    for block_path in block_files:
        nloci = sum(
            1
            for line in block_path.open()
            if re.match(r"^[0-9]+\s+[0-9]+", line)
        )
        base = block_path.stem

        for chain in range(1, n_chains + 1):
            job_name = f"{base}_chain{chain}"
            ctl_path = contig_out / f"{job_name}.ctl"
            ctl_path.write_text(
                f"""seed = -1

seqfile = {block_path.name}
Imapfile = {imap_path.name}

speciesdelimitation = 0
speciestree = 1
speciesmodelprior = 1

nloci = {nloci}

species&tree = {nspecies} {" ".join(species_list)}
 {count_line}
 {starting_tree}

phase = {phase_flags}

usedata = 1
cleandata = 0
clock = 1

thetaprior = {thetaprior}
tauprior = {tauprior}

finetune = 1
burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}

jobname = {job_name}
Threads = {threads} 1 1
print = 1 0 0 0
"""
            )


keep_rows = select_keep_samples()
keep_ids = {sample_id for sample_id, _ in keep_rows}

for contig_dir in contig_dirs():
    contig_name = contig_dir.name
    contig_out = output_dir / contig_name
    contig_out.mkdir(parents=True, exist_ok=True)

    block_files = []
    buffer = []
    block_index = 1

    for locus_path in sorted(contig_dir.glob("*.bpp.txt")):
        block = parse_locus(locus_path, keep_ids)
        if block is None:
            continue

        buffer.append(block)
        if len(buffer) == block_size:
            block_files.append(write_block(contig_out, contig_name, block_index, buffer))
            buffer = []
            block_index += 1

    if buffer:
        block_files.append(write_block(contig_out, contig_name, block_index, buffer))

    if block_files:
        write_controls(contig_out, block_files, keep_rows)
PY
