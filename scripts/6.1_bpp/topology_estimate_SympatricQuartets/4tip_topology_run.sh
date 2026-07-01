#!/usr/bin/env bash
set -euo pipefail

# Purpose: prepare BPP A01 four-tip topology runs from raw BPP loci.
# Input: raw CDS/nCDS BPP loci, sample-population TSVs, and focal population list.
# Output: 100-locus concat files, Imap.tsv, and two A01 control files per block.
# Software: bash, Python, BPP

NCDS_BPP_ROOT="<path/to/ncds_bpp_locus_root>"
CDS_BPP_ROOT="<path/to/cds_bpp_locus_root>"
AUTO_SAMPLE_POP_TSV="<path/to/autosome_sample_population_labels.tsv>"
SEX_SAMPLE_POP_TSV="<path/to/sex_chromosome_sample_population_labels.tsv>"
POPULATION_LIST="<path/to/focal_population_list.txt>"
OUTPUT_DIR="<path/to/four_tip_topology_output_dir>"

DATASET_LABEL="<dataset_label>"
SEX_CONTIG="<sex_chromosome_contig_label>"

BLOCK_SIZE=100
N_CHAINS=2
THREADS=16

BURNIN=100000
SAMPFREQ=10
NSAMPLE=100000

THETAPRIOR="2.1 0.03"
TAUPRIOR="2.1 0.03"

SPECIES_LIST="<outgroup_label> <taxon_A_region_1> <taxon_B_region_1> <taxon_A_region_2> <taxon_B_region_2>"
STARTING_TREE="(<outgroup_label>,((<taxon_A_region_1>,<taxon_B_region_1>),(<taxon_A_region_2>,<taxon_B_region_2>)));"

if [[ $# -gt 0 ]]; then
  source "$1"
fi

mkdir -p "$OUTPUT_DIR"

export NCDS_BPP_ROOT CDS_BPP_ROOT AUTO_SAMPLE_POP_TSV SEX_SAMPLE_POP_TSV
export POPULATION_LIST OUTPUT_DIR DATASET_LABEL SEX_CONTIG
export BLOCK_SIZE N_CHAINS THREADS BURNIN SAMPFREQ NSAMPLE
export THETAPRIOR TAUPRIOR SPECIES_LIST STARTING_TREE

python3 - <<'PY'
from pathlib import Path
import os
import re

ncds_root = Path(os.environ["NCDS_BPP_ROOT"])
cds_root = Path(os.environ["CDS_BPP_ROOT"])
auto_sample_pop_tsv = Path(os.environ["AUTO_SAMPLE_POP_TSV"])
sex_sample_pop_tsv = Path(os.environ["SEX_SAMPLE_POP_TSV"])
population_list = Path(os.environ["POPULATION_LIST"])
output_dir = Path(os.environ["OUTPUT_DIR"])
dataset_label = os.environ["DATASET_LABEL"]
sex_contig = os.environ["SEX_CONTIG"]

block_size = int(os.environ["BLOCK_SIZE"])
n_chains = int(os.environ["N_CHAINS"])
threads = int(os.environ["THREADS"])

burnin = os.environ["BURNIN"]
sampfreq = os.environ["SAMPFREQ"]
nsample = os.environ["NSAMPLE"]
thetaprior = os.environ["THETAPRIOR"]
tauprior = os.environ["TAUPRIOR"]
species_list = os.environ["SPECIES_LIST"].split()
starting_tree = os.environ["STARTING_TREE"]


def read_populations(path):
    return [
        line.strip()
        for line in path.open()
        if line.strip() and not line.startswith("#")
    ]


def read_sample_pop(path):
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


def keep_ids_for_map(map_path, allowed_pops):
    return {
        sample_id
        for sample_id, pop in read_sample_pop(map_path)
        if pop in allowed_pops
    }


def contig_dirs(root):
    dirs = [path for path in sorted(root.iterdir()) if path.is_dir()]
    return dirs if dirs else [root]


def subset_locus(path, keep_ids):
    with path.open() as handle:
        header = handle.readline().strip()
        parts = header.split()
        seq_len = int(parts[1])

        seqs = []
        for line in handle:
            if not line.startswith("^"):
                continue
            match = re.match(r"^\^(\S+)\s+(.*)$", line.rstrip("\n"))
            if match is None:
                continue
            sample_id = match.group(1)
            if sample_id in keep_ids:
                seqs.append((sample_id, match.group(2).replace(" ", "")))

    if not seqs:
        return None

    lines = [f"{len(seqs)} {seq_len}"]
    lines.extend(f"^{sample_id}  {seq}" for sample_id, seq in seqs)
    return "\n".join(lines) + "\n\n"


def collect_ids(block_text):
    ids = set()
    for line in block_text.splitlines():
        if line.startswith("^"):
            ids.add(line[1:].split(None, 1)[0])
    return ids


def write_imap(contig_dir, contig, present_ids, allowed_pops):
    map_path = sex_sample_pop_tsv if contig == sex_contig else auto_sample_pop_tsv
    id_to_pop = dict(read_sample_pop(map_path))

    imap_path = contig_dir / "Imap.tsv"
    with imap_path.open("w") as out:
        for sample_id in sorted(present_ids):
            pop = id_to_pop.get(sample_id)
            if pop in allowed_pops:
                out.write(f"{sample_id}\t{pop}\n")

    return imap_path


def count_line(imap_path):
    counts = {pop: 0 for pop in species_list}
    for sample_id, pop in read_sample_pop(imap_path):
        counts[pop] = counts.get(pop, 0) + 1
    return " ".join(str(counts.get(pop, 0)) for pop in species_list)


def count_loci(path):
    return sum(
        1
        for line in path.open()
        if re.match(r"^[0-9]+\s+[0-9]+", line)
    )


def write_ctl(contig_dir, concat_path, nloci, counts, chain):
    base = concat_path.name.replace(".concat.txt", "")
    job_name = f"{base}_chain{chain}"
    ctl_path = contig_dir / f"{job_name}.ctl"
    nspecies = len(species_list)
    phase_flags = " ".join(["1"] * nspecies)

    ctl_path.write_text(
        f"""seed = -1

seqfile = {concat_path.name}
Imapfile = Imap.tsv

speciesdelimitation = 0
speciestree = 1
speciesmodelprior = 1

nloci = {nloci}

species&tree = {nspecies} {" ".join(species_list)}
 {counts}
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


def process_contig(coding_type, contig_dir, allowed_pops):
    contig = contig_dir.name
    map_path = sex_sample_pop_tsv if contig == sex_contig else auto_sample_pop_tsv
    keep_ids = keep_ids_for_map(map_path, allowed_pops)
    out_dir = output_dir / f"{dataset_label}_{coding_type}_{contig}"
    out_dir.mkdir(parents=True, exist_ok=True)

    buffer = []
    present_ids = set()
    block_index = 1
    concat_files = []

    for locus_path in sorted(contig_dir.glob("*.bpp.txt")):
        block = subset_locus(locus_path, keep_ids)
        if block is None:
            continue

        buffer.append(block)
        present_ids |= collect_ids(block)

        if len(buffer) == block_size:
            concat_path = out_dir / f"{dataset_label}_{coding_type}_{contig}_block{block_index}.concat.txt"
            concat_path.write_text("".join(buffer))
            concat_files.append(concat_path)
            buffer = []
            block_index += 1

    if not concat_files:
        return

    imap_path = write_imap(out_dir, contig, present_ids, allowed_pops)
    counts = count_line(imap_path)

    for concat_path in concat_files:
        nloci = count_loci(concat_path)
        for chain in range(1, n_chains + 1):
            write_ctl(out_dir, concat_path, nloci, counts, chain)


allowed_pops = set(read_populations(population_list))

for root, coding_type in [(ncds_root, "ncds"), (cds_root, "cds")]:
    for contig_dir in contig_dirs(root):
        process_contig(coding_type, contig_dir, allowed_pops)
PY
