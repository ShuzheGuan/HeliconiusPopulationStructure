#!/usr/bin/env bash
set -euo pipefail

# Purpose: prepare and optionally run final BPP MSC-I parameter-estimate models.
# Input: BPP locus files, focal-population list, and sample-population TSVs.
# Output: trimmed loci, concatenated seq files, Imap files, and BPP control files.
# Software: bash, Python, BPP

RAW_BPP_ROOT="<path/to/raw_bpp_locus_root>"
FOCAL_POP_LIST="<path/to/focal_population_list.txt>"
TRIM_MAP_AUTO="<path/to/autosome_sample_population_labels.tsv>"
TRIM_MAP_SEX="<path/to/sex_chromosome_sample_population_labels.tsv>"
IMAP_MAP_AUTO="<path/to/autosome_imap_sample_population_labels.tsv>"
IMAP_MAP_SEX="<path/to/sex_chromosome_imap_sample_population_labels.tsv>"
OUTPUT_ROOT="<path/to/final_model_output_root>"

SEX_CONTIG="<sex_chromosome_contig_label>"
CHR1_CONTIG="<chromosome_1_contig_label>"
TARGET_RANDOM_LOCI=4000
RANDOM_SEED=20260504

MODEL_POPS=("<population_1>" "<population_2>" "<population_3>" "<population_4>")
MODEL_PREFIX=""
MODEL_TREE="(<final_msci_network>);"

DATASETS=(sex_chr chr1 random4k)
NCHAINS=3
THREADS=16
RUN_BPP=false

BURNIN=100000
SAMPFREQ=10
NSAMPLE=50000

THETAPRIOR="2.1 0.03 e"
TAUPRIOR="2.1 0.03"
PHIPRIOR="1 1"

if [[ $# -gt 0 ]]; then
  source "$1"
fi

join_by() {
  local IFS="$1"
  shift
  echo "$*"
}

if [[ -z "$MODEL_PREFIX" ]]; then
  MODEL_PREFIX="$(join_by _vs_ "${MODEL_POPS[@]}")"
fi

MODEL_POPS_JOINED="$(IFS=$'\t'; echo "${MODEL_POPS[*]}")"
DATASETS_JOINED="$(IFS=$'\t'; echo "${DATASETS[*]}")"

export RAW_BPP_ROOT FOCAL_POP_LIST TRIM_MAP_AUTO TRIM_MAP_SEX
export IMAP_MAP_AUTO IMAP_MAP_SEX OUTPUT_ROOT SEX_CONTIG CHR1_CONTIG
export TARGET_RANDOM_LOCI RANDOM_SEED MODEL_PREFIX MODEL_TREE
export MODEL_POPS_JOINED DATASETS_JOINED NCHAINS THREADS
export BURNIN SAMPFREQ NSAMPLE THETAPRIOR TAUPRIOR PHIPRIOR

python3 - <<'PY'
from pathlib import Path
import os
import random
import re

raw_root = Path(os.environ["RAW_BPP_ROOT"])
focal_pop_list = Path(os.environ["FOCAL_POP_LIST"])
trim_map_auto = Path(os.environ["TRIM_MAP_AUTO"])
trim_map_sex = Path(os.environ["TRIM_MAP_SEX"])
imap_map_auto = Path(os.environ["IMAP_MAP_AUTO"])
imap_map_sex = Path(os.environ["IMAP_MAP_SEX"])
output_root = Path(os.environ["OUTPUT_ROOT"])

sex_contig = os.environ["SEX_CONTIG"]
chr1_contig = os.environ["CHR1_CONTIG"]
target_random_loci = int(os.environ["TARGET_RANDOM_LOCI"])
random_seed = int(os.environ["RANDOM_SEED"])

model_prefix = os.environ["MODEL_PREFIX"]
model_tree = os.environ["MODEL_TREE"]
model_pops = os.environ["MODEL_POPS_JOINED"].split("\t")
datasets = os.environ["DATASETS_JOINED"].split("\t")

nchains = int(os.environ["NCHAINS"])
threads = os.environ["THREADS"]
burnin = os.environ["BURNIN"]
sampfreq = os.environ["SAMPFREQ"]
nsample = os.environ["NSAMPLE"]
thetaprior = os.environ["THETAPRIOR"]
tauprior = os.environ["TAUPRIOR"]
phiprior = os.environ["PHIPRIOR"]


def read_populations(path):
    return {line.strip() for line in path.open() if line.strip() and not line.startswith("#")}


def read_sample_pop(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] not in {"sample", "sample_id", "id"}:
                rows.append((parts[0], parts[1]))
    return rows


def map_for_dataset(dataset, auto_map, sex_map):
    return sex_map if dataset == "sex_chr" else auto_map


def selected_loci(dataset, all_loci):
    if dataset == "sex_chr":
        return [path for path in all_loci if sex_contig in str(path)]
    if dataset == "chr1":
        return [path for path in all_loci if chr1_contig in str(path)]

    candidates = [path for path in all_loci if sex_contig not in str(path)]
    if len(candidates) <= target_random_loci:
        return candidates
    return sorted(random.Random(random_seed).sample(candidates, target_random_loci))


def trim_locus(in_path, out_path, keep_ids):
    with in_path.open() as handle:
        header = handle.readline().strip()
        parts = header.split()
        seq_len = int(parts[1])

        kept = []
        for line in handle:
            if not line.startswith("^"):
                continue
            match = re.match(r"^\^(\S+)\s+(.*)$", line.rstrip("\n"))
            if match and match.group(1) in keep_ids:
                kept.append((match.group(1), match.group(2).replace(" ", "")))

    if not kept:
        return

    with out_path.open("w") as out:
        out.write(f"{len(kept)} {seq_len}\n")
        for sample_id, seq in kept:
            out.write(f"^{sample_id}  {seq}\n")


def write_imap(path, map_path):
    rows = read_sample_pop(map_path)
    counts = {pop: 0 for pop in model_pops}

    with path.open("w") as out:
        for pop in model_pops:
            for sample_id, sample_pop in rows:
                if sample_pop == pop:
                    out.write(f"{sample_id}\t{pop}\n")
                    counts[pop] += 1

    return counts


def count_loci(path):
    return sum(1 for line in path.open() if re.match(r"^[0-9]+\s+[0-9]+", line))


def write_ctl(pair_dir, dataset, concat_path, imap_path, counts, nloci, chain):
    job_name = f"{model_prefix}_{dataset}_msci_chain{chain}"
    ctl_path = pair_dir / f"{job_name}.ctl"
    count_line = " ".join(str(counts[pop]) for pop in model_pops)
    phase_line = " ".join(["1"] * len(model_pops))

    ctl_path.write_text(
        f"""seed = -1

seqfile = {concat_path.name}
Imapfile = {imap_path.name}

speciesdelimitation = 0
speciestree = 0
usedata = 1
cleandata = 0
clock = 1

nloci = {nloci}

species&tree = {len(model_pops)} {" ".join(model_pops)}
 {count_line}
{model_tree}

phase = {phase_line}

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


focal_pops = read_populations(focal_pop_list)
all_loci = [
    path
    for path in sorted(raw_root.rglob("*.bpp.txt"))
    if not path.name.endswith(".trimmed.bpp.txt")
]

for dataset in datasets:
    trim_map = map_for_dataset(dataset, trim_map_auto, trim_map_sex)
    imap_map = map_for_dataset(dataset, imap_map_auto, imap_map_sex)
    keep_ids = {sample_id for sample_id, pop in read_sample_pop(trim_map) if pop in focal_pops}

    pair_dir = output_root / dataset
    raw_dir = pair_dir / "9_raw_data"
    raw_dir.mkdir(parents=True, exist_ok=True)

    for in_path in selected_loci(dataset, all_loci):
        out_path = raw_dir / f"{in_path.name.removesuffix('.bpp.txt')}.trimmed.bpp.txt"
        trim_locus(in_path, out_path, keep_ids)

    imap_path = pair_dir / f"{model_prefix}.imap"
    counts = write_imap(imap_path, imap_map)

    concat_path = pair_dir / f"{model_prefix}_{dataset}_concat_allloci.txt"
    with concat_path.open("w") as out:
        for path in sorted(raw_dir.glob("*.trimmed.bpp.txt")):
            out.write(path.read_text())

    nloci = count_loci(concat_path)
    for chain in range(1, nchains + 1):
        write_ctl(pair_dir, dataset, concat_path, imap_path, counts, nloci, chain)
PY

if [[ "$RUN_BPP" == "true" ]]; then
  find "$OUTPUT_ROOT" -name "*.ctl" -print | sort | while read -r ctl; do
    (
      cd "$(dirname "$ctl")"
      bpp --no-pin --cfile "$(basename "$ctl")"
    )
  done
fi
