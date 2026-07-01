#!/usr/bin/env bash
set -euo pipefail

# Purpose: rarefy two groups to equal sample size and count private/shared derived alleles.
# Input: polarized PLINK files and a two-column sample-group TSV.
# Output: one derived-allele sharing SFS table per replicate.
# Software: PLINK, Python, awk

PLINK_PREFIX="<path/to/polarized_plink_prefix_without_extension>"
GROUP_LABEL_TSV="<path/to/sample_group_labels.tsv>"
OUTPUT_DIR="<path/to/allele_sharing_rarefaction_output_dir>"
DATASET_LABEL="<dataset_label>"
CHR_SET="<species_chr_set>"
REP_START=1
REP_END=100

FOCAL_GROUP="<focal_group_label>"
COMPARISON_GROUP="<comparison_group_label>"

mkdir -p "$OUTPUT_DIR"

TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p "$TMP_DIR"

FAM_MAP="${TMP_DIR}/fam_iid_to_fid.tsv"
FAM_IIDS="${TMP_DIR}/fam_iids.txt"
FOCAL_IIDS="${TMP_DIR}/${FOCAL_GROUP}.iids"
COMPARISON_IIDS="${TMP_DIR}/${COMPARISON_GROUP}.iids"

awk 'BEGIN {OFS="\t"} {print $2, $1}' "${PLINK_PREFIX}.fam" > "$FAM_MAP"
cut -f1 "$FAM_MAP" > "$FAM_IIDS"

awk -F'\t' -v group="$FOCAL_GROUP" '$2 == group {print $1}' "$GROUP_LABEL_TSV" \
| sort -u \
| awk 'NR==FNR {x[$1]=1; next} $1 in x' "$FAM_IIDS" - \
> "$FOCAL_IIDS"

awk -F'\t' -v group="$COMPARISON_GROUP" '$2 == group {print $1}' "$GROUP_LABEL_TSV" \
| sort -u \
| awk 'NR==FNR {x[$1]=1; next} $1 in x' "$FAM_IIDS" - \
> "$COMPARISON_IIDS"

N_FOCAL="$(wc -l < "$FOCAL_IIDS" | awk '{print $1}')"
N_COMPARISON="$(wc -l < "$COMPARISON_IIDS" | awk '{print $1}')"
N_MIN="$N_FOCAL"
[[ "$N_COMPARISON" -lt "$N_MIN" ]] && N_MIN="$N_COMPARISON"

sample_ids() {
  local infile="$1"
  local n="$2"
  local seed="$3"
  local outfile="$4"

  INFILE="$infile" N="$n" SEED="$seed" OUTFILE="$outfile" python3 - <<'PY'
import os, random
ids = [x.strip() for x in open(os.environ["INFILE"]) if x.strip()]
random.seed(os.environ["SEED"])
chosen = random.sample(ids, int(os.environ["N"]))
open(os.environ["OUTFILE"], "w").write("\n".join(chosen) + "\n")
PY
}

make_keep2() {
  local iid_file="$1"
  local keep2_file="$2"

  awk '
    NR == FNR {fid[$1] = $2; next}
    $1 in fid {print fid[$1], $1}
  ' "$FAM_MAP" "$iid_file" > "$keep2_file"
}

freq_counts() {
  local group="$1"
  local keep2="$2"
  local prefix="$3"
  local out_c1="$4"

  plink \
    --bfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --chr-set "$CHR_SET" \
    --keep "$keep2" \
    --freq counts \
    --out "$prefix" >/dev/null

  awk 'BEGIN {FS="[[:space:]]+"; OFS="\t"} NR > 1 {print $2, $5}' "${prefix}.frq.counts" \
  | sort -t $'\t' -k1,1 \
  > "$out_c1"
}

for REP in $(seq "$REP_START" "$REP_END"); do
  REP_DIR="${OUTPUT_DIR}/rep${REP}"
  mkdir -p "$REP_DIR"

  FOCAL_KEEP_IIDS="${TMP_DIR}/${FOCAL_GROUP}.rep${REP}.iids"
  COMPARISON_KEEP_IIDS="${TMP_DIR}/${COMPARISON_GROUP}.rep${REP}.iids"
  FOCAL_KEEP2="${TMP_DIR}/${FOCAL_GROUP}.rep${REP}.keep2"
  COMPARISON_KEEP2="${TMP_DIR}/${COMPARISON_GROUP}.rep${REP}.keep2"
  FOCAL_C1="${TMP_DIR}/${FOCAL_GROUP}.rep${REP}.C1.tsv"
  COMPARISON_C1="${TMP_DIR}/${COMPARISON_GROUP}.rep${REP}.C1.tsv"
  JOINED="${TMP_DIR}/rep${REP}.joined.tsv"
  FINAL="${OUTPUT_DIR}/derived_allele_sharing_SFS_rep${REP}.tsv"

  sample_ids "$FOCAL_IIDS" "$N_MIN" "${DATASET_LABEL}:${FOCAL_GROUP}:${REP}" "$FOCAL_KEEP_IIDS"
  sample_ids "$COMPARISON_IIDS" "$N_MIN" "${DATASET_LABEL}:${COMPARISON_GROUP}:${REP}" "$COMPARISON_KEEP_IIDS"
  make_keep2 "$FOCAL_KEEP_IIDS" "$FOCAL_KEEP2"
  make_keep2 "$COMPARISON_KEEP_IIDS" "$COMPARISON_KEEP2"

  freq_counts "$FOCAL_GROUP" "$FOCAL_KEEP2" "${TMP_DIR}/${FOCAL_GROUP}.rep${REP}" "$FOCAL_C1"
  freq_counts "$COMPARISON_GROUP" "$COMPARISON_KEEP2" "${TMP_DIR}/${COMPARISON_GROUP}.rep${REP}" "$COMPARISON_C1"

  join -t $'\t' -a 1 -a 2 -e 0 -o 0,1.2,2.2 "$FOCAL_C1" "$COMPARISON_C1" > "$JOINED"

  printf 'dataset\trep\tgroup\tac_derived\tcount_sites\n' > "$FINAL"
  awk -v dataset="$DATASET_LABEL" -v rep="$REP" -v focal="$FOCAL_GROUP" -v comparison="$COMPARISON_GROUP" '
    BEGIN {FS=OFS="\t"}
    {
      a = $2 + 0
      b = $3 + 0
      if (a > 0 && b == 0) sfs[focal "_uniq", a]++
      else if (b > 0 && a == 0) sfs[comparison "_uniq", b]++
      else if (a > 0 && b > 0) sfs["Shared", a + b]++
    }
    END {
      for (key in sfs) {
        split(key, x, SUBSEP)
        print dataset, rep, x[1], x[2], sfs[key]
      }
    }
  ' "$JOINED" | sort -t $'\t' -k3,3 -k4,4n >> "$FINAL"

  printf 'dataset\trep\t%s_n\t%s_n\tdownsample_n\n' "$FOCAL_GROUP" "$COMPARISON_GROUP" > "${REP_DIR}/group_sizes.tsv"
  printf '%s\t%s\t%s\t%s\t%s\n' "$DATASET_LABEL" "$REP" "$N_FOCAL" "$N_COMPARISON" "$N_MIN" >> "${REP_DIR}/group_sizes.tsv"
done
