#!/usr/bin/env bash
set -euo pipefail

# Purpose: assemble mitochondrial genomes de novo from trimmed paired-end reads.
# Input: sample list, sample-to-read mapping file, trimmed paired-end FASTQ files, and mitochondrial seed FASTA.
# Output: GetOrganelle assembly directories and renamed mitochondrial FASTA files.
# Software: GetOrganelle

SAMPLE_LIST="<path/to/sample_list.txt>"
MAPFILE="<path/to/sample_to_read_prefix_map.txt>"
READ_DIR="<path/to/trimmed_reads_dir>"
SEED="<path/to/mitochondrial_seed.fasta>"
OUT_DIR="<path/to/getorganelle_output_dir>"
SEQUENCES_DIR="<path/to/mitochondrial_fasta_dir>"
THREADS=16
ORGANELLE_TYPE="animal_mt"

mkdir -p "$OUT_DIR" "$SEQUENCES_DIR"

while read -r SAMPLE_ID; do
  [[ -n "$SAMPLE_ID" ]] || continue

  MERGED_NAME="$(awk -v id="$SAMPLE_ID" '$1 == id {print $2; exit}' "$MAPFILE")"
  [[ -n "$MERGED_NAME" ]] || continue

  READ_PREFIX="$MERGED_NAME"

  READ1="${READ_DIR}/${READ_PREFIX}_1_trimmed.fastq.gz"
  READ2="${READ_DIR}/${READ_PREFIX}_2_trimmed.fastq.gz"
  SAMPLE_OUT="${OUT_DIR}/${SAMPLE_ID}"

  [[ -f "$READ1" && -f "$READ2" ]] || continue

  rm -rf "$SAMPLE_OUT"

  get_organelle_from_reads.py \
    -1 "$READ1" \
    -2 "$READ2" \
    -o "$SAMPLE_OUT" \
    -F "$ORGANELLE_TYPE" \
    -s "$SEED" \
    -t "$THREADS"
done < "$SAMPLE_LIST"

for path in "$OUT_DIR"/*/*graph1.1.path_sequence.fasta; do
  [[ -e "$path" ]] || continue

  sample="$(basename "$(dirname "$path")")"

  if [[ "$path" == *complete.graph1.1.path_sequence.fasta ]]; then
    new_name="${sample}_complete_mito.fasta"
  else
    new_name="${sample}_partial_mito.fasta"
  fi

  cp "$path" "${SEQUENCES_DIR}/${new_name}"
done

for file in "$SEQUENCES_DIR"/*.fasta; do
  [[ -e "$file" ]] || continue

  sample="$(basename "$file" .fasta)"
  tmp="${file}.tmp"

  if [[ "$file" == *_complete_mito.fasta ]]; then
    printf '>%s\n' "$sample" > "$tmp"
    awk 'NR > 1' "$file" >> "$tmp"
    mv "$tmp" "$file"
  elif [[ "$file" == *_partial_mito.fasta ]]; then
    awk -v prefix="$sample" '
      /^>/ {
        i++
        print ">" prefix "_partial_" i
        next
      }
      { print }
    ' "$file" > "$tmp"
    mv "$tmp" "$file"
  fi
done
