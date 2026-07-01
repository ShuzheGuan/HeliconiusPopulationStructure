#!/usr/bin/env bash
set -euo pipefail

# Purpose: align paired-end reads to the reference genome with BWA-MEM2 and sort/index BAM files.
# Input: trimmed paired-end FASTQ files, reference genome FASTA, sample ID list, and read-group metadata table.
# Output: sorted BAM files with read-group tags and BAM index files.
# Software: BWA-MEM2 and samtools

FASTQ_DIR="<path/to/trimmed_reads_dir>"
GENOME="<path/to/reference_genome.fasta>"
OUTPUT_DIR="<path/to/aligned_sorted_bam_dir>"
RG_FILE="<path/to/read_group_metadata.csv>"
FILE_LIST="<path/to/sample_id_list.txt>"
THREADS=32
SORT_THREADS=4
SORT_MEMORY="4G"

mkdir -p "$OUTPUT_DIR"

bwa-mem2 index "$GENOME"

while read -r PREFIX; do
  [[ -n "$PREFIX" ]] || continue

  READ1="${FASTQ_DIR}/${PREFIX}_1_trimmed.fastq.gz"
  READ2="${FASTQ_DIR}/${PREFIX}_2_trimmed.fastq.gz"
  OUTPUT_BAM="${OUTPUT_DIR}/${PREFIX}_sorted_RG.bam"

  IFS=',' read -r id run biosample platform flowcell < <(grep "^.*,${PREFIX}," "$RG_FILE")

  RGLB="${biosample}_${run}"
  RGPU="$flowcell"

  bwa-mem2 mem -v 2 -M -t "$THREADS" \
    -R "@RG\tID:${run}\tLB:${RGLB}\tPL:${platform}\tPU:${RGPU}\tSM:${biosample}" \
    "$GENOME" "$READ1" "$READ2" \
  | samtools sort \
      -m "$SORT_MEMORY" \
      -@ "$SORT_THREADS" \
      -O BAM \
      -o "$OUTPUT_BAM"

  samtools quickcheck "$OUTPUT_BAM"
  samtools index "$OUTPUT_BAM"
done < "$FILE_LIST"
