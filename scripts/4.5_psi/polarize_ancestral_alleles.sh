#!/usr/bin/env bash
set -euo pipefail

# Purpose: infer ancestral state from fixed outgroup genotypes and polarize SNPs so A1 is the derived allele.
# Input: merged ingroup/outgroup biallelic VCF and a one-column outgroup sample list.
# Output: derived-allele map, polarized-site VCF, polarized PLINK files, and PLINK .raw genotype file.
# Software: bcftools, PLINK

INPUT_VCF="<path/to/ingroup_outgroup_biallelic.vcf.gz>"
OUTGROUP_SAMPLE_LIST="<path/to/outgroup_sample_ids.txt>"
OUTPUT_DIR="<path/to/polarized_plink_output_dir>"
RAW_OUTPUT_DIR="<path/to/polarized_raw_output_dir>"
OUTPUT_PREFIX="<dataset_label>"

mkdir -p "$OUTPUT_DIR" "$RAW_OUTPUT_DIR"

[[ -s "$INPUT_VCF" ]] || exit 1
[[ -s "$OUTGROUP_SAMPLE_LIST" ]] || exit 1

DERIVED_ALLELE_MAP="${OUTPUT_DIR}/${OUTPUT_PREFIX}_derived_alleles.tsv"
TARGETS="${OUTPUT_DIR}/${OUTPUT_PREFIX}_polarized_targets.tsv"
FILTERED_VCF="${OUTPUT_DIR}/${OUTPUT_PREFIX}_polarized_sites.vcf.gz"
POLARIZED_ID_VCF="${OUTPUT_DIR}/${OUTPUT_PREFIX}_polarized_sites_with_ids.vcf.gz"
POLARIZED_PREFIX="${OUTPUT_DIR}/${OUTPUT_PREFIX}_polarized"
RAW_PREFIX="${RAW_OUTPUT_DIR}/${OUTPUT_PREFIX}_polarized"

bcftools view -S "$OUTGROUP_SAMPLE_LIST" -m2 -M2 -v snps -Ou "$INPUT_VCF" \
| bcftools query -f '%CHROM:%POS:%REF:%ALT\t%REF\t%ALT[\t%GT]\n' \
| awk '
  BEGIN { OFS = "\t" }
  {
    id = $1
    ref = $2
    alt = $3
    n00 = 0
    n11 = 0
    nhet = 0
    nmiss = 0
    ntotal = 0

    for (i = 4; i <= NF; i++) {
      g = $i
      ntotal++
      if (g == "./.") { nmiss++; continue }
      if (g == "0/0") { n00++; continue }
      if (g == "1/1") { n11++; continue }
      nhet++
    }

    nonmiss = ntotal - nmiss
    if (nonmiss == 0) next
    if (nhet > 0) next
    if (n00 > 0 && n11 == 0) {
      print id, alt
    } else if (n11 > 0 && n00 == 0) {
      print id, ref
    }
  }
' > "$DERIVED_ALLELE_MAP"

[[ -s "$DERIVED_ALLELE_MAP" ]] || exit 1

cut -f 1 "$DERIVED_ALLELE_MAP" \
| tr ':' '\t' \
| awk 'BEGIN { OFS = "\t" } { print $1, $2 }' \
> "$TARGETS"

bcftools view -R "$TARGETS" -Oz -o "$FILTERED_VCF" "$INPUT_VCF"
bcftools index -f "$FILTERED_VCF"

bcftools annotate -Oz -o "$POLARIZED_ID_VCF" -I +'%CHROM:%POS:%REF:%ALT' "$FILTERED_VCF"
bcftools index -f "$POLARIZED_ID_VCF"

plink \
  --vcf "$POLARIZED_ID_VCF" \
  --a1-allele "$DERIVED_ALLELE_MAP" 2 1 \
  --allow-extra-chr \
  --make-bed \
  --out "$POLARIZED_PREFIX"

plink \
  --bfile "$POLARIZED_PREFIX" \
  --allow-extra-chr \
  --recode A \
  --out "$RAW_PREFIX"
