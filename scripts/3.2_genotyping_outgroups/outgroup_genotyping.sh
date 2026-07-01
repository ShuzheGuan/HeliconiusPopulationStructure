#!/usr/bin/env bash
set -euo pipefail

# Purpose: genotype outgroups at known ingroup sites and merge with ingroup calls.
# Input: outgroup BAMs, reference, ingroup biallelic-sites VCF, and ingroup VCF.
# Output: filtered merged ingroup/outgroup biallelic VCF.
# Software: GATK, bcftools, tabix

BAM_DIR="<path/to/outgroup_deduplicated_bam_dir>"
GVCF_DIR="<path/to/outgroup_gvcf_dir>"
OUT_DIR="<path/to/outgroup_genotyping_output_dir>"
REFERENCE="<path/to/reference_genome.fasta>"
SITES_VCF="<path/to/ingroup_biallelic_sites.vcf.gz>"
INGROUP_VCF="<path/to/ingroup_biallelic.vcf.gz>"
THREADS=8
INTERVAL_PADDING=100
GQ=30
AD=3

mkdir -p "$GVCF_DIR" "$OUT_DIR"

for BAM in "$BAM_DIR"/*_deduplicated.bam; do
  SAMPLE="$(basename "${BAM%_deduplicated.bam}")"
  gatk HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$BAM" \
    -O "${GVCF_DIR}/${SAMPLE}.g.vcf.gz" \
    --emit-ref-confidence GVCF \
    --output-mode EMIT_ALL_CONFIDENT_SITES \
    --alleles "$SITES_VCF" \
    -L "$SITES_VCF" \
    -ip "$INTERVAL_PADDING" \
    --native-pair-hmm-threads "$THREADS"
  tabix -p vcf "${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
done

GVCF_ARGS=()
for GVCF in "$GVCF_DIR"/*.g.vcf.gz; do
  GVCF_ARGS+=(-V "$GVCF")
done

COMBINED_GVCF="${OUT_DIR}/outgroup_combined.g.vcf.gz"
OUTGROUP_VCF="${OUT_DIR}/outgroup_genotyped_exact_sites.vcf.gz"
OUTGROUP_BIALLELIC="${OUT_DIR}/outgroup_genotyped_exact_sites_biallelic.vcf.gz"
MERGED_VCF="${OUT_DIR}/ingroup_outgroup_merged.vcf.gz"
MERGED_BIALLELIC="${OUT_DIR}/ingroup_outgroup_merged_biallelic_only.vcf.gz"
FILTERED_VCF="${OUT_DIR}/ingroup_outgroup_merged_biallelic_only_GATKfiltering_GQ${GQ}_AD${AD}.vcf.gz"
REGIONS_BED="${OUT_DIR}/outgroup_sites.bed"

gatk CombineGVCFs -R "$REFERENCE" "${GVCF_ARGS[@]}" -O "$COMBINED_GVCF"
gatk GenotypeGVCFs -R "$REFERENCE" -V "$COMBINED_GVCF" -L "$SITES_VCF" -O "$OUTGROUP_VCF"

bcftools view --threads "$THREADS" -v snps -m2 -M2 "$OUTGROUP_VCF" -Oz -o "$OUTGROUP_BIALLELIC"
tabix -p vcf "$OUTGROUP_BIALLELIC"

bcftools query -f'%CHROM\t%POS\n' "$OUTGROUP_BIALLELIC" \
| awk 'BEGIN {OFS="\t"} {print $1, $2 - 1, $2}' \
| sort -k1,1 -k2,2n -u \
> "$REGIONS_BED"

bcftools merge --threads "$THREADS" -R "$REGIONS_BED" "$INGROUP_VCF" "$OUTGROUP_BIALLELIC" -Oz -o "$MERGED_VCF"
tabix -p vcf "$MERGED_VCF"

bcftools view --threads "$THREADS" -v snps -m2 -M2 "$MERGED_VCF" -Oz -o "$MERGED_BIALLELIC"
tabix -p vcf "$MERGED_BIALLELIC"

gatk VariantFiltration \
  -V "$MERGED_BIALLELIC" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  --genotype-filter-expression "GQ < ${GQ}.0" --genotype-filter-name "lowGQ" \
  --genotype-filter-expression "isHet == 1 && (AD[0] < ${AD} || AD[1] < ${AD})" --genotype-filter-name "lowSupport_Het" \
  --genotype-filter-expression "isHomRef == 1 && AD[0] < ${AD}" --genotype-filter-name "lowSupport_HomRef" \
  --genotype-filter-expression "isHomVar == 1 && AD[1] < ${AD}" --genotype-filter-name "lowSupport_HomAlt" \
  --set-filtered-genotype-to-no-call \
  -O "$FILTERED_VCF"
tabix -p vcf "$FILTERED_VCF"
