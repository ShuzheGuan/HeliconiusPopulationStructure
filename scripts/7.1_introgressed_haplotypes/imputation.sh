#!/usr/bin/env bash
set -euo pipefail

# Purpose: filter variants by missingness, then phase/impute with Beagle.
# Input: one VCF, a genetic map, and the Beagle jar.
# Output: missingness-filtered VCF and Beagle-imputed VCF.
# Software: bcftools, tabix, Beagle, Java

INPUT_VCF="<path/to/input.vcf.gz>"
FILTERED_VCF="<path/to/output_miss10.vcf.gz>"
IMPUTED_PREFIX="<path/to/output_miss10_impute_prefix>"
GENETIC_MAP="<path/to/genetic_map.txt>"
BEAGLE_JAR="<path/to/beagle.jar>"
MAX_MISSING=0.1
THREADS=32
JAVA_MEM="95g"

bcftools +fill-tags "$INPUT_VCF" -Ou -- -t F_MISSING \
  | bcftools view -Oz -o "$FILTERED_VCF" -i "F_MISSING<=${MAX_MISSING}"

tabix -p vcf "$FILTERED_VCF"

java -Xmx"$JAVA_MEM" -jar "$BEAGLE_JAR" \
  gt="$FILTERED_VCF" \
  map="$GENETIC_MAP" \
  out="$IMPUTED_PREFIX" \
  nthreads="$THREADS"

tabix -p vcf "${IMPUTED_PREFIX}.vcf.gz"
