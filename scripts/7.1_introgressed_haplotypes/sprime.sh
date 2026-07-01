#!/usr/bin/env bash
set -euo pipefail

# Purpose: identify introgressed haplotypes with SPrime from an imputed VCF.
# Input: imputed VCF, outgroup sample list, genetic map, and SPrime jar.
# Output: SPrime score and log files.
# Software: SPrime, Java

IMPUTED_VCF="<path/to/imputed.vcf.gz>"
OUTGROUP_LIST="<path/to/outgroup_samples.list>"
GENETIC_MAP="<path/to/genetic_map.txt>"
SPRIME_JAR="<path/to/sprime.jar>"
OUTPUT_PREFIX="<path/to/sprime_output_prefix>"
MUTATION_RATE="2.9e-9"
JAVA_MEM="95g"

java -Xmx"$JAVA_MEM" -jar "$SPRIME_JAR" \
  gt="$IMPUTED_VCF" \
  outgroup="$OUTGROUP_LIST" \
  map="$GENETIC_MAP" \
  out="$OUTPUT_PREFIX" \
  mu="$MUTATION_RATE"
