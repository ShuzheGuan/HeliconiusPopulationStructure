#!/usr/bin/env bash
set -euo pipefail

# Purpose: summarize posterior BEAST trees into a consensus/MCC tree.
# Input: BEAST posterior .trees file after burn-in removal/thinning.
# Output: annotated BEAST consensus tree.
# Software: TreeAnnotator

TREEANNOTATOR_BIN="<path/to/treeannotator>"
INPUT_TREES="<path/to/thinned_posterior.trees>"
OUTPUT_CONSENSUS_TREE="<path/to/beast_consensus.tre>"
JAVA_MEMORY="30g"
BURNIN_PERCENT=0

mkdir -p "$(dirname "$OUTPUT_CONSENSUS_TREE")"

[[ -x "$TREEANNOTATOR_BIN" ]] || exit 1
[[ -s "$INPUT_TREES" ]] || exit 1

export JAVA_TOOL_OPTIONS="-Xmx${JAVA_MEMORY}"

"$TREEANNOTATOR_BIN" \
  -burnin "$BURNIN_PERCENT" \
  -file "$INPUT_TREES" \
  "$OUTPUT_CONSENSUS_TREE"
