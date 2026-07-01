#!/usr/bin/env bash
set -euo pipefail

# Purpose: run BEAST on a prepared mitochondrial XML file.
# Input: BEAST XML file containing the substitution model, clock model, tree prior, MCMC settings, and output paths.
# Output: BEAST log and posterior tree files defined by the XML.
# Software: BEAST

BEAST_BIN="<path/to/beast>"
XML_FILE="<path/to/mitochondrial_beast.xml>"
THREADS=32

[[ -x "$BEAST_BIN" ]] || exit 1
[[ -s "$XML_FILE" ]] || exit 1

"$BEAST_BIN" \
  -threads "$THREADS" \
  -overwrite \
  "$XML_FILE"
