#!/usr/bin/env bash
# SINEderella-route vs COSEG on a hand-labelled copy set. Run on KIT.
#   run_compare.sh SPECIES [WORKDIR]       SPECIES = saq | ccr | teu | dmo  (POS__SPECIES__* in aln_c)
# Steps: labelled copies -> reference = consensus of ALL copies (no subfamily is favoured) ->
# COSEG input (drop = Price's rule, keep = all copies) -> COSEG over -m sweep ->
# SubFam chunks on the same copies -> score everything against the owner's groups.
set -euo pipefail
SP=${1:?species}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
W=${2:-/data/V/toki/coseg_cmp/$SP}
ALN=${ALN:-/data/W/toki/SINE_disc/aln_c}
COSEG=${COSEG:-/data/V/toki/coseg/src}
SUBFAM=${SUBFAM:-$HERE/../../SubFam}
PY=${PYTHON:-python3.12}
THREADS=${THREADS:-16}
mkdir -p "$W"
cd "$W" || exit 1

$PY "$HERE/prep_labeled.py" "$ALN" "$SP" copies
exec bash "$HERE/run_core.sh" "$W"
