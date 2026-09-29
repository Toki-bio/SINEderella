#!/usr/bin/env bash
# fs_all.sh RUN_DIR OUT_DIR [THREADS=16] [MINK=50] [FLANK=1000]
# All four flankscan stages on one SINEderella run, timed; writes OUT_DIR/DONE (or FAILED) at the end.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN=${1:?RUN_DIR}; OUT=${2:?OUT_DIR}; T=${3:-16}; MINK=${4:-50}; F=${5:-1000}
mkdir -p "$OUT"; rm -f "$OUT/DONE" "$OUT/FAILED"
step() { local t0=$SECONDS; echo "== $1 $(date +%T)" >&2; shift
         if "$@"; then echo "   ok, $((SECONDS - t0)) s" >&2; else echo "   FAILED" >&2; touch "$OUT/FAILED"; exit 1; fi; }
step fs1 bash "$HERE/fs1_extract.sh"   "$RUN" "$OUT" "$F"
step fs2 bash "$HERE/fs2_trf.sh"       "$OUT" "$T"
step fs3 bash "$HERE/fs3_partners.sh"  "$OUT" "$RUN/consensuses.clean.fa" "$T"
step fs4 bash "$HERE/fs4_junctions.sh" "$OUT" "$RUN" "$MINK"
touch "$OUT/DONE"
