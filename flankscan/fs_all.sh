#!/usr/bin/env bash
# fs_all.sh RUN_DIR OUT_DIR [THREADS=16] [MINK=50] [FLANK=1000] [PARTNERS.fa]
# MINKC (default 25): minimum copies carrying the flank unit for a chain extension (stage 6c); ends are checked on the 60 best copies only.
# (PARTNERS.fa: optional partner library for stage 3, from fs0_partnerlib.sh)
# All flankscan stages (1-7, 6b) on one SINEderella run, timed; writes OUT_DIR/DONE (or FAILED) at the end.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN=${1:?RUN_DIR}; OUT=${2:?OUT_DIR}; T=${3:-16}; MINK=${4:-50}; F=${5:-1000}; PART=${6:-}
mkdir -p "$OUT"; rm -f "$OUT/DONE" "$OUT/FAILED"
step() { local t0=$SECONDS; echo "== $1 $(date +%T)" >&2; shift
         if "$@"; then echo "   ok, $((SECONDS - t0)) s" >&2; else echo "   FAILED" >&2; touch "$OUT/FAILED"; exit 1; fi; }
step fs1 bash "$HERE/fs1_extract.sh"   "$RUN" "$OUT" "$F"
step fs2 bash "$HERE/fs2_trf.sh"       "$OUT" "$T"
step fs3 bash "$HERE/fs3_partners.sh"  "$OUT" "$RUN/consensuses.clean.fa" "$T" $PART
step fs4 bash "$HERE/fs4_junctions.sh" "$OUT" "$RUN" "$MINK"
step fs5 bash "$HERE/fs5_build.sh"     "$OUT" "$RUN" 60 "$MINK" "$T"
step fs6b bash "$HERE/fs6b_ends.sh"    "$OUT" "$RUN" "$T"
step fs6c bash "$HERE/fs6c_chain.sh"   "$OUT" "$RUN" "$T" "${MINKC:-25}"
step fs6 bash "$HERE/fs6_reassign.sh"  "$OUT" "$RUN" "$T"
step fs7 bash "$HERE/fs7_hierarchy.sh" "$OUT"
touch "$OUT/DONE"
