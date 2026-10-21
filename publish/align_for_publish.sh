#!/usr/bin/env bash
# Publish-quality alignments for the HTML report (step7 → border loop → step8a → DISC).
#
# Usage: align_for_publish.sh <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]
#
# Environment:
#   SINEDERELLA_BIN — repo root with step7/step8a (default: parent of publish/)
#   DISC            — SINE-discriminator tree with boundary_justify.py, etc.
#   PEEL_FLAGS, PEEL_ALN_DIR, SKIP_STEP7, SKIP_BORDER_LOOP, SKIP_BORDER_SCAN
#   TRIM_DISPLAY_MODE — occupancy (default) or hybrid
set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]}"
OUT_DIR="${3:-$RUN_ROOT/results/alignments}"

PUBLISH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINEDERELLA_BIN="${SINEDERELLA_BIN:-$(dirname "$PUBLISH_DIR")}"
DISC="${DISC:-}"
PEEL_FLAGS="${PEEL_FLAGS:-}"
PEEL_ALN_DIR="${PEEL_ALN_DIR:-}"

if [[ -z "$DISC" ]]; then
  # disc/ in this repository first (2026-10-21); an explicit DISC or SINE_DISC still wins
  for d in "${SINE_DISC:-}" "$SINEDERELLA_BIN/disc" "/staging/tmp/sinedisc" "$HOME/SINE_discriminator/site" \
           "$HOME/SINE_discriminator"; do
    [[ -n "$d" && -f "$d/boundary_justify.py" ]] && DISC="$d" && break
  done
fi
[[ -n "$DISC" && -f "$DISC/boundary_justify.py" ]] || {
  echo "ERROR: set DISC to SINE-discriminator root (need boundary_justify.py)" >&2
  exit 1
}
# the border scan (needs_border_loop.py, a child python process) imports boundary.py from $DISC; without the export it fell back to
# /staging/tmp/sinedisc, failed with "No module named boundary" and the scan reported nothing flagged (found 2026-10-02)
export DISC

# DRAGEN toolchain only where it exists. On other hosts (therioserver) the unconditional
# prepend put /usr/bin ahead of the active conda env: /usr/bin/python3 has no numpy, so the
# border scan failed with a WARN and every subfamily was reported "none flagged".
[[ -d /staging/conda/envs/bioinfo/bin ]] && \
  export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

# The verdict chain (step 7 -> border loop -> flank verdicts; docs/VERDICTS.md) runs in every run since 2026-10-21
# (publish/verdict_chain.sh, called by the orchestrator). Reused here when it finished after the assignment and after stage 9;
# otherwise (a run made before, or VERDICT_RERUN=1) it runs now.
S2V="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
DONE_V="$S2V/verdict_chain.done"
if [[ "${VERDICT_RERUN:-0}" != "1" && -f "$DONE_V" && "$DONE_V" -nt "$S2V/assignment_full.tsv" \
      && ( ! -e "$RUN_ROOT/results/flank_twins.tsv" || "$DONE_V" -nt "$RUN_ROOT/results/flank_twins.tsv" ) ]]; then
  log "verdicts: from the run ($(tr '\n' ' ' < "$DONE_V"))"
else
  log "verdicts: running publish/verdict_chain.sh"
  SINEDERELLA_BIN="$SINEDERELLA_BIN" bash "$PUBLISH_DIR/verdict_chain.sh" "$RUN_ROOT" "$SPECIES" || {
    echo "ERROR: verdict_chain.sh failed" >&2; exit 1; }
fi

for need in step7_boundary_refine.sh step8a_extract_alignments.sh; do
  [[ -x "$SINEDERELLA_BIN/$need" ]] || {
    echo "ERROR: missing $SINEDERELLA_BIN/$need" >&2; exit 1; }
done

log "orient consensuses.publish.fa (simple-repeat tail at 3 prime, pre-step8a)"
ORIENT_PY="$SINEDERELLA_BIN/tools/orient_consensus_bank.py"
if [[ "${SKIP_ORIENT_BANK:-0}" != "1" && -f "$ORIENT_PY" ]]; then
  python3 "$ORIENT_PY" "$RUN_ROOT/consensuses.publish.fa" \
    --report "$RUN_ROOT/consensuses.publish.orient.tsv"
else
  log "SKIP orient consensuses.publish.fa (SKIP_ORIENT_BANK=1 or tool missing)"
fi

log "step8a: extract publish alignments"
STEP8_OUT="$RUN_ROOT/results/alignments"
mkdir -p "$STEP8_OUT"
# step8a reads the searched bank and assignments; the publish versions are swapped in for it and swapped back
# afterwards. A backup left by an interrupted earlier run holds the searched files: restore it, never overwrite it.
SWAP=("consensuses.clean.fa" "step2/step2_output/assigned.fasta")
for f in "${SWAP[@]}"; do
  if [[ -f "$RUN_ROOT/$f.pre_publish.bak" ]]; then
    log "restoring $f from an interrupted publish"; cp -f "$RUN_ROOT/$f.pre_publish.bak" "$RUN_ROOT/$f"
  fi
  cp -f "$RUN_ROOT/$f" "$RUN_ROOT/$f.pre_publish.bak"
done
restore_swap() { for f in "${SWAP[@]}"; do cp -f "$RUN_ROOT/$f.pre_publish.bak" "$RUN_ROOT/$f"; done; }
trap restore_swap EXIT
cp -f "$RUN_ROOT/consensuses.publish.fa" "$RUN_ROOT/consensuses.clean.fa"
cp -f "$RUN_ROOT/step2/step2_output/assigned.publish.fasta" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta"
DISC="$DISC" "$SINEDERELLA_BIN/step8a_extract_alignments.sh" "$RUN_ROOT" "$SPECIES"
restore_swap; trap - EXIT
for f in "${SWAP[@]}"; do rm -f "$RUN_ROOT/$f.pre_publish.bak"; done

if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  mkdir -p "$OUT_DIR"
  cp -a "$STEP8_OUT"/*.aln.fa "$OUT_DIR/" 2>/dev/null || true
fi

log "DISC: boundary justify (must run before rebuild_consensus_row -- it is what"
log "  establishes the element=upper/flank=lower case convention rebuild_consensus_row"
log "  depends on; before this runs, row case is just the raw genome's native"
log "  soft-masking, unrelated to element boundaries -- confirmed root cause of the"
log "  oma_SINE10 top100 row-0 collapsing to 5bp, 2026-09-09)"
shopt -s nullglob
ORIENT_ALN="$SINEDERELLA_BIN/tools/orient_publish_aln.py"
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/boundary_justify.py" "$f"
done
if [[ "${SKIP_ORIENT_MSA:-0}" != "1" && -f "$ORIENT_ALN" ]]; then
  log "orient published MSAs (tail check on row 1, post-justify)"
  python3 "$ORIENT_ALN" "$STEP8_OUT"/*.aln.fa || true
fi
log "DISC: rebuild consensus row (now runs on correctly-cased rows)"
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/rebuild_consensus_row.py" "$f"
  python3 "$DISC/trim_display_flanks.py" "$f" \
    --mode "${TRIM_DISPLAY_MODE:-occupancy}" || true
done

# After rebuild. Do not run rebuild_consensus_row again: it would rewrite the
# consensus from the copies and put a reverse-complement or a short edge back.
[[ -f "$DISC/correct_published_aln.py" ]] || {
  echo "ERROR: missing $DISC/correct_published_aln.py" >&2; exit 1; }
BANK="$RUN_ROOT/consensuses.publish.fa"
[[ -f "$BANK" ]] || {
  echo "ERROR: missing forward seed bank $BANK" >&2; exit 1; }
log "orient to the forward seed and extend shared flanks into the consensus"
for f in "$STEP8_OUT"/*_top100.aln.fa "$STEP8_OUT"/*_rand100.aln.fa "$STEP8_OUT"/*_subfam.aln.fa; do
  [[ -f "$f" ]] || continue
  python3 "$DISC/correct_published_aln.py" "$BANK" "$f" --write
done
# Last: the consensus exactly as the genome was searched goes back on as row 2
# (<subfamily>_seed_as_searched). Row 1 is rebuilt above and may differ; the seed
# row lets both be read together. The DISC tools skip rows with that name, so a
# later re-score or re-correction measures the copies only.
if [[ "${SKIP_SEED_ROW:-0}" != "1" && -f "$SINEDERELLA_BIN/tools/add_seed_row.py" ]]; then
  log "add the consensus as searched as row 2"
  python3 "$SINEDERELLA_BIN/tools/add_seed_row.py" "$RUN_ROOT/consensuses.clean.fa" "$SPECIES"     "$STEP8_OUT"/*_top100.aln.fa "$STEP8_OUT"/*_rand100.aln.fa "$STEP8_OUT"/*_subfam.aln.fa     --threads "${THREADS:-4}" || { log "ERROR: add_seed_row failed: the plates would go out without the consensus as searched"; exit 1; }   # a WARNING until 2026-10-19 (GLM audit 2026-10-19)
fi
if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  cp -a "$STEP8_OUT"/*_top100.aln.fa "$STEP8_OUT"/*_rand100.aln.fa \
    "$STEP8_OUT"/*_subfam.aln.fa "$OUT_DIR/" 2>/dev/null || true
fi

# Consensus alignments, a regular output (owner 2026-10-18): every consensus of the searched bank, and separately the
# consensuses of subfamilies whose flanks were shown unique (flank_verdict.tsv: interspersed). MAFFT L-INS-i with
# --threadit 0 (same input -> same alignment). Written after the DISC tools above, which must not touch them.
CONS_BANK="$RUN_ROOT/consensuses.clean.fa"
if [[ -s "$CONS_BANK" ]]; then
  log "consensus alignments"
  mafft --quiet --localpair --maxiterate 1000 --threadit 0 --thread "${THREADS:-4}" "$CONS_BANK" \
    > "$STEP8_OUT/${SPECIES}_consensuses.aln.fa" || log "WARNING: consensus alignment failed"
  FV="$RUN_ROOT/step2/step2_output/flank_verdict.tsv"
  if [[ -s "$FV" ]]; then
    awk -F'\t' 'FILENAME == ARGV[1] { if (FNR > 1 && $2 == "interspersed") ok[$1] = 1; next }
                /^>/ { n = substr($1, 2); p = (n in ok) } p' "$FV" "$CONS_BANK" > "$STEP8_OUT/.sine_cons.fa"
    if [[ $(grep -c '^>' "$STEP8_OUT/.sine_cons.fa") -ge 2 ]]; then
      mafft --quiet --localpair --maxiterate 1000 --threadit 0 --thread "${THREADS:-4}" "$STEP8_OUT/.sine_cons.fa" \
        > "$STEP8_OUT/${SPECIES}_consensuses_sine.aln.fa" || log "WARNING: SINE consensus alignment failed"
    fi
    rm -f "$STEP8_OUT/.sine_cons.fa"
  fi
  if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
    cp -a "$STEP8_OUT/${SPECIES}"_consensuses*.aln.fa "$OUT_DIR/" 2>/dev/null || true
  fi
fi

log "Done: $STEP8_OUT"
