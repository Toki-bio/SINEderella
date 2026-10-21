#!/usr/bin/env bash
# verdict_chain.sh RUN_ROOT [SPECIES_CODE] - the verdict of every copy (docs/VERDICTS.md): step 7 subfamily flank test ->
# border loop for subfamilies step 7 did not confirm -> tools/flank_verdict.py (loci_verdict.tsv, flank_verdict.tsv).
#
# Run by the orchestrator in every run after flankscan stage 9 (docs/AIM.md: the verdicted loci list is an output of the
# workflow itself, not of the publish step), and by publish/align_for_publish.sh when the run has no fresh verdicts. Until
# 2026-10-21 it lived inside align_for_publish.sh, so a run without --publish had no loci_verdict.tsv.
# Writes step2_output/verdict_chain.done when it finishes.
#
# Environment: SINEDERELLA_BIN, DISC (default: disc/ of this repository), PEEL_FLAGS, PEEL_ALN_DIR, SKIP_STEP7,
#   SKIP_BORDER_LOOP, BORDER_SCAN / SKIP_BORDER_SCAN (as before in align_for_publish.sh)
set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> [SPECIES_CODE]}"
SPECIES="${2:-}"
PUBLISH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINEDERELLA_BIN="${SINEDERELLA_BIN:-$(dirname "$PUBLISH_DIR")}"
PEEL_FLAGS="${PEEL_FLAGS:-}"
PEEL_ALN_DIR="${PEEL_ALN_DIR:-}"
# needs_border_loop.py imports boundary.py from DISC
DISC="${DISC:-${SINE_DISC:-$SINEDERELLA_BIN/disc}}"
[[ -f "$DISC/boundary.py" ]] || { echo "ERROR: verdict_chain: no boundary.py in DISC=$DISC" >&2; exit 1; }
export DISC
S2="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
[[ -n "$S2" && -d "$S2" ]] || { echo "ERROR: verdict_chain: no step2 output in $RUN_ROOT" >&2; exit 1; }
rm -f "$S2/verdict_chain.done"

# DRAGEN toolchain only where it exists. On other hosts (therioserver) the unconditional
# prepend put /usr/bin ahead of the active conda env: /usr/bin/python3 has no numpy, so the
# border scan failed with a WARN and every subfamily was reported "none flagged".
[[ -d /staging/conda/envs/bioinfo/bin ]] && \
  export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

if [[ "${SKIP_STEP7:-0}" != "1" ]]; then
  log "step7: boundary refinement"
  "$SINEDERELLA_BIN/step7_boundary_refine.sh" "$RUN_ROOT" 50 50 1000
else
  log "step7: skipped (SKIP_STEP7=1)"
fi

if [[ "${SKIP_BORDER_LOOP:-0}" == "1" ]]; then
  # SKIP_BORDER_LOOP means "don't RE-RUN the expensive scan", not "discard
  # whatever border-loop widening already exists". This used to unconditionally
  # cp the raw, un-widened assigned.fasta/consensuses.clean.fa over the
  # *.publish.* files, silently reverting every subfamily's already-computed
  # border-loop extension back to zero on every skip-flagged rerun -- caught
  # 2026-09-10 on oma_SINE16: a real, correctly-computed +4/+75bp extension
  # from an earlier run vanished after two SKIP_BORDER_LOOP=1 reruns done for
  # an unrelated SINE10 fix. Only fall back to the raw copy when no prior
  # publish file exists at all.
  log "border loop: skipped (SKIP_BORDER_LOOP=1) -- preserving existing *.publish.* if present"
  if [[ ! -f "$RUN_ROOT/step2/step2_output/assigned.publish.fasta" ]]; then
    cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
      "$RUN_ROOT/step2/step2_output/assigned.publish.fasta"
  fi
  if [[ ! -f "$RUN_ROOT/consensuses.publish.fa" ]]; then
    cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.publish.fa"
  fi
else
  log "detecting subfamilies needing border loop"
  NEED_ARGS=(python3 "$PUBLISH_DIR/needs_border_loop.py" "$RUN_ROOT")
  [[ -n "$PEEL_FLAGS" && -f "$PEEL_FLAGS" ]] && NEED_ARGS+=(--peel-flags "$PEEL_FLAGS")
  # the column-walk border scan is opt-in (BORDER_SCAN=1): on unaligned extracts it reads flank noise below ~20 copies (random
  # DNA, 8 copies: 50 of 50 flagged) and the loop then replaced the searched consensus (audit 2026-10-14)
  [[ "${BORDER_SCAN:-0}" == "1" && "${SKIP_BORDER_SCAN:-0}" != "1" ]] && NEED_ARGS+=(--scan)
  # a crash of the detection must stop the publish, not publish as "none flagged" (GLM audit 2026-10-19)
  "${NEED_ARGS[@]}" > "$RUN_ROOT/.needs_border_loop.out" || {
    echo "ERROR: needs_border_loop.py failed - the border loop cannot be decided" >&2; exit 1; }
  mapfile -t BORDER_SFS < "$RUN_ROOT/.needs_border_loop.out"
  rm -f "$RUN_ROOT/.needs_border_loop.out"

  # results of an earlier publish are not re-applied: apply_border_to_assigned / merge_rebuilt_consensus read every file there
  rm -rf "$RUN_ROOT/border_loop" "$RUN_ROOT/rebuilt_consensus"
  if ((${#BORDER_SFS[@]})); then
    log "border loop: ${#BORDER_SFS[@]} subfamilies"
    for sf in "${BORDER_SFS[@]}"; do
      PEEL_ALN=""
      if [[ -n "$PEEL_ALN_DIR" ]]; then
        [[ -f "$PEEL_ALN_DIR/${sf}.aln" ]] && PEEL_ALN="$PEEL_ALN_DIR/${sf}.aln"
        [[ -z "$PEEL_ALN" && -f "$PEEL_ALN_DIR/${sf#${SPECIES}_}.aln" ]] && \
          PEEL_ALN="$PEEL_ALN_DIR/${sf#${SPECIES}_}.aln"
      fi
      python3 "$PUBLISH_DIR/border_loop_subfam.py" "$RUN_ROOT" "$sf" \
        ${PEEL_ALN:+--peel-aln "$PEEL_ALN"} \
        --sine-script "$SINEDERELLA_BIN/sine_consensus.sh" || {
          log "WARN border loop failed for $sf — continuing"; }
    done
    python3 "$PUBLISH_DIR/apply_border_to_assigned.py" "$RUN_ROOT"
    python3 "$PUBLISH_DIR/merge_rebuilt_consensus.py" "$RUN_ROOT"
  else
    log "border loop: none flagged"
    cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
      "$RUN_ROOT/step2/step2_output/assigned.publish.fasta"
    cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.publish.fa"
  fi
fi

# flank verdicts: the expansion's end state decides (docs/AIM.md). A subfamily whose flanks are still shared at the limit,
# or could not be tested, is not interspersed: its copies are not SINE loci in loci_verdict.tsv and step6 leaves them out
# of the SINE totals. Until 2026-10-17 nothing read the end state of step7 or of the border loop.
log "flank verdicts"
python3 "$SINEDERELLA_BIN/tools/flank_verdict.py" "$RUN_ROOT" >&2 || {
  echo "ERROR: tools/flank_verdict.py failed" >&2; exit 1; }

printf 'finished\t%s\ncode\t%s\n' "$(date '+%F %T')" "$(git -C "$SINEDERELLA_BIN" log -1 --format=%h 2>/dev/null || echo unknown)" \
  > "$S2/verdict_chain.done"
