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
  for d in "${SINE_DISC:-}" "/staging/tmp/sinedisc" "$HOME/SINE_discriminator/site" \
           "$HOME/SINE_discriminator"; do
    [[ -n "$d" && -f "$d/boundary_justify.py" ]] && DISC="$d" && break
  done
fi
[[ -n "$DISC" && -f "$DISC/boundary_justify.py" ]] || {
  echo "ERROR: set DISC to SINE-discriminator root (need boundary_justify.py)" >&2
  exit 1
}

export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

if [[ "${SKIP_STEP7:-0}" != "1" ]]; then
  log "step7: boundary refinement"
  "$SINEDERELLA_BIN/step7_boundary_refine.sh" "$RUN_ROOT" 50 50 1000
else
  log "step7: skipped (SKIP_STEP7=1)"
fi

if [[ "${SKIP_BORDER_LOOP:-0}" == "1" ]]; then
  log "border loop: skipped (SKIP_BORDER_LOOP=1)"
  cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
    "$RUN_ROOT/step2/step2_output/assigned.publish.fasta"
  cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.publish.fa"
else
  log "detecting subfamilies needing border loop"
  NEED_ARGS=(python3 "$PUBLISH_DIR/needs_border_loop.py" "$RUN_ROOT")
  [[ -n "$PEEL_FLAGS" && -f "$PEEL_FLAGS" ]] && NEED_ARGS+=(--peel-flags "$PEEL_FLAGS")
  [[ "${SKIP_BORDER_SCAN:-0}" != "1" ]] && NEED_ARGS+=(--scan)
  mapfile -t BORDER_SFS < <("${NEED_ARGS[@]}" || true)

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

for need in step7_boundary_refine.sh step8a_extract_alignments.sh; do
  [[ -x "$SINEDERELLA_BIN/$need" ]] || {
    echo "ERROR: missing $SINEDERELLA_BIN/$need" >&2; exit 1; }
done

log "step8a: extract publish alignments"
STEP8_OUT="$RUN_ROOT/results/alignments"
mkdir -p "$STEP8_OUT"
cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.clean.fa.pre_publish.bak"
cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta.pre_publish.bak"
cp -f "$RUN_ROOT/consensuses.publish.fa" "$RUN_ROOT/consensuses.clean.fa"
cp -f "$RUN_ROOT/step2/step2_output/assigned.publish.fasta" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta"
"$SINEDERELLA_BIN/step8a_extract_alignments.sh" "$RUN_ROOT" "$SPECIES"
cp -f "$RUN_ROOT/consensuses.clean.fa.pre_publish.bak" "$RUN_ROOT/consensuses.clean.fa"
cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta.pre_publish.bak" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta"

if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  mkdir -p "$OUT_DIR"
  cp -a "$STEP8_OUT"/*.aln.fa "$OUT_DIR/" 2>/dev/null || true
fi

log "DISC: rebuild consensus row + justify + trim"
shopt -s nullglob
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/rebuild_consensus_row.py" "$f"
  python3 "$DISC/boundary_justify.py" "$f"
  python3 "$DISC/trim_display_flanks.py" "$f" \
    --mode "${TRIM_DISPLAY_MODE:-occupancy}" || true
done

log "Done: $STEP8_OUT"
