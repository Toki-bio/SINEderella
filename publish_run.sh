#!/usr/bin/env bash
# End-to-end publish: step4 (if needed) → publish alignments → step6 HTML → discriminator overlay.
#
# Usage:
#   publish_run.sh <RUN_ROOT> <SPECIES_CODE>
#
# Environment:
#   SINEDERELLA_BIN  — SINEderella repo root (default: dir containing this script)
#   DISC             — SINE-discriminator root (auto-detected if unset)
#   RAW_ALN_BASE     — published raw URL prefix for MSA viewer links, e.g.
#                      https://raw.githubusercontent.com/org/repo/main/mysp/alignments/
#   PAGES_INDEX      — optional URL for multi-species “All species” nav link
#   USE_DISC         — 1 (default) inject verdict columns when DISC found
#   THREADS          — passed to step4
#   SKIP_ALIGN       — 1 to skip align_for_publish (alignments already built)
#   SKIP_STEP4       — 1 to skip step4 even if plots missing
#
# Output: <RUN_ROOT>/results/report.html
set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE>}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE>}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINEDERELLA_BIN="${SINEDERELLA_BIN:-$SCRIPT_DIR}"
PUBLISH_DIR="$SINEDERELLA_BIN/publish"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ log "ERROR: $*"; exit 1; }

[[ -d "$RUN_ROOT" ]] || die "RUN_ROOT not found: $RUN_ROOT"

# Resolve DISC (SINE-discriminator)
if [[ -z "${DISC:-}" ]]; then
  for d in "${SINE_DISC:-}" "/staging/tmp/sinedisc" \
           "$HOME/SINE_discriminator/site" "$HOME/SINE_discriminator" \
           "C:/work/SINE_discriminator/site" "C:/work/SINE_discriminator"; do
    [[ -n "$d" && -f "$d/verdict.py" ]] && DISC="$d" && break
  done
fi
if [[ -n "${DISC:-}" && -f "$DISC/verdict.py" ]]; then
  log "DISC=$DISC"
else
  log "WARNING: DISC not found — publish alignments need boundary_justify.py; verdict overlay skipped"
  DISC=""
fi

# Locate step2 output
S2="$(find "$RUN_ROOT/step2" -maxdepth 1 -type d -name 'step2_output*' 2>/dev/null \
      | sort | tail -1)"
[[ -n "$S2" && -d "$S2" ]] || die "No step2_output under $RUN_ROOT"

# ── 1. Step4 plots (pctid TSVs + PNGs) ─────────────────────────────────────
if [[ "${SKIP_STEP4:-0}" != "1" ]]; then
  need4=0
  [[ -d "$S2/plots" ]] || need4=1
  [[ "$need4" -eq 0 ]] && \
    ls "$S2/plots"/*_pctid.tsv >/dev/null 2>&1 || need4=1
  if [[ "$need4" -eq 1 ]]; then
    [[ -x "$RUN_ROOT/step4_plots.sh" ]] || cp -f "$SINEDERELLA_BIN/step4_plots.sh" "$RUN_ROOT/"
    [[ -f "$RUN_ROOT/plot_subfamily.py" ]] || \
      cp -f "$SINEDERELLA_BIN/plot_subfamily.py" "$RUN_ROOT/"
    log "step4: generating plots + pctid TSVs"
    (cd "$RUN_ROOT" && ./step4_plots.sh "$(pwd -P)" "$THREADS")
  else
    log "step4: plots present, skipping"
  fi
fi

# ── 2. Publish alignments (step7/8 + DISC MSA tools) ───────────────────────
if [[ "${SKIP_ALIGN:-0}" != "1" ]]; then
  log "publish: align_for_publish"
  bash "$PUBLISH_DIR/align_for_publish.sh" "$RUN_ROOT" "$SPECIES"
else
  log "publish alignments skipped (SKIP_ALIGN=1)"
fi

# ── 3. Step6 HTML report ───────────────────────────────────────────────────
[[ -f "$RUN_ROOT/step6_report.py" ]] || cp -f "$SINEDERELLA_BIN/step6_report.py" "$RUN_ROOT/"
mkdir -p "$RUN_ROOT/results"
OUT_HTML="$RUN_ROOT/results/report.html"

STEP6_ARGS=( "$RUN_ROOT" "--out" "$OUT_HTML" "--no-sineplot" )
[[ -n "${RAW_ALN_BASE:-}" ]] && \
  STEP6_ARGS+=( "--aln-base" "$RAW_ALN_BASE" )
[[ -n "$SPECIES" ]] && \
  STEP6_ARGS+=( "--species-code" "$SPECIES" )
[[ -n "${PAGES_INDEX:-}" ]] && \
  STEP6_ARGS+=( "--pages-index" "$PAGES_INDEX" )

log "step6: building HTML → $OUT_HTML"
python3 "$RUN_ROOT/step6_report.py" "${STEP6_ARGS[@]}"

# ── 4. SINE-discriminator overlay (verdict columns) ────────────────────────
if [[ "${USE_DISC:-1}" == "1" && -n "$DISC" ]]; then
  INJECT="$DISC/inject_disc_report.py"
  [[ -f "$INJECT" ]] || INJECT="$DISC/site/inject_disc_report.py"
  if [[ -f "$INJECT" ]]; then
    log "discriminator: inject alignment verdicts"
    INJ_ARGS=( "$OUT_HTML" "$RUN_ROOT/results/alignments" "$SPECIES" )
    [[ -n "${RAW_ALN_BASE:-}" ]] && INJ_ARGS+=( "--raw-base" "$RAW_ALN_BASE" )
    INJ_ARGS+=( "--summary-tsv" "$S2/summary.by_subfam.tsv" )
    python3 "$INJECT" "${INJ_ARGS[@]}"
  else
    log "WARNING: inject_disc_report.py not found under DISC"
  fi
fi

log "Published report: $OUT_HTML"
