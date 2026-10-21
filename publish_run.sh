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
  # disc/ = the SINE-discriminator modules SINEderella uses, kept in this repository (2026-10-21; before that only
  # unversioned copies on each machine existed). An explicit DISC or SINE_DISC still wins.
  for d in "${SINE_DISC:-}" "$SINEDERELLA_BIN/disc" "/staging/tmp/sinedisc" \
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
# Always the installed step6: publishing is presentation, and a copy frozen in the run dir kept old
# page bugs alive (the "Genome: ?" header fix and the seed label never reached the 25 bat runs,
# 2026-09-28). The previous copy is kept as step6_report.py.prev.
if [[ -f "$RUN_ROOT/step6_report.py" ]] && ! cmp -s "$SINEDERELLA_BIN/step6_report.py" "$RUN_ROOT/step6_report.py"; then
    cp -f "$RUN_ROOT/step6_report.py" "$RUN_ROOT/step6_report.py.prev"
fi
cp -f "$SINEDERELLA_BIN/step6_report.py" "$RUN_ROOT/"
# its helper modules too (report_hierarchy, report_profile, report_flank_uniqueness, ...): step6 runs from
# the run dir, and a helper it cannot import silently drops its panel (rsi_v3: no hierarchy section)
cp -f "$SINEDERELLA_BIN"/report_*.py "$RUN_ROOT/" 2>/dev/null || true
# report_blocks imports tools/consensus_blocks.py from the run dir; without it the whole similarity section (with the
# length-version, array, satellite and audit tables) was dropped with a one-line warning (toy run 2026-10-05)
mkdir -p "$RUN_ROOT/tools" && cp -f "$SINEDERELLA_BIN/tools/consensus_blocks.py" "$RUN_ROOT/tools/" 2>/dev/null || true
mkdir -p "$RUN_ROOT/results"
OUT_HTML="$RUN_ROOT/results/report.html"
# composites of the bank as units, for the Element hierarchy (grey = stretches no unit explains; 2026-10-18)
python3 "$SINEDERELLA_BIN/tools/composite_units.py" "$RUN_ROOT" "${THREADS:-4}" \
  || log "WARNING: tools/composite_units.py failed; the hierarchy is drawn without unit composition"
# tRNA / 5S / 7SL origin of every consensus and composite part (needs TRNA_REF; RFAM_CM and RNATOOLS_BIN optional)
bash "$SINEDERELLA_BIN/tools/trna_heads.sh" "$RUN_ROOT" "${THREADS:-4}" \
  || log "WARNING: tools/trna_heads.sh failed; the report has no tRNA section"

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
    # flank verdicts (tools/flank_verdict.py): copies column = SINE loci, and a row for every subfamily, also those
    # without SINE loci and therefore without plates (the row says why). The overlay ignores the option if it is older.
    if [[ -s "$S2/flank_verdict.tsv" ]] && grep -q -- "--verdict-tsv" "$INJECT"; then
      INJ_ARGS+=( "--verdict-tsv" "$S2/flank_verdict.tsv" )
    fi
    # consensus length column, and rows grouped by family when the run has a family table (--families; owner 2026-10-21)
    if grep -q -- "--bank" "$INJECT"; then
      INJ_ARGS+=( "--bank" "$RUN_ROOT/consensuses.clean.fa" )
      [[ -s "$RUN_ROOT/families/families.clean.tsv" ]] && INJ_ARGS+=( "--families" "$RUN_ROOT/families/families.clean.tsv" )
    fi
    python3 "$INJECT" "${INJ_ARGS[@]}"
  else
    log "WARNING: inject_disc_report.py not found under DISC"
  fi
fi

# ── 5. The consensus exactly as searched, as row 2 of every plate ─────────────
# Row 1 is rebuilt from the copies and its edge moved; the seed the genome was searched with would
# otherwise be gone. Added last, so the report's verdict columns are computed without it.
# align_for_publish.sh already adds it; a second pass upper-cased the trims the first one had restored and dropped the trim
# proposals from the plate (audit 2026-10-14). Only for plates built without align_for_publish: SKIP_ALIGN=1 LEGACY_SEED_ROW=1.
if [[ "${LEGACY_SEED_ROW:-0}" == "1" && "${SKIP_SEED_ROW:-0}" != "1" && -f "$SINEDERELLA_BIN/tools/add_seed_row.py" ]]; then
  log "plates: add the consensus as searched (row 2)"
  python3 "$SINEDERELLA_BIN/tools/add_seed_row.py" "$RUN_ROOT/consensuses.clean.fa" "$SPECIES"     "$RUN_ROOT"/results/alignments/*_{top100,rand100,subfam}.aln.fa --threads "${THREADS:-4}" ||     die "add_seed_row.py failed: the plates would go out without the consensus as searched"   # a WARNING until 2026-10-19 (GLM audit 2026-10-19)
fi

log "Published report: $OUT_HTML"
