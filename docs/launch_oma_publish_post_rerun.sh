#!/usr/bin/env bash
# Publish oma alignments + report from consensus-repair rerun (merged bank, step2 rerun dir).
#
# Wires step2_output -> step2_output_rerun_* so hardcoded publish paths see new assignments.
# Does NOT touch step2_output_pre_repair backup.
#
# Usage: bash launch_oma_publish_post_rerun.sh [RUN_ROOT]
set -euo pipefail

RUN="${1:-/staging/tmp/scorpions/oma/run_oma}"
SINEDERELLA="${SINEDERELLA:-/staging/tmp/SINEderella}"
DISC="${DISC:-/staging/tmp/sinedisc}"
RERUN_DIR="${RERUN_DIR:-}"
LOG="${RUN}/publish_post_rerun.log"
THREADS="${THREADS:-32}"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }
die(){ log "ERROR: $*"; exit 1; }

[[ -d "$RUN" ]] || die "RUN not found: $RUN"
[[ -d "$SINEDERELLA" ]] || die "SINEDERELLA not found: $SINEDERELLA"
[[ -f "$DISC/rebuild_consensus_row.py" ]] || die "DISC not found: $DISC"

exec >>"$LOG" 2>&1
log "=== oma publish post-rerun ==="
log "RUN=$RUN DISC=$DISC THREADS=$THREADS"

export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
export PEEL_FLAGS="${PEEL_FLAGS:-$SINEDERELLA/oma_peel_border_flags.tsv}"
export PEEL_ALN_DIR="${PEEL_ALN_DIR:-/staging/tmp/scorpions/oma/sd/rebuild}"
export TRIM_DISPLAY_MODE="${TRIM_DISPLAY_MODE:-occupancy}"

cd "$SINEDERELLA"
git fetch origin main 2>/dev/null || true
git reset --hard origin/main 2>/dev/null || true
log "SINEderella at $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
chmod +x step7_boundary_refine.sh step8a_extract_alignments.sh run_publish_alignments.sh \
  step2_asSINEment.sh step3_postprocess.sh step4_plots.sh 2>/dev/null || true

cd "$RUN"

# ── pick rerun step2 dir ───────────────────────────────────────────────────
if [[ -z "$RERUN_DIR" ]]; then
  RERUN_DIR="$(ls -dt step2/step2_output_rerun_* 2>/dev/null | head -1 || true)"
fi
[[ -n "$RERUN_DIR" && -d "$RERUN_DIR" ]] || die "no step2_output_rerun_* under $RUN/step2"
log "rerun step2: $RERUN_DIR"

# ── point step2_output at rerun (backup real dir once) ───────────────────────
cd step2
if [[ -L step2_output ]]; then
  rm -f step2_output
elif [[ -d step2_output && ! -e step2_output_pre_repair ]]; then
  mv step2_output step2_output_pre_repair
  log "backed up step2_output -> step2_output_pre_repair"
fi
ln -sfn "$(basename "$RERUN_DIR")" step2_output
log "step2_output -> $(readlink -f step2_output)"
cd "$RUN"

[[ -s consensuses.clean.fa ]] || die "missing consensuses.clean.fa"
[[ -s consensuses.rebuilt.fa ]] || {
  log "rebuild consensus bank..."
  cp -f "$SINEDERELLA/rebuild_consensus_bank.py" .
  python3 rebuild_consensus_bank.py "$(pwd -P)"
}

# ── publish alignments (step7 + border loop + step8a + DISC) ────────────────
log "publish alignments starting..."
export DISC
bash "$SINEDERELLA/run_publish_alignments.sh" "$RUN" oma "$RUN/alignments"
log "publish alignments done"

# ── drop obsolete merged-away subfamilies from alignments/ ──────────────────
for dead in oma_big76 oma_group34 oma_sub515; do
  rm -f alignments/${dead}_*.aln.fa alignments/oma_${dead}_*.aln.fa 2>/dev/null || true
done
log "removed obsolete big76/group34/sub515 alignment symlinks/files if any"

# ── copy rerun plots into site-facing tree ─────────────────────────────────
PLOTS_SRC="$RERUN_DIR/plots"
PLOTS_DST="$RUN/alignments/oma/plots"
if [[ -d "$PLOTS_SRC" ]]; then
  mkdir -p "$PLOTS_DST"
  cp -f "$PLOTS_SRC"/*_pctid.tsv "$PLOTS_DST"/ 2>/dev/null || true
  cp -f "$PLOTS_SRC"/*_divergence.png "$PLOTS_DST"/ 2>/dev/null || true
  cp -f "$PLOTS_SRC"/*_nucfreq.png "$PLOTS_DST"/ 2>/dev/null || true
  log "plots copied to $PLOTS_DST"
fi

# ── step6 HTML (embed images is default; no --embed-images flag) ───────────
cp -f "$SINEDERELLA/step6_report.py" .
mkdir -p results
log "step6 report..."
python3 step6_report.py "$(pwd -P)" \
  --out results/report_post_rerun.html \
  --no-sineplot \
  --species-code oma \
  --aln-base "https://raw.githubusercontent.com/Toki-bio/SINE-discriminator/main/alignments/"

if [[ -f "$DISC/inject_disc_report.py" ]]; then
  log "inject discriminator verdicts..."
  python3 "$DISC/inject_disc_report.py" \
    results/report_post_rerun.html \
    "$RUN/alignments" \
    --out results/report_post_rerun.html || log "WARN inject failed"
fi

log "=== PUBLISH POST-RERUN COMPLETE ==="
log "alignments: $RUN/alignments/"
log "report:     $RUN/results/report_post_rerun.html"
log "next: sync alignments/ + plots/ to Toki-bio/SINE-discriminator and run patch_oma_report.py"
