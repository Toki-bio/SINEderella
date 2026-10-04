#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# step4_diagnostic.sh â€” SINEderella diagnostic analysis wrapper
#
# Runs step4_diagnostic.py which performs:
#   1. PCA on subfamily bitscore matrix
#   2. Diagnostic position weight computation (MI, RF, KL)
#   3. Per-copy diagnostic state extraction
#   4. Discordance detection (fp_risk, fn_risk, chimera)
#   5. Diagnostic heatmap, volcano plots, diagnostic vs bitscore scatter
#   6. SINEplot input generation (all-vs-all ssearch36 matrix)
#
# Usage:
#   bash step4_diagnostic.sh RUN_ROOT [THREADS] [TOP_K]
#
# Dependencies (Python):
#   pip install numpy scipy pandas matplotlib scikit-learn
# Dependencies (system):
#   mafft, ssearch36
# =============================================================================

log(){ printf '[%s] [diag] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ printf '[%s] [diag] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

RUN_ROOT="${1:-}"
[[ -n "$RUN_ROOT" ]] || die "usage: $0 RUN_ROOT [THREADS] [TOP_K]"
RUN_ROOT="$(readlink -f "$RUN_ROOT")"
[[ -d "$RUN_ROOT" ]] || die "Run directory not found: $RUN_ROOT"

THREADS="${2:-$(nproc 2>/dev/null || echo 4)}"
TOP_K="${3:-20}"

# Locate this script's directory (where step4_diagnostic.py lives)
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
DIAG_PY="$SCRIPT_DIR/step4_diagnostic.py"

[[ -f "$DIAG_PY" ]] || die "step4_diagnostic.py not found in: $SCRIPT_DIR"

# Check Python dependencies
log "Checking Python dependencies ..."
python3 -c "import numpy, scipy, pandas, matplotlib, sklearn" 2>/dev/null || {
    log "Missing Python packages. Install with:"
    log "  pip install numpy scipy pandas matplotlib scikit-learn"
    die "Python dependencies not satisfied"
}

# Check system dependencies
for tool in mafft ssearch36; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool not found in PATH"
done

log "Run root: $RUN_ROOT"
log "Threads: $THREADS"
log "Top-K diagnostic positions: $TOP_K"

# Run the Python analysis
python3 "$DIAG_PY" "$RUN_ROOT" --threads "$THREADS" --top-k "$TOP_K"

# â”€â”€ Generate SINEplot HTML if input was created â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
SINEPLOT_INPUT="$RUN_ROOT/step2/step2_output/diagnostic/sineplot_input.txt"
SINEPLOT_HTML="$RUN_ROOT/step2/step2_output/diagnostic/sineplot.html"

if command -v python3 >/dev/null 2>&1 && [[ -f "$SINEPLOT_INPUT" ]]; then
    SINEPLOT_PY="$SCRIPT_DIR/SINEplot.py"
    # Also check common locations
    [[ -f "$SINEPLOT_PY" ]] || SINEPLOT_PY="$HOME/SINEplot/SINEplot.py"
    [[ -f "$SINEPLOT_PY" ]] || SINEPLOT_PY="$HOME/git/SINEplot/SINEplot.py"

    if [[ -f "$SINEPLOT_PY" ]]; then
        log "Generating SINEplot interactive HTML ..."
        python3 "$SINEPLOT_PY" "$SINEPLOT_INPUT" -o "$SINEPLOT_HTML" \
            -t "SINE Subfamily Diagnostic View" 2>&1 | while IFS= read -r line; do
            log "[SINEplot] $line"
        done
        if [[ -f "$SINEPLOT_HTML" ]]; then
            log "SINEplot HTML: $SINEPLOT_HTML"
        fi
    else
        log "SINEplot.py not found â€” skip HTML generation"
        log "To generate: python SINEplot.py $SINEPLOT_INPUT -o $SINEPLOT_HTML"
    fi
fi

# â”€â”€ Symlink outputs to results/ â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
RESULTS_DIR="$RUN_ROOT/results"
if [[ -d "$RESULTS_DIR" ]]; then
    DIAG_DIR="$RUN_ROOT/step2/step2_output/diagnostic"
    if [[ -d "$DIAG_DIR" ]]; then
        ln -sfn "../step2/step2_output/diagnostic" "$RESULTS_DIR/diagnostic" 2>/dev/null || true
        log "Symlinked diagnostic/ into results/"
    fi
fi

log "Done."
