#!/usr/bin/env bash
# ==========================================================
# run_step5_wrapper.sh — Universal wrapper for step5
#
# Usage:
#   run_step5_wrapper.sh <run_dir> [step5_script_path]
#
# Arguments:
#   <run_dir>              Path to SINEderella run directory
#   [step5_script_path]    Path to step5_align_subfamilies.sh (auto-detected if omitted)
#
# Works from anywhere: uses absolute paths internally
# ==========================================================
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

# ── Parse arguments ──────────────────────────────────────────────────────────
RUN_DIR="${1:-}"
STEP5_SCRIPT="${2:-}"

[[ -n "$RUN_DIR" ]] || die "Usage: $0 <run_dir> [step5_script_path]"
[[ -d "$RUN_DIR" ]] || die "Run directory not found: $RUN_DIR"

RUN_DIR="$(readlink -f "$RUN_DIR")"
CONS_BANK="$RUN_DIR/consensuses.clean.fa"
[[ -f "$CONS_BANK" ]] || die "Consensus bank not found: $CONS_BANK"

# Find or validate step5 script
if [[ -z "$STEP5_SCRIPT" ]]; then
  # Try to find step5_align_subfamilies.sh in current directory or PATH
  if [[ -f "step5_align_subfamilies.sh" ]]; then
    STEP5_SCRIPT="$(readlink -f "step5_align_subfamilies.sh")"
  elif command -v step5_align_subfamilies.sh >/dev/null 2>&1; then
    STEP5_SCRIPT="$(command -v step5_align_subfamilies.sh)"
  else
    die "step5_align_subfamilies.sh not found in current directory or PATH. Provide explicit path as argument 2."
  fi
else
  [[ -f "$STEP5_SCRIPT" ]] || die "step5_align_subfamilies.sh not found: $STEP5_SCRIPT"
  STEP5_SCRIPT="$(readlink -f "$STEP5_SCRIPT")"
fi

# ── Extract species name ─────────────────────────────────────────────────────
SPECIES="$(basename "$(dirname "$RUN_DIR")")"
[[ -n "$SPECIES" ]] || SPECIES="genome"

log "Run directory: $RUN_DIR"
log "Species:       $SPECIES"
log "Consensus:     $CONS_BANK"
log "Step5 script:  $STEP5_SCRIPT"

# ── Build a minimal single-species layout for step5 in a temp dir ─────────────
# step5 expects: <multi_dir>/<Species>/run_<timestamp>/step2/...
# We create that structure with a symlink, run step5, then move results
# into the run directory and clean up.
TMPDIR_MULTI="$(mktemp -d "${TMPDIR:-/tmp}/step5_wrapper.XXXXXX")"
trap 'rm -rf "$TMPDIR_MULTI"' EXIT

mkdir -p "$TMPDIR_MULTI/$SPECIES"
ln -s "$RUN_DIR" "$TMPDIR_MULTI/$SPECIES/$(basename "$RUN_DIR")"

# ── Run step5 ────────────────────────────────────────────────────────────────
log "Running step5..."
"$STEP5_SCRIPT" "$CONS_BANK" "$TMPDIR_MULTI"

# ── Move results into the run directory ──────────────────────────────────────
SRC_DIR="$TMPDIR_MULTI/alignments/$SPECIES"
DEST_DIR="$RUN_DIR/alignments/$SPECIES"

if [[ -d "$SRC_DIR" ]] && [[ -n "$(find "$SRC_DIR" -type f -name '*.aln.fa' 2>/dev/null)" ]]; then
  mkdir -p "$DEST_DIR"
  mv "$SRC_DIR"/*.aln.fa "$DEST_DIR"/
  log "✓ Success! Alignments in: $DEST_DIR"
else
  log "⚠ No alignments produced. Check step5 output."
fi
