#!/usr/bin/env bash
# step6_report.sh -- thin wrapper around step6_report.py
# Usage:  step6_report.sh <RUN_ROOT> [extra args forwarded to python script]
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

RUN_ROOT="${1:-}"
[[ -n "$RUN_ROOT" ]] || die "usage: $0 RUN_ROOT [extra args]"
shift || true

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PY="$SCRIPT_DIR/step6_report.py"
[[ -f "$PY" ]] || die "missing helper: $PY"

# Pick the newest python on PATH that is >= 3.8 (SINEplot needs it).
pick_python() {
  local cand
  if [[ -n "${SINEDERELLA_PYTHON:-}" ]] && command -v "$SINEDERELLA_PYTHON" >/dev/null 2>&1; then
    echo "$SINEDERELLA_PYTHON"; return
  fi
  for cand in python3.12 python3.11 python3.10 python3.9 python3.8 python3; do
    if command -v "$cand" >/dev/null 2>&1; then
      if "$cand" -c 'import sys; sys.exit(0 if sys.version_info >= (3,8) else 1)' 2>/dev/null; then
        echo "$cand"; return
      fi
    fi
  done
  echo ""
}
PYBIN="$(pick_python)"
[[ -n "$PYBIN" ]] || die "no python >= 3.8 found in PATH (set SINEDERELLA_PYTHON to override)"

log "Building HTML report for $RUN_ROOT  (python: $PYBIN)"
"$PYBIN" "$PY" "$RUN_ROOT" "$@"
log "Done. Open: $RUN_ROOT/results/report.html"
