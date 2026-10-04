#!/usr/bin/env bash
# Compare conse vs sine_consensus.sh on one top100 alignment.
# Usage: flank_border_consensus_test.sh <top100.aln.fa> <out_dir> [peel_loci.aln]
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:/usr/bin:$PATH

ALN="${1:?alignment required}"
OUT="${2:?output dir required}"
PEEL_ALN="${3:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$SCRIPT_DIR/flank_border_consensus_test.py" \
  "$ALN" "$OUT" "$PEEL_ALN" "$SCRIPT_DIR/sine_consensus.sh"
