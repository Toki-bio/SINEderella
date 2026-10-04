#!/bin/bash
# step1b_cluster_subfamilies_assist.sh — First-pass subfamily clustering ASSIST
# Usage: bash step1b_cluster_subfamilies_assist.sh input.clw output_dir [--recurse]
#
# Optional step between SubFam (which produces input.clw, a collection of
# per-chunk consensuses, NOT a converged alignment — see MANUAL.md §6.1.1) and
# the manual subfamily-clustering step (§6.1.2).
#
# Pipeline:
#   1. Degap input.clw and realign properly with mafft (§6.1.1's fix).
#   2. Run cluster_assist.js (Node, vendors MSA-viewer's SINEClusterer) as a
#      first-pass automated clustering pass.
#
# THIS IS AN ASSIST, NOT A REPLACEMENT for manual review — validated (2026-08-21)
# to reliably resolve only the coarsest structure (e.g. a true 2-subfamily split)
# and to under-split real multi-subfamily datasets (tested on a 9-subfamily
# ground-truth case: collapsed to 2 blobs + 1 tiny cluster, and recursive
# re-clustering of the blobs did not recover further structure). Always inspect
# summary.tsv cluster sizes against your prior expectation, and fall back to
# manual clustering (§6.1.2) whenever the split looks too coarse.
#
# Requires: seqkit, mafft, node (>=14)
#
# Output (in output_dir):
#   input.realigned.fasta      — properly realigned input (see MANUAL.md 6.1.1)
#   clusters/summary.tsv       — cluster_id, size
#   clusters/cluster_N.fasta, cluster_N_consensus.fasta, unassigned.fasta

set -e

INPUT_CLW=${1:?Usage: bash step1b_cluster_subfamilies_assist.sh input.clw output_dir [--recurse]}
OUT_DIR=${2:?Usage: bash step1b_cluster_subfamilies_assist.sh input.clw output_dir [--recurse]}
shift 2 || true
EXTRA_ARGS="$@"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT_DIR"

echo "[1/2] Degapping and realigning $INPUT_CLW ..."
seqkit seq -w 0 "$INPUT_CLW" | seqkit seq -g > "$OUT_DIR/input.degapped.fasta"
mafft --auto --quiet "$OUT_DIR/input.degapped.fasta" > "$OUT_DIR/input.realigned.fasta"

echo "[2/2] Running cluster_assist.js (first-pass ASSIST, review before trusting) ..."
node "$SCRIPT_DIR/cluster_assist.js" "$OUT_DIR/input.realigned.fasta" "$OUT_DIR/clusters" $EXTRA_ARGS

echo ""
echo "Done. Inspect $OUT_DIR/clusters/summary.tsv — if cluster count/sizes look"
echo "too coarse relative to what you expect, fall back to manual clustering"
echo "per MANUAL.md section 6.1.2."
