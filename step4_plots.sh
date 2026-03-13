#!/usr/bin/env bash
# ==========================================================
# step4_plots.sh - Generate per-subfamily diagnostic plots
# ==========================================================
#   Plot 1: Histogram of % divergence from consensus per copy
#   Plot 2: Per-position nucleotide frequency (stacked bar chart)
#
# Usage:  step4_plots.sh <RUN_ROOT> [THREADS] [MAX_SEQS]
#
# Output: <step2_output>/plots/
#           <subfamily>_divergence.png
#           <subfamily>_divergence.pdf
#           <subfamily>_nucfreq.png
#           <subfamily>_nucfreq.pdf
# ==========================================================
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

RUN_ROOT="${1:-}"
[[ -n "$RUN_ROOT" ]] || die "usage: $0 RUN_ROOT [THREADS] [MAX_SEQS]"
RUN_ROOT="$(readlink -f "$RUN_ROOT")"

THREADS="${2:-$(nproc 2>/dev/null || echo 1)}"
MAX_SEQS="${3:-10000}"

# Locate step2 output
OUT="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
[[ -n "$OUT" && -d "$OUT" ]] || die "cannot find $RUN_ROOT/step2/step2_output*"

CONS="$RUN_ROOT/consensuses.clean.fa"
SUBFAM_DIR="$OUT/subfamilies"
PLOTS_DIR="$OUT/plots"

[[ -s "$CONS" ]]        || die "missing/empty: $CONS"
[[ -d "$SUBFAM_DIR" ]]  || die "missing directory: $SUBFAM_DIR"
command -v ssearch36 >/dev/null 2>&1 || die "ssearch36 not in PATH"
command -v mafft     >/dev/null 2>&1 || die "mafft not in PATH"
command -v seqkit    >/dev/null 2>&1 || die "seqkit not in PATH"
command -v python3   >/dev/null 2>&1 || die "python3 not in PATH"
python3 -c "import matplotlib" 2>/dev/null || die "python3 matplotlib not installed"

# Log tool versions for diagnostics
SSEARCH_VER=$(ssearch36 -h 2>&1 | grep -i 'version\|ssearch' | head -1 || echo "unknown")
log "ssearch36: $(command -v ssearch36) [$SSEARCH_VER]"

mkdir -p "$PLOTS_DIR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT_PY="$SCRIPT_DIR/plot_subfamily.py"
[[ -f "$PLOT_PY" ]] || die "missing Python helper: $PLOT_PY"

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/plots_tmp.XXXXXX")
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

count_fasta(){
  awk 'BEGIN{c=0} /^>/{c++} END{printf("%d\n", c)}' "$1"
}

###############################################################################
# Process each subfamily
###############################################################################
mapfile -t subfam_files < <(find "$SUBFAM_DIR" -maxdepth 1 -name "*.fasta" -size +0c | sort)
nsf=${#subfam_files[@]}
[[ $nsf -gt 0 ]] || die "No non-empty subfamily FASTA files in $SUBFAM_DIR"

log "Processing $nsf subfamilies (max_seqs=$MAX_SEQS, threads=$THREADS)"
log "  RUN_ROOT=$RUN_ROOT"
log "  OUT=$OUT"
log "  CONS=$CONS"
log "  SUBFAM_DIR=$SUBFAM_DIR"

isf=0
for sf_fasta in "${subfam_files[@]}"; do
  sf_name="$(basename "$sf_fasta" .fasta)"
  isf=$((isf+1))

  nseqs=$(count_fasta "$sf_fasta")
  [[ $nseqs -ge 1 ]] || { log "  Skipping $sf_name (0 sequences)"; continue; }

  log "[$isf/$nsf] $sf_name ($nseqs copies)"

  # --- Extract single consensus for this subfamily ---
  awk -v sf="$sf_name" '
    /^>/{name=$1; sub(/^>/,"",name); found=(name==sf)}
    found{print}
  ' "$CONS" > "$tmpdir/cons_${sf_name}.fa"

  [[ -s "$tmpdir/cons_${sf_name}.fa" ]] || {
    log "  WARNING: consensus not found for $sf_name, skipping"
    continue
  }

  # --- Subsample if needed ---
  copies_fa="$sf_fasta"
  if [[ $nseqs -gt $MAX_SEQS ]]; then
    log "  Subsampling $nseqs -> $MAX_SEQS copies"
    seqkit sample -n "$MAX_SEQS" -s 42 "$sf_fasta" > "$tmpdir/sampled_${sf_name}.fa" 2>/dev/null
    copies_fa="$tmpdir/sampled_${sf_name}.fa"
    nseqs=$MAX_SEQS
  fi

  # --- Strip pipe-delimited fields from copy headers for clean IDs ---
  awk '/^>/{
    h=$0; sub(/^>/,"",h)
    split(h,p,"|"); print ">"p[1]
    next
  } {print}' "$copies_fa" > "$tmpdir/copies_clean_${sf_name}.fa"

  [[ -s "$tmpdir/copies_clean_${sf_name}.fa" ]] || {
    log "  WARNING: copies_clean is empty for $sf_name, skipping"
    continue
  }

  # ===================================================================
  # Plot 1: Divergence histogram from ssearch36 %identity
  # ===================================================================
  log "  Running ssearch36 for %identity..."
  ssearch36 -Q -n -z 11 -E 100 -T "$THREADS" -m 8 \
    "$tmpdir/cons_${sf_name}.fa" "$tmpdir/copies_clean_${sf_name}.fa" \
    > "$tmpdir/sim_${sf_name}.m8" 2>"$tmpdir/ssearch_err_${sf_name}.txt" || true

  # Best hit per copy (highest bitscore): extract %identity (col3)
  awk -F'\t' '{
    seq=$2; pctid=$3+0; bs=$12+0
    if(!(seq in best) || bs>best[seq]){ best[seq]=bs; pid[seq]=pctid }
  } END {
    for(s in pid) printf "%s\t%.2f\n", s, pid[s]
  }' "$tmpdir/sim_${sf_name}.m8" > "$tmpdir/pctid_${sf_name}.tsv"

  npctid=$(wc -l < "$tmpdir/pctid_${sf_name}.tsv")
  if [[ $npctid -eq 0 ]]; then
    log "  WARNING: no ssearch36 hits for $sf_name"
    # --- Diagnostic info ---
    cons_sz=$(wc -c < "$tmpdir/cons_${sf_name}.fa" 2>/dev/null || echo 0)
    copies_sz=$(wc -c < "$tmpdir/copies_clean_${sf_name}.fa" 2>/dev/null || echo 0)
    m8_sz=$(wc -c < "$tmpdir/sim_${sf_name}.m8" 2>/dev/null || echo 0)
    err_sz=$(wc -c < "$tmpdir/ssearch_err_${sf_name}.txt" 2>/dev/null || echo 0)
    log "    consensus file:  $cons_sz bytes  copies file: $copies_sz bytes"
    log "    m8 output: $m8_sz bytes  stderr: $err_sz bytes"
    if [[ -s "$tmpdir/ssearch_err_${sf_name}.txt" ]]; then
      log "    ssearch36 stderr: $(head -3 "$tmpdir/ssearch_err_${sf_name}.txt" | tr '\n' ' ')"
    fi
    log "    consensus hdr: $(head -1 "$tmpdir/cons_${sf_name}.fa")"
    log "    1st copy hdr:  $(head -1 "$tmpdir/copies_clean_${sf_name}.fa")"
    continue
  fi

  # ===================================================================
  # Plot 2: Per-position nucleotide frequency from pairwise alignment
  # ===================================================================
  log "  Running MAFFT pairwise alignments for nucleotide profile..."

  # --addfragments: align each copy INDIVIDUALLY against the consensus
  # (pairwise, not a full MSA).  The consensus is the reference.
  mafft --addfragments "$tmpdir/copies_clean_${sf_name}.fa" \
    --adjustdirection --thread "$THREADS" --quiet \
    "$tmpdir/cons_${sf_name}.fa" \
    > "$tmpdir/msa_${sf_name}.fa" 2>/dev/null || true

  [[ -s "$tmpdir/msa_${sf_name}.fa" ]] || {
    log "  WARNING: MAFFT produced empty output for $sf_name"
    continue
  }

  # ===================================================================
  # Generate both plots via Python
  # ===================================================================
  log "  Generating plots..."
  python3 "$PLOT_PY" \
    --subfamily "$sf_name" \
    --pctid "$tmpdir/pctid_${sf_name}.tsv" \
    --msa "$tmpdir/msa_${sf_name}.fa" \
    --outdir "$PLOTS_DIR" \
    --consensus-name "$sf_name" \
  || log "  WARNING: plot generation failed for $sf_name"

done

log "All plots written to: $PLOTS_DIR"
log "Done."
