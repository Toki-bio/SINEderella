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

# Divergence is measured against the consensus the copies were searched and assigned with. Until 2026-10-21 a bank
# rebuilt from the copies (consensuses.rebuilt.fa) was used when present; it could differ from the searched one (rsi
# r1_r2_r3_r3_C11: 889 bp rebuilt against 608 searched, 78 % identical), and the copies showed 23 % divergence instead
# of 6 %. The copies' own consensus is checked separately by the consensus audit (tools/consensus_audit.py).
CONS="$RUN_ROOT/consensuses.clean.fa"
SUBFAM_DIR="$OUT/subfamilies"
PLOTS_DIR="$OUT/plots"

[[ -s "$CONS" ]]        || die "missing/empty: $CONS"
[[ -d "$SUBFAM_DIR" ]]  || die "missing directory: $SUBFAM_DIR"
command -v ssearch36 >/dev/null 2>&1 || die "ssearch36 not in PATH"
command -v mafft     >/dev/null 2>&1 || die "mafft not in PATH"
command -v seqkit    >/dev/null 2>&1 || die "seqkit not in PATH"
# Prefer a newer python3.X with matplotlib already importable over a bare
# "python3" that may resolve to a system-ancient interpreter (e.g. 3.5) with
# no usable matplotlib wheel available for it. Falls back to plain python3
# unchanged if no better candidate is found, so this is a no-op on systems
# where python3 already has matplotlib.
PYTHON3="python3"
if ! python3 -c "import matplotlib" 2>/dev/null; then
  for cand in python3.13 python3.12 python3.11 python3.10 python3.9; do
    if command -v "$cand" >/dev/null 2>&1 && "$cand" -c "import matplotlib" 2>/dev/null; then
      PYTHON3="$cand"
      break
    fi
  done
fi
command -v "$PYTHON3" >/dev/null 2>&1 || die "$PYTHON3 not in PATH"
"$PYTHON3" -c "import matplotlib" 2>/dev/null || die "$PYTHON3 matplotlib not installed"

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

###############################################################################
# One subfamily: subsample, %identity (ssearch36), consensus profile (MAFFT --addfragments), plots
###############################################################################
# Runs in its own subshell with errexit on, as the old serial loop did: what stopped the old step (a failing seqkit,
# awk, sort or cp) stops this subfamily with a non-zero exit, and the parent then stops the step non-zero. What the old
# loop only warned about (no ssearch36 hits, empty MAFFT output, plot failure) is still only a warning.
process_sf(){
  local sf_fasta="$1" isf="$2"
  local sf_name nseqs copies_fa npctid cons_sz copies_sz m8_sz err_sz
  sf_name="$(basename "$sf_fasta" .fasta)"

  nseqs=$(count_fasta "$sf_fasta")
  [[ $nseqs -ge 1 ]] || { log "  Skipping $sf_name (0 sequences)"; return 0; }

  log "[$isf/$nsf] $sf_name ($nseqs copies)"

  # --- Extract single consensus for this subfamily ---
  awk -v sf="$sf_name" '
    /^>/{name=$1; sub(/^>/,"",name); found=(name==sf)}
    found{print}
  ' "$CONS" > "$tmpdir/cons_${sf_name}.fa"

  [[ -s "$tmpdir/cons_${sf_name}.fa" ]] || {
    log "  WARNING: consensus not found for $sf_name, skipping"
    return 0
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
    return 0
  }

  # ===================================================================
  # Plot 1: Divergence histogram from ssearch36 %identity
  # ===================================================================
  log "  Running ssearch36 for %identity (both strands, matches step2)..."
  # -z 3 (2026-11-21, docs/REPRODUCIBILITY.md; was -z 11): fixed statistics; the lines are sorted below because several
  # threads write them in another order from run to run and the first of two equal-score alignments was taken
  # bare file names inside $tmpdir (2026-11-23): ssearch36 cannot open a path longer than ~120 characters (FASTA36
  # file-name limit) - it truncates the name and reports "cannot open library ...". $tmpdir plus a long subfamily name
  # crosses that. Same pattern as tools/continuation_check.py and the satellite tools.
  ( cd "$tmpdir" && ssearch36 -g -3 -Q -n -z 3 -E 100 -T "$SF_THREADS" -m 8 \
      "cons_${sf_name}.fa" "copies_clean_${sf_name}.fa" ) \
    > "$tmpdir/sim_${sf_name}.m8" 2>"$tmpdir/ssearch_err_${sf_name}.txt" || true

  # Best hit per copy (highest bitscore): extract %identity (col3)
  awk -F'\t' '{
    seq=$2; pctid=$3+0; bs=$12+0
    if(!(seq in best) || bs>best[seq]){ best[seq]=bs; pid[seq]=pctid }
  } END {
    for(s in pid) printf "%s\t%.2f\n", s, pid[s]
  }' <(LC_ALL=C sort "$tmpdir/sim_${sf_name}.m8") | LC_ALL=C sort > "$tmpdir/pctid_${sf_name}.tsv"

  cp -f "$tmpdir/pctid_${sf_name}.tsv" "$PLOTS_DIR/${sf_name}_pctid.tsv"

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
    return 0
  fi

  # ===================================================================
  # Plot 2: Per-position nucleotide frequency from pairwise alignment
  # ===================================================================
  log "  Running MAFFT pairwise alignments for nucleotide profile..."

  # --addfragments: align each copy INDIVIDUALLY against the consensus
  # (pairwise, not a full MSA).  The consensus is the reference.
  mafft --addfragments "$tmpdir/copies_clean_${sf_name}.fa" \
    --adjustdirection --thread "$SF_THREADS" --quiet \
    "$tmpdir/cons_${sf_name}.fa" \
    > "$tmpdir/msa_${sf_name}.fa" 2>/dev/null || true

  [[ -s "$tmpdir/msa_${sf_name}.fa" ]] || {
    log "  WARNING: MAFFT produced empty output for $sf_name"
    return 0
  }

  # ===================================================================
  # Generate both plots via Python
  # ===================================================================
  log "  Generating plots..."
  "$PYTHON3" "$PLOT_PY" \
    --subfamily "$sf_name" \
    --pctid "$tmpdir/pctid_${sf_name}.tsv" \
    --msa "$tmpdir/msa_${sf_name}.fa" \
    --outdir "$PLOTS_DIR" \
    --consensus-name "$sf_name" \
  || log "  WARNING: plot generation failed for $sf_name"
}

###############################################################################
# Run the subfamilies JOBS at a time (2026-11-27). Each subfamily runs the old commands with the same arguments; only
# the thread count of ssearch36 and MAFFT changes (SF_THREADS = THREADS / JOBS, default 4). Neither output depends on
# it: the ssearch36 lines are sorted before the best hit is taken, and mafft --addfragments wrote the same alignment at
# --thread 1, 4, 8 and 24 (Ata e2-4, 10,000 copies). MAFFT is the slow part and hardly gains past 4 threads (e2-4: 245 s
# at 1, 150 s at 4, 145 s at 24), so several subfamilies at 4 threads each use the cores better than one at THREADS.
# STEP4_JOBS overrides the number of subfamilies at a time. Log lines of a subfamily stay together, in subfamily order.
###############################################################################
JOBS="${STEP4_JOBS:-$(( THREADS / 4 ))}"
(( JOBS >= 1 )) || JOBS=1
(( JOBS <= nsf )) || JOBS=$nsf
SF_THREADS=$(( THREADS / JOBS ))
(( SF_THREADS >= 1 )) || SF_THREADS=1
log "  $JOBS subfamilies at a time, $SF_THREADS threads each"

# print the finished subfamilies in order; a subfamily that exited non-zero is an error of the step
next_out=1
failed=0
flush_done(){
  local rc
  while (( next_out <= nsf )) && [[ -e "$tmpdir/job_${next_out}.rc" ]]; do
    cat "$tmpdir/job_${next_out}.out"
    cat "$tmpdir/job_${next_out}.err" >&2
    rc=$(< "$tmpdir/job_${next_out}.rc")
    if [[ "$rc" != 0 ]]; then
      log "ERROR: subfamily $(basename "${subfam_files[next_out-1]}" .fasta) failed (exit $rc)"
      failed=1
    fi
    next_out=$((next_out+1))
  done
}

trap 'kill $(jobs -pr) 2>/dev/null; wait; rm -rf "$tmpdir"' EXIT
isf=0
for sf_fasta in "${subfam_files[@]}"; do
  isf=$((isf+1))
  while (( $(jobs -pr | wc -l) >= JOBS )); do wait -n || true; flush_done; done
  (( failed == 0 )) || break
  (
    set +e
    ( set -euo pipefail; process_sf "$sf_fasta" "$isf" ) > "$tmpdir/job_${isf}.out" 2> "$tmpdir/job_${isf}.err"
    echo "$?" > "$tmpdir/job_${isf}.rc.part"
    mv "$tmpdir/job_${isf}.rc.part" "$tmpdir/job_${isf}.rc"
  ) &
done
wait
flush_done
(( failed == 0 )) || die "step 4 stopped: a subfamily failed (see the ERROR line above)"
(( next_out > nsf )) || die "step 4 stopped: only $((next_out-1)) of $nsf subfamilies finished"

log "All plots written to: $PLOTS_DIR"
log "Done."
