#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# step1c_deplete.sh — deplete-and-resample SubFam round for rare subfamilies
#
# Usage:
#   step1c_deplete.sh <RUN_ROOT> [PERCENTILE=25] [BIN_SIZE=20] [SAMPLE=30000] [THREADS]
#
# Standalone/modular (like step7/8): runs against a completed run after step2
# and step3, changes nothing in it, and writes to RUN_ROOT/step1c/.
#
# Why: step1 samples SAMPLE copies (default 30,000) from extracted.fasta for
# SubFam. A subfamily with a few hundred genomic copies is a few dozen in that
# sample and cannot fill a SubFam chunk, so it never gets a consensus of its own
# and step2 assigns its copies to an abundant sister. This step removes the
# copies the current bank explains well and runs SubFam on what is left:
#
#   residual = unassigned copies (step2)
#            + assigned copies whose single-alignment similarity to their
#              subfamily consensus (step3 sim_ratio) is BELOW the PERCENTILE-th
#              percentile of that subfamily's copies
#
# i.e. only the close top (100-PERCENTILE) % of every subfamily is removed. Rare
# sister subfamilies sit only a few % below the sister's own copies, so a hard
# cut (5th percentile) throws most of them away; 25 keeps ~2/3 of them while
# shrinking the pool 4x. BIN_SIZE 20 (not 50): the residual is small and its
# rare members are few.
#
# Output (RUN_ROOT/step1c/):
#   residual.fasta        copies not explained by the bank
#   residual.tsv          seqID, subfamily-or-unassigned, sim_ratio, reason
#   subfam_input/input.fasta        the (sampled) residual given to SubFam
#   subfam_input/input.clw          SubFam chunk consensuses (unaligned, as in step1)
#   subfam_input/input.clw.al       the consensuses aligned (MAFFT L-INS-i)
#
# What to do with it: MANUAL.md §6.1 — review input.clw.al like step1's, pick the
# chunk consensuses that form a coherent new subfamily, add its consensus to the
# bank, re-run step2/step3, and run this step again until a round adds nothing.
# Most rows will be tails of known subfamilies, truncated or chimeric copies;
# a new subfamily shows as a block of several near-identical chunk consensuses.
#
# Requires: SubFam (same directory as this script, or on PATH), mafft, seqkit.
###############################################################################

RUN_ROOT="${1:-}"
PERCENTILE="${2:-25}"
BIN_SIZE="${3:-20}"
SAMPLE="${4:-30000}"
THREADS="${5:-$(nproc 2>/dev/null || echo 1)}"

if [[ -z "$RUN_ROOT" || "$RUN_ROOT" == "-h" || "$RUN_ROOT" == "--help" ]]; then
  sed -n '/^# Usage:/,/^# Requires/p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 1
fi
RUN_ROOT="$(readlink -f "$RUN_ROOT")"
[[ "$PERCENTILE" =~ ^[0-9]+$ && "$PERCENTILE" -ge 0 && "$PERCENTILE" -le 100 ]] || { echo "ERROR: PERCENTILE must be 0-100" >&2; exit 1; }
[[ "$BIN_SIZE" =~ ^[0-9]+$ && "$BIN_SIZE" -ge 2 ]] || { echo "ERROR: BIN_SIZE must be an integer >= 2" >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ log "ERROR: $*"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBFAM="$SCRIPT_DIR/SubFam"
[[ -x "$SUBFAM" ]] || SUBFAM="$(command -v SubFam || true)"
[[ -n "$SUBFAM" ]] || die "SubFam not found (expected next to this script or on PATH)"
for t in mafft seqkit; do command -v "$t" >/dev/null 2>&1 || die "$t not found on PATH"; done

EXTRACTED="$RUN_ROOT/genome.clean_step1/extracted.fasta"
STEP2_OUT="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
[[ -n "$STEP2_OUT" ]] || die "no step2 output under $RUN_ROOT/step2"
ASSIGN="$STEP2_OUT/assignment_full.tsv"     # seqID subfam thr votes status bitscore
SIM="$STEP2_OUT/sim_scores.tsv"             # seqID sim_bs self_bs sim_ratio   (step3)
UNASSIGNED="$STEP2_OUT/unassigned.fasta"
for f in "$EXTRACTED" "$ASSIGN" "$SIM" "$UNASSIGNED"; do
  [[ -f "$f" ]] || die "missing: $f (run step2 and step3 first)"
done

OUT="$RUN_ROOT/step1c"
mkdir -p "$OUT/subfam_input"
TMP="$(mktemp -d "$RUN_ROOT/.step1c_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

log "RUN_ROOT=$RUN_ROOT  percentile=$PERCENTILE  bin_size=$BIN_SIZE  sample=$SAMPLE  threads=$THREADS"
log "step2 output: $STEP2_OUT"

###############################################################################
# 1) Per-subfamily cut-off on sim_ratio, then the residual list
###############################################################################
# percentile over the assigned copies of each subfamily (nearest-rank)
awk -F'\t' -v OFS='\t' -v pct="$PERCENTILE" '
  FILENAME==ARGV[1] { if ($5=="assigned") sf[$1]=$2; next }          # ASSIGN
  FILENAME==ARGV[2] { if ($1 in sf && $4!="." && $4!="") { r[sf[$1]]=r[sf[$1]] " " $4 } ; next }   # SIM
  END {
    for (s in r) {
      n=split(substr(r[s],2), v, " ")
      # insertion sort is fine: a subfamily has at most a few hundred thousand copies, but
      # to stay cheap on big ones, bucket to 0.0001 resolution instead
      delete b; for (i=1;i<=n;i++) { k=int(v[i]*10000); b[k]++ }
      need=(pct*n)/100; if (need<1) need=1; acc=0; cut=""
      for (k=0; k<=10000 && cut==""; k++) { if (k in b) { acc+=b[k]; if (acc>=need) cut=k/10000 } }
      if (cut=="") cut=1
      print s, n, cut
    }
  }' "$ASSIGN" "$SIM" | sort -k1,1 > "$OUT/cutoffs.tsv"
log "per-subfamily cut-offs (subfam, n_assigned, sim_ratio at P$PERCENTILE):"
awk '{printf "  %-24s n=%-8s cut=%s\n", $1, $2, $3}' "$OUT/cutoffs.tsv" >&2

awk -F'\t' -v OFS='\t' '
  FILENAME==ARGV[1] { cut[$1]=$3; next }                              # cutoffs
  FILENAME==ARGV[2] { if ($5=="assigned") sf[$1]=$2; else status[$1]=$5; next }   # ASSIGN
  FILENAME==ARGV[3] { ratio[$1]=$4; next }                            # SIM
  /^>/ {                                                              # EXTRACTED headers
    id=substr($1,2)
    if (id in sf) {
      s=sf[id]; r=(id in ratio ? ratio[id] : ".")
      if (r=="." ) { print id, s, r, "no_similarity_score"; next }
      if (r+0 < cut[s]+0) print id, s, r, "below_P"
    } else {
      print id, (id in status ? status[id] : "not_in_step2"), ".", "unassigned"
    }
  }' "$OUT/cutoffs.tsv" "$ASSIGN" "$SIM" "$EXTRACTED" > "$OUT/residual.tsv"

cut -f1 "$OUT/residual.tsv" > "$TMP/residual.ids"
seqkit grep -f "$TMP/residual.ids" "$EXTRACTED" 2>/dev/null > "$OUT/residual.fasta"
n_all=$(grep -c '^>' "$EXTRACTED" || true)
n_res=$(grep -c '^>' "$OUT/residual.fasta" || true)
log "residual: $n_res of $n_all copies ($(awk -F'\t' '$4=="unassigned"' "$OUT/residual.tsv" | wc -l) unassigned, $(awk -F'\t' '$4=="below_P"' "$OUT/residual.tsv" | wc -l) below the cut)"
[[ "$n_res" -ge "$BIN_SIZE" ]] || die "residual too small for SubFam ($n_res < $BIN_SIZE)"

###############################################################################
# 2) SubFam on the residual (sampled to SAMPLE), as in step1
###############################################################################
if [[ "$n_res" -le "$SAMPLE" ]]; then
  cp "$OUT/residual.fasta" "$OUT/subfam_input/input.fasta"
else
  log "sampling $SAMPLE of $n_res residual copies"
  seqkit sample -n "$SAMPLE" -s 11 "$OUT/residual.fasta" > "$OUT/subfam_input/input.fasta" 2>/dev/null
fi

log "running SubFam (bin size $BIN_SIZE) on $(grep -c '^>' "$OUT/subfam_input/input.fasta") copies"
( cd "$OUT/subfam_input" && "$SUBFAM" input.fasta "$BIN_SIZE" > subfam_log.txt 2>&1 ) || die "SubFam failed, see $OUT/subfam_input/subfam_log.txt"
[[ -s "$OUT/subfam_input/input.clw" ]] || die "SubFam produced no input.clw"

n_cons=$(grep -c '^>' "$OUT/subfam_input/input.clw" || true)
log "aligning $n_cons chunk consensuses (MAFFT L-INS-i)"
mafft --thread "$THREADS" --threadit 0 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet \
  "$OUT/subfam_input/input.clw" > "$OUT/subfam_input/input.clw.al"
[[ -s "$OUT/subfam_input/input.clw.al" ]] || die "MAFFT produced empty output"

log "done: $OUT/subfam_input/input.clw.al ($n_cons rows). Review per MANUAL.md §6.1; a new subfamily is a block of several near-identical rows."
