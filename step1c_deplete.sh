#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# step1c_deplete.sh — deplete-and-resample SubFam round for rare subfamilies
#
# Usage:
#   step1c_deplete.sh <RUN_ROOT> [PERCENTILE=25] [BIN_SIZE=20] [SAMPLE=30000] [THREADS] [MIN_LEN_FRAC=0.8]
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
# Length filter: sim_ratio is a bitscore ratio, so a truncated copy scores low
# even when every base it has matches. Only copies at least MIN_LEN_FRAC of their
# consensus's length can enter the residual; shorter ones are set aside (listed
# in set_aside.tsv, assignment untouched): a fragment cannot define a subfamily.
# The consensus length comes from RUN_ROOT/consensuses.clean.fa (an unassigned
# copy uses its best-vote consensus from step2, or the bank median).
#
# Output (RUN_ROOT/step1c/):
#   residual.fasta        copies not explained by the bank
#   residual.tsv          seqID, subfamily-or-unassigned, sim_ratio, reason
#   set_aside.tsv         copies kept out of the residual for being too short (same columns)
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
MIN_LEN_FRAC="${6:-0.8}"

if [[ -z "$RUN_ROOT" || "$RUN_ROOT" == "-h" || "$RUN_ROOT" == "--help" ]]; then
  sed -n '/^# Usage:/,/^# Requires/p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 1
fi
RUN_ROOT="$(readlink -f "$RUN_ROOT")"
[[ "$PERCENTILE" =~ ^[0-9]+$ && "$PERCENTILE" -ge 0 && "$PERCENTILE" -le 100 ]] || { echo "ERROR: PERCENTILE must be 0-100" >&2; exit 1; }
[[ "$BIN_SIZE" =~ ^[0-9]+$ && "$BIN_SIZE" -ge 2 ]] || { echo "ERROR: BIN_SIZE must be an integer >= 2" >&2; exit 1; }
awk -v f="$MIN_LEN_FRAC" 'BEGIN { exit !(f >= 0 && f <= 1) }' || { echo "ERROR: MIN_LEN_FRAC must be in [0, 1]" >&2; exit 1; }

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
BANK="$RUN_ROOT/consensuses.clean.fa"
for f in "$EXTRACTED" "$ASSIGN" "$SIM" "$UNASSIGNED" "$BANK"; do
  [[ -f "$f" ]] || die "missing: $f (run step2 and step3 first)"
done

OUT="$RUN_ROOT/step1c"
mkdir -p "$OUT/subfam_input"
TMP="$(mktemp -d "$RUN_ROOT/.step1c_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

log "RUN_ROOT=$RUN_ROOT  percentile=$PERCENTILE  bin_size=$BIN_SIZE  sample=$SAMPLE  threads=$THREADS  min_len_frac=$MIN_LEN_FRAC"
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

# consensus lengths (bank) and copy lengths (extracted), for the length filter
awk '/^>/ { if (n) print n "\t" l; n=substr($1,2); l=0; next } { l+=length($0) } END { if (n) print n "\t" l }' "$BANK" > "$TMP/cons_len.tsv"
awk '/^>/ { if (n) print n "\t" l; n=substr($1,2); l=0; next } { l+=length($0) } END { if (n) print n "\t" l }' "$EXTRACTED" > "$TMP/copy_len.tsv"

awk -F'\t' -v OFS='\t' -v frac="$MIN_LEN_FRAC" -v aside="$OUT/set_aside.tsv" '
  FILENAME==ARGV[1] { cut[$1]=$3; next }                              # cutoffs
  FILENAME==ARGV[2] { if ($5=="assigned") sf[$1]=$2; else { status[$1]=$5; best[$1]=$2 }; next }   # ASSIGN
  FILENAME==ARGV[3] { ratio[$1]=$4; next }                            # SIM
  FILENAME==ARGV[4] { clen[$1]=$2; cl[++nc]=$2; next }                # bank consensus lengths
  FILENAME==ARGV[5] { len[$1]=$2; next }                              # copy lengths
  function median(   i, j, t) { for (i=2;i<=nc;i++) { t=cl[i]; for (j=i-1; j>0 && cl[j]>t; j--) cl[j+1]=cl[j]; cl[j+1]=t } return cl[int((nc+1)/2)] }
  function longenough(id, s,   need) {
    need = (s in clen ? clen[s] : med) * frac
    return len[id] + 0 >= need
  }
  BEGIN { med = "" }
  /^>/ {                                                              # EXTRACTED headers
    if (med == "") med = (nc ? median() : 0)
    id=substr($1,2)
    if (id in sf) {
      s=sf[id]; r=(id in ratio ? ratio[id] : ".")
      if (r==".") { reason="no_similarity_score" }
      else if (r+0 < cut[s]+0) { reason="below_P" }
      else next
      if (longenough(id, s)) print id, s, r, reason; else print id, s, r, reason "_short" > aside
    } else {
      s=(id in best ? best[id] : "NA")
      if (longenough(id, s)) print id, (id in status ? status[id] : "not_in_step2"), ".", "unassigned"
      else print id, (id in status ? status[id] : "not_in_step2"), ".", "unassigned_short" > aside
    }
  }' "$OUT/cutoffs.tsv" "$ASSIGN" "$SIM" "$TMP/cons_len.tsv" "$TMP/copy_len.tsv" "$EXTRACTED" > "$OUT/residual.tsv"
[[ -f "$OUT/set_aside.tsv" ]] || : > "$OUT/set_aside.tsv"

cut -f1 "$OUT/residual.tsv" > "$TMP/residual.ids"
seqkit grep -f "$TMP/residual.ids" "$EXTRACTED" 2>/dev/null > "$OUT/residual.fasta"
n_all=$(grep -c '^>' "$EXTRACTED" || true)
n_res=$(grep -c '^>' "$OUT/residual.fasta" || true)
log "residual: $n_res of $n_all copies ($(awk -F'\t' '$4=="unassigned"' "$OUT/residual.tsv" | wc -l) unassigned, $(awk -F'\t' '$4=="below_P"' "$OUT/residual.tsv" | wc -l) below the cut); $(wc -l < "$OUT/set_aside.tsv") set aside as shorter than $MIN_LEN_FRAC x their consensus"
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
