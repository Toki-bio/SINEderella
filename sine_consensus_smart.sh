#!/bin/bash
# sine_consensus_smart.sh — Enhanced version with variance detection and early stopping
#
# Additions to default sine_consensus.sh:
#   - Detects high-variance families (garbage/FP-contaminated)
#   - Implements early stopping if variance doesn't decrease
#   - Flags low-confidence output
#
# Usage:  bash sine_consensus_smart.sh [options] input.fasta [basename]

set -euo pipefail

# ── Defaults ───────────────────────────────────────────────────────────────
SUBSAMPLE_SIZE=100
MAX_ITERS=50
STABILITY_THRESH=0.005
CONS_THRESH=50       # percent
MIN_COVERAGE=30      # percent
KEEP=0
MIN_ITERS=10
VARIANCE_WINDOW=3    # detect variance over last N iterations
VARIANCE_LIMIT=0.15  # if avg change in window > 15%, likely bad data
EARLY_STOP=1         # enable early stopping on high variance

usage() {
    sed -n '2,/^$/s/^# //p' "$0"
    exit 0
}

while getopts "n:m:s:t:c:kh" opt; do
    case $opt in
        n) SUBSAMPLE_SIZE=$OPTARG ;;
        m) MAX_ITERS=$OPTARG ;;
        s) STABILITY_THRESH=$OPTARG ;;
        t) CONS_THRESH=$OPTARG ;;
        c) MIN_COVERAGE=$OPTARG ;;
        k) KEEP=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUT_FASTA=${1:?Error: input FASTA required}
[ -f "$INPUT_FASTA" ] || { echo "Error: $INPUT_FASTA not found" >&2; exit 1; }

BASENAME=$(basename "$INPUT_FASTA" .fasta)
WORKDIR=$(mktemp -d "${BASENAME}_cons_XXXXXX")
LOG="${BASENAME}_consensus.log"

cleanup() {
    rm -rf "$WORKDIR"
    [ "$KEEP" -eq 0 ] && rm -f "$LOG"
}
trap cleanup EXIT

progress() {
    printf '\r\033[K%s' "$*" >&2
}

random_subsample() {
    local INFILE=$1 N=$2 OUTFILE=$3
    awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
         { printf("%s",$0) }
         END { printf("\n") }' "$INFILE" \
        | shuf > "$WORKDIR/shuffled.tmp"
    head -n "$N" "$WORKDIR/shuffled.tmp" \
        | awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' > "$OUTFILE"
    rm -f "$WORKDIR/shuffled.tmp"
}

compute_consensus() {
    local ALN=$1 CONS=$2
    awk -v thresh="$CONS_THRESH" -v mincov="$MIN_COVERAGE" '
    /^>/ {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i,c]++
            }
        }
        seq = ""
        next
    }
    { seq = seq $0 }
    END {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i,c]++
            }
        }
        thr = thresh / 100.0
        cov = mincov / 100.0
        printf ">consensus\n"
        for (i = 1; i <= maxcol; i++) {
            total = nseq
            gap = (count[i,"-"] + 0) + (count[i,"."] + 0)
            nongap = total - gap
            if (nongap < total * cov) { printf "-"; continue }
            split("A C G T", bases, " ")
            max_count = 0
            best = "-"
            for (b = 1; b <= 4; b++) {
                base = bases[b]
                c = count[i,base] + 0
                if (c > max_count) {
                    max_count = c
                    best = base
                } else if (c == max_count && c > 0) {
                    if (base < best) best = base
                }
            }
            freq = max_count / total
            if (freq >= thr) {
                printf best
            } else {
                printf "-"
            }
        }
        print ""
    }
    ' "$ALN" > "$CONS"
}

approx_distance() {
    awk '
    NR==FNR { if (/^>/) next; seq1 = seq1 $0; next }
    /^>/ { next }
    { seq2 = seq2 $0 }
    END {
        len1 = length(seq1); len2 = length(seq2)
        if (len1 == 0 || len2 == 0) { print 1.0; exit }
        lendiff = (len1 > len2 ? len1 - len2 : len2 - len1)
        if (lendiff > len1 * 0.1) { print 1.0; exit }
        minlen = (len1 < len2 ? len1 : len2)
        diff = 0
        for (i = 1; i <= minlen; i++)
            if (substr(seq1, i, 1) != substr(seq2, i, 1)) diff++
        print diff / len1
    }
    ' "$1" "$2"
}

# ── Init ───────────────────────────────────────────────────────────────────
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA")
QUALITY_FLAG=""

{
    echo "=== sine_consensus_smart.sh ==="
    echo "Input: $INPUT_FASTA  ($SEQ_COUNT sequences)"
    echo "Params: subsample=$SUBSAMPLE_SIZE  max_iters=$MAX_ITERS  thresh=$CONS_THRESH%  min_cov=$MIN_COVERAGE%"
    echo "Variance detection: window=$VARIANCE_WINDOW iters, limit=$VARIANCE_LIMIT"
    echo "---"
} > "$LOG"

if [ "$SEQ_COUNT" -lt 3 ]; then
    echo "Error: need at least 3 sequences, got $SEQ_COUNT" >&2
    exit 1
fi

[ "$SUBSAMPLE_SIZE" -gt "$SEQ_COUNT" ] && SUBSAMPLE_SIZE=$SEQ_COUNT

MASTER_RAW="$WORKDIR/master_raw.fasta"
MASTER_ALN="$WORKDIR/master.aln.fasta"
PREV_CONS="$WORKDIR/prev_cons.fasta"
CURR_CONS="$WORKDIR/curr_cons.fasta"

# ── Iteration 1 ───────────────────────────────────────────────────────────
progress "Iter 1/$MAX_ITERS: subsampling..."
random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" "$WORKDIR/sub_1.fasta"
mafft --auto --quiet "$WORKDIR/sub_1.fasta" > "$WORKDIR/aln_1.fasta"
compute_consensus "$WORKDIR/aln_1.fasta" "$WORKDIR/mini_1.fasta"

cp "$WORKDIR/mini_1.fasta" "$MASTER_RAW"
cp "$WORKDIR/mini_1.fasta" "$MASTER_ALN"
cp "$WORKDIR/mini_1.fasta" "$PREV_CONS"
echo "Iter  1: init" >> "$LOG"

HAS_ALISTAT=0
command -v esl-alistat >/dev/null 2>&1 && HAS_ALISTAT=1

# ── Main loop with variance detection ───────────────────────────────────────
CONVERGED=0
CHANGES=()  # Track recent changes for variance

for ITER in $(seq 2 "$MAX_ITERS"); do
    progress "Iter $ITER/$MAX_ITERS: subsampling..."
    random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" "$WORKDIR/sub_${ITER}.fasta"

    progress "Iter $ITER/$MAX_ITERS: aligning subsample..."
    mafft --auto --quiet "$WORKDIR/sub_${ITER}.fasta" > "$WORKDIR/aln_${ITER}.fasta"

    compute_consensus "$WORKDIR/aln_${ITER}.fasta" "$WORKDIR/mini_${ITER}.fasta"
    cat "$WORKDIR/mini_${ITER}.fasta" >> "$MASTER_RAW"

    progress "Iter $ITER/$MAX_ITERS: aligning master..."
    mafft --auto --quiet "$MASTER_RAW" > "$MASTER_ALN"

    compute_consensus "$MASTER_ALN" "$CURR_CONS"

    # Distance
    CHANGE=$(approx_distance "$PREV_CONS" "$CURR_CONS")
    CHANGES+=("$CHANGE")

    LOG_LINE="Iter $(printf '%2d' "$ITER"): change=$CHANGE"
    if [ "$HAS_ALISTAT" -eq 1 ]; then
        AVGID=$(esl-alistat "$MASTER_ALN" 2>/dev/null | awk '/Average identity/ {print $NF}' || echo "?")
        LOG_LINE="$LOG_LINE  avgid=$AVGID"
    fi
    echo "$LOG_LINE" >> "$LOG"

    progress "Iter $ITER/$MAX_ITERS: change=$CHANGE"

    # ── Convergence check (standard threshold) ──────────────────────────────
    if [ "$ITER" -ge "$MIN_ITERS" ] && [ "$(echo "$CHANGE < $STABILITY_THRESH" | bc -l)" -eq 1 ]; then
        echo "Converged at iter $ITER (change=$CHANGE < $STABILITY_THRESH)" >> "$LOG"
        CONVERGED=1
        break
    fi

    # ── Variance detection (early stopping) ─────────────────────────────────
    if [ "$EARLY_STOP" -eq 1 ] && [ "$ITER" -ge $((MIN_ITERS + VARIANCE_WINDOW)) ]; then
        # Calculate average change over last VARIANCE_WINDOW iterations
        window_avg=0
        count=0
        for (( i=${#CHANGES[@]}-VARIANCE_WINDOW; i<${#CHANGES[@]}; i++ )); do
            if [ $i -ge 0 ]; then
                window_avg=$(echo "$window_avg + ${CHANGES[$i]}" | bc -l)
                count=$((count + 1))
            fi
        done
        
        if [ "$count" -gt 0 ]; then
            window_avg=$(echo "$window_avg / $count" | bc -l)
            if [ "$(echo "$window_avg > $VARIANCE_LIMIT" | bc -l)" -eq 1 ]; then
                echo "WARNING: High variance detected (avg change=$window_avg > $VARIANCE_LIMIT)" >> "$LOG"
                echo "Early stopping — likely garbage/FP-contaminated data" >> "$LOG"
                QUALITY_FLAG="LOW_CONFIDENCE"
                break
            fi
        fi
    fi

    cp "$CURR_CONS" "$PREV_CONS"
    rm -f "$WORKDIR/sub_${ITER}.fasta" "$WORKDIR/aln_${ITER}.fasta"
done

progress ""
echo "" >&2

[ "$CONVERGED" -eq 0 ] && echo "Warning: did not converge within $MAX_ITERS iterations" >> "$LOG"
[ -n "$QUALITY_FLAG" ] && echo "QUALITY: $QUALITY_FLAG" >> "$LOG"
[ ! -f "$CURR_CONS" ] && cp "$PREV_CONS" "$CURR_CONS"

# ── Output ─────────────────────────────────────────────────────────────────
OUTFILE="${BASENAME}_consensus.fasta"
printf ">%s\n" "$BASENAME" > "$OUTFILE"
grep -v '^>' "$CURR_CONS" | tr -d '-' >> "$OUTFILE"
CONS_LEN=$(grep -v '^>' "$OUTFILE" | tr -d '\n' | wc -c)
echo "Output: $OUTFILE ($CONS_LEN bp)" >> "$LOG"

# ── Trace ──────────────────────────────────────────────────────────────────
if [ "$KEEP" -eq 1 ]; then
    TRACE="${BASENAME}_trace.aln.fasta"
    > "$TRACE"
    for f in $(ls "$WORKDIR"/mini_*.fasta 2>/dev/null | sort -t_ -k2 -n); do
        ITER_NUM=$(basename "$f" | grep -oE '[0-9]+')
        awk -v n="$ITER_NUM" '/^>/{print ">iter_"n; next}{print}' "$f" >> "$TRACE"
    done
    awk -v name="$BASENAME" '/^>/{print ">FINAL_"name; next}{print}' "$CURR_CONS" >> "$TRACE"
    echo "Trace: $TRACE" >> "$LOG"
fi

# ── Summary to stderr ──────────────────────────────────────────────────────
if [ -n "$QUALITY_FLAG" ]; then
    echo "$BASENAME: $QUALITY_FLAG, $CONS_LEN bp (stopped at iter $ITER)" >&2
elif [ "$CONVERGED" -eq 1 ]; then
    echo "$BASENAME: converged at iter $ITER, $CONS_LEN bp" >&2
else
    echo "$BASENAME: NOT converged after $MAX_ITERS iters, $CONS_LEN bp" >&2
fi
