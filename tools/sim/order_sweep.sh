#!/usr/bin/env bash
# Chunk purity of SubFam against the true subfamilies for different orderings, on finished simulation dirs.
#   order_sweep.sh OUT.tsv JOBS CHUNK_N RUNDIR [RUNDIR ...]
# Orderings: k-mer tree with k = 4 5 6 8 10 (default k is 6), and the MAFFT guide tree (-m).
# OUT.tsv: rundir, ordering, chunks, purity, homogeneity, completeness, V, ARI (every copy is placed in a chunk).
set -uo pipefail
OUT=${1:?out.tsv}; JOBS=${2:?jobs}; N=${3:?chunk n}; shift 3
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBFAM=${SUBFAM:?set SUBFAM to the SubFam.sh to test}
PY=${PYTHON:-python3.12}
: > "$OUT.jobs"
for d in "$@"; do for o in k4 k5 k6 k8 k10 mafft; do echo "$d $o" >> "$OUT.jobs"; done; done
one() {
    d=$1; o=$2; w="$d/order_$o"; rm -rf "$w"; mkdir -p "$w"
    case $o in mafft) opt="-m" ;; k*) opt="-k ${o#k}" ;; esac
    bash "$SUBFAM" -n "$N" -t 4 $opt -x copies -o "$w" "$d/copies.fa" > "$w/log.txt" 2>&1 || { echo "FAILED $d $o" >&2; return; }
    row=$($PY "$HERE/../coseg_compare/score.py" "$d/copies.labels" "$w/copies.chunks.tsv:x" 2> /dev/null | tail -1)
    nch=$(cut -f2 "$w/copies.chunks.tsv" | sort -u | wc -l)
    pur=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[1]}'); h=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[2]}')
    c=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[3]}'); v=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[4]}')
    a=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[5]}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$d" "$o" "$nch" "$pur" "$h" "$c" "$v" "$a"
}
export -f one
export HERE SUBFAM PY N
xargs -a "$OUT.jobs" -P "$JOBS" -L 1 bash -c 'one "$@"' _ > "$OUT.rows" 2> "$OUT.err"
{ printf 'rundir\tordering\tchunks\tpurity\thomog\tcompl\tV\tARI\n'; sort "$OUT.rows"; } > "$OUT"
echo "wrote $OUT"
