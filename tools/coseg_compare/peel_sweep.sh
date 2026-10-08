#!/usr/bin/env bash
# Sweep the peel's settings on the SubFam chunks of finished run directories (no tool is re-run except the peel).
#   peel_sweep.sh OUT.tsv JOBS RUNDIR [RUNDIR ...]
# For every run dir with sf_n<N>/copies.clw (SubFam chunk consensuses) and every combination of
#   GRID_G  PEEL_GLOBAL_CONS  (columns where one base holds more than this share of the alignment are skipped; default 0.80)
#   GRID_M  PEEL_MIN_SET      (a feature must mark at least this many chunks; default 5)
#   GRID_B  PEEL_MIN_BLOCK    (a block needs this many co-occurring features; default 3)
#   GRID_J  PEEL_FEAT_JACCARD (features join a block above this overlap; default 0.45)
# runs peel_features.py and scores the copy-level partition against copies.labels.
# OUT.tsv: rundir, n, G, M, B, J, chunks, groups, placed, copies, ARI(placed), ARI(all copies), V(placed)
# The peel is by the owner (SINE-discriminator): its settings are environment variables so that they can be swept.
set -uo pipefail
OUT=${1:?out.tsv}; JOBS=${2:?jobs}; shift 2
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PEELPY=${PEELPY:?set PEELPY to peel_features.py}
PY=${PYTHON:-python3.12}
GRID_G=${GRID_G:-"0.80 0.90 0.95 0.99"}
GRID_M=${GRID_M:-"3 5 8"}
GRID_B=${GRID_B:-"3"}
GRID_J=${GRID_J:-"0.45"}
: > "$OUT.jobs"
for d in "$@"; do
    for c in "$d"/sf_n*/copies.clw; do
        [ -f "$c" ] || continue
        n=$(basename "$(dirname "$c")" | sed 's/sf_n//')
        for g in $GRID_G; do for m in $GRID_M; do for b in $GRID_B; do for j in $GRID_J; do
            echo "$d $n $g $m $b $j" >> "$OUT.jobs"
        done; done; done; done
    done
done
one() {
    d=$1; n=$2; g=$3; m=$4; b=$5; j=$6
    w="$d/sweep/n${n}_g${g}_m${m}_b${b}_j${j}"; mkdir -p "$w"
    chunks="$d/sf_n$n/copies.chunks.tsv"
    [ -f "$d/sf_n$n/chunk_truth.json" ] || $PY "$HERE/peel_labels.py" truth "$chunks" "$d/copies.labels" "$d/sf_n$n/chunk_truth.json"
    nch=$(cut -f2 "$chunks" | sort -u | wc -l)
    PEEL_GLOBAL_CONS=$g PEEL_MIN_SET=$m PEEL_MIN_BLOCK=$b PEEL_FEAT_JACCARD=$j $PY "$PEELPY" "$d/sf_n$n/copies.clw" "$w" "$d/sf_n$n/chunk_truth.json" > "$w/peel.log" 2>&1
    ncopies=$(wc -l < "$d/copies.labels")
    if [ -f "$w/peel_features.json" ]; then
        $PY "$HERE/peel_labels.py" groups "$chunks" "$w/peel_features.json" "$w/copy_groups.tsv" 2> /dev/null
    else : > "$w/copy_groups.tsv"; fi
    if [ -s "$w/copy_groups.tsv" ]; then
        row=$($PY "$HERE/score.py" "$d/copies.labels" "$w/copy_groups.tsv:peel" 2> /dev/null | tail -1)
        placed=$(echo "$row" | awk '{print $2}'); groups=$(echo "$row" | awk '{print $3}')
        v=$(echo "$row" | awk -F'|' '{split($2,x," "); print x[4]}'); ari=$(echo "$row" | awk -F'|' '{split($2,x," "); print x[5]}')
        ariall=$(echo "$row" | awk -F'|' '{split($4,x," "); print x[1]}')
    else placed=0; groups=0; v=0; ari=0; ariall=0; fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$d" "$n" "$g" "$m" "$b" "$j" "$nch" "$groups" "$placed" "$ncopies" "$ari" "$ariall" "$v"
}
export -f one
export HERE PEELPY PY
xargs -a "$OUT.jobs" -P "$JOBS" -L 1 bash -c 'one "$@"' _ > "$OUT.rows" 2> "$OUT.err"
{ printf 'rundir\tn\tG\tM\tB\tJ\tchunks\tgroups\tplaced\tcopies\tARI_placed\tARI_all\tV_placed\n'; sort "$OUT.rows"; } > "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT.rows") settings)"
