#!/usr/bin/env bash
# Sweep the peel's settings on the SubFam chunks of finished run directories (no tool is re-run except the peel).
#   peel_sweep.sh OUT.tsv JOBS RUNDIR [RUNDIR ...]
# For every run dir with sf_n<N>/copies.clw (SubFam chunk consensuses) and every setting of PEEL_GLOBAL_CONS x PEEL_MIN_SET
# (override with GRID_G, GRID_M), runs peel_features.py and scores the copy-level partition against copies.labels.
# OUT.tsv: rundir, n, global_cons, min_set, chunks, groups, placed, copies, ARI(placed), ARI(all copies), V(placed)
# The peel is by the owner (SINE-discriminator): its settings are environment variables precisely so that they can be swept.
set -uo pipefail
OUT=${1:?out.tsv}; JOBS=${2:?jobs}; shift 2
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PEELPY=${PEELPY:?set PEELPY to peel_features.py}
PY=${PYTHON:-python3.12}
GRID_G=${GRID_G:-"0.80 0.90 0.95 0.99"}
GRID_M=${GRID_M:-"3 5 8"}
: > "$OUT.jobs"
for d in "$@"; do
    for c in "$d"/sf_n*/copies.clw; do
        [ -f "$c" ] || continue
        n=$(basename "$(dirname "$c")" | sed 's/sf_n//')
        for g in $GRID_G; do for m in $GRID_M; do echo "$d $n $g $m" >> "$OUT.jobs"; done; done
    done
done
one() {
    d=$1; n=$2; g=$3; m=$4
    w="$d/sweep/n${n}_g${g}_m${m}"; mkdir -p "$w"
    chunks="$d/sf_n$n/copies.chunks.tsv"
    [ -f "$d/sf_n$n/chunk_truth.json" ] || $PY "$HERE/peel_labels.py" truth "$chunks" "$d/copies.labels" "$d/sf_n$n/chunk_truth.json"
    nch=$(cut -f2 "$chunks" | sort -u | wc -l)
    PEEL_GLOBAL_CONS=$g PEEL_MIN_SET=$m $PY "$PEELPY" "$d/sf_n$n/copies.clw" "$w" "$d/sf_n$n/chunk_truth.json" > "$w/peel.log" 2>&1
    ncopies=$(wc -l < "$d/copies.labels")
    if [ -f "$w/peel_features.json" ]; then
        $PY "$HERE/peel_labels.py" groups "$chunks" "$w/peel_features.json" "$w/copy_groups.tsv" 2> /dev/null
    else : > "$w/copy_groups.tsv"; fi
    row=$($PY "$HERE/score.py" "$d/copies.labels" "$w/copy_groups.tsv:peel" 2> /dev/null | tail -1)
    if [ -s "$w/copy_groups.tsv" ]; then
        # row: peel placed groups | purity homog compl V ARI | ARIcom Vcom grpcom | ARIall
        placed=$(echo "$row" | awk '{print $2}'); groups=$(echo "$row" | awk '{print $3}')
        v=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[4]}'); ari=$(echo "$row" | awk -F'|' '{split($2,b," "); print b[5]}')
        ariall=$(echo "$row" | awk -F'|' '{split($4,d," "); print d[1]}')
    else placed=0; groups=0; v=0; ari=0; ariall=0; fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$d" "$n" "$g" "$m" "$nch" "$groups" "$placed" "$ncopies" "$ari" "$ariall" "$v"
}
export -f one
export HERE PEELPY PY
xargs -a "$OUT.jobs" -P "$JOBS" -L 1 bash -c 'one "$@"' _ > "$OUT.rows" 2> "$OUT.err"
{ printf 'rundir\tn\tglobal_cons\tmin_set\tchunks\tgroups\tplaced\tcopies\tARI_placed\tARI_all\tV_placed\n'; sort "$OUT.rows"; } > "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT.rows") settings)"
