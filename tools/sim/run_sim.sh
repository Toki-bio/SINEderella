#!/usr/bin/env bash
# Simulated benchmark with a known subfamily tree. Run on KIT.
#   run_sim.sh OUTDIR [SEEDS="1 2 3"] [JOBS=4] [SCENARIOS=scenarios.tsv]
# Each (scenario, seed) gets OUTDIR/<scenario>_s<seed>/ with copies.fa, copies.labels, the tree, and the
# comparison of tools/coseg_compare/run_core.sh (COSEG, SubFam chunks, SubFam+peel). scores are collected by
# aggregate.py into OUTDIR/summary.tsv.
set -euo pipefail
OUT=${1:?outdir}
SEEDS=${2:-"1 2 3"}
JOBS=${3:-4}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SC=${4:-$HERE/scenarios.tsv}
PY=${PYTHON:-python3.12}
mkdir -p "$OUT" || exit 1
cd "$OUT" || exit 1
: > jobs.txt
while IFS=$'\t' read -r name args; do
    [ -n "$name" ] || continue
    for s in $SEEDS; do echo "$name $s $args" >> jobs.txt; done
done < "$SC"
one() {
    name=$1; seed=$2; shift 2
    d="$OUT/${name}_s$seed"
    [ -f "$d/scores.txt" ] && { echo "skip $d (done)"; return 0; }
    rm -rf "$d"; mkdir -p "$d"
    $PY "$HERE/sim_family.py" "$d/sim" --seed "$seed" "$@" 2> "$d/sim.log"
    cp "$d/sim.fa" "$d/copies.fa"; cp "$d/sim.labels" "$d/copies.labels"
    THREADS=${THREADS:-8} bash "$HERE/../coseg_compare/run_core.sh" "$d" > "$d/run.log" 2>&1 || echo "FAILED $d (see run.log)"
    echo "done $d"
}
export -f one
export OUT HERE PY
xargs -a jobs.txt -P "$JOBS" -L 1 bash -c 'one "$@"' _
$PY "$HERE/aggregate.py" "$OUT" > "$OUT/summary.tsv"
echo "summary: $OUT/summary.tsv"
