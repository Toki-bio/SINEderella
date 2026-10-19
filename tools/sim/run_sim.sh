#!/usr/bin/env bash
# Simulated benchmark with a known subfamily tree. Run on KIT.
#   run_sim.sh OUTDIR [SEEDS="1 2 3"] [JOBS=4] [SCENARIOS=scenarios.tsv]
# scenarios.tsv: name <TAB> simulator arguments [<TAB> seeds for this scenario only]
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
OUT=$(cd "$OUT" && pwd)     # relative OUTDIR gave OUTDIR/OUTDIR paths after the cd below, until 2026-10-19 (GLM audit 2026-10-19)
cd "$OUT" || exit 1
: > jobs.txt
while IFS=$'\t' read -r name args seeds; do
    [ -n "$name" ] || continue
    for s in ${seeds:-$SEEDS}; do echo "$name $s $args" >> jobs.txt; done
done < "$SC"
one() {
    name=$1; seed=$2; shift 2
    d="$OUT/${name}_s$seed"
    [ -f "$d/scores.txt" ] && { echo "skip $d (done)"; return 0; }
    rm -rf "$d"; mkdir -p "$d"
    # a failed simulation or comparison is recorded in FAILED (read by aggregate.py) and in the exit status; until 2026-10-19
    # run_core.sh ran on whatever sim_family.py left and the seed vanished from the summary (GLM audit 2026-10-19)
    $PY "$HERE/sim_family.py" "$d/sim" --seed "$seed" "$@" 2> "$d/sim.log" || { echo "sim_family.py failed" > "$d/FAILED"; echo "FAILED $d (see sim.log)"; return 1; }
    cp "$d/sim.fa" "$d/copies.fa"; cp "$d/sim.labels" "$d/copies.labels"
    THREADS=${THREADS:-8} bash "$HERE/../coseg_compare/run_core.sh" "$d" > "$d/run.log" 2>&1 || { echo "run_core.sh failed" > "$d/FAILED"; echo "FAILED $d (see run.log)"; return 1; }
    echo "done $d"
}
export -f one
export OUT HERE PY
xargs -a jobs.txt -P "$JOBS" -L 1 bash -c 'one "$@"' _; rc=$?
$PY "$HERE/aggregate.py" "$OUT" > "$OUT/summary.tsv"
echo "summary: $OUT/summary.tsv"
[ $rc -eq 0 ] || { echo "some runs FAILED (see */FAILED)"; exit 1; }
