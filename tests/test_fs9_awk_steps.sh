#!/usr/bin/env bash
# test_fs9_awk_steps.sh [WORKDIR] - the awk versions of two stage-9 steps give the same rows as the Python ones:
#   flankscan/fs9_repmask.sh rep|mask  vs  flankscan/fs9_repmask.py rep|mask
#   flankscan/fs9_seeds.sh             vs  flankscan/fs9_seeds.py   (| sort -u, as fs9_twins.sh runs them)
# Inputs: the two toys of tests/test_fs9_engines.sh (run first; its default-engine results are read from WORKDIR/toy1/new and
# WORKDIR/toy2/new), plus a KCOUNT variant with lower-case k-mers, a duplicated k-mer line and a one-field line.
# Every comparison runs with mawk and gawk (when on PATH), 1 and 4 workers, and FS9_Q 16 and 3 (residue classes of fs9_seeds.sh).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"; F="$SD/flankscan"
W=${1:-${TMPDIR:-/tmp}/fs9_engines}
export LC_ALL=C
bash "$HERE/test_fs9_engines.sh" "$W" > "$W.engines.log" 2>&1 || { cat "$W.engines.log"; echo "FAIL: test_fs9_engines.sh"; exit 1; }
AWKS=(); for a in mawk gawk; do command -v $a > /dev/null && AWKS+=("$(command -v $a)"); done
FAILS=0; N=0
same() {   # NAME FILE_PY FILE_AWK
    N=$((N + 1)); if cmp -s "$2" "$3"; then :; else echo "FAIL $1 ($(wc -l < "$2") vs $(wc -l < "$3") rows)"; FAILS=$((FAILS + 1)); fi
}
for toy in toy1 toy2; do
    D=$W/$toy/new; O=$W/$toy/awksteps; mkdir -p "$O"
    gawk -F'\t' -v OFS='\t' '{ print $1, $2, toupper($3) }' "$D/seqs.raw.tsv" > "$O/nomask.tsv"     # make_seqs 0
    CAPT=20; [[ $toy == toy2 ]] && CAPT=6
    if [[ $toy == toy2 ]]; then
        KC=$W/toy2/kcount.txt
        { cat "$KC"; gawk 'NR % 3 == 1 { print tolower($1), $2 } NR == 2 { print; print $1 }' "$KC"; } > "$O/kcount2.txt"
        for kc in "$KC" "$O/kcount2.txt"; do
            python3 "$F/fs9_repmask.py" rep "$O/nomask.tsv" "$kc" 20 > "$O/rep.py"
            for a in "${AWKS[@]}"; do for p in 1 4; do
                FS9_REP_AWK=$a bash "$F/fs9_repmask.sh" rep "$O/nomask.tsv" "$kc" 20 $p > "$O/rep.awk"
                same "$toy rep $(basename "$kc") $(basename "$a") P=$p" "$O/rep.py" "$O/rep.awk"; done; done
        done
        python3 "$F/fs9_repmask.py" rep "$O/nomask.tsv" "$KC" 20 > "$O/rep20.tsv"
        same "$toy rep vs fs9_twins.sh rep20.tsv" "$D/rep20.tsv" "$O/rep20.tsv"
    else
        cp "$D/rep20.tsv" "$O/rep20.tsv"                          # toy1: no KCOUNT, rep20.tsv counted from the genome
    fi
    python3 "$F/fs9_repmask.py" mask "$O/nomask.tsv" "$O/rep20.tsv" > "$O/masked.py"
    for a in "${AWKS[@]}"; do for p in 1 4; do
        FS9_REP_AWK=$a bash "$F/fs9_repmask.sh" mask "$O/nomask.tsv" "$O/rep20.tsv" $p > "$O/masked.awk"
        same "$toy mask $(basename "$a") P=$p" "$O/masked.py" "$O/masked.awk"; done; done
    for sq in nomask masked; do
        S=$O/$sq.tsv; [[ $sq == masked ]] && S=$O/masked.py
        for pm in "111111111111 6" "111010010100110111 3"; do
            set -- $pm
            python3 "$F/fs9_seeds.py" "$S" $1 $2 --cap $CAPT | sort -u > "$O/cand.py"
            [[ -s $O/cand.py ]] || echo "note: $toy $sq $1: no candidates"
            for a in "${AWKS[@]}"; do for p in 1 4; do for q in 16 3; do
                FS9_SEEDS_AWK=$a FS9_Q=$q bash "$F/fs9_seeds.sh" "$S" $1 $2 --cap $CAPT --threads $p | sort -u > "$O/cand.awk"
                same "$toy seeds $sq $1 $(basename "$a") P=$p Q=$q" "$O/cand.py" "$O/cand.awk"; done; done; done
        done
    done
done
(( FAILS == 0 )) && echo "test_fs9_awk_steps: $N comparisons, all identical" || { echo "test_fs9_awk_steps: $FAILS of $N comparisons differ"; exit 1; }
