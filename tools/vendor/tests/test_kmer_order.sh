#!/usr/bin/env bash
# Byte-for-byte test of kmer_order (C) against the Python/numpy oracle that SubFam.sh carries
# (KMER_ORDER_PY, extracted from SubFam.sh itself so the oracle cannot drift from the shipped code).
#
#   test_kmer_order.sh WORKDIR [REAL_FASTA ...]
#
# Needs: cc, python3 >= 3.6 with numpy (set PYTHON=). Extra arguments are real FASTA files, each
# tested with k=6 canonical+orient and k=8 plain. Exits 1 on the first difference.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
VENDOR=$(dirname "$HERE")
W=${1:?usage: test_kmer_order.sh WORKDIR [REAL_FASTA ...]}
shift
PYTHON=${PYTHON:-python3}
mkdir -p "$W" || exit 1
cd "$W" || exit 1

# the oracle, cut out of SubFam.sh
awk '/^KMER_ORDER_PY=\$\(cat <<.PY.$/ {f=1; next} f && /^PY$/ {exit} f' "$VENDOR/SubFam.sh" > oracle.py
[ -s oracle.py ] || { echo "could not extract the oracle from SubFam.sh"; exit 1; }
cc -O2 -ffp-contract=off -fopenmp -o kmer_order "$VENDOR/kmer_order.c" -lm 2>/dev/null ||
    cc -O2 -ffp-contract=off -o kmer_order "$VENDOR/kmer_order.c" -lm || { echo "compile failed"; exit 1; }
"$PYTHON" "$HERE/gen_kmer_cases.py" cases "${SCALE:-1}" || exit 1

pass=0; fail=0
check() {   # file k canonical orient
    local f=$1 k=$2 c=$3 o=$4 b
    b=$(basename "$f" .fa)
    "$PYTHON" oracle.py "$f" "py.out" "$k" "$c" "$o" >py.log 2>&1; local rp=$?
    ./kmer_order "$f" "c.out" "$k" "$c" "$o" >c.log 2>&1; local rc=$?
    if [ $rp -ne 0 ] && [ $rc -ne 0 ]; then
        echo "both fail   $b k=$k c=$c o=$o (python: $(tail -1 py.log))"; pass=$((pass+1)); return
    fi
    if [ $rp -ne $rc ]; then
        echo "EXIT DIFFERS $b k=$k c=$c o=$o python=$rp c=$rc ($(tail -1 py.log) | $(tail -1 c.log))"; fail=$((fail+1)); return
    fi
    if cmp -s py.out c.out; then pass=$((pass+1)); echo "ok   $b k=$k c=$c o=$o ($(grep -c '^>' c.out) seqs)"
    else fail=$((fail+1)); echo "DIFF $b k=$k c=$c o=$o"; cp py.out "fail_${b}_$k$c$o.py.out"; cp c.out "fail_${b}_$k$c$o.c.out"; fi
}

for f in cases/*.fa; do
    case $(basename "$f") in big*) continue ;; esac
    for k in 3 6 8 12; do
        for co in "0 0" "1 0" "0 1" "1 1"; do check "$f" "$k" ${co% *} ${co#* }; done
    done
done
for f in cases/big*.fa "$@"; do
    [ -f "$f" ] || continue
    check "$f" 6 1 1
    check "$f" 8 0 0
done
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
