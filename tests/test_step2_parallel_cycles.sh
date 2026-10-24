#!/usr/bin/env bash
# test_step2_parallel_cycles.sh [WORKDIR] - the ten step-2 voting cycles can run in parallel (VOTE_PAR > 1, run_vote_cycles, 2026-10-24):
#  1. 600 copies (one part of >= 500): a run with VOTE_PAR=4 (4 cycles at a time, -T 1 each) and a default run (one cycle
#     after another) both give 10 votes per copy to the clear copies and the same assignment for them; 75 copies (a part under
#     500 copies, where ssearch36's statistics depend on -T): the cycles run one after another at -T 4, as before;
#  2. the cycles are independent stochastic estimates: the raw ssearch36 outputs of the 10 cycles (seen through a wrapper on
#     PATH) are not all the same (E-values of -z 11 come from clock-seeded shuffles); identical outputs would mean the cycles
#     shared a seed and the vote lost its resampling (-z 11 IS the quality filter);
#  3. a failing ssearch36 stops step 2 (exit != 0, "ssearch36 failed (cycle"); so does one that prints *** and exits 0.
# Needs ssearch36, seqkit.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD=${SD:-$(dirname "$HERE")}
W=${1:-${TMPDIR:-/tmp}/step2_parallel_toy}; rm -rf "$W"; mkdir -p "$W/bin" "$W/calls"; cd "$W" || exit 1
REAL=$(command -v ssearch36) || { echo "test_step2_parallel_cycles: no ssearch36"; exit 1; }
python3 - <<'PY'
import random
random.seed(7)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
mut = lambda s, p: "".join(random.choice("ACGT".replace(c, "")) if random.random() < p else c for c in s)
cons = {"TA": rnd(180), "TB": rnd(180), "TC": rnd(180)}
open("cons.fa", "w").write("".join(">%s\n%s\n" % kv for kv in cons.items()))
with open("ext.fa", "w") as f:
    i = 0
    for k, v in cons.items():
        for _ in range(200):
            i += 1; f.write(">ctg%d:1000-1180(+)\n%s\n" % (i, mut(v, 0.06)))
PY
cat > bin/ssearch36 <<EOF
#!/usr/bin/env bash
# test wrapper: MODE=fail exits 1, MODE=stars prints *** and exits 0, otherwise the real ssearch36 with its output kept
case "\${MODE:-}" in fail) exit 1 ;; stars) echo "*** cannot open library" >&2; exit 0 ;; esac
out=\$(mktemp "$W/calls/call_XXXXXX")
"$REAL" "\$@" | tee "\$out"
EOF
chmod +x bin/ssearch36
FAILS=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; FAILS=$((FAILS + 1)); fi; }
VOTE_PAR=4 PATH="$W/bin:$PATH" bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 4 par > par.log 2>&1; RC1=$?
mv calls calls_par; mkdir calls
PATH="$W/bin:$PATH" bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 4 seq > seq.log 2>&1; RC2=$?
chk "both runs exit 0" '[ $RC1 -eq 0 ] && [ $RC2 -eq 0 ]'
chk "parallel run logged 4 cycles at a time" 'grep -q "4 at a time, 1 ssearch36 threads each" par.log'
chk "default (VOTE_PAR unset): one cycle after another at -T 4, as before" 'grep -q "1 at a time, 4 ssearch36 threads each" seq.log'
chk "10 ssearch36 calls in each run (one part)" '[ $(ls calls_par | wc -l) -eq 10 ] && [ $(ls calls | wc -l) -eq 10 ]'
chk "all 600 copies 10/10 in both runs" '[ $(awk -F"\t" "\$3 == 10" par/all_votes.tsv | wc -l) -eq 600 ] && [ $(awk -F"\t" "\$3 == 10" seq/all_votes.tsv | wc -l) -eq 600 ]'
chk "same assignment in both runs" 'cmp -s <(cut -f1,2,5 par/assignment_full.tsv | sort) <(cut -f1,2,5 seq/assignment_full.tsv | sort)'
NDIST=$(for f in calls_par/*; do md5sum < "$f"; done | sort -u | wc -l)
chk "parallel cycles are independent estimates ($NDIST of 10 raw outputs distinct, want >= 9)" '[ $NDIST -ge 9 ]'
head -150 ext.fa > small.fa
VOTE_PAR=4 PATH="$W/bin:$PATH" bash "$SD/step2_asSINEment.sh" cons.fa small.fa 4 small > small.log 2>&1; RC5=$?
chk "75 copies (a part under 500): cycles one after another at -T 4" '[ $RC5 -eq 0 ] && grep -q "1 at a time, 4 ssearch36 threads each" small.log'
MODE=fail VOTE_PAR=4 PATH="$W/bin:$PATH" bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 4 fail > fail.log 2>&1; RC3=$?
chk "a failing ssearch36 stops step 2" '[ $RC3 -ne 0 ] && grep -q "ssearch36 failed (cycle 1" fail.log'
MODE=stars VOTE_PAR=4 PATH="$W/bin:$PATH" bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 4 stars > stars.log 2>&1; RC4=$?
chk "an ssearch36 that prints *** stops step 2" '[ $RC4 -ne 0 ] && grep -q "ssearch36 error in cycle 1" stars.log'
(( FAILS == 0 )) && echo "test_step2_parallel_cycles: passed" || { echo "test_step2_parallel_cycles: $FAILS failures"; exit 1; }
