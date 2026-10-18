#!/usr/bin/env bash
# test_step2_assigned.sh [WORKDIR] - step 2 writes ONLY copies with status "assigned" into assigned.fasta and the per-subfamily
# FASTA files (until 2026-10-19 rejected_low_bitscore copies went in too: rsi 1,731 of 65,500; GLM audit finding).
# Toy: two unrelated 150 bp consensuses, 30 copies each at 8 % divergence and 6 each at 40 % (votes unanimous, bitscore
# under the threshold -> rejected_low_bitscore). Needs ssearch36, seqkit (as step 2 does).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"
W=${1:-${TMPDIR:-/tmp}/step2_assigned_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(11)
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))
mut = lambda s, p: "".join(random.choice(B.replace(c, "")) if random.random() < p else c for c in s)
cons = {"TA": rnd(150), "TB": rnd(150)}
with open("cons.fa", "w") as f:
    for k, v in cons.items():
        f.write(">%s\n%s\n" % (k, v))
with open("ext.fa", "w") as f:
    i = 0
    for k, v in cons.items():
        for p, n in ((0.08, 30), (0.40, 6)):
            for _ in range(n):
                i += 1
                f.write(">ctg%d:%d-%d(+)\n%s\n" % (i, 1000, 1150, mut(v, p)))
PY
bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 2 out > step2.log 2>&1 || { tail -20 step2.log; echo "test_step2_assigned: step 2 failed"; exit 1; }
O=$(ls -d out* | head -1)
A=$(awk -F'\t' 'NR > 1 && $5 == "assigned"' "$O/assignment_full.tsv" | wc -l)
REJ=$(awk -F'\t' 'NR > 1 && $5 == "rejected_low_bitscore"' "$O/assignment_full.tsv" | wc -l)
F=$(grep -c '^>' "$O/assigned.fasta")
S=$(cat "$O"/subfamilies/*.fasta 2> /dev/null | grep -c '^>')
awk -F'\t' 'NR > 1 && $5 == "rejected_low_bitscore" { print $1 }' "$O/assignment_full.tsv" | sort > rej.ids
grep '^>' "$O/assigned.fasta" | sed 's/^>//; s/|.*//' | sort > af.ids
LEAK=$(comm -12 rej.ids af.ids | wc -l)
echo "assigned=$A rejected=$REJ assigned.fasta=$F subfamily files=$S rejected in assigned.fasta=$LEAK"
FAILS=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; FAILS=$((FAILS + 1)); fi; }
chk "the toy produces rejected_low_bitscore copies"            '[ "$REJ" -gt 0 ]'
chk "assigned.fasta holds exactly the assigned copies"          '[ "$F" -eq "$A" ]'
chk "no rejected copy in assigned.fasta"                        '[ "$LEAK" -eq 0 ]'
chk "per-subfamily files hold exactly the assigned copies"      '[ "$S" -eq "$A" ]'
[ "$FAILS" -eq 0 ] && echo "test_step2_assigned: passed" || { echo "test_step2_assigned: FAILED"; exit 1; }
