#!/usr/bin/env bash
# test_step7_small_families.sh [WORKDIR] - step7_boundary_refine.sh tests the flanks of families from 7 copies.
# Until 2026-10-21 families with fewer than 20 copies got "insufficient_data" and all their copies were grey, only for their
# size (rsi MEG-TR: 16 copies grey in one run, 20 copies tested in the next).
# Planted genome, 3 Mb random DNA, a 200 bp element per family:
#   A 10 copies, independent insertions                        -> both sides confirmed
#   B  8 copies, one pair duplicated with its flanks (twins)   -> both sides confirmed (1 of 28 pairs, under background + 5 %)
#   C  7 copies, each inside one identical 3 kb duplicated block -> both sides undetermined (flanks shared past 1000 bp)
#   D  4 copies, independent                                   -> insufficient_data (below 7, one twin pair could decide)
#   E  7 copies, independent (the smallest family tested)      -> both sides confirmed
# Needs samtools, python3. STEP7 overrides the script tested.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"
S7=${STEP7:-$SD/step7_boundary_refine.sh}
W=${1:-${TMPDIR:-/tmp}/step7_small_toy}; rm -rf "$W"; mkdir -p "$W/run/step2/step2_output"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(3)
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))
mut = lambda s, p: "".join(random.choice(B.replace(c, "")) if random.random() < p else c for c in s)
g = list(rnd(3000000))
heads, slot = [], iter(range(20000, 2980000, 20000))
def put(seq):
    p = next(slot)
    g[p:p + len(seq)] = list(seq)
    return p
def copy(fam, el):
    p = put(mut(el, 0.05))
    heads.append("chr1:%d-%d(+)|%s|100" % (p, p + len(el), fam))
for fam, n in (("A", 10), ("D", 4), ("E", 7)):
    el = rnd(200)
    for _ in range(n):
        copy(fam, el)
el = rnd(200)
for _ in range(6):
    copy("B", el)
twin = rnd(1400) + mut(el, 0.05) + rnd(1400)                 # one copy with 1.4 kb flanks, placed twice
for _ in range(2):
    p = put(twin)
    heads.append("chr1:%d-%d(+)|B|100" % (p + 1400, p + 1600))
el = rnd(200)
block = rnd(1400) + el + rnd(1400)                           # every copy of C inside the same 3 kb block
for _ in range(7):
    p = put(mut(block, 0.01))
    heads.append("chr1:%d-%d(+)|C|100" % (p + 1400, p + 1600))
open("run/genome.clean.fa", "w").write(">chr1\n" + "".join(g) + "\n")
with open("run/step2/step2_output/assigned.fasta", "w") as f:
    for h in heads:
        f.write(">%s\nACGT\n" % h)
PY
samtools faidx run/genome.clean.fa
bash "$S7" "$W/run" > s7.log 2>&1 || { tail -20 s7.log; echo "test_step7_small_families: FAIL (step7 failed)"; exit 1; }
T=run/step2/step2_output/boundary_refinement.tsv
st(){ awk -F'\t' -v f="$1" '$1 == f && $2 == "general" { print $5 }' "$T" | sort -u | tr '\n' ' ' | sed 's/ $//'; }
fail=0
for e in "A confirmed" "B confirmed" "C undetermined" "D insufficient_data" "E confirmed"; do
    set -- $e
    got=$(st "$1")
    if [[ "$got" == "$2" ]]; then echo "  $1: $got"; else echo "  $1: expected $2, got '$got'"; fail=1; fi
done
if (( fail )); then cat "$T"; echo "test_step7_small_families: FAIL"; exit 1; fi
echo "test_step7_small_families: PASS"
