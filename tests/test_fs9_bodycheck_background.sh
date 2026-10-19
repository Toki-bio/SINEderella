#!/usr/bin/env bash
# test_fs9_bodycheck_background.sh [WORKDIR] - fs9_bodycheck.py draws its random background only from copies in no flank pair.
# Until 2026-10-19 it drew it from all copies, so in a family made mostly of duplicated copies the background was the
# duplication identity itself and the duplicated pairs fell back to "unique" (GLM audit finding 45).
# Toy: 20 copies with one identical body in 10 flank pairs (T1 twins) + 12 unpaired copies at 25 % divergence.
# Expected: all 20 paired copies stay twin1. Needs bedtools, ssearch36. BODYCHECK_PY overrides the script tested.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"
PY=${BODYCHECK_PY:-$SD/flankscan/fs9_bodycheck.py}
W=${1:-${TMPDIR:-/tmp}/fs9_bodycheck_bg_toy}; rm -rf "$W"; mkdir -p "$W/out"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(5)
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))
mut = lambda s, p: "".join(random.choice(B.replace(c, "")) if random.random() < p else c for c in s)
body = rnd(200)
g, bed, names = [], [], []
pos = 0
for i in range(32):
    s = body if i < 20 else mut(body, 0.25)
    pad = rnd(500)
    g.append(pad + s); st = pos + 500; en = st + 200; pos = en
    n = "ctg:%d-%d(+)" % (st, en); names.append(n); bed.append("ctg\t%d\t%d\t%s\t0\t+\n" % (st, en, n))
g.append(rnd(500))
open("g.fa", "w").write(">ctg\n" + "".join(g) + "\n")
open("copies.bed", "w").write("".join(bed))
with open("out/flank_groups.tsv", "w") as f:
    f.write("copy_a\tcopy_b\tside\tidentity\tlength\ttier\tclass\n")
    for k in range(10):
        f.write("%s\t%s\t5\t100.0\t100\tT1\ttwin\n" % (names[2 * k], names[2 * k + 1]))
with open("out/copy_status.tsv", "w") as f:
    for i, n in enumerate(names):
        f.write("%s\t%s\n" % (n, "twin1\t1\tG%d" % (i // 2) if i < 20 else "unique\t0\t-"))
PY
python3 -I "$PY" g.fa copies.bed out > bc.log 2>&1 || { cat bc.log; echo "test_fs9_bodycheck_background: fs9_bodycheck.py failed"; exit 1; }
T=$(awk -F'\t' '$2 == "twin1"' out/copy_status.tsv | wc -l)
if [[ $T -eq 20 ]]; then echo "test_fs9_bodycheck_background: PASS (20 of 20 duplicated copies stay twin1)"
else cat bc.log; echo "test_fs9_bodycheck_background: FAIL ($T of 20 duplicated copies stay twin1)"; exit 1; fi
