#!/usr/bin/env bash
# test_step2_dupnames.sh [WORKDIR] - step 2 stops on two copies with one name: its votes are kept per name, so the two copies
# would share one verdict and one subfamily (until 2026-10-19 it went on; GLM audit finding 113). Needs ssearch36, seqkit.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD=${SD:-$(dirname "$HERE")}
W=${1:-${TMPDIR:-/tmp}/step2_dupnames_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(2)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
mut = lambda s, p: "".join(random.choice("ACGT".replace(c, "")) if random.random() < p else c for c in s)
cons = {"TA": rnd(150), "TB": rnd(150)}
open("cons.fa", "w").write("".join(">%s\n%s\n" % kv for kv in cons.items()))
with open("ext.fa", "w") as f:
    for i in range(20):
        f.write(">ctg%d:1000-1150(+)\n%s\n" % (i, mut(cons["TA"], 0.08)))
    f.write(">ctg0:1000-1150(+)\n%s\n" % mut(cons["TB"], 0.08))      # same name, a TB copy
PY
bash "$SD/step2_asSINEment.sh" cons.fa ext.fa 2 out > step2.log 2>&1; RC=$?
if [[ $RC -ne 0 ]] && grep -q "duplicate copy names" step2.log; then echo "test_step2_dupnames: PASS (step 2 stops on ctg0:1000-1150(+))"
else echo "test_step2_dupnames: FAIL (exit $RC, no duplicate-name stop)"; exit 1; fi
