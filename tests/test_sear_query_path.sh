#!/usr/bin/env bash
# test_sear_query_path.sh [WORKDIR] - sear names its outputs after the query FILE name, also when the query is given with a
# directory (extract_alignments.sh tier D passes /tmp/.../consensus.fa). Until 2026-10-19 the directory went into the output
# name, the hit file could not be written and tier D evidence was always "No hits found" (GLM audit findings 4, 217).
# Toy: 6 copies of a 200 bp element in a 30 kb genome. Needs what sear needs (ssearch36, bedtools, samtools, seqkit).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SEAR=${SEAR:-$(dirname "$HERE")/sear}
W=${1:-${TMPDIR:-/tmp}/sear_query_path_toy}; rm -rf "$W"; mkdir -p "$W/q" "$W/run"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(7)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
mut = lambda s, p: "".join(random.choice("ACGT".replace(c, "")) if random.random() < p else c for c in s)
e = rnd(200)
open("q/consensus.fa", "w").write(">cons\n%s\n" % e)
open("run/gen.fa", "w").write(">c1\n%s\n" % "".join(rnd(4000) + mut(e, 0.05) for _ in range(6)) + rnd(4000) + "\n")
PY
cd run || exit 1
bash "$SEAR" "$W/q/consensus.fa" gen.fa 0.5 30 20 > sear.log 2>&1
N=$( [[ -f gen-consensus.bed ]] && wc -l < gen-consensus.bed || echo none)
if [[ "$N" != none && "$N" -ge 6 ]]; then echo "test_sear_query_path: PASS (gen-consensus.bed, $N hits)"
else tail -5 sear.log; echo "test_sear_query_path: FAIL (gen-consensus.bed: $N)"; exit 1; fi
