#!/usr/bin/env bash
# test_sine_consensus_label.sh [WORKDIR] - sine_consensus.sh labels the last trace record FINAL_not_converged when it did not
# converge (until 2026-10-19 always FINAL_converged; GLM audit finding 65). max_iters 2 cannot converge (needs >= 5 rounds).
# Needs mafft, bc. SC overrides the script tested.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SC=${SC:-$(dirname "$HERE")/sine_consensus.sh}
W=${1:-${TMPDIR:-/tmp}/sine_consensus_label_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
python3 -c '
import random
random.seed(4); e = "".join(random.choice("ACGT") for _ in range(150))
mut = lambda s: "".join(random.choice("ACGT") if random.random() < 0.1 else c for c in s)
open("copies.fasta", "w").write("".join(">c%d\n%s\n" % (i, mut(e)) for i in range(30)))'
bash "$SC" copies.fasta 10 2 0.01 > sc.log 2>&1
L=$(grep '^>FINAL' convergence_trace.fasta 2> /dev/null)
if [[ "$L" == ">FINAL_not_converged" ]]; then echo "test_sine_consensus_label: PASS"
else echo "test_sine_consensus_label: FAIL (last trace record '$L')"; exit 1; fi
