#!/usr/bin/env bash
# test_reproducible_run.sh GENOME.fa BANK.fa [WORKDIR] [SPECIES_CODE]
# The whole workflow twice on the same input (with publish: verdicts and plates), each run in its own directory, then
# tests/compare_runs.py: loci, hits, assignments, verdicts, published consensuses and every plate must be identical.
# The two runs use different numbers of threads (THREADS_A, THREADS_B; default 8 and 12): the result must not depend on it.
# Minutes on a small genome (11 Mb, 21 consensuses: about 15 min per run); run it on a server, not in pytest.
# Exit 0 = identical. docs/REPRODUCIBILITY.md.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD="$(dirname "$HERE")"
G=$(readlink -f "${1:?genome fasta}"); K=$(readlink -f "${2:?consensus bank fasta}")
W=${3:-${TMPDIR:-/tmp}/sinederella_reproducible_run}; SP=${4:-rep}
mkdir -p "$W"; cd "$W" || exit 1
pids=()
for x in A B; do
  rm -rf "$W/$x"; mkdir -p "$W/$x" "$W/tmp_$x"
  cp "$G" "$W/$x/genome.fa"; cp "$K" "$W/$x/bank.fa"
  if [ "$x" = A ]; then T=${THREADS_A:-8}; else T=${THREADS_B:-12}; fi
  ( cd "$W/$x" && THREADS=$T PUBLISH=1 SPECIES_CODE="$SP" TMPDIR="$W/tmp_$x" "$SD/SINEderella" genome.fa bank.fa > run.log 2>&1 ) &
  pids+=($!)
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ "$rc" -eq 0 ] || { tail -n 5 "$W"/A/run.log "$W"/B/run.log; echo "test_reproducible_run: a run failed (see run.log in $W/A and $W/B)"; exit 1; }
RA=$(ls -d "$W"/A/run_* | head -n 1); RB=$(ls -d "$W"/B/run_* | head -n 1)
python3 "$HERE/compare_runs.py" "$RA" "$RB" --loci "$W/loci_A_B.tsv" | tee "$W/compare_A_B.txt"
[ "${PIPESTATUS[0]}" -eq 0 ] && echo "test_reproducible_run: passed" || { echo "test_reproducible_run: FAILED (see $W/compare_A_B.txt, $W/loci_A_B.tsv)"; exit 1; }
