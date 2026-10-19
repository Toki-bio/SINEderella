#!/usr/bin/env bash
# test_research_failures.sh [WORKDIR] - research scripts report a failed step instead of a normal result (GLM audit 2026-10-19):
#   28/29/91 peel_sweep.sh: a failed truth step / peel is a FAILED row, not a row of zeros
#   37/107   run_sim.sh: a failed sim_family.py marks the run FAILED and the exit status is not 0
#   30       rescore.sh: a failed score.py keeps the existing scores.txt
# Every Python call goes to a stub that fails (PYTHON). Needs only bash and coreutils. TOOLS overrides the tools directory.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; TOOLS=${TOOLS:-$(dirname "$HERE")/tools}
W=${1:-${TMPDIR:-/tmp}/research_failures_toy}; rm -rf "$W"; mkdir -p "$W/run/sf_n5"; cd "$W" || exit 1
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }
printf '#!/usr/bin/env bash\necho "stub python fails" >&2; exit 1\n' > fakepy; chmod +x fakepy
export PYTHON=$W/fakepy
printf ">x\nACGT\n" > run/sf_n5/copies.clw; printf "c1\tx_1\n" > run/sf_n5/copies.chunks.tsv; printf "c1\tA\n" > run/copies.labels
PEELPY=$W/fakepy GRID_G=0.8 GRID_M=5 bash "$TOOLS/coseg_compare/peel_sweep.sh" sweep.tsv 1 "$W/run" > ps.log 2>&1
chk "28/29/91 peel_sweep: failed step is a FAILED row" 'grep -q FAILED sweep.tsv'
printf "scen\t--x 1\t1\n" > sc.tsv
bash "$TOOLS/sim/run_sim.sh" "$W/sim" 1 1 "$W/sc.tsv" > rs.log 2>&1; RC=$?
chk "37/107 run_sim: failed simulation, exit $RC, FAILED marker" '[[ $RC -ne 0 && -f sim/scen_s1/FAILED ]]'
mkdir -p old; echo "kept" > old/scores.txt
bash "$TOOLS/coseg_compare/rescore.sh" "$W/old" > rsc.log 2>&1
chk "30 rescore: failed score.py keeps scores.txt" 'grep -q kept old/scores.txt'
[[ $fails -eq 0 ]] && echo "test_research_failures: PASS" || { echo "test_research_failures: $fails FAIL"; exit 1; }
