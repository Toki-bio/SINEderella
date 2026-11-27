#!/usr/bin/env bash
# step4_plots.sh runs subfamilies in parallel (2026-11-27). Checks on a toy run (seconds):
#   1. STEP4_JOBS=1 and STEP4_JOBS=3 give the same tables and the same PNG bytes, and the log blocks in subfamily order
#   2. a subfamily whose command fails (stub seqkit that exits 1 for one subfamily) stops the step non-zero,
#      with an ERROR line naming it
# Needs ssearch36, mafft, seqkit and python3 with matplotlib (conda env sinederella).
# Usage: bash tests/test_step4_parallel.sh [SINEderella dir]
set -uo pipefail
[ -f ~/miniforge3/etc/profile.d/conda.sh ] && { source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella 2>/dev/null || true; }
SD="$(cd "${1:-$(dirname "${BASH_SOURCE[0]}")/..}" && pwd)"
T="$(mktemp -d "${TMPDIR:-/tmp}/test_step4.XXXXXX")"
trap 'rm -rf "$T"' EXIT
FAILS=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; FAILS=$((FAILS + 1)); fi; }

# toy: three consensuses, copies = consensus with a few substitutions (fixed, no randomness)
make_run() {
  local R="$1"
  mkdir -p "$R/step2/step2_output/subfamilies"
  awk 'BEGIN{
    srand(7); b="ACGT"
    for (k = 1; k <= 3; k++) {
      name = "SF" k; s = ""
      for (i = 1; i <= 150; i++) s = s substr(b, int(rand()*4)+1, 1)
      print ">" name > (R "/consensuses.clean.fa")
      print s > (R "/consensuses.clean.fa")
      f = R "/step2/step2_output/subfamilies/" name ".fasta"
      for (c = 1; c <= 12; c++) {
        t = s
        for (m = 1; m <= c; m++) { p = int(rand()*150)+1; t = substr(t, 1, p-1) substr(b, int(rand()*4)+1, 1) substr(t, p+1) }
        printf ">%s_c%d|chr1:%d-%d|+\n%s\n", name, c, c*1000, c*1000+150, t > f
      }
      close(f)
    }
  }' R="$R"
}

run_step4() {  # RUN_ROOT JOBS MAX_SEQS [PATH prefix]
  ( [ -n "${4:-}" ] && export PATH="$4:$PATH"
    STEP4_JOBS="$2" bash "$SD/step4_plots.sh" "$1" 3 "$3" > "$1/step4.out" 2> "$1/step4.err" )
}

for j in 1 3; do make_run "$T/j$j"; run_step4 "$T/j$j" $j 8; echo "jobs=$j exit=$?" > "$T/j$j/rc"; done
chk "jobs=1 exit 0"                     'grep -q "exit=0" "$T/j1/rc"'
chk "jobs=3 exit 0"                     'grep -q "exit=0" "$T/j3/rc"'
chk "jobs=3: 3 pctid tables and 12 plots" '[ $(ls "$T/j3/step2/step2_output/plots" | grep -c "_pctid.tsv$") -eq 3 ] && [ $(ls "$T/j3/step2/step2_output/plots" | grep -c "\.png$\|\.pdf$") -eq 12 ]'
chk "jobs=1 and jobs=3: same tables and PNG bytes" \
    '[ "$(cd "$T/j1/step2/step2_output/plots" && md5sum *_pctid.tsv *.png)" = "$(cd "$T/j3/step2/step2_output/plots" && md5sum *_pctid.tsv *.png)" ]'
order() { awk '$3 ~ /^\[[0-9]+\/3\]$/ {printf "%s ", $3}' "$1"; }
chk "jobs=3: log blocks in subfamily order" '[ "$(order "$T/j3/step4.err")" = "[1/3] [2/3] [3/3] " ]'

# failing child: seqkit stub fails for SF2 only (MAX_SEQS=8 < 12 copies, so every subfamily is subsampled)
mkdir -p "$T/stub"
REAL_SEQKIT="$(command -v seqkit)"
cat > "$T/stub/seqkit" <<EOF
#!/usr/bin/env bash
case "\$*" in *SF2.fasta*) echo "stub seqkit: failing on purpose" >&2; exit 1;; esac
exec "$REAL_SEQKIT" "\$@"
EOF
chmod +x "$T/stub/seqkit"
make_run "$T/fail"; run_step4 "$T/fail" 3 8 "$T/stub"; rc=$?
chk "failing subfamily: step exits non-zero (got $rc)" '[ $rc -ne 0 ]'
chk "failing subfamily: ERROR line names SF2"         'grep -q "ERROR: subfamily SF2 failed" "$T/fail/step4.err"'
chk "failing subfamily: no Done line"                 '! grep -q "Done\." "$T/fail/step4.err"'

echo "$FAILS failed"
exit $(( FAILS > 0 ))
