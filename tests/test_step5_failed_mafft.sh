#!/usr/bin/env bash
# test_step5_failed_mafft.sh [WORKDIR] - the legacy step 5 scripts (run by hand; GLM audit 2026-10-19, findings 26, 62, 63, 79):
#   a failed mafft leaves no .aln.fa behind and is not counted (step5_direct.sh, step5_align_subfamilies.sh, run_step5_wrapper.sh);
#   run_step5_wrapper.sh works on a run directory whose name does not start with run_.
# mafft is a stub on PATH (MAFFT_FAIL=1: writes a partial record and exits 1). Needs only bash and coreutils.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; SD=${SD:-$(dirname "$HERE")}
W=${1:-${TMPDIR:-/tmp}/step5_failed_mafft_toy}; rm -rf "$W"; mkdir -p "$W/bin" "$W/sp/myrun/step2/step2_output/subfamilies"; cd "$W" || exit 1
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }
cat > bin/mafft <<'EOF'
#!/usr/bin/env bash
[[ -n "${MAFFT_FAIL:-}" ]] && { echo ">partial"; exit 1; }
for a in "$@"; do f=$a; done; cat "$f"
EOF
chmod +x bin/mafft; export PATH="$W/bin:$PATH"
R=$W/sp/myrun
printf ">cons\nACGTACGTAC\n" > "$R/consensuses.clean.fa"
for i in 1 2 3; do printf ">c%d\nACGTACGTAC\n" $i; done > "$R/step2/step2_output/subfamilies/A.fasta"

MAFFT_FAIL=1 bash "$SD/step5_direct.sh" "$R" > d1.log 2>&1
chk "79 step5_direct: failed mafft leaves no A.aln.fa" '[[ ! -e $R/alignments/A.aln.fa ]]'
chk "79 step5_direct: failed mafft is not counted" 'grep -q "Done! 0 alignments" d1.log'
rm -rf "$R/alignments"

cp -r "$R" "$W/sp/run_x"
MAFFT_FAIL=1 bash "$SD/run_step5_wrapper.sh" "$W/sp/run_x" "$SD/step5_align_subfamilies.sh" > w1.log 2>&1
chk "26/62 wrapper + step5: failed mafft, no alignment shipped" '[[ -z "$(find $W/sp/run_x/alignments -name "*.aln.fa" 2> /dev/null)" ]]'
bash "$SD/run_step5_wrapper.sh" "$R" "$SD/step5_align_subfamilies.sh" > w2.log 2>&1
chk "63 wrapper: run directory not named run_* gets its alignment" '[[ -s $R/alignments/sp/A.aln.fa ]]'

[[ $fails -eq 0 ]] && echo "test_step5_failed_mafft: PASS" || { echo "test_step5_failed_mafft: $fails FAIL"; exit 1; }
