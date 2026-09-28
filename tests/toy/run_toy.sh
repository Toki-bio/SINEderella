#!/usr/bin/env bash
# Toy end-to-end test of step8a and align_for_publish (seconds). Run BEFORE any real publish on changed code.
# Usage: [SD=<SINEderella dir>] [DISCD=<SINE_discriminator/site>] bash tests/toy/run_toy.sh
# Unit test of tandem-array selection at realistic density: python3 tests/toy/test_array_order.py tools/array_order.py
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$HOME/SINEderella:$HOME/SINEderella/bin:$HOME/SubFam:$PATH" TMPDIR=$HOME/tmp
SD=${SD:-$HOME/SINEderella}; DISCD=${DISCD:-$HOME/SINE_discriminator/site}
T=~/tmp/toy8a; rm -rf $T; mkdir -p $T
python3 "$HERE/make_toy.py" $T
SECONDS=0
DISC=$DISCD bash $SD/step8a_extract_alignments.sh $T toy > $T/step8a.log 2>&1; echo "step8a exit=$? in ${SECONDS}s"
grep -E "re-extract|WARN|ERROR|rror" $T/step8a.log
A=$T/results/alignments
echo "--- plates"; for f in $A/*.aln.fa; do printf "%-24s rows %3d soft %2d array %2d | row1 %s row2 %s row3 %s\n" $(basename $f) $(grep -c "^>" $f) $(grep -c "\[soft\]" $f) $(grep -c "\[array\]" $f) $(grep "^>" $f | head -3 | cut -c2-40 | tr -d ' '); done
echo "--- checks"
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fi; }
chk "TOYS top100: 5 array units marked"                 '[ $(grep -c "\[array\]" $A/toy_TOYS_top100.aln.fa) -eq 5 ]'
# selection of one unit per array into top100 is covered by test_array_order.py (plate row order is MAFFT's)
chk "TOYS: 4 soft rows"                                 '[ $(grep -c "\[soft\]" $A/toy_TOYS_top100.aln.fa) -eq 4 ]'
chk "TOYB/TOYC/TOYL: no array marks"                    '[ $(cat $A/toy_TOYB_*.aln.fa $A/toy_TOYC_*.aln.fa $A/toy_TOYL_*.aln.fa | grep -c "\[array\]") -eq 0 ]'
chk "TOYC soft-only plates exist, 3 soft rows"          '[ $(grep -c "\[soft\]" $A/toy_TOYC_top100.aln.fa) -eq 3 ]'
chk "every plate: row2 = original, row1 = _extended"    '! for f in $A/*.aln.fa; do s=$(basename $f .aln.fa | cut -d_ -f2); [ "$(grep "^>" $f | sed -n 1p)" = ">${s}_extended" ] && [ "$(grep "^>" $f | sed -n 2p)" = ">$s" ] || echo bad; done | grep -q bad'
chk "TOYA not re-extracted (array units skipped)"   '! grep -q "toy_TOYA_.*re-extracting" $T/step8a.log'
chk "TOYL re-extracted (continuation loop ran)"         'grep -q "toy_TOYL_top100.*re-extracting" $T/step8a.log'
for f in $A/toy_TOYL_top100.aln.fa $A/toy_TOYB_top100.aln.fa; do echo "   $(basename $f): $(python3 $DISCD/continuation.py $f | tr '\n' ' ')"; done

echo "=== align_for_publish on a fresh toy (as repub_soft: no step7 / border loop)"
P=~/tmp/toypub; rm -rf $P; mkdir -p $P; python3 "$HERE/make_toy.py" $P > /dev/null
SECONDS=0
THREADS=4 SKIP_STEP7=1 SKIP_BORDER_LOOP=1 DISC=$DISCD SINEDERELLA_BIN=$SD bash $SD/publish/align_for_publish.sh $P toy > $P/afp.log 2>&1; echo "align_for_publish exit=$? in ${SECONDS}s"
grep -iE "error|traceback|warn|fail" $P/afp.log | head -15
A=$P/results/alignments
ls $A | head -30
echo "--- proposals.tsv"; column -t -s$'\t' $A/proposals.tsv 2>/dev/null | cut -c1-200
echo "--- continuation.tsv"; cat $A/continuation.tsv 2>/dev/null
for f in $A/toy_TOYL_top100.aln.fa $A/toy_TOYS_top100.aln.fa; do echo "--- $f"; grep "^>" $f | head -4; python3 - "$f" <<'PY'
import sys; sys.path.insert(0, "/home/toki/SINE_discriminator/site")
from fix_alignments import read_fa
n, s = read_fa(sys.argv[1])
u = [j for j, c in enumerate(s[0]) if c.isupper()]
print("  row1 upper span %d-%d, letters %d (lower %d)" % (u[0], u[-1], sum(c not in "-." for c in s[0]), sum(c.islower() for c in s[0])))
for x, y in list(zip(n, s))[:5]: print("  %-30s %s" % (x[:30], y[max(0, u[-1] - 15):u[-1] + 60]))
PY
done
echo "=== full publish_run.sh on the toy (step6 may stop on missing step3 files - reported, not hidden)"
SECONDS=0
THREADS=4 SKIP_STEP7=1 SKIP_BORDER_LOOP=1 SKIP_STEP4=1 SKIP_ALIGN=1 DISC=$DISCD SINEDERELLA_BIN=$SD bash $SD/publish_run.sh $P toy > $P/pr.log 2>&1; echo "publish_run exit=$? in ${SECONDS}s"; tail -8 $P/pr.log
