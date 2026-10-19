#!/usr/bin/env bash
# test_flankscan_fixes.sh [WORKDIR] - flankscan stages 1, 6b, 7 and 8 (fs_all.sh, run by hand; GLM audit 2026-10-19):
#   44  fs_endcheck: a 5' flank cut short by a contig end, made entirely of a bank unit up to the junction, is an open end
#   58  fs7_hierarchy.sh stops when ends.tsv (stage 6b) is missing instead of passing every end as closed
#   116 fs1_extract.sh rebuilds a genome index older than the genome
#   118 fs8_singletons.sh stops when the family is not exactly one record of cons.masked.fa
# Needs bedtools, seqkit, samtools, ssearch36, gawk. FS overrides the flankscan directory tested.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; FS=${FS:-$(dirname "$HERE")/flankscan}
W=${1:-${TMPDIR:-/tmp}/flankscan_fixes_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }
python3 - <<'PY'
import random
random.seed(3)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
unit, elem = rnd(100), rnd(300)
open("cons.masked.fa", "w").write(">r1\n%s\n>r1b\n%s\n" % (unit, rnd(100)))
# two contigs: 100 bp of unit r1 right up to the element (the 5' flank is cut to 100 bp by the contig start)
open("genome.clean.fa", "w").write("".join(">c%d\n%s\n" % (i, unit + elem + rnd(600)) for i in (1, 2)))
open("cand.bed", "w").write("c1\t100\t400\tx1\t0\t+\nc2\t100\t400\tx2\t0\t+\n")
PY
samtools faidx genome.clean.fa

# 44
mkdir -p e/cand; cp cons.masked.fa e/; cp cand.bed e/cand/E1.bed
OPEN=$(cd e && RUN=$W T=2 bash -c "source '$FS/fs_lib.sh'; fs_endcheck E1" 2> e.log | cut -f8)
chk "44 fs_endcheck: unit up to a contig-cut 5' junction is open on 5' (got '$OPEN')" '[[ "$OPEN" == *5* ]]'

# 58
mkdir -p h; for f in cons.tsv peaks.tsv reassign.tsv candidates.tsv; do printf "h\n" > h/$f; done; : > h/candidates.fa
bash "$FS/fs7_hierarchy.sh" h > h.log 2>&1; RC=$?
chk "58 fs7 without ends.tsv stops (exit $RC)" '[[ $RC -ne 0 ]] && grep -q ends.tsv h.log'

# 116: index of an older, shorter genome; the copy is 600 bp from the contig end in the new one
mkdir -p r/results; printf ">c1\n%s\n" "$(head -c 500 /dev/zero | tr '\0' A)" > r/genome.clean.fa; samtools faidx r/genome.clean.fa
sleep 1; cp genome.clean.fa r/genome.clean.fa
printf "locus\tfam\tx\ty\tstatus\nc1:100-400(+)\tF\t.\t.\tassigned\n" > r/results/assignment_full.tsv
bash "$FS/fs1_extract.sh" r r/fs 500 > r.log 2>&1
C3=$(gawk -F'\t' 'NR == 2 { print $12 }' r/fs/loci.tsv 2> /dev/null)
chk "116 fs1 with a stale index: clamp3 = 0 (got '$C3')" '[[ "$C3" == 0 ]]'

# 118: family absent from cons.masked.fa
mkdir -p s; cp cons.masked.fa s/
REF=$(cd s && seqkit grep -p "r1x" cons.masked.fa | grep -c '^>')
R1=$(cd s && seqkit grep -p "r1" cons.masked.fa | grep -c '^>')      # exact ID match: r1 does not pick r1b
grep -q 'not exactly one record' "$FS/fs8_singletons.sh"; HAS=$?
chk "118 fs8 checks the reference record (absent family matched $REF records, r1 matched $R1)" '[[ $HAS -eq 0 && $REF -eq 0 && $R1 -eq 1 ]]'

[[ $fails -eq 0 ]] && echo "test_flankscan_fixes: PASS" || { echo "test_flankscan_fixes: $fails FAIL"; exit 1; }
