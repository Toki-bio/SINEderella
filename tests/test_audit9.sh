#!/usr/bin/env bash
# test_audit9.sh [WORKDIR] - the low findings of docs/AUDIT_2026-11-10.md (flankscan), one block each:
#   a  fs5_build.sh: the distinct elements of a peak, found through 1 kb bins, are the rows a comparison of every pair gives
#   b  fs8_singletons.sh: a group whose alignment has no conserved window does not stop stage 8
#   f  fs8_singletons.sh: the reference row keeps a dust-masked stretch inside the consensus (only the A tail is cut)
#   c  fs9_dupregion.sh: a failed minimap2 stops it
#   d  fs3_partners.sh: a missing dustmasker is named
#   e  fs9_twins.sh, fs9_dupregion.sh, fs9_run.sh: a genome index older than the genome is made again
#   i  fs9_run.sh: a "(+,-)" locus is among the copies of stage 9
#   g  fs6_reassign.sh: the no-candidate header has the 8 columns of the regular one
#   h  fs6c_chain.sh: a second run on a directory that already holds chains stops and leaves chains.tsv alone
# Needs gawk, bedtools, samtools, seqkit, ssearch36, mafft, python3 + numpy. FS overrides the flankscan directory tested.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; FS=${FS:-$(dirname "$HERE")/flankscan}
W=${1:-${TMPDIR:-/tmp}/audit9_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }

# ---------- a: the element step of fs5_build.sh itself (cut out of the script) against every-pair comparison in python
mkdir -p a/cand
python3 - <<'PY'
import random
random.seed(9)
rows = []
for i in range(4000):                      # 2 peaks, 3 contigs, spans 150-2500 bp, many overlaps, some end-to-end
    c = "c%d" % random.randint(1, 3); s = random.randint(0, 400) * 500 + random.choice([0, 0, 137]); e = s + random.choice([150, 300, 500, 1000, 2500])
    rows.append((random.choice(["P1", "P2"]), random.randint(100, 900) + random.random(), c, s, e, "w%d_3" % i, random.choice("+-")))
rows.sort(key=lambda r: (r[0], -r[1]))
with open("a/cand/spans.tsv", "w") as f:
    for r in rows: f.write("%s\t%.3f\t%s\t%d\t%d\t%s\t%s\n" % r)
kept = []
for r in rows:
    if r[0] != "P1": continue
    if any(k[2] == r[2] and r[3] < k[4] and r[4] > k[3] for k in kept): continue
    kept.append(r)
with open("a/expect.bed", "w") as f:
    for r in kept[:60]: f.write("%s\t%d\t%d\t%s\t%.3f\t%s\n" % (r[2], r[3], r[4], r[5], r[1], r[6]))
open("a/expect.n", "w").write("%d\n" % len(kept))
PY
( cd a && P=P1 NBEST=60 NAME=x && source <(sed -n '/# overlapping spans = one element/,/cand\/spans.tsv > "cand\/\$NAME.bed"/p' "$FS/fs5_build.sh") )
chk "a fs5: distinct elements of a peak = every-pair comparison ($(cat a/cand/n.tmp 2> /dev/null) elements, expected $(cat a/expect.n))" \
    "cmp -s a/cand/x.bed a/expect.bed && cmp -s a/cand/n.tmp a/expect.n"

# ---------- b, f: stage 8 on a family whose 30 "single" copies share nothing (no conserved window)
mkdir -p s8
python3 - <<'PY'
import random
random.seed(4)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
cons = rnd(80) + "AT" * 10 + rnd(80)                     # 180 bp, a 20 bp dust-masked stretch inside, then the A tail
open("s8/lib.in.fa", "w").write(">F\n%s\n" % (cons + "A" * 12))
open("s8/cons.masked.fa", "w").write(">F\n%s\n" % (cons[:80] + "N" * 20 + cons[100:] + "N" * 12))
open("s8/cons.tsv", "w").write("F\t192\t181\t20\n")
cp = open("s8/copies.tsv", "w"); cp.write("wid\tlocus\tfamily\tclass\tside5\tside3\n")
un = open("s8/units.tsv", "w"); un.write("wid\tfamily\tk\trole\tunit\tstrand\tcons_s\tcons_e\ttag\twin_s\twin_e\tbits\tpid\n")
jn = open("s8/junctions.tsv", "w"); jn.write("\t".join(["h"] * 26) + "\n")
wf = open("s8/windows.fa", "w")
for i in range(1, 31):
    w = "w%d" % i
    wf.write(">%s\n%s\n" % (w, rnd(700)))                # unrelated sequence: unit 301-480 of a 700 bp window
    cp.write("%s\tc1:%d-%d(+)\tF\tsingle\t-\t-\n" % (w, i * 5000, i * 5000 + 180))
    un.write("%s\tF\t1\tmain\tF\t+\t1\t180\the\t301\t480\t200\t90\n" % w)
    jn.write("\t".join([w, "F", "3"] + ["."] * 20 + ["0", ".", "."]) + "\n")
PY
bash "$FS/fs8_singletons.sh" s8 F "" 2 > s8.log 2>&1; RC8=$?
chk "b fs8: a group with no conserved window does not stop stage 8 (exit $RC8)" \
    "[ $RC8 -eq 0 ] && [ -s s8/singletons/F/summary.tsv ] && gawk -F'\t' '\$1 == \"S\" && \$4 == \"-\" { f = 1 } END { exit !f }' s8/singletons/F/ends.tsv && grep -q 'no conserved window' s8.log"
RL=$(seqkit head -n 1 s8/singletons/F/S.aln.fa 2> /dev/null | seqkit seq -s -g -w 0 2> /dev/null | tr -d '\n' | wc -c)
chk "f fs8: reference row = the consensus up to its A tail, masked stretch kept ($RL bp, expected 180)" "[ $RL -eq 180 ]"

# ---------- toy genome for c, e, i
python3 "$FS/tests/make_toy_flank.py" "$W/toy9" 1 > /dev/null
[[ -s toy9/genome.fa.fai ]] || samtools faidx toy9/genome.fa
GL=$(gawk '{ s += $2 } END { print s }' toy9/genome.fa.fai)
stale() {   # stale DIR/NAME: the genome with the index of a 500 bp genome made before it
    mkdir -p "$1"; printf ">old\n%s\n" "$(head -c 500 /dev/zero | tr '\0' A)" > "$1/$2"; samtools faidx "$1/$2"
    touch -d '2020-01-01' "$1/$2.fai"; cp toy9/genome.fa "$1/$2"
}
fresh() { [[ $(gawk '{ s += $2 } END { print s + 0 }' "$1" 2> /dev/null) -eq $GL ]]; }

# ---------- c: minimap2 that fails
mkdir -p nomm; printf '#!/bin/sh\necho "stub minimap2: failed" >&2\nexit 1\n' > nomm/minimap2; chmod +x nomm/minimap2
( PATH="$W/nomm:$PATH" bash "$FS/fs9_dupregion.sh" toy9/genome.fa toy9/copies.bed "$W/dupfail" 2 > dupfail.log 2>&1 ); RCD=$?
chk "c fs9_dupregion: a failed minimap2 stops it (exit $RCD)" "[ $RCD -ne 0 ] && grep -q 'minimap2 failed' dupfail.log"

# ---------- d: no dustmasker on PATH
mkdir -p nodust d3; for X in bash readlink gawk cat seqkit ssearch36 dirname; do ln -sf "$(command -v $X)" nodust/$X; done
printf ">F\nACGTACGTAC\n" > d3/cons.fa; : > d3/windows.masked.fa
( PATH="$W/nodust" "$(command -v bash)" "$FS/fs3_partners.sh" d3 d3/cons.fa 1 > d3.log 2>&1 ); RC3=$?
chk "d fs3: a missing dustmasker is named (exit $RC3)" "[ $RC3 -ne 0 ] && grep -q 'dustmasker not found' d3.log"

# ---------- e: stale index
stale e1 genome.fa
BODYCHECK=0 bash "$FS/fs9_twins.sh" e1/genome.fa toy9/copies.bed "$W/e1out" > e1.log 2>&1
chk "e fs9_twins: an index older than the genome is made again" "fresh e1/genome.fa.fai"
if command -v minimap2 > /dev/null; then
    stale e2 genome.fa
    bash "$FS/fs9_dupregion.sh" e2/genome.fa toy9/copies.bed "$W/e2out" 2 > e2.log 2>&1
    chk "e fs9_dupregion: an index older than the genome is made again" "fresh e2/genome.fa.fai"
else echo "SKIP e fs9_dupregion (no minimap2 on PATH)"; fi
# ---------- e, i: fs9_run; the first locus carries the strand "+,-"
mkdir -p run9/results; stale run9 genome.clean.fa
{ printf "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\n"
  gawk -F'\t' '{ print $1 ":" $2 "-" $3 "(" (NR == 1 ? "+,-" : $6) ")\tFAM\t500\t10\tassigned\t100" }' toy9/copies.bed; } > run9/results/assignment_full.tsv
MIX=$(gawk -F'\t' 'NR == 2 { print $1 }' run9/results/assignment_full.tsv)
BODYCHECK=0 bash "$FS/fs9_run.sh" run9 2 > run9.log 2>&1
chk "e fs9_run: an index older than the genome is made again" "fresh run9/genome.clean.fa.fai"
chk "i fs9_run: the (+,-) locus has a stage-9 status ($MIX)" "gawk -F'\t' -v L='$MIX' '\$1 == L { f = 1 } END { exit !f }' run9/results/flank_twins/FAM/copy_status.tsv"

# ---------- g: stage 6 without candidates
mkdir -p g6; : > g6/candidates.fa
bash "$FS/fs6_reassign.sh" g6 g6 1 > g6.log 2>&1
chk "g fs6: no-candidate header of reassign.tsv has 8 columns ($(head -1 g6/reassign.tsv | gawk -F'\t' '{ print NF }'))" "[ \"\$(head -1 g6/reassign.tsv | gawk -F'\t' '{ print NF }')\" -eq 8 ]"

# ---------- h: stage 6c a second time
mkdir -p h6/cand; : > h6/genome.clean.fa
printf "name\tpeak\ttype\tupstream\tdownstream\tn_peak\tmirrors\telements\tn_used\tlen_p10\tlen_med\tlen_p90\tcons_len\tstatus\n" > h6/candidates.tsv
printf "r1__r3_P1\tP1\tcomposite\tr1\tr3\t60\t-\t60\t60\t300\t310\t320\t310\textended:chain_C1\n" >> h6/candidates.tsv
printf "chain_C1\tC1\tchain\tr1\tr3\t40\t-\t40\t40\t500\t510\t520\t510\tkept\n" >> h6/candidates.tsv
printf "name\tside5_share\tside5_unit\tside5_second\tside3_share\tside3_unit\tside3_second\topen\nchain_C1\t0.00\t-\t-\t0.00\t-\t-\t-\n" > h6/ends.tsv
printf "name\tparent\tpeak\tside5_unit\ts5_qs\ts5_qe\ts5_gap\tside3_unit\ts3_qs\ts3_qe\ts3_gap\tcopies\nchain_C1\tr1__r3_P1\tC1\t-\t-\t-\t-\tr3\t1\t200\t0\t40\n" > h6/chains.tsv
: > h6/members.tsv; : > h6/candidates.fa
bash "$FS/fs6c_chain.sh" h6 h6 1 > h6.log 2>&1; RCH=$?
chk "h fs6c: a second run stops and chains.tsv keeps its chain row (exit $RCH, $(wc -l < h6/chains.tsv) lines)" "[ $RCH -ne 0 ] && [ \$(wc -l < h6/chains.tsv) -eq 2 ] && grep -q 'already holds chains' h6.log"

echo "fails=$fails"; [[ $fails -eq 0 ]] && echo "test_audit9: PASS"; exit $((fails > 0))
