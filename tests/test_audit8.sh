#!/usr/bin/env bash
# test_audit8.sh [WORKDIR] - the four findings of docs/AUDIT_2026-11-10.md (flankscan), one block each:
#   1  fs3_partners.sh: a host hit that runs across an inserted copy is cut at the gap; 30 planted nested copies
#      (15 on each strand, 212 bp inside a 312 bp host) are all called nested by stage 4
#   2  fs_fold (fs_lib.sh): a candidate that lacks 148 bp of a longer one is not folded into it
#   3  fs2_trf.sh: a trf that ends a good run with exit 253 (TRF 4.09.1) is accepted; a trf killed by a signal stops stage 2
#   4  fs9_run.sh: the 20-mer count is made once with jellyfish for a small genome too, and fs9_twins.sh gives the same
#      copy_status.tsv with that count as with its own awk count
# Needs gawk, bedtools, samtools, seqkit, ssearch36, dustmasker, python3 + numpy, jellyfish (block 4), no trf (stubs).
# FS overrides the flankscan directory tested (old code: FS=<copy of main>/flankscan).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; FS=${FS:-$(dirname "$HERE")/flankscan}
W=${1:-${TMPDIR:-/tmp}/audit8_toy}; rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails + 1)); fi; }
# trf stubs: no repeat found in any window ("@name" lines only), exit code as given / killed
mkdir -p stub253 stubkill
printf '#!/bin/sh\ngrep "^>" "$1" | sed "s/^>/@/"\nexit 253\n' > stub253/trf
printf '#!/bin/sh\nkill -9 $$\n' > stubkill/trf
chmod +x stub253/trf stubkill/trf

# ---------- 1 and 3: stages 1-4 on 30 copies of TI nested in a host H
mkdir -p run/results
gawk 'function rnd(n,   s, i) { s = ""; for (i = 1; i <= n; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1); return s }
      function mut(s, p,   o, i, c) { o = ""; for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); if (rand() < p) c = substr("ACGT", int(rand() * 4) + 1, 1); o = o c } return o }
      function rc(s,   o, i, c) { o = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1); o = o (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : "A") } return o }
      BEGIN { srand(11); H = rnd(300); I = rnd(200); A = "AAAAAAAAAAAA"
          print ">H\n" H A > "run/consensuses.clean.fa"; print ">TI\n" I A > "run/consensuses.clean.fa"
          print "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold" > "run/results/assignment_full.tsv"
          g = ""
          for (k = 1; k <= 30; k++) { g = g rnd(3000); pos = length(g)
              el = mut(substr(H, 1, 150), 0.04) mut(I A, 0.04) mut(substr(H, 151) A, 0.04)       # 150 + 212 + 162
              if (k % 2) { s = pos + 150; st = "+" } else { el = rc(el); s = pos + 162; st = "-" }
              g = g el
              print "chr1:" s "-" s + 212 "(" st ")\tTI\t500\t10\tassigned\t100" > "run/results/assignment_full.tsv" }
          g = g rnd(3000); print ">chr1"; for (i = 1; i <= length(g); i += 80) print substr(g, i, 80) }' > run/genome.clean.fa
samtools faidx run/genome.clean.fa
bash "$FS/fs1_extract.sh" run out 1000 2> /dev/null
( PATH="$W/stub253:$PATH" bash "$FS/fs2_trf.sh" out 4 2> fs2.err ); RC2=$?
chk "3 fs2: a trf that exits 253 on success is accepted (stage 2 exit $RC2, trf.tsv written)" "[ $RC2 -eq 0 ] && [ -s out/trf.tsv ] && [ -s out/windows.masked.fa ]"
if [[ ! -s out/windows.masked.fa ]]; then    # old code: stage 2 stopped; stage 3 still needs its output (no flank repeat = unmasked)
    printf "wid\tfamily\tclass\tstart\tend\tperiod\tcopies\tpct_match\tscore\tmotif\n" > out/trf.tsv
    seqkit seq -w 0 out/windows.fa > out/windows.masked.fa 2> /dev/null
fi
mkdir -p outk; cp out/loci.tsv out/windows.fa outk/
( PATH="$W/stubkill:$PATH" bash "$FS/fs2_trf.sh" outk 4 2> fs2k.err ); RCK=$?
chk "3 fs2: a trf killed by a signal stops stage 2 (exit $RCK)" "[ $RCK -ne 0 ] && grep -q 'trf failed' fs2k.err"
bash "$FS/fs3_partners.sh" out run/consensuses.clean.fa 4 2> /dev/null
bash "$FS/fs4_junctions.sh" out run 20 2> /dev/null
NN=$(gawk -F'\t' 'NR > 1 && $4 == "nested"' out/copies.tsv | wc -l)
chk "1 fs3 + fs4: nested copies called nested ($NN of 30, both strands)" "[ $NN -eq 30 ]"
# the host ends facing the copy, in H coordinates: 150 on one side, 151 on the other (+- 5)
NJ=$(gawk -F'\t' 'NR > 1 && $9 == "H" && (($3 == 5 && $11 >= 145 && $11 <= 155) || ($3 == 3 && $11 >= 146 && $11 <= 156))' out/junctions.tsv | wc -l)
chk "1 fs3: host junction coordinates 150 | 151 +- 5 on both sides ($NJ of 60)" "[ $NJ -ge 57 ]"

# ---------- 2: fs_fold on a 381 bp candidate and a 233 bp one that lacks its middle 148 bp
mkdir -p fold/cand; (
    cd fold || exit 1
    gawk 'function rnd(n,   s, i) { s = ""; for (i = 1; i <= n; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1); return s }
          BEGIN { srand(5); a = rnd(100); b = rnd(148); c = rnd(133)
                  print ">long\n" a b c > "cand/long.fa"; print ">short\n" a c > "cand/short.fa"; print ">u1\n" rnd(200) > "cons.masked.fa" }'
    printf "long\t381\t-\nshort\t233\t-\n" > list.tsv
    source "$FS/fs_lib.sh"; fs_fold list.tsv status.tsv 2> fold.err )
chk "2 fs_fold: the 233 bp candidate that lacks 148 bp of the 381 bp one is kept ($(tr '\t\n' ' ;' < fold/status.tsv))" \
    "[ \"\$(gawk -F'\t' '\$2 == \"kept\"' fold/status.tsv | wc -l)\" -eq 2 ]"

# ---------- 4: stage 9 with the jellyfish count
if command -v jellyfish > /dev/null; then
    python3 "$FS/tests/make_toy_flank.py" "$W/toy9" 1 > /dev/null
    mkdir -p run9/results; cp toy9/genome.fa run9/genome.clean.fa; samtools faidx run9/genome.clean.fa
    { printf "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\n"
      gawk -F'\t' '{ print $1 ":" $2 "-" $3 "(" $6 ")\tFAM\t500\t10\tassigned\t100" }' toy9/copies.bed; } > run9/results/assignment_full.tsv
    bash "$FS/fs9_run.sh" run9 4 > fs9_run.log 2>&1; RC9=$?
    chk "4 fs9_run: 20-mer count made with jellyfish for a small genome (exit $RC9)" "[ $RC9 -eq 0 ] && [ -s run9/genome.clean.fa.k20gt20.tsv ] && [ -s run9/results/flank_twins/FAM/copy_status.tsv ]"
    # the same family through fs9_twins.sh with no KCOUNT (its own awk count): the same copy status and pairs
    bash "$FS/fs9_twins.sh" run9/genome.clean.fa run9/results/flank_twins/FAM/copies.bed "$W/awk9" > awk9.log 2>&1
    chk "4 fs9_twins: copy_status.tsv and flank_groups.tsv identical with the jellyfish count and with the awk count ($(wc -l < awk9/copy_status.tsv) copies)" \
        "cmp -s awk9/copy_status.tsv run9/results/flank_twins/FAM/copy_status.tsv && cmp -s awk9/flank_groups.tsv run9/results/flank_twins/FAM/flank_groups.tsv"
else
    echo "SKIP 4 (no jellyfish on PATH)"
fi
echo "fails=$fails"; [[ $fails -eq 0 ]] && echo "test_audit8: PASS"; exit $((fails > 0))
