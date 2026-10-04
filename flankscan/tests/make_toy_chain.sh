#!/usr/bin/env bash
# make_toy_chain.sh OUT_DIR - a SINEderella-like run dir with ONE case: a chain of three units.
# Consensuses TD (150 bp), TE (160 bp), TF (170 bp + 12 A). Planted, 8 % point mutations, random strand,
# >= 3 kb apart on a 2 Mb contig:
#   chain   60 x TD + 30 bp linker + TE + 40 bp linker + TF  (each unit assigned as its own locus, as
#           SINEderella splits a composite)
#   single  30 x each of TD, TE, TF alone
# The pairs TD-TE and TE-TF are each found by flankscan stages 3-5; each candidate built alone has an OPEN
# end (the third unit), which stage 6c must close by building the whole chain.
# Writes OUT/genome.clean.fa, consensuses.clean.fa, results/assignment_full.tsv, planted.fa, truth.tsv.
set -euo pipefail
OUT=$1; mkdir -p "$OUT/results"
gawk -v OUT="$OUT" '
function rnd(n,   s,i){ s=""; for(i=0;i<n;i++) s=s substr("ACGT",int(rand()*4)+1,1); return s }
function mut(s,   o,i,c){ o=""; for(i=1;i<=length(s);i++){ c=substr(s,i,1); if(rand()<0.08) c=substr("ACGT",int(rand()*4)+1,1); o=o c }; return o }
function rc(s,   o,i,c){ o=""; for(i=length(s);i>=1;i--){ c=substr(s,i,1); o=o (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A") }; return o }
function reps(u,n,   o,i){ o=""; for(i=0;i<n;i++) o=o u; return o }
function loc(base, x, y, st, L) { return "chr1:" base+(st=="+" ? x-1 : L-y) "-" base+(st=="+" ? y : L-x+1) "(" st ")" }
function put3(a, l1, b, l2, c, f1, f2, f3, cs,   e, st, s, L, ha, hb, a2, b2, a3, base, g){
    e=a l1 b l2 c; L=length(e); ha=length(a); hb=length(b)
    a2=ha+length(l1)+1; b2=a2+hb-1; a3=b2+length(l2)+1
    st=(rand()<0.5)?"+":"-"; s=(st=="-")?rc(e):e; g=rnd(3000+int(rand()*1500)); g1=g1 g
    base=length(g1); g1=g1 s
    id=loc(base,1,ha,st,L);  print id "\t" f1 "\t500\t10\tassigned\t100" >> ASG; print id "\t" cs "_1\t" f1 >> TR
    id=loc(base,a2,b2,st,L); print id "\t" f2 "\t500\t10\tassigned\t100" >> ASG; print id "\t" cs "_2\t" f2 >> TR
    id=loc(base,a3,L,st,L);  print id "\t" f3 "\t500\t10\tassigned\t100" >> ASG; print id "\t" cs "_3\t" f3 >> TR
}
function put1(e, f, cs,   st, s, L, base, g, id){
    L=length(e); st=(rand()<0.5)?"+":"-"; s=(st=="-")?rc(e):e; g=rnd(3000+int(rand()*1500)); g1=g1 g
    base=length(g1); g1=g1 s
    id=loc(base,1,L,st,L); print id "\t" f "\t500\t10\tassigned\t100" >> ASG; print id "\t" cs "\t" f >> TR
}
BEGIN{
    srand(11); ASG=OUT "/results/assignment_full.tsv"; TR=OUT "/truth.tsv"
    print "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold" > ASG
    print "locus\tcase\tfamily" > TR
    TD=rnd(150); TE=rnd(160); TF=rnd(170) reps("A",12); L1=rnd(30); L2=rnd(40)
    printf(">TD\n%s\n>TE\n%s\n>TF\n%s\n", TD, TE, TF) > (OUT "/consensuses.clean.fa")
    printf(">chain\n%s\n", TD L1 TE L2 TF) > (OUT "/planted.fa")
    g1=rnd(2000)
    for(i=0;i<60;i++) put3(mut(TD), L1, mut(TE), L2, mut(TF), "TD", "TE", "TF", "chain")
    for(i=0;i<30;i++){ put1(mut(TD), "TD", "single"); put1(mut(TE), "TE", "single"); put1(mut(TF), "TF", "single") }
    g1=g1 rnd(3000)
    printf(">chr1\n") > (OUT "/genome.clean.fa")
    for(i=1;i<=length(g1);i+=80) print substr(g1,i,80) >> (OUT "/genome.clean.fa")
}'
samtools faidx "$OUT/genome.clean.fa"
echo "toy chain: $(($(wc -l < "$OUT/truth.tsv") - 1)) planted loci in $OUT"
