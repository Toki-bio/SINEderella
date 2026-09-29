#!/usr/bin/env bash
# make_toy.sh OUT_DIR - a SINEderella-like run dir with every case flankscan must recognise, planted
# at known places. Pure bash + gawk. Deterministic (srand(7)).
#
# Consensuses (random sequence, fixed seed):
#   TA  180 bp + 12 A tail   the "head" family
#   TB  200 bp + 12 A tail   partner of the structural dimer
#   TC  150 bp + 15 A tail   host / chance partner
# Planted, 8 % point mutations per copy, random strand, >= 3 kb apart on a 2 Mb contig + a 20 kb one:
#   single     40 x TA alone
#   dimer      40 x TA[1-130] + fixed 39 bp linker + TB   (structural composite; locus = TB part,
#              assigned TB; the TA head is a separate assigned TA locus, as SINEderella splits it)
#   chance     20 x TA[1-(90..170)] + random gap 0-25 + TC  (locus = TC)
#   atail      15 x TC with its A tail + A run of 3-20 + TA   (TA inserted into TC's tail; locus = TA)
#   nested     15 x TC[1-80] + TA + TC[81-165]               (TA inside TC; locus = TA)
#   satellite  1 array of 6 units, unit = TA + 300 bp, on the 20 kb contig (loci = the 6 TA copies)
#   tatail     20 x TB whose A tail is followed by (TA)15    (a microsatellite SINE tail; locus = TB)
#   homodimer  30 x TB + fixed 20 bp + TB   (loci = both TB parts)
#   piecewise  30 x TA[1-100] + TB[80-212], no gap  (loci = the TA part and the TB part, which
#              starts mid-consensus)
#   contigend  1 x TA 200 bp from the start of the 20 kb contig (clamp test)
# Writes OUT/genome.clean.fa, OUT/consensuses.clean.fa, OUT/results/assignment_full.tsv (the
# SINEderella layout) and OUT/truth.tsv (locus id, case, family).
set -euo pipefail
OUT=$1; mkdir -p "$OUT/results"
gawk -v OUT="$OUT" '
function rnd(n,   s,i){ s=""; for(i=0;i<n;i++) s=s substr("ACGT",int(rand()*4)+1,1); return s }
function mut(s,   o,i,c){ o=""; for(i=1;i<=length(s);i++){ c=substr(s,i,1); if(rand()<0.08) c=substr("ACGT",int(rand()*4)+1,1); o=o c }; return o }
function rc(s,   o,i,c){ o=""; for(i=length(s);i>=1;i--){ c=substr(s,i,1); o=o (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A") }; return o }
function reps(u,n,   o,i){ o=""; for(i=0;i<n;i++) o=o u; return o }
# place: insert element text at position p of contig 1 (as a list of pieces, assembled at the end)
function put(e, cas, fam, coreFrom, coreTo,   st, s, a, b){
    # coreFrom/coreTo: 1-based span of the ASSIGNED core inside e (forward orientation of e)
    st = (rand()<0.5) ? "+" : "-"
    if (st=="-") { s=rc(e); a=length(e)-coreTo+1; b=length(e)-coreFrom+1 } else { s=e; a=coreFrom; b=coreTo }
    gap = rnd(3000 + int(rand()*1500))
    g1 = g1 gap
    # core genomic 0-based start / end on contig 1
    cs = length(g1) + a - 1; ce = length(g1) + b
    g1 = g1 s
    id = "chr1:" cs "-" ce "(" st ")"
    print id "\t" fam "\t500\t10\tassigned\t100" >> ASG
    print id "\t" cas "\t" fam >> TR
}
# put2: a two-part element h + link + t; SINEderella assigns its two parts as two loci - plant both
function put2(h, link, t, f1, f2, c1, c2,   e, st, s, L, hl, bs, a1, b1, a2, b2, base){
    e=h link t
    st=(rand()<0.5)?"+":"-"; s=(st=="-")?rc(e):e; gap=rnd(3000+int(rand()*1500)); g1=g1 gap
    L=length(e); hl=length(h); bs=hl+length(link)+1
    if(st=="+"){ a1=1; b1=hl; a2=bs; b2=L } else { a1=L-hl+1; b1=L; a2=1; b2=L-bs+1 }
    base=length(g1); g1=g1 s
    print "chr1:" base+a1-1 "-" base+b1 "(" st ")\t" f1 "\t500\t10\tassigned\t100" >> ASG
    print "chr1:" base+a1-1 "-" base+b1 "(" st ")\t" c1 "\t" f1 >> TR
    print "chr1:" base+a2-1 "-" base+b2 "(" st ")\t" f2 "\t500\t10\tassigned\t100" >> ASG
    print "chr1:" base+a2-1 "-" base+b2 "(" st ")\t" c2 "\t" f2 >> TR
}
BEGIN{
    srand(7); ASG=OUT "/results/assignment_full.tsv"; TR=OUT "/truth.tsv"
    print "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold" > ASG
    print "locus\tcase\tfamily" > TR
    TA=rnd(180) reps("A",12); TB=rnd(200) reps("A",12); TC=rnd(150) reps("A",15)
    LINK=rnd(39)
    printf(">TA\n%s\n>TB\n%s\n>TC\n%s\n", TA, TB, TC) > (OUT "/consensuses.clean.fa")
    g1=rnd(2000)
    for(i=0;i<40;i++) put(mut(TA), "single", "TA", 1, length(TA))
    for(i=0;i<40;i++){ h=mut(substr(TA,1,130)); put2(h, LINK, mut(TB), "TA", "TB", "dimer_left", "dimer_right") }
    for(i=0;i<20;i++){ k=90+int(rand()*81); gp=int(rand()*26); h=mut(substr(TA,1,k)) rnd(gp); c=mut(TC)
        put(h c, "chance_right", "TC", length(h)+1, length(h)+length(c)) }
    for(i=0;i<15;i++){ c=mut(TC) reps("A",3+int(rand()*18)); a=mut(TA)
        put(c a, "atail_right", "TA", length(c)+1, length(c)+length(a)) }
    for(i=0;i<15;i++){ c=mut(TC); a=mut(TA); e=substr(c,1,80) a substr(c,81,85)
        put(e, "nested", "TA", 81, 80+length(a)) }
    for(i=0;i<20;i++){ b=mut(TB) reps("TA",15); put(b, "tatail", "TB", 1, length(TB)) }
    # homodimer: TB + fixed 20 bp + TB; piecewise: TA[1-100] joined directly to TB[80-] (the
    # downstream part starts mid-consensus: one element that two consensuses each cover in part)
    LINK2=rnd(20)
    for(i=0;i<30;i++){ h=mut(TB); put2(h, LINK2, mut(TB), "TB", "TB", "homo_left", "homo_right") }
    for(i=0;i<30;i++){ h=mut(substr(TA,1,100)); put2(h, "", mut(substr(TB,80)), "TA", "TB", "piece_left", "piece_right") }
    # contig 2: satellite array (6 x (TA + 300 bp)) and one TA near the contig start
    unit=rnd(300); g2=rnd(200)
    e=mut(TA); cs=length(g2); g2=g2 e
    print "chr2:" cs "-" cs+length(e) "(+)\tTA\t500\t10\tassigned\t100" >> ASG
    print "chr2:" cs "-" cs+length(e) "(+)\tcontigend\tTA" >> TR
    g2=g2 rnd(4000)
    for(i=0;i<6;i++){ e=mut(TA); cs=length(g2); g2=g2 e mut(unit)
        print "chr2:" cs "-" cs+length(e) "(+)\tTA\t500\t10\tassigned\t100" >> ASG
        print "chr2:" cs "-" cs+length(e) "(+)\tsatellite\tTA" >> TR }
    g2=g2 rnd(4000)
    g1=g1 rnd(3000)
    printf(">chr1\n") > (OUT "/genome.clean.fa")
    for(i=1;i<=length(g1);i+=80) print substr(g1,i,80) >> (OUT "/genome.clean.fa")
    printf(">chr2\n") >> (OUT "/genome.clean.fa")
    for(i=1;i<=length(g2);i+=80) print substr(g2,i,80) >> (OUT "/genome.clean.fa")
}'
samtools faidx "$OUT/genome.clean.fa"
echo "toy: $(($(wc -l < "$OUT/truth.tsv") - 1)) planted loci in $OUT"
