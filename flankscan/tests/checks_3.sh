# checks_3.sh - stage 3 (partners) against the planted truth; sourced by run_tests.sh (T, O, FS, chk set)
bash "$FS/fs3_partners.sh" "$O" "$T/consensuses.clean.fa" 4

gawk -F'\t' 'FNR==NR { if (FNR>1) cas[$1]=$2; next } FNR>1 { print $1 "\t" cas[$2] }' \
    "$T/truth.tsv" "$O/loci.tsv" > "$O/wid_case.tsv"
# J CASE 'gawk condition on a junctions.tsv row' -> number of distinct windows of CASE with such a row
J() {
    gawk -F'\t' -v C="$1" 'FNR==NR { if ($2==C) want[$1]; next }
        FNR>1 && ($1 in want) && ('"$2"') && !($1 in hit) { hit[$1]; n++ } END { print n+0 }' \
        "$O/wid_case.tsv" "$O/junctions.tsv"
}

chk "fs3: main unit = assigned family (>= 195/197)" \
    "[ \$(gawk -F'\t' 'NR>1 && \$3==3 && \$4==\$2' $O/junctions.tsv | wc -l) -ge 195 ]"
# dimer: TA[1-130] + 39 bp linker + TB. However the local alignments split the linker, the main
# overrun past TA 130, the gap and the partner start past TB 1 add up to 39.
chk "fs3: dimer_left: TB on side 3, same strand, linker 39 +- 3 (>= 36/40)" \
    "[ \$(J dimer_left '\$3==3 && \$9==\"TB\" && \$10==\"same\" && (d=(\$8-130)+\$16+(\$11-1)) >= 36 && d <= 42') -ge 36 ]"
chk "fs3: dimer_right: TA on side 5, same strand, linker 39 +- 3 (>= 36/40)" \
    "[ \$(J dimer_right '\$3==5 && \$9==\"TA\" && \$10==\"same\" && (d=(\$11-130)+\$16+(\$8-1)) >= 36 && d <= 42') -ge 36 ]"
chk "fs3: chance: TA head on side 5 within 30 bp, tag h (>= 16/20)" \
    "[ \$(J chance_right '\$3==5 && \$9==\"TA\" && \$10==\"same\" && \$16<=30 && \$13 ~ /h/') -ge 16 ]"
chk "fs3: A-tail insertion: TC on side 5 ending at its tail (tag e), gap >= 70 % A (>= 13/15)" \
    "[ \$(J atail_right '\$3==5 && \$9==\"TC\" && \$10==\"same\" && \$13 ~ /e/ && \$20>=0.7') -ge 13 ]"
# nested: TC[1-80] + TA + TC[81-165]: TC on both sides, and its two facing ends are consecutive
chk "fs3: nested: TC split around the copy, facing ends 80|81 +- 10 (>= 13/15)" \
    "[ \$(gawk -F'\t' 'FNR==NR{if(\$2==\"nested\")w[\$1];next}
         FNR>1 && (\$1 in w) && \$9==\"TC\" && \$10==\"same\" { if (\$3==5) a[\$1]=\$11; else b[\$1]=\$11 }
         END { for (x in a) if ((x in b) && a[x]>=70 && a[x]<=90 && b[x]-a[x]>=-9 && b[x]-a[x]<=11) n++; print n+0 }' \
         $O/wid_case.tsv $O/junctions.tsv) -ge 13 ]"
chk "fs3: single copies: no partner within 200 bp (<= 1 of 40)" \
    "[ \$(J single '\$9!=\"-\" && \$9!=\"nomain\" && \$16<=200') -le 1 ]"
chk "fs3: satellite: next TA unit on side 3 after the 300 bp spacer (>= 4/6)" \
    "[ \$(J satellite '\$3==3 && \$9==\"TA\" && \$10==\"same\" && \$16>=285 && \$16<=325') -ge 4 ]"
chk "fs3: (TA)n-tailed TB: no partner within 100 bp on side 3 (0/20)" \
    "[ \$(J tatail '\$3==3 && \$9!=\"-\" && \$16<=100') -eq 0 ]"
chk "fs3: contig-end copy: side 5 reports clamp 800" "[ \$(J contigend '\$3==5 && \$19==800') -eq 1 ]"
# planted tails 12 / 12 / 15 A; the walker may add 1-3 A-rich random bases before them
chk "fs3: consensus A tails masked (planted 12/12/15 A, +0..3 bp)" \
    "gawk -F'\t' '{t[\$1]=\$2-\$3+1} END{exit !(t[\"TA\"]>=12 && t[\"TA\"]<=15 && t[\"TB\"]>=12 && t[\"TB\"]<=15 && t[\"TC\"]>=15 && t[\"TC\"]<=18)}' $O/cons.tsv"
