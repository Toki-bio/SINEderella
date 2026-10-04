# checks_2.sh - stage 2 (TRF) against the planted truth; sourced by run_tests.sh (T, O, FS, chk set)
bash "$FS/fs2_trf.sh" "$O"

# case of every window: wid -> case (truth.tsv locus = loci.tsv column 2)
gawk -F'\t' 'FNR==NR { if (FNR>1) cas[$1]=$2; next } FNR>1 { print $1 "\t" cas[$2] }' \
    "$T/truth.tsv" "$O/loci.tsv" > "$O/wid_case.tsv"
# windows of case C that have at least one repeat of class K (motif filter optional): count
n_with() {  # n_with CASE CLASS [MOTIF_REGEX]
    gawk -F'\t' -v C="$1" -v K="$2" -v M="${3:-.}" '
        FNR==NR { if ($2==C) want[$1]; next }
        FNR>1 && ($1 in want) && $3==K && $10 ~ M && !($1 in hit) { hit[$1]; n++ }
        END { print n+0 }' "$O/wid_case.tsv" "$O/trf.tsv"
}

chk "fs2: (TA)n-tailed TB copies get a tail repeat with motif TA/AT (>= 18/20)" \
    "[ \$(n_with tatail tail '^(TA|AT)+\$') -ge 18 ]"
chk "fs2: satellite array units are classed satellite (>= 4/6)" \
    "[ \$(n_with satellite satellite) -ge 4 ]"
chk "fs2: single TA copies: no satellite (0/40)" "[ \$(n_with single satellite) -eq 0 ]"
# a non-A tail on a plain copy would be a false boundary extension
chk "fs2: single TA copies: no tail other than (A)n (0/40)" \
    "[ \$(gawk -F'\t' 'FNR==NR{if(\$2==\"single\")w[\$1];next} FNR>1 && (\$1 in w) && \$3==\"tail\" && \$10 !~ /^A+\$/' $O/wid_case.tsv $O/trf.tsv | wc -l) -eq 0 ]"
for c in dimer_left dimer_right nested atail_right chance_right; do
    chk "fs2: $c copies: no satellite" "[ \$(n_with $c satellite) -eq 0 ]"
done
# masking touches flanks only: same lengths, identical cores, and every masked base lies in a flank5/flank3 repeat
chk "fs2: masked windows keep length and core exactly (all)" \
    "gawk -F'\t' 'FNR==NR { if (FNR>1) { cs[\$1]=\$8; ce[\$1]=\$9 }; next }
                 /^>/ { w=substr(\$0,2); next }
                 { if (w in s) { if (length(\$0)!=length(s[w]) || substr(\$0,cs[w],ce[w]-cs[w]+1)!=substr(s[w],cs[w],ce[w]-cs[w]+1)) bad++ } else s[w]=\$0 }
                 END { exit bad>0 }' $O/loci.tsv <(seqkit seq -w 0 $O/windows.fa) <(seqkit seq -w 0 $O/windows.masked.fa)"
chk "fs2: every N in a masked window lies inside a flank5/flank3 repeat" \
    "gawk -F'\t' 'FNR==NR { if (FNR>1 && (\$3==\"flank5\"||\$3==\"flank3\")) iv[\$1]=iv[\$1] \" \" \$4 \"-\" \$5; next }
                 /^>/ { w=substr(\$0,2); next }
                 { n=split(\$0,c,\"\"); for (i=1;i<=n;i++) if (c[i]==\"N\") { ok=0; k=split(iv[w],q,\" \")
                       for (j=1;j<=k;j++) { split(q[j],ab,\"-\"); if (i>=ab[1] && i<=ab[2]) { ok=1; break } }
                       if (!ok) bad++ } }
                 END { exit bad>0 }' $O/trf.tsv <(seqkit seq -w 0 $O/windows.masked.fa)"
# the element's own A tail (inside the core, ending at its 3' end) is class tail, never core
chk "fs2: no (A)n ending within 15 bp of the core 3' end is classed core (all)" \
    "[ \$(gawk -F'\t' 'FNR==NR{if(FNR>1)ce[\$1]=\$9;next} FNR>1 && \$3==\"core\" && \$10 ~ /^A+\$/ && \$5>=ce[\$1]-15' $O/loci.tsv $O/trf.tsv | wc -l) -eq 0 ]"
chk "fs2: single TA copies: own A tail found as tail (>= 25/40)" "[ \$(n_with single tail '^A+\$') -ge 25 ]"
# the host's A tail right before a copy inserted into it is an (A)n head
chk "fs2: A-tail insertions: (A)n head (>= 13/15)" "[ \$(n_with atail_right head '^A+\$') -ge 13 ]"
chk "fs2: masking is not vacuous (> 0 N)" "[ \$(grep -v '^>' $O/windows.masked.fa | tr -cd N | wc -c) -gt 0 ]"
chk "fs2: trf_summary.tsv header is the first line" "head -1 $O/trf_summary.tsv | grep -q '^family'"
