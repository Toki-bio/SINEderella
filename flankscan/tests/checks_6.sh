# checks_6.sh - stage 3 with a partner library (fs3 4th argument): a part missing from the bank
# gets a name, and never becomes a copy's main unit. Sourced by run_tests.sh after stage 3.
P6=$O/px; rm -rf "$P6"; mkdir -p "$P6"
cp "$O/loci.tsv" "$O/windows.fa" "$O/windows.masked.fa" "$P6/"
bash "$FS/fs3_partners.sh" "$P6" "$T/consensuses.clean.fa" 4 "$T/partnerlib.fa"
[[ -s "$O/wid_case.tsv" ]] || gawk -F'\t' 'FNR==NR { if (FNR>1) c[$1]=$2; next } FNR>1 { print $1 "\t" c[$2] }' \
    "$T/truth.tsv" "$O/loci.tsv" > "$O/wid_case.tsv"
# S DIR CASE 'condition on the side-5 junction row' -> copies of that case
S() { gawk -F'\t' -v C="$2" 'FNR==NR { if ($2==C) w[$1]; next }
          FNR>1 && ($1 in w) && $3==5 && ('"$3"') { n++ } END { print n+0 }' "$O/wid_case.tsv" "$1/junctions.tsv"; }
chk "fs3+lib: without the library, TA after a LINE end has no partner on side 5 (>= 18/20)" \
    "[ \$(S $O after_line '\$9==\"-\"') -ge 18 ]"
chk "fs3+lib: with it, the partner on side 5 is x.LE, same strand, gap <= 20 = its A10 tail + slack (>= 18/20)" \
    "[ \$(S $P6 after_line '\$9==\"x.LE\" && \$10==\"same\" && \$16<=20') -ge 18 ]"
chk "fs3+lib: no copy has a partner-library unit as its main unit" \
    "[ \$(gawk -F'\t' 'FNR>1 && \$4==\"main\" && \$5 ~ /^x\./' $P6/units.tsv | wc -l) -eq 0 ]"
chk "fs3+lib: main unit family unchanged by the library (all copies)" \
    "cmp -s <(gawk -F'\t' 'FNR>1 && \$4==\"main\" {print \$1, \$5}' $O/units.tsv | sort) <(gawk -F'\t' 'FNR>1 && \$4==\"main\" {print \$1, \$5}' $P6/units.tsv | sort)"
chk "fs3+lib: the dimer linker is still 39 +- 3 on side 3 of dimer_left (>= 36/40)" \
    "[ \$(gawk -F'\t' 'FNR==NR { if (\$2==\"dimer_left\") w[\$1]; next } FNR>1 && (\$1 in w) && \$3==3 && \$9==\"TB\" && (d=(\$8-130)+\$16+(\$11-1))>=36 && d<=42 {n++} END {print n+0}' $O/wid_case.tsv $P6/junctions.tsv) -ge 36 ]"
# stage 4 on the library run: a library partner has no copy count (crashed with a division by zero
# on rsi before E was guarded); no x.* row in the family summary
cp "$O/trf.tsv" "$P6/"
chk "fs4+lib: runs, and family_summary.tsv has no partner-library row" \
    "bash $FS/fs4_junctions.sh $P6 $T 20 && ! grep -q '^x\.' $P6/family_summary.tsv"
