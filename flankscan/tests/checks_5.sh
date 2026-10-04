# checks_5.sh - stages 5 (candidate consensuses) and 6 (re-assign once) against the planted truth;
# sourced by run_tests.sh after stage 4 (MINK 20 on the toy, as in checks_4.sh)
bash "$FS/fs5_build.sh" "$O" "$T" 60 20 4
bash "$FS/fs6_reassign.sh" "$O" "$T" 4

# K 'condition on a candidates.tsv row' -> number of such rows
K() { gawk -F'\t' 'NR>1 && ('"$1"') { n++ } END { print n+0 }' "$O/candidates.tsv"; }
chk "fs5: 6 peaks -> 3 built (mirrors folded), all kept" "[ \$(K '1') -eq 3 ] && [ \$(K '\$14==\"kept\"') -eq 3 ]"
chk "fs5: one composite TA__TB, one homodimer TB__TB, one piecewise TA__TB" \
    "[ \$(K '\$3==\"composite\" && \$4==\"TA\" && \$5==\"TB\"') -eq 1 ] && [ \$(K '\$3==\"homodimer\" && \$4==\"TB\"') -eq 1 ] && [ \$(K '\$3==\"piecewise\" && \$4==\"TA\" && \$5==\"TB\"') -eq 1 ]"
chk "fs5: distinct elements: composite 40, homodimer 30, piecewise 30 (>= 90 %)" \
    "[ \$(K '\$3==\"composite\" && \$8>=36') -eq 1 ] && [ \$(K '\$3==\"homodimer\" && \$8>=27') -eq 1 ] && [ \$(K '\$3==\"piecewise\" && \$8>=27') -eq 1 ]"

# every candidate against the planted elements: the right one, >= 95 % identity over >= 95 % of it
ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 "$O/candidates.fa" "$T/planted.fa" 2> /dev/null > "$O/cand_vs_planted.m8"
gawk '/^>/ { n = substr($1, 2); next } { L[n] += length($0) }
      END { for (n in L) print n "\t" L[n] }' "$T/planted.fa" > "$O/planted.len"
V() { gawk -F'\t' -v C="$1" -v P="$2" 'FNR==NR { L[$1]=$2; next }
          $1 ~ C && !seen[$1]++ { ok = ($2 == P && $3 >= 95 && $4 >= 0.95 * L[P]) } END { exit !ok }' \
          "$O/planted.len" "$O/cand_vs_planted.m8"; }
chk "fs5: composite consensus = planted TA[1-130] + 39 bp + TB (best hit, >= 95 % id, >= 95 % length)" "V \"^\$(gawk -F'\t' '\$3==\"composite\"{print \$1}' $O/candidates.tsv)\$\" dimer"
chk "fs5: homodimer consensus = planted TB + 20 bp + TB" "V '^TB__TB_' homodimer"
chk "fs5: piecewise consensus = planted TA[1-100] + TB[80-]" "V \"^\$(gawk -F'\t' '\$3==\"piecewise\"{print \$1}' $O/candidates.tsv)\$\" piecewise"
chk "fs5: consensus lengths within 5 % of the planted ones" \
    "gawk -F'\t' 'FNR==NR{L[\$1]=\$2;next} FNR>1{p=(\$3==\"composite\"?\"dimer\":\$3); d=\$13-L[p]; if (d<0) d=-d; if (d > 0.05*L[p]) bad++} END{exit bad>0}' $O/planted.len $O/candidates.tsv"

# stage 6
R() { gawk -F'\t' -v T="$1" 'FNR==NR { if ($3==T) c=$1; next } $1==c { print $'"$2"' }' "$O/candidates.tsv" "$O/reassign.tsv"; }
chk "fs6: every candidate accepted (>= 70 % of its elements one full unit)" \
    "[ \$(gawk -F'\t' 'NR>1 && \$6==\"accept\"' $O/reassign.tsv | wc -l) -eq 3 ]"
chk "fs6: composite: >= 36 of 40 planted elements read as one full unit" "[ \$(R composite 4) -ge 36 ]"
chk "fs6: homodimer and piecewise: >= 27 of 30 each" "[ \$(R homodimer 4) -ge 27 ] && [ \$(R piecewise 4) -ge 27 ]"
# no false moves: single TA copies and (TA)n-tailed TB copies keep their family
M() { gawk -F'\t' -v C="$1" -v F="$2" 'FNR==NR { if ($2==C) w[$1]; next }
          FNR>1 && $4=="main" && ($1 in w) && $5==F { n++ } END { print n+0 }' "$O/wid_case.tsv" "$O/reassign/units.tsv"; }
chk "fs6: single TA copies stay TA (>= 38/40)" "[ \$(M single TA) -ge 38 ]"
# a single TB copy may score as well against TB__TB as against TB (main = TB__TB, in part); what
# must not happen is a single copy reading as a whole candidate
MF() { gawk -F'\t' -v C="$1" 'FNR==NR { if ($2==C) w[$1]; next }
          FNR>1 && $4=="main" && ($1 in w) && $5 ~ /__/ && $9=="he" { n++ } END { print n+0 }' "$O/wid_case.tsv" "$O/reassign/units.tsv"; }
chk "fs6: no single TA / (TA)n-tailed TB copy reads as a whole candidate (0/60)" "[ \$(MF single) -eq 0 ] && [ \$(MF tatail) -eq 0 ]"
chk "fs6: both halves of a dimer counted as ONE element (composite elements <= 40)" "[ \$(R composite 2) -le 40 ] && [ \$(R composite 2) -ge 36 ]"
