# checks_4.sh - stage 4 (peaks + per-copy classes) against the planted truth; sourced by run_tests.sh
# The toy has 30-40 copies per planted structure, so MINK = 20 here (real data: 50).
bash "$FS/fs4_junctions.sh" "$O" "$T" 20

# C CASE CLASS -> number of copies of that planted case given that class
C() { gawk -F'\t' -v C="$1" -v K="$2" 'FNR==NR { if ($2==C) want[$1]; next }
          FNR>1 && ($1 in want) && $4==K { n++ } END { print n+0 }' "$O/wid_case.tsv" "$O/copies.tsv"; }
# P 'condition on a peaks.tsv row' -> number of such peaks
P() { gawk -F'\t' 'NR>1 && ('"$1"') { n++ } END { print n+0 }' "$O/peaks.tsv"; }

chk "fs4: peak TA side 3 TB composite, linker 39 +- 3 (whole-TA linker dimer, >= 36)" \
    "[ \$(P '\$1==\"TA\" && \$2==3 && \$3==\"TB\" && \$5==\"composite\" && \$6>=36 && (d=(\$9-130)+\$11+(\$10-1))>=36 && d<=42') -eq 1 ]"
chk "fs4: peak TB side 5 TA composite (>= 36)" \
    "[ \$(P '\$1==\"TB\" && \$2==5 && \$3==\"TA\" && \$5==\"composite\" && \$6>=36') -eq 1 ]"
chk "fs4: peak TB side 5 TA piecewise, TB starts ~80 after TA ~100 (>= 27)" \
    "[ \$(P '\$1==\"TB\" && \$2==5 && \$3==\"TA\" && \$5==\"piecewise\" && \$6>=27 && \$9>=70 && \$9<=90 && \$10>=90 && \$10<=110') -eq 1 ]"
chk "fs4: peak TA side 3 TB piecewise (>= 27)" \
    "[ \$(P '\$1==\"TA\" && \$2==3 && \$3==\"TB\" && \$5==\"piecewise\" && \$6>=27') -eq 1 ]"
chk "fs4: peaks TB side 3 and side 5 TB homodimer, gap 20 +- 3 (>= 27 each)" \
    "[ \$(P '\$1==\"TB\" && \$3==\"TB\" && \$5==\"homodimer\" && \$6>=27 && \$11>=17 && \$11<=23') -eq 2 ]"
chk "fs4: the dimer linker consensus is recovered (>= 90 % id)" \
    "[ \$(P '\$1==\"TA\" && \$2==3 && \$5==\"composite\" && \$14!=\".\" && \$13>=0.9') -eq 1 ]"
chk "fs4: no peak for the chance TA+TC tandems (variable cut and gap)" "[ \$(P '\$1==\"TC\"') -eq 0 ]"
chk "fs4: no other peaks (exactly 6)" "[ \$(P '1') -eq 6 ]"

chk "fs4: dimer_left -> composite (>= 36/40)"  "[ \$(C dimer_left composite) -ge 36 ]"
chk "fs4: dimer_right -> composite (>= 36/40)" "[ \$(C dimer_right composite) -ge 36 ]"
chk "fs4: homo_left / homo_right -> homodimer (>= 27/30 each)" \
    "[ \$(C homo_left homodimer) -ge 27 ] && [ \$(C homo_right homodimer) -ge 27 ]"
chk "fs4: piece_left / piece_right -> piecewise (>= 27/30 each)" \
    "[ \$(C piece_left piecewise) -ge 27 ] && [ \$(C piece_right piecewise) -ge 27 ]"
chk "fs4: nested -> nested (>= 13/15)"          "[ \$(C nested nested) -ge 13 ]"
chk "fs4: atail_right -> atail (>= 13/15)"      "[ \$(C atail_right atail) -ge 13 ]"
chk "fs4: chance_right -> chance (>= 16/20)"    "[ \$(C chance_right chance) -ge 16 ]"
chk "fs4: satellite -> satellite (>= 4/6)"      "[ \$(C satellite satellite) -ge 4 ]"
chk "fs4: single -> single (>= 38/40)"          "[ \$(C single single) -ge 38 ]"
chk "fs4: tatail -> single (>= 18/20)"          "[ \$(C tatail single) -ge 18 ]"
chk "fs4: family_summary.tsv header first, TA TB TC listed" \
    "head -1 $O/family_summary.tsv | grep -q '^family' && [ \$(tail -n +2 $O/family_summary.tsv | wc -l) -eq 3 ]"
