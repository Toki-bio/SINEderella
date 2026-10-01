#!/usr/bin/env bash
# fs6c_chain.sh OUT_DIR RUN_DIR [THREADS=8] [MINK=25] [MAXROUND=3]
#
# Stage 6c of flankscan: chains. A candidate whose end is OPEN (stage 6b: >= 50 % of its copies continue into
# a known bank unit at that end) is the middle of a longer chain - stages 3-5 pair every copy with its nearest
# neighbour only, so r1 + r3 + r3 is found as the pairs r1-r3 and r3-r3. Here every open candidate is
# EXTENDED over the unit found in the flank:
#   1) per copy, the flank unit hit of stage 6b gives the extension (the element grows to the outer end of
#      the hit; the hit's consensus span and the gap to the element are kept for the hierarchy);
#   2) only the copies that carry the commonest unit on every open side are kept (both sides open needs
#      >= MINK such copies, else the side with the larger share is extended alone); fewer than MINK
#      copies: the candidate stays open (chains.skip);
#   3) the extended copies are built exactly like a stage-5 candidate (fs_build_cons: flanked alignment,
#      consensus extended into the flanks) as a NEW candidate chain_Ck with its own peak id Ck, whose
#      member copies are registered in members.tsv (stage 6 counts them); the parent gets the status
#      extended:chain_Ck;
#   4) the new candidates are merged with the kept ones by best hit (fs_fold), and their own ends are checked;
#      a chain of four or more units is reached over MAXROUND rounds; ends still open after that stay open.
#
# In : OUT/candidates.tsv members.tsv ends.tsv cand/ (stages 4-6b), RUN_DIR/genome.clean.fa
# Out: OUT/candidates.tsv (+ chain rows, parents extended:...), candidates.fa, members.tsv, ends.tsv (+ chain
#      rows), cand/chain_Ck.{bed,fa,aln.fa}
#      OUT/chains.tsv   name parent peak side5_unit s5_qs s5_qe s5_gap side3_unit s3_qs s3_qe s3_gap copies
#                       (the unit added at each end: its consensus span and the gap to the element)
#      OUT/chains.skip  open candidates with too few copies carrying the flank unit
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/fs_lib.sh"
OUT=${1:?OUT_DIR}; RUN=$(readlink -f "${2:?RUN_DIR}"); T=${3:-8}; MINK=${4:-25}; MAXROUND=${5:-3}
FL=100; SHARE=0.6; SKIP=0.3; EFL=250
cd "$OUT"
[[ -s ends.tsv ]] || { echo "fs6c: no ends.tsv - run fs6b first" >&2; exit 1; }
printf "name\tparent\tpeak\tside5_unit\ts5_qs\ts5_qe\ts5_gap\tside3_unit\ts3_qs\ts3_qe\ts3_gap\tcopies\n" > chains.tsv
: > chains.skip
: > cand_families.tsv      # chain TAB comma-separated families of ALL its units (stage 6 re-assigns their copies)
K=$(gawk -F'\t' 'NR > 1 && $2 ~ /^C[0-9]+$/ { n = substr($2, 2) + 0; if (n > m) m = n } END { print m + 0 }' candidates.tsv)
new_total=0

for ((round = 1; round <= MAXROUND; round++)); do
    # open candidates that are still kept and not skipped
    gawk -F'\t' 'FILENAME == ARGV[1] { if (FNR > 1 && $14 == "kept") kept[$1]; next }
                 FILENAME == ARGV[2] { skip[$1]; next }
                 FNR > 1 && $8 != "-" && ($1 in kept) && !($1 in skip) { print }' candidates.tsv chains.skip ends.tsv > "cand/todo.tsv"
    [[ -s cand/todo.tsv ]] || break
    made=0; : > cand/status_changes.tsv
    while IFS=$'\t' read -r N S5 U5 _ S3 U3 _ OP <&3; do
        # the copies carrying the commonest unit on each open side (per-copy junction hits of stage 6b)
        : > cand/sel5.tsv; : > cand/sel3.tsv
        [[ "$OP" == *5* ]] && gawk -F'\t' -v U="$U5" '$2 == U' "cand/ends/$N.h5" > cand/sel5.tsv
        [[ "$OP" == *3* ]] && gawk -F'\t' -v U="$U3" '$2 == U' "cand/ends/$N.h3" > cand/sel3.tsv
        SIDES=$OP
        if [[ "$OP" == 53 ]]; then
            NB=$(gawk -F'\t' 'FNR == NR { a[$1]; next } ($1 in a)' cand/sel5.tsv cand/sel3.tsv | wc -l)
            if (( NB < MINK )); then
                if [[ $(gawk -v a="$S5" -v b="$S3" 'BEGIN { print (a >= b) }') == 1 ]]; then SIDES=5; : > cand/sel3.tsv; else SIDES=3; : > cand/sel5.tsv; fi
            fi
        fi
        # extended bed: the genomic window grows by the hit's outer end per open side
        gawk -F'\t' -v OFS='\t' -v SIDES="$SIDES" -v EFL=$EFL '
            FILENAME == ARGV[1] { len[$1] = $2; next }
            FILENAME == ARGV[2] { if ($1 != "") { e5[$1] = EFL - $5 + 1; h5[$1] = 1 }; next }     # sel5: flank lo
            FILENAME == ARGV[3] { e3[$1] = $6; h3[$1] = 1; next }                                 # sel3: flank hi
            { n = $4; if ((SIDES ~ /5/ && !(n in h5)) || (SIDES ~ /3/ && !(n in h3))) next
              a = (SIDES ~ /5/) ? e5[n] : 0; b = (SIDES ~ /3/) ? e3[n] : 0
              if ($6 == "-") { s = $2 - b; e = $3 + a } else { s = $2 - a; e = $3 + b }
              if (s < 0) s = 0; if (e > len[$1]) e = len[$1]
              print $1, s, e, $4, $5, $6 }' "$RUN/genome.clean.fa.fai" cand/sel5.tsv cand/sel3.tsv "cand/$N.bed" > cand/ext.bed
        NSEL=$(wc -l < cand/ext.bed)
        if (( NSEL < MINK )); then echo "$N" >> chains.skip; continue; fi
        K=$((K + 1)); CN="chain_C$K"; PK="C$K"
        cp cand/ext.bed "cand/$CN.bed"
        fs_build_cons "$CN"
        # what was added at each end: median consensus span of the unit and median gap to the element
        part() {   # part SIDE -> "unit qs qe gap"
            local F="cand/sel$1.tsv"; [[ "$SIDES" == *$1* ]] || { echo "- - - -"; return; }
            gawk -F'\t' -v S=$1 -v EFL=$EFL 'function med(a, n,   i, j, t) { for (i = 2; i <= n; i++) { t = a[i]; for (j = i - 1; j >= 1 && a[j] > t; j--) a[j + 1] = a[j]; a[j + 1] = t } return a[int((n + 1) / 2)] }
                { u = $2; qs[NR] = $3; qe[NR] = $4; g[NR] = (S == 5) ? EFL - $6 : $5 - 1; if (g[NR] < 0) g[NR] = 0 }
                END { print u, med(qs, NR), med(qe, NR), med(g, NR) }' OFS=' ' "$F"
        }
        read -r P5U P5S P5E P5G <<< "$(part 5)"; read -r P3U P3S P3E P3G <<< "$(part 3)"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$CN" "$N" "$PK" "$P5U" "$P5S" "$P5E" "$P5G" "$P3U" "$P3S" "$P3E" "$P3G" "$NSEL" >> chains.tsv
        gawk -F'\t' -v OFS='\t' -v P="$PK" '{ n = $4; sub(/_.*/, "", n); print P, n, 0 }' "cand/$CN.bed" >> members.tsv
        LENS=$(gawk -F'\t' '{ print $3 - $2 }' "cand/$CN.bed" | sort -n | gawk -v OFS='\t' '{ v[NR] = $1 } END { print v[int(NR * 0.1) + 1], v[int((NR + 1) / 2)], v[(int(NR * 0.9) > 0) ? int(NR * 0.9) : 1] }')
        CL=$(tail -1 "cand/$CN.fa" | tr -d '\n' | wc -c)
        # first / last unit family (the full bank names) for stage 6
        UP=$(gawk -F'\t' -v N="$N" -v U="$P5U" '$1 == N { print (U != "-" ? U : $4) }' candidates.tsv)
        DN=$(gawk -F'\t' -v N="$N" -v U="$P3U" '$1 == N { print (U != "-" ? U : $5) }' candidates.tsv)
        printf "%s\t%s\tchain\t%s\t%s\t%s\t-\t%s\t%s\t%s\t%s\tkept\n" "$CN" "$PK" "$UP" "$DN" "$NSEL" "$NSEL" "$NSEL" "$LENS" "$CL" >> candidates.tsv
        # every unit family of the chain (the parent's units + the added ones): stage 6 re-assigns the copies of ALL of them
        PF=$(gawk -F'\t' -v N="$N" 'FILENAME == ARGV[1] { f[$1] = $2; next } $1 == N { print (N in f) ? f[N] : $4 "," $5 }' cand_families.tsv candidates.tsv | head -1)
        printf "%s\t%s\n" "$CN" "$(echo "$PF,$P5U,$P3U" | tr ',' '\n' | grep -v '^-$' | grep -v '^$' | sort -u | paste -sd, -)" >> cand_families.tsv
        printf "%s\textended:%s\n" "$N" "$CN" >> cand/status_changes.tsv
        made=$((made + 1)); new_total=$((new_total + 1))
        echo "fs6c: $N (open $SIDES) -> $CN: $NSEL copies, $CL bp, added ${P5U}(5') ${P3U}(3')" >&2
    done 3< cand/todo.tsv
    # parents extended; merge the new chains with the kept candidates by best hit
    { head -1 candidates.tsv; gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { st[$1] = $2; next } FNR > 1 { if ($1 in st) $14 = st[$1]; print }' \
        cand/status_changes.tsv candidates.tsv; } > cand/cand.new
    mv cand/cand.new candidates.tsv
    if (( made > 0 )); then
        gawk -F'\t' 'NR > 1 && $14 == "kept" { print $6 "\t" $1 "\t" $13 }' candidates.tsv | sort -t$'\t' -k1,1nr | cut -f2,3 > cand/list.tsv
        fs_fold cand/list.tsv cand/status.tsv
        gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { st[$1] = $2; next } FNR == 1 { print; next } { if ($1 in st) $14 = st[$1]; print }' cand/status.tsv candidates.tsv > cand/cand.new
        mv cand/cand.new candidates.tsv
        for CN in $(gawk -F'\t' 'NR > 1 && $3 == "chain" && $14 == "kept" { print $1 }' candidates.tsv); do
            grep -q "^$CN	" ends.tsv || fs_endcheck "$CN" >> ends.tsv
        done
    fi
    gawk -F'\t' 'NR > 1 && $14 == "kept" { print $1 }' candidates.tsv | while read -r X; do cat "cand/$X.fa"; done > candidates.fa
    rm -f cand/todo.tsv cand/sel5.tsv cand/sel3.tsv cand/ext.bed cand/list.tsv cand/status.tsv cand/status_changes.tsv
    (( made > 0 )) || break
done
echo "fs6c: $new_total chain candidates in $((MAXROUND)) rounds max; $(gawk -F'\t' 'NR > 1 && $14 == "kept"' candidates.tsv | wc -l) kept; open and skipped: $(wc -l < chains.skip) -> $OUT/chains.tsv" >&2
