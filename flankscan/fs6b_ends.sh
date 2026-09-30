#!/usr/bin/env bash
# fs6b_ends.sh OUT_DIR RUN_DIR [THREADS=8]
#
# Stage 6b of flankscan: are the ENDS of each kept candidate closed? A candidate is built from a pair of
# units (stages 3-5 pair every copy with its nearest neighbour per side), so in a chain of three or more
# units (rsi r1 + r3 + r3) the pairs are found separately, and a pair built alone is cut short: its copies
# still continue into a known unit beyond the end, yet stage 6 accepts it because its copies read as full
# units of it. Found on rsi P12 (r3 + r3, 153 copies, "accept" 88 %): its 5' flank is r1 in 57 / 60 copies.
#
# Per kept candidate and side, over the copies the candidate was built from (cand/NAME.bed): the flank
# (FL bp beyond the element end, strand-aware) is searched against the bank (masked consensuses, stage 3;
# ssearch36 E <= 1e-3); a copy's end is OPEN when its best hit lies within NEAR bp of the junction and is
# >= MINLEN bp. The side is open when >= OPENFRAC of the copies are open; the unit reported is the commonest.
#
# In : OUT/candidates.tsv, OUT/cand/NAME.bed, OUT/cons.masked.fa ; RUN_DIR/genome.clean.fa(.fai)
# Out: OUT/ends.tsv   name  side5_share side5_unit side5_second side3_share side3_unit side3_second open
#                      (share = fraction of the candidate's copies open on that side; open = 5 / 3 / 53 / -)
set -euo pipefail
OUT=${1:?OUT_DIR}; RUN=$(readlink -f "${2:?RUN_DIR}"); T=${3:-8}
FL=250; NEAR=30; MINLEN=40; OPENFRAC=0.5
cd "$OUT"; mkdir -p cand/ends; G="$RUN/genome.clean.fa"
printf "name\tside5_share\tside5_unit\tside5_second\tside3_share\tside3_unit\tside3_second\topen\n" > ends.tsv
gawk -F'\t' 'NR > 1 && $14 == "kept" { print $1 }' candidates.tsv | while read -r N; do
    B=cand/$N.bed; [[ -s "$B" ]] || continue
    gawk -F'\t' -v OFS='\t' '{ print $1, $2, $3, $4, 0, $6 }' "$B" > "cand/ends/$N.b6"
    bedtools flank -i "cand/ends/$N.b6" -g "$G.fai" -l $FL -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null \
        | seqkit seq -w 0 > "cand/ends/$N.f5.fa"
    bedtools flank -i "cand/ends/$N.b6" -g "$G.fai" -l 0 -r $FL -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null \
        | seqkit seq -w 0 > "cand/ends/$N.f3.fa"
    NC=$(grep -c '>' "cand/ends/$N.f5.fa" || true)
    RES=""
    for S in 5 3; do
        # the junction is at the END of a 5' flank (position = its length) and at the START of a 3' flank
        ssearch36 -m 8 -E 1e-3 -Z 1000 -z 11 -T "$T" cons.masked.fa "cand/ends/$N.f$S.fa" 2> /dev/null \
        | gawk -F'\t' -v S=$S -v FL=$FL -v NEAR=$NEAR -v ML=$MINLEN -v NC=$NC -v OFS='\t' '
            { lo = ($9 < $10 ? $9 : $10); hi = ($9 < $10 ? $10 : $9); if (hi - lo + 1 < ML) next
              near = (S == 5) ? (FL - hi <= NEAR) : (lo <= NEAR + 1)
              if (!near) next
              if (!($2 in be) || $12 > be[$2]) { be[$2] = $12; bu[$2] = $1 } }
            END { for (k in bu) c[bu[k]]++; n = 0; for (k in bu) n++
                  t = ""; tn = 0; s2 = ""; sn = 0
                  for (u in c) { if (c[u] > tn) { s2 = t; sn = tn; t = u; tn = c[u] } else if (c[u] > sn) { s2 = u; sn = c[u] } }
                  printf "%.2f\t%s\t%s\n", (NC ? n / NC : 0), (t == "" ? "-" : t), (s2 == "" ? "-" : s2 ":" sn) }' > "cand/ends/$N.s$S"
    done
    paste <(printf "%s\n" "$N") "cand/ends/$N.s5" "cand/ends/$N.s3" \
    | gawk -F'\t' -v OFS='\t' -v F=$OPENFRAC '{ o = ($2 >= F ? "5" : "") ($5 >= F ? "3" : ""); print $1, $2, $3, $4, $5, $6, $7, (o == "" ? "-" : o) }' >> ends.tsv
done
echo "fs6b: $(($(wc -l < ends.tsv) - 1)) candidates checked, $(gawk -F'\t' 'NR > 1 && $8 != "-"' ends.tsv | wc -l) with an open end -> $OUT/ends.tsv" >&2
