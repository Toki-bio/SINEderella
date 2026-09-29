#!/usr/bin/env bash
# fs0_partnerlib.sh OUT.fa [TAXID=9397] [LINE3=400] [TRNA_FA=~/refs/smallrna/hg38-tRNAs.fa] [RFAM_CM=~/refs/smallrna/Rfam.cm]
#
# A PARTNER library for stage 3: the pieces a SINE is commonly built from or sits next to, that a
# SINE consensus bank does not hold. A part missing from the bank shows only as spacer / flank;
# with these added it gets a name (a tRNA head, a 7SL or 5S part, a LINE 3' end).
# Every name starts with "x." - stage 3 uses them as partners only, never as the copy's main unit.
#   x.tRNA-<aa>-<anticodon>   one human tRNA per anticodon (the first listed), from TRNA_FA
#   x.7SL, x.5S               Rfam consensus (cmemit -c, U -> T) of RF00017 Metazoa_SRP and RF00001 5S_rRNA
#   x.<LINE>_3end             the last LINE3 bp of every Dfam LINE family of the clade TAXID and its
#                             ancestors (default 9397 Chiroptera: the Mammalia-wide L1M / L2 / CR1_Mam
#                             families plus bat ones); the 5' end / ORF fragments Dfam lists
#                             separately (_5end, _orf1, _orf2) are skipped - a SINE sits by a 3' end
# Needs: conda env rnatools (cmfetch, cmemit), curl (dfam.org API). Writes OUT.fa and OUT.fa.README.
set -euo pipefail
OUT=${1:?OUT.fa}; TAX=${2:-9397}; LINE3=${3:-400}
TRNA=${4:-$HOME/refs/smallrna/hg38-tRNAs.fa}; CM=${5:-$HOME/refs/smallrna/Rfam.cm}
TMP=$(mktemp -d "$HOME/tmp/fs0.XXXX"); trap 'rm -rf "$TMP"' EXIT

# tRNAs: header ">Homo_sapiens_tRNA-Ala-AGC-1-1 (...)" -> the first gene of each anticodon; Undet / Sup skipped
seqkit seq -w 0 "$TRNA" | gawk '
    /^>/ { keep = 0
           if (match($1, /tRNA-([A-Za-z]+)-([ACGTN]{3})-[0-9]+-[0-9]+$/, m) && m[1] != "Undet" && m[1] != "Sup") {
               k = m[1] "-" m[2]; if (!(k in seen)) { seen[k] = 1; keep = 1; print ">x.tRNA-" k } }
           next }
    keep { print toupper($0) }' > "$TMP/trna.fa"

# 7SL and 5S: the Rfam model consensus (capital = conserved, lower case kept as bases)
RNA="$(conda info --base)/envs/rnatools/bin"                     # Infernal lives in env rnatools
for pair in RF00017:7SL RF00001:5S; do
    acc=${pair%%:*}; nm=${pair#*:}
    "$RNA/cmfetch" "$CM" "$acc" > "$TMP/$nm.cm"
    "$RNA/cmemit" -c "$TMP/$nm.cm" | seqkit seq -w 0 | gawk -v n="$nm" '/^>/ { print ">x." n; next } { s = toupper($0); gsub(/U/, "T", s); print s }'
done > "$TMP/rna.fa"

# Dfam LINEs of the clade and its ancestors; the 3' LINE3 bp of each. The API answers
# {"content_type": ..., "body": "<fasta, newlines escaped as \n>"}; headers ">DF000000008.4 L1M5_orf2"
curl -s -m 120 "https://dfam.org/api/families?clade=$TAX&clade_relatives=ancestors&type=LINE&format=fasta&limit=2000" \
| sed -e 's/^.*"body":"//' -e 's/"}[[:space:]]*$//' -e 's/\\n/\n/g' > "$TMP/dfam_line.fa"
seqkit seq -w 0 "$TMP/dfam_line.fa" | gawk -v L=$LINE3 '
    /^>/ { n = $2; skip = (n ~ /_(5end|orf1|orf2)$/); next }
    !skip && n != "" { s = toupper($0); if (length(s) > L) s = substr(s, length(s) - L + 1)
                       sub(/_3end$/, "", n); print ">x." n "_3end"; print s }' > "$TMP/line.fa"
[[ -s "$TMP/line.fa" ]] || { echo "fs0: no LINE from Dfam (network?) - stop" >&2; exit 1; }

cat "$TMP/trna.fa" "$TMP/rna.fa" "$TMP/line.fa" > "$OUT"
{ echo "fs0_partnerlib.sh $(date +%F): taxid $TAX, LINE 3' $LINE3 bp"
  echo "tRNA $(grep -c '^>' "$TMP/trna.fa") from $TRNA"
  echo "Rfam $(grep -c '^>' "$TMP/rna.fa") (RF00017 7SL, RF00001 5S) from $CM"
  echo "Dfam LINE $(grep -c '^>' "$TMP/line.fa") of $(grep -c '^>' "$TMP/dfam_line.fa") families (clade $TAX + ancestors; _5end/_orf skipped)"
  grep '^>' "$TMP/line.fa" | tr '\n' ' '; echo; } > "$OUT.README"
cat "$OUT.README" >&2
