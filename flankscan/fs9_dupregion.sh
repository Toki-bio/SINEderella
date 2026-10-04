#!/usr/bin/env bash
# fs9_dupregion.sh GENOME.fa COPIES.bed OUT_DIR [THREADS=8]
#
# Stage 9b of flankscan: is the NEIGHBOURHOOD of a copy duplicated elsewhere in the genome (class R, "dup-region")?
# A copy inserted independently into a segmentally duplicated region has unique insertion sites, so it is not a twin
# (fs9_twins.sh), but its flanks are not unique either. Each flank (W2 bp outward from the junction, in element
# orientation) is mapped to the genome; a hit counts when it is >= MINAL bp long, >= IDR % identical (1 - NM / alignment
# length), on a place at least ARRKB away from the copy (a hit closer than that is the copy itself or its array
# neighbours), and made of unmasked sequence: lower-case (soft-masked) bases in the query are cut away and only
# stretches of >= MINAL bp are mapped, so repeats the genome already knows are never evidence.
#
# In : GENOME.fa (soft-masked if possible), COPIES.bed (bed6)
# Out: OUT/dupregion.tsv  copy side nhits best_identity best_length best_target     (only copies with a hit)
#      OUT/dupregion.copies  one copy name per line
# Env: W2=200 MINAL=100 IDR=85 ARRKB=2000; needs bedtools, minimap2, gawk
set -euo pipefail
G=${1:?GENOME.fa}; B=${2:?COPIES.bed}; OUT=${3:?OUT_DIR}; T=${4:-8}
W2=${W2:-200}; MINAL=${MINAL:-100}; IDR=${IDR:-85}; ARRKB=${ARRKB:-2000}
mkdir -p "$OUT"; G=$(readlink -f "$G"); B=$(readlink -f "$B"); cd "$OUT"
[[ -s "$G.fai" ]] || samtools faidx "$G"
export LC_ALL=C
bedtools flank -i "$B" -g "$G.fai" -l $W2 -r 0 -s | gawk -F'\t' -v OFS='\t' '{ $4 = $4 "|5"; print }' > fl5.bed
bedtools flank -i "$B" -g "$G.fai" -l 0 -r $W2 -s | gawk -F'\t' -v OFS='\t' '{ $4 = $4 "|3"; print }' > fl3.bed
cat fl5.bed fl3.bed > fl.bed
bedtools getfasta -fi "$G" -bed fl.bed -s -nameOnly 2> /dev/null | seqkit seq -w 0 | sed '/^>/s/([+-])$//' > fl.fa
# unmasked stretches only (upper case runs of >= MINAL bp); the stretch start is kept in the name
gawk -v M=$MINAL '/^>/ { n = substr($0, 2); next } { s = $0; i = 1; L = length(s)
    while (i <= L) { c = substr(s, i, 1); if (c ~ /[ACGT]/) { j = i; while (j <= L && substr(s, j, 1) ~ /[ACGT]/) j++; if (j - i >= M) print ">" n "@" i "\n" substr(s, i, j - i); i = j } else i++ } }' fl.fa > q.fa
[[ -s q.fa ]] || { : > dupregion.tsv; : > dupregion.copies; echo "fs9b: no unmasked flank stretches"; exit 0; }
minimap2 -x asm20 -k12 -w4 -N 10 -p 0.2 --secondary=yes -t "$T" "$G" q.fa 2> /dev/null > q.paf || true
# own coordinates of every flank (bed) to exclude self and neighbours
gawk -F'\t' -v OFS='\t' -v M=$MINAL -v IDR=$IDR -v KB=$ARRKB '
    FILENAME == ARGV[1] { fc[$4] = $1; fs[$4] = $2; fe[$4] = $3; next }
    { split($1, a, "@"); f = a[1]; split(f, nm, "|"); copy = nm[1]; side = nm[2]
      al = $11; if (al < M) next; nmm = 0; for (t = 13; t <= NF; t++) if ($t ~ /^NM:i:/) nmm = substr($t, 6); id = 100 * (al - nmm) / al; if (id < IDR) next
      if ($6 == fc[f] && $8 < fe[f] + KB && $9 > fs[f] - KB) next          # the copy itself or a neighbour within ARRKB
      n[copy SUBSEP side]++
      if (id * al > bi[copy SUBSEP side]) { bi[copy SUBSEP side] = id * al; bid[copy SUBSEP side] = id; bal[copy SUBSEP side] = al; bt[copy SUBSEP side] = $6 ":" $8 "-" $9 } }
    END { for (k in n) { split(k, x, SUBSEP); printf "%s\t%s\t%d\t%.1f\t%d\t%s\n", x[1], x[2], n[k], bid[k], bal[k], bt[k] } }' fl.bed q.paf | sort > dupregion.tsv
cut -f1 dupregion.tsv | sort -u > dupregion.copies
echo "fs9b: $(wc -l < dupregion.copies) of $(wc -l < "$B") copies have a flank duplicated elsewhere"
