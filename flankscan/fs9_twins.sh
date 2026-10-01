#!/usr/bin/env bash
# fs9_twins.sh GENOME.fa COPIES.bed OUT_DIR [THREADS=8]
#
# Stage 9 of flankscan: flank uniqueness over ALL copies of one family. Independent insertions have unrelated flanks;
# copies that were multiplied together with their neighbourhood (segmental duplication, array unit, carried by another
# element) have flanks that are identical from the junction outward. A pair of copies is a TWIN pair when the flank
# sequences next to the element (the proximal W bp, read outward from the junction, in the copy's own orientation:
# 5' with 5', 3' with 3', so inverted duplicates are found too) are similar from the junction on, colinear, over
# >= MINLEN bp.
#
#   tier 1  identity >= ID1 (85 %)      twin
#   tier 2  identity >= ID2 (70 %)      possible old duplicate (reported separately)
#   array   twin pair on one contig less than ARRKB (2 kb) apart: reported as array, not as twin
#   untestable  a flank shorter than MINFL (contig end): never counted as unique
#
# Method (readable on purpose, all awk): (1) flanks with bedtools; masked (lower case) bases become N, so repeats the
# genome already knows (LINE fragments, simple repeats) are never evidence; (2) every K-mer of the first IW bp with its
# offset from the junction; k-mers that occur more than GCAP times in the genome or in more than CAP flanks are dropped
# (hot spots and repeats); (3) candidate pairs = >= MINKEYS shared k-mers on one diagonal (band +-BAND bp), first one
# within START bp of the junction in BOTH copies; (4) each candidate is confirmed by scoring the diagonal (match +1,
# mismatch or N -1, best segment, start <= START bp from the junction, a few diagonals around the seed to cross indels).
# Pass 1 uses K=12 (tier 1); pass 2 uses K=9 and a lower identity floor and only reports pairs pass 1 did not (tier 2).
#
# In : GENOME.fa (+ .fai made if missing; soft-masked if possible), COPIES.bed (bed6, strand = element orientation)
# Out: OUT/flank_groups.tsv  copy_a copy_b side identity length tier class     (one row per pair and side)
#      OUT/copy_status.tsv   copy status (twin1 twin2 array untestable unique) n_partners group
#      OUT/summary.txt       counts and shares
# Env: W=100 IW=80 MINFL=50 MINLEN=50 START=20 BAND=3 CAP=20 GCAP=20 ID1=85 ID2=70 ARRKB=2000 KCOUNT=<k-mer count file
#      "kmer TAB count" for the genome, e.g. from jellyfish dump; counted here for genomes up to 300 Mb>
set -euo pipefail
G=${1:?GENOME.fa}; B=${2:?COPIES.bed}; OUT=${3:?OUT_DIR}; T=${4:-8}
W=${W:-100}; IW=${IW:-80}; MINFL=${MINFL:-50}; MINLEN=${MINLEN:-50}; START=${START:-20}; BAND=${BAND:-3}
CAP=${CAP:-20}; GCAP=${GCAP:-20}; ID1=${ID1:-85}; ID2=${ID2:-70}; ARRKB=${ARRKB:-2000}
mkdir -p "$OUT"; cd "$OUT"; G=$(readlink -f "$G"); B=$(readlink -f "$B")
[[ -s "$G.fai" ]] || samtools faidx "$G"
export LC_ALL=C

# 1) flanks, outward from the junction (5' flank reversed), masked bases -> N
bedtools flank -i "$B" -g "$G.fai" -l $W -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f5.fa
bedtools flank -i "$B" -g "$G.fai" -l 0 -r $W -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f3.fa
for S in 5 3; do
    gawk -v S=$S 'function rev(s,   r, i) { r = ""; for (i = length(s); i >= 1; i--) r = r substr(s, i, 1); return r }
        /^>/ { n = substr($0, 2); sub(/\([+-]\)$/, "", n); next }
        { x = $0; gsub(/[acgtn]/, "N", x); x = toupper(x); print n "\t" S "\t" (S == 5 ? rev(x) : x) }' "f$S.fa"
done > seqs.tsv
# untestable: a side that is missing or has fewer than MINFL unmasked-or-masked bases
gawk -F'\t' -v M=$MINFL 'FILENAME == ARGV[1] { c[$4]; next } { if (length($3) >= M) ok[$1 SUBSEP $2] = 1; seen[$1] }
    END { for (n in c) if (!((n SUBSEP 5) in ok) || !((n SUBSEP 3) in ok)) print n }' "$B" seqs.tsv | sort > untestable.txt

# genome k-mer counts (for the repeat filter), per K
kcount() {   # kcount K -> file "kmer TAB count" (only k-mers seen > GCAP times)
    local K=$1
    if [[ -n "${KCOUNT:-}" && -s "${KCOUNT:-}" ]]; then gawk -F'\t' -v C=$GCAP '$2 > C' "$KCOUNT"; return; fi
    local SZ; SZ=$(gawk '{ s += $2 } END { print s + 0 }' "$G.fai")
    (( SZ <= 300000000 )) || { echo "fs9: genome > 300 Mb: give KCOUNT (jellyfish dump of k=$K)" >&2; exit 1; }
    gawk -v K=$K -v C=$GCAP 'function rc(s,   r, i, c) { r = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1); r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N") } return r }
        /^>/ { next } { s = toupper($0); L = length(s); seq = seq s }
        END { n = length(seq); for (i = 1; i + K - 1 <= n; i++) { m = substr(seq, i, K); if (m ~ /N/) continue; r = rc(m); if (r < m) m = r; c[m]++ }
              for (m in c) if (c[m] > C) print m "\t" c[m] }' "$G"
}

# 2-4) one pass: candidate pairs from shared keys, then confirmation. pass_run K MINKEYS IDMIN TAG
pass_run() {
    local K=$1 MK=$2 IDM=$3 TAG=$4
    kcount $K > "rep.$TAG.tsv"
    # keys: k-mer, offset from the junction, copy, side; low-complexity and masked k-mers skipped; repeat k-mers skipped
    gawk -F'\t' -v K=$K -v IW=$IW 'function rc(s,   r, i, c) { r = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1); r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N") } return r }
        FILENAME == ARGV[1] { rep[$1]; next }
        { s = substr($3, 1, IW + K - 1); n = length(s)
          for (i = 1; i + K - 1 <= n && i <= IW; i++) { m = substr(s, i, K); if (m ~ /N/) continue
              split(m, ch, ""); nd = 0; for (j = 1; j <= K; j++) if (!(ch[j] in seen)) { seen[ch[j]] = 1; nd++ } delete seen; if (nd <= 2) continue
              cm = m; r = rc(m); if (r < cm) cm = r
              if (cm in rep) continue
              print m "\t" i "\t" $1 "\t" $2 } }' "rep.$TAG.tsv" seqs.tsv | sort -t$'\t' -k1,1 -T . > "keys.$TAG.tsv"
    # buckets of equal k-mer (same orientation: both flanks are in element orientation, so no reverse complement)
    gawk -F'\t' -v CAP=$CAP -v OFS='\t' 'function flush(   i, j) { if (n > 1 && n <= CAP) for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++)
            if (sd[i] == sd[j] && cp[i] != cp[j]) { if (cp[i] < cp[j]) print sd[i], cp[i], cp[j], of[i] - of[j], of[i], of[j]; else print sd[i], cp[j], cp[i], of[j] - of[i], of[j], of[i] } n = 0 }
        { if ($1 != last) { flush(); last = $1 } n++; of[n] = $2; cp[n] = $3; sd[n] = $4 } END { flush() }' "keys.$TAG.tsv" > "hits.$TAG.tsv"
    # candidate pairs: >= MK keys within one band of delta, first key within START bp of the junction in both copies
    gawk -F'\t' -v OFS='\t' -v MK=$MK -v BAND=$BAND -v START=$START '
        { p = $1 SUBSEP $2 SUBSEP $3; c[p, $4]++; if (!(p in seen)) { seen[p] = 1; order[++np] = p }
          dl[p] = dl[p] " " $4; if (!((p, "a") in mina) || $5 < mina[p, "a"]) mina[p, "a"] = $5; if (!((p, "b") in minb) || $6 < minb[p, "b"]) minb[p, "b"] = $6 }
        END { for (q = 1; q <= np; q++) { p = order[q]; nd = split(dl[p], D, " "); best = 0; bd = 0
              for (a = 1; a <= nd; a++) { s = 0; for (d = D[a] - BAND; d <= D[a] + BAND; d++) s += c[p, d] + 0; if (s > best) { best = s; bd = D[a] } }
              if (best >= MK) { split(p, f, SUBSEP); print f[1], f[2], f[3], bd, best } } }' "hits.$TAG.tsv" | sort -u > "cand.$TAG.tsv"
    # confirmation: score the diagonal (and a few around it): match +1, mismatch or N -1; best segment starting <= START
    gawk -F'\t' -v OFS='\t' -v IDM=$IDM -v ML=$MINLEN -v START=$START -v BAND=$BAND -v TAG=$TAG '
        FILENAME == ARGV[1] { sq[$1 SUBSEP $2] = $3; next }
        { A = sq[$2 SUBSEP $1]; Bq = sq[$3 SUBSEP $1]; la = length(A); lb = length(Bq); bestid = 0; bestl = 0
          for (dd = $4 - BAND; dd <= $4 + BAND; dd++) {           # offset in A = offset in B + dd
              sc = 0; st = 0; mt = 0; ln = 0; bs = 0; bm = 0; bl = 0
              for (j = 1; j <= lb; j++) { i = j + dd; if (i < 1 || i > la) continue
                  x = substr(A, i, 1); y = substr(Bq, j, 1)
                  if (sc <= 0) { sc = 0; st = (i <= START && j <= START) ? 1 : 0; mt = 0; ln = 0 }
                  if (!st) continue
                  ln++; if (x == y && x != "N") { sc++; mt++ } else sc--
                  if (sc > bs) { bs = sc; bm = mt; bl = ln } }
              if (bl >= ML && 100 * bm / bl >= IDM && (bl > bestl)) { bestl = bl; bestid = 100 * bm / bl } }
          if (bestl >= ML) printf "%s\t%s\t%s\t%.1f\t%d\t%s\n", $2, $3, $1, bestid, bestl, TAG }' seqs.tsv "cand.$TAG.tsv" > "conf.$TAG.tsv"
}

pass_run 12 8 $ID1 p1
pass_run 9 6 $ID2 p2
# pairs found in pass 2 only are tier 2; a pair found in both keeps its pass 1 row (tier 1 needs identity >= ID1)
gawk -F'\t' -v OFS='\t' -v ID1=$ID1 -v KB=$ARRKB 'FILENAME == ARGV[1] { ctg[$4] = $1; pos[$4] = $2; next }
    { key = $1 SUBSEP $2 SUBSEP $3; if (key in done) next; done[key] = 1
      tier = ($4 >= ID1) ? "T1" : "T2"
      d = pos[$1] - pos[$2]; if (d < 0) d = -d
      cl = (ctg[$1] == ctg[$2] && d < KB) ? "array" : "twin"
      print $1, $2, $3, $4, $5, tier, cl }' "$B" conf.p1.tsv conf.p2.tsv | sort -k1,1 -k2,2 -k3,3 > flank_groups.body
{ printf "copy_a\tcopy_b\tside\tidentity\tlength\ttier\tclass\n"; cat flank_groups.body; } > flank_groups.tsv

# per-copy status and groups (union-find over twin pairs); precedence twin1 > twin2 > array > untestable > unique
gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { allc[$4]; next } FILENAME == ARGV[2] { unt[$1]; next }
    function root(x) { while (up[x] != x) x = up[x]; return x }
    { a = $1; b = $2; if (!(a in up)) up[a] = a; if (!(b in up)) up[b] = b
      st = ($7 == "array") ? "array" : ($6 == "T1" ? "twin1" : "twin2")
      for (v = 1; v <= 2; v++) { c = (v == 1) ? a : b; rank = (st == "twin1") ? 4 : (st == "twin2") ? 3 : 2
          if (rank > R[c]) { R[c] = rank; S[c] = st }; np[c]++ }
      if ($7 == "twin") { ra = root(a); rb = root(b); if (ra != rb) up[ra] = rb } }
    END { for (c in allc) { st = (c in S) ? S[c] : ((c in unt) ? "untestable" : "unique")
              if ((c in unt) && !(c in S)) st = "untestable"
              print c, st, np[c] + 0, ((c in up) && ((c in S) && (S[c] ~ /twin/)) ? "G" root(c) : "-") } }' "$B" untestable.txt flank_groups.body | sort > copy_status.tsv
gawk -F'\t' '{ n[$2]++; t++ } END { printf "copies\t%d\n", t; for (k in n) printf "%s\t%d\t%.1f%%\n", k, n[k], 100 * n[k] / t }' copy_status.tsv | sort > summary.txt
cat summary.txt
