#!/usr/bin/env bash
# fs9_twins.sh GENOME.fa COPIES.bed OUT_DIR
#
# Stage 9 of flankscan: flank uniqueness over ALL copies of one family. Independent insertions have unrelated flanks;
# copies that were multiplied together with their neighbourhood (segmental duplication, array unit, carried by another
# element) have flanks that are similar from the junction outward. A pair of copies is a TWIN pair when the flank
# sequences next to the element (the proximal W bp, read outward from the junction, in the copy's own orientation:
# 5' with 5', 3' with 3', so inverted duplicates are found too) are similar from the junction on, colinear, over
# >= MINLEN bp.
#
#   tier 1  identity >= ID1 (85 %)      twin
#   tier 2  identity >= ID2 (70 %)      possible old duplicate (reported separately)
#   array   twin pairs of copies on one contig whose neighbours lie less than ARRKB (2 kb) apart: reported as array
#   untestable  a flank shorter than MINFL (contig end): never counted as unique
#
# Method (all awk, readable on purpose):
#  1) flanks with bedtools; masked (lower case) bases become N, and so does every base covered by a 20-mer that occurs
#     more than GCAP times in the genome: repeats the genome already knows (LINE fragments, simple repeats) are never
#     evidence of a shared flank
#  2) seeds of the first IW bp, each with its offset from the junction. Pass 1: contiguous 12-mers (tier 1). Pass 2: the
#     spaced seed 111010010100110111 (weight 11, finds old duplicates). A seed that occurs in more than CAP flanks at
#     the same offset is dropped (hot spots)
#  3) candidate pair = >= MK shared seeds on one diagonal (band +-BAND bp) with the first seed within START bp of the
#     junction in BOTH copies (MK1 = 6 in pass 1, MK2 = 3 in pass 2)
#  4) confirmation: score the diagonal (and BAND diagonals around it, to cross small indels): match +1, mismatch -1, masked bases
#     skipped, best segment, start <= START bp from the junction, length >= MINLEN, identity >= ID1 (tier 1) or ID2 (tier 2)
#
# In : GENOME.fa (+ .fai made if missing; soft-masked if possible), COPIES.bed (bed6, strand = element orientation)
# Out: OUT/flank_groups.tsv  copy_a copy_b side identity length tier class     (one row per pair and side)
#      OUT/copy_status.tsv   copy status (twin1 twin2 array untestable masked unique) n_partners group
#                            (masked = fewer than MINFL testable bases in a flank after the repeat filter: nothing to compare)
#      OUT/summary.txt       counts and shares
# Env: MK1=6 MK2=3 (shared seeds needed in pass 1 / 2) W=100 IW=80 MINFL=50 MINLEN=50 START=20 BAND=3 CAP=20 GCAP=20 ID1=85 ID2=70 ARRKB=2000
#      KCOUNT=<file "20-mer TAB count" of the genome, canonical, e.g. jellyfish dump -c>; counted here up to 300 Mb
set -euo pipefail
G=${1:?GENOME.fa}; B=${2:?COPIES.bed}; OUT=${3:?OUT_DIR}
W=${W:-100}; IW=${IW:-80}; MINFL=${MINFL:-50}; MINLEN=${MINLEN:-50}; START=${START:-20}; BAND=${BAND:-3}
CAP=${CAP:-20}; GCAP=${GCAP:-20}; ID1=${ID1:-85}; ID2=${ID2:-70}; ARRKB=${ARRKB:-2000}; MK1=${MK1:-6}; MK2=${MK2:-3}
mkdir -p "$OUT"; G=$(readlink -f "$G"); B=$(readlink -f "$B"); cd "$OUT"
[[ -s "$G.fai" ]] || samtools faidx "$G"
export LC_ALL=C
RC='function rc(s,   r, i, c) { r = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1); r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N") } return r }'

# 1) flanks, outward from the junction (5' flank reversed). SOFTMASK=1 turns soft-masked (lower-case) bases into N as the first
#    version did; the default keeps them as bases, because an assembly masked by RepeatMasker is lower case over most of its
#    length (rsi: 79 %), which made the flanks of copies that sit in duplicated repeats untestable and reported them as "unique"
#    (rsi MEG-RS 2026-10-05: 11 copies with 96-100 % identical flanks on other contigs, all-N in both flanks). The k-mer filter
#    below (GCAP) still masks the hot spots: a 20-mer seen more than GCAP times in the genome is never evidence of a shared flank.
SOFTMASK=${SOFTMASK:-0}
bedtools flank -i "$B" -g "$G.fai" -l $W -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f5.fa
bedtools flank -i "$B" -g "$G.fai" -l 0 -r $W -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f3.fa
for S in 5 3; do
    gawk -v S=$S -v SM=$SOFTMASK 'function rev(s,   r, i) { r = ""; for (i = length(s); i >= 1; i--) r = r substr(s, i, 1); return r }
        /^>/ { n = substr($0, 2); sub(/\([+-]\)$/, "", n); next }
        { x = $0; if (SM == 1) gsub(/[acgtn]/, "N", x); x = toupper(x); gsub(/[^ACGT]/, "N", x); print n "\t" S "\t" (S == 5 ? rev(x) : x) }' "f$S.fa"
done > seqs.tsv
# untestable: a copy with a missing or short flank (contig end) on either side
gawk -F'\t' -v M=$MINFL 'FILENAME == ARGV[1] { c[$4]; next } { if (length($3) >= M) ok[$1 SUBSEP $2] = 1 }
    END { for (n in c) if (!((n SUBSEP 5) in ok) || !((n SUBSEP 3) in ok)) print n }' "$B" seqs.tsv | sort > untestable.txt

# repeat filter: 20-mers seen more than GCAP times in the genome mask the flank bases they cover
GK=20
if [[ -n "${KCOUNT:-}" && -s "${KCOUNT:-}" ]]; then
    # jellyfish dump -c writes "KMER COUNT" separated by a SPACE; the first version split on tabs only, so $2 was empty,
    # rep20.tsv came out empty and the repeat filter was silently off in every run with a KCOUNT (rsi, 2026-10-01)
    gawk -v C=$GCAP 'NF >= 2 && $2 + 0 > C { print toupper($1) }' "$KCOUNT" > rep20.tsv
    [[ -s rep20.tsv ]] || { echo "fs9: KCOUNT $KCOUNT gave no 20-mers above $GCAP: wrong file or format (expect 'KMER COUNT' lines)" >&2; exit 1; }
else
    SZ=$(gawk '{ s += $2 } END { print s + 0 }' "$G.fai")
    (( SZ <= 300000000 )) || { echo "fs9: genome > 300 Mb: give KCOUNT (jellyfish dump -c, k=20, canonical)" >&2; exit 1; }
    gawk -v K=$GK -v C=$GCAP "$RC"'
        function scan(   n, i, m, r) { n = length(seq); for (i = 1; i + K - 1 <= n; i++) { m = substr(seq, i, K); if (m ~ /N/) continue; r = rc(m); if (r < m) m = r; c[m]++ } }
        /^>/ { if (seq != "") scan(); seq = ""; next } { seq = seq toupper($0) }
        END { scan(); for (m in c) if (c[m] > C) print m }' "$G" > rep20.tsv
fi
gawk -F'\t' -v K=$GK "$RC"'
    FILENAME == ARGV[1] { rep[$1]; next }
    { s = $3; n = length(s); for (i = 1; i <= n; i++) mk[i] = 0
      for (i = 1; i + K - 1 <= n; i++) { m = substr(s, i, K); if (m ~ /N/) continue; r = rc(m); if (r < m) m = r; if (m in rep) for (j = i; j < i + K; j++) mk[j] = 1 }
      o = ""; for (i = 1; i <= n; i++) o = o (mk[i] ? "N" : substr(s, i, 1)); print $1 "\t" $2 "\t" o }' rep20.tsv seqs.tsv > seqs.masked.tsv
mv seqs.masked.tsv seqs.tsv
# masked: a copy with fewer than MINFL testable (non-N) bases in a flank after masking has nothing to compare on that side; it is
# reported as "masked", never as "unique" (the first version called such copies unique)
gawk -F'\t' -v M=$MINFL '{ s = $3; nn = gsub(/N/, "", s); if (length(s) < M) print $1 }' seqs.tsv | sort -u > masked.txt

# 2-4) one pass: pass_run PATTERN MINKEYS IDMIN TAG
pass_run() {
    local PAT=$1 MK=$2 IDM=$3 TAG=$4
    # seeds: the pattern's 1-positions of each window, with the offset of the window start from the junction; windows with N
    # or with at most 2 distinct bases are skipped
    gawk -F'\t' -v PAT="$PAT" -v IW=$IW 'BEGIN { L = length(PAT); for (j = 1; j <= L; j++) if (substr(PAT, j, 1) == "1") pos[++np] = j }
        { s = substr($3, 1, IW + L - 1); n = length(s)
          for (i = 1; i + L - 1 <= n && i <= IW; i++) { w = substr(s, i, L); if (w ~ /N/) continue
              m = ""; for (j = 1; j <= np; j++) m = m substr(w, pos[j], 1)
              nd = 0; for (j = 1; j <= np; j++) { c = substr(m, j, 1); if (!(c in seen)) { seen[c] = 1; nd++ } } delete seen
              if (nd <= 2) continue
              print m "\t" i "\t" $1 "\t" $2 } }' seqs.tsv | sort -t$'\t' -k1,1 -T . > "keys.$TAG.tsv"
    # buckets of one seed (same orientation: both flanks are in element orientation): pairs, with the delta of the offsets
    gawk -F'\t' -v CAP=$CAP -v OFS='\t' '
        function flush(   i, j) { if (n > 1) { split("", cnt); for (i = 1; i <= n; i++) cnt[of[i]]++
              for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++) if (sd[i] == sd[j] && cp[i] != cp[j] && cnt[of[i]] <= CAP && cnt[of[j]] <= CAP) {
                  if (cp[i] < cp[j]) print sd[i], cp[i], cp[j], of[i] - of[j], of[i], of[j]; else print sd[i], cp[j], cp[i], of[j] - of[i], of[j], of[i] } } n = 0 }
        { if ($1 != last) { flush(); last = $1 } n++; of[n] = $2; cp[n] = $3; sd[n] = $4 } END { flush() }' "keys.$TAG.tsv" > "hits.$TAG.tsv"
    # candidate pairs: >= MK seeds within one band of delta, first seed within START bp of the junction in both copies.
    # Streaming: the hits are sorted (stable) by side and copy pair, so that one pair's seeds are consecutive and the awk holds one pair at
    # a time. The first version kept every pair of the whole family in memory and was killed by the system on 538 486 copies (Sicista
    # B1) and 1 185 121 copies (DIP), 2026-10-06. Same output: the order of one pair's seeds is kept (sort -s), and the final sort -u
    # does not depend on the order of the pairs.
    sort -s -t$'\t' -k1,1 -k2,2 -k3,3 -T . "hits.$TAG.tsv" | gawk -F'\t' -v OFS='\t' -v MK=$MK -v BAND=$BAND -v START=$START -v SLACK=25 '
        function flush(   nd, a, d, s, best, bd, f) {
            if (cur == "") return
            if (ma > START + SLACK || mb > START + SLACK) return     # masked bases at the junction push the first seed outward; the confirmation applies START
            nd = split(dl, D, " "); best = 0; bd = 0
            for (a = 1; a <= nd; a++) { s = 0; for (d = D[a] - BAND; d <= D[a] + BAND; d++) s += c[d] + 0; if (s > best) { best = s; bd = D[a] } }
            if (best >= MK) { split(cur, f, SUBSEP); print f[1], f[2], f[3], bd, best } }
        { p = $1 SUBSEP $2 SUBSEP $3
          if (p != cur) { flush(); cur = p; delete c; dl = ""; ma = 1e9; mb = 1e9 }
          c[$4]++; dl = dl " " $4; if ($5 < ma) ma = $5; if ($6 < mb) mb = $6 }
        END { flush() }' | sort -u > "cand.$TAG.tsv"
    # confirmation: score the diagonal and BAND diagonals around it (offset in copy a = offset in copy b + delta)
    gawk -F'\t' -v OFS='\t' -v IDM=$IDM -v ML=$MINLEN -v START=$START -v BAND=$BAND -v TAG=$TAG '
        FILENAME == ARGV[1] { sq[$1 SUBSEP $2] = $3; next }
        { A = sq[$2 SUBSEP $1]; Bq = sq[$3 SUBSEP $1]; la = length(A); lb = length(Bq); bestid = 0; bestl = 0
          for (dd = $4 - BAND; dd <= $4 + BAND; dd++) {
              sc = 0; st = 0; mt = 0; ln = 0; bs = 0; bm = 0; bl = 0
              for (j = 1; j <= lb; j++) { i = j + dd; if (i < 1 || i > la) continue
                  x = substr(A, i, 1); y = substr(Bq, j, 1)
                  if (sc <= 0) { sc = 0; st = (i <= START && j <= START) ? 1 : 0; mt = 0; ln = 0 }
                  if (!st) continue
                  if (x == "N" || y == "N") continue                    # a masked base is neither a match nor a mismatch
                  ln++; if (x == y) { sc++; mt++ } else sc--
                  if (sc > bs) { bs = sc; bm = mt; bl = ln } }
              if (bl >= ML && 100 * bm / bl >= IDM && bl > bestl) { bestl = bl; bestid = 100 * bm / bl } }
          if (bestl >= ML) printf "%s\t%s\t%s\t%.1f\t%d\t%s\n", $2, $3, $1, bestid, bestl, TAG }' seqs.tsv "cand.$TAG.tsv" > "conf.$TAG.tsv"
}

pass_run 111111111111 $MK1 $ID1 p1
pass_run 111010010100110111 $MK2 $ID2 p2
# a pair found in both passes keeps its pass 1 row; the tier follows the identity
gawk -F'\t' -v OFS='\t' -v ID1=$ID1 '{ key = $1 SUBSEP $2 SUBSEP $3; if (key in done) next; done[key] = 1
      print $1, $2, $3, $4, $5, ($4 >= ID1 ? "T1" : "T2"), "twin" }' conf.p1.tsv conf.p2.tsv | sort -k1,1 -k2,2 -k3,3 > pairs.body
# arrays: copies on one contig closer than ARRKB to a twin partner form one array; every pair inside it is "array", not "twin"
gawk -F'\t' -v OFS='\t' -v KB=$ARRKB 'FILENAME == ARGV[1] { ctg[$4] = $1; pos[$4] = $2; next }
    function root(x) { while (up[x] != x) x = up[x]; return x }
    FILENAME == ARGV[2] { a[++n] = $0; if (!($1 in up)) up[$1] = $1; if (!($2 in up)) up[$2] = $2
        d = pos[$1] - pos[$2]; if (d < 0) d = -d
        if (ctg[$1] == ctg[$2] && d < KB) { ra = root($1); rb = root($2); if (ra != rb) up[ra] = rb; arr[$1]; arr[$2] }; next }
    END { for (i = 1; i <= n; i++) { split(a[i], f, "\t"); cl = (f[1] in arr && f[2] in arr && root(f[1]) == root(f[2])) ? "array" : "twin"; print f[1], f[2], f[3], f[4], f[5], f[6], cl } }' "$B" pairs.body > flank_groups.body
{ printf "copy_a\tcopy_b\tside\tidentity\tlength\ttier\tclass\n"; cat flank_groups.body; } > flank_groups.tsv

# per-copy status; precedence twin1 > twin2 > array > untestable > masked > unique; groups = twin pairs joined
gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { allc[$4]; next } FILENAME == ARGV[2] { unt[$1]; next } FILENAME == ARGV[3] { msk[$1]; next }
    function root(x) { while (up[x] != x) x = up[x]; return x }
    { a = $1; b = $2; if (!(a in up)) up[a] = a; if (!(b in up)) up[b] = b
      st = ($7 == "array") ? "array" : ($6 == "T1" ? "twin1" : "twin2"); rank = (st == "twin1") ? 4 : (st == "twin2") ? 3 : 2
      for (v = 1; v <= 2; v++) { c = (v == 1) ? a : b; if (rank > R[c]) { R[c] = rank; S[c] = st }; np[c]++ }
      if ($7 == "twin") { ra = root(a); rb = root(b); if (ra != rb) up[ra] = rb } }
    END { for (c in allc) { st = (c in S) ? S[c] : ((c in unt) ? "untestable" : ((c in msk) ? "masked" : "unique"))
              print c, st, np[c] + 0, ((c in S) && S[c] ~ /twin/ ? "G" root(c) : "-") } }' "$B" untestable.txt masked.txt flank_groups.body | sort > copy_status.tsv
gawk -F'\t' '{ n[$2]++; t++ } END { printf "copies\t%d\n", t; for (k in n) printf "%s\t%d\t%.1f%%\n", k, n[k], 100 * n[k] / t }' copy_status.tsv | sort > summary.txt
cat summary.txt
