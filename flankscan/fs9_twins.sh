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
# Soft-masked (lower-case) flank bases and the body check (2026-10-14, flankscan/tests/real_twins/):
#   SOFTMASK=0 (default) masked bases are compared like any base; 1 = they become N (misses twins inside repeats, rsi MEG-RS);
#   both = (a) SOFTMASK=1 at the tiers above united with (b) masked bases kept at identity >= IDMASK (95 %).
#   BODYCHECK=1 (default, fs9_bodycheck.py): a flank pair stays a twin only if the two SINE BODIES are more similar than random
#   copy pairs of the family (quantile 1 - 0.05/k, k = partners of the busier copy). Hot spots (independent insertions into two
#   copies of a young repeat or another SINE) fail it: Python Squam3C 930 of 995 SOFTMASK=0 pairs, Anolis Squam1 143,029 of
#   154,360; copies duplicated with their neighbourhood pass (body identity tracks flank identity). Toy: 0 of 40 hot-spot
#   copies, 4 of 4 soft-masked duplications, 12 of 16 duplicated copies (the 78 % and 73 % pairs are at the family background).
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
# Env: SOFTMASK=both|0|1 IDMASK=95 BODYCHECK=1|0 MK1=6 MK2=3 (shared seeds needed in pass 1 / 2) W=100 IW=80 MINFL=50 MINLEN=50 START=20 BAND=3 CAP=20 GCAP=20 ID1=85 ID2=70 ARRKB=2000
#      KCOUNT=<file "20-mer TAB count" of the genome, canonical, e.g. jellyfish dump -c>; counted here up to 300 Mb
#      FS9_THREADS (default BODYCHECK_THREADS, else 1): workers of fs9_seeds.py and of the confirmation step
#      FS9_AWK=1: the repeat filter and steps 2-3 in awk, the reference implementation (it writes keys/hits files of tens of GB
#      on families of 10^5 copies); default fs9_repmask.py and fs9_seeds.py: the same rows, nothing large on disk
set -euo pipefail
G=${1:?GENOME.fa}; B=${2:?COPIES.bed}; OUT=${3:?OUT_DIR}
W=${W:-100}; IW=${IW:-80}; MINFL=${MINFL:-50}; MINLEN=${MINLEN:-50}; START=${START:-20}; BAND=${BAND:-3}
CAP=${CAP:-20}; GCAP=${GCAP:-20}; ID1=${ID1:-85}; ID2=${ID2:-70}; ARRKB=${ARRKB:-2000}; MK1=${MK1:-6}; MK2=${MK2:-3}
SOFTMASK=${SOFTMASK:-0}; IDMASK=${IDMASK:-95}
FS9_AWK=${FS9_AWK:-0}; FS9_THREADS=${FS9_THREADS:-${BODYCHECK_THREADS:-1}}
HERE9="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
case $SOFTMASK in both|0|1) ;; *) echo "fs9: SOFTMASK must be both, 0 or 1" >&2; exit 1 ;; esac
mkdir -p "$OUT"; G=$(readlink -f "$G"); B=$(readlink -f "$B"); cd "$OUT"
[[ -s "$G.fai" ]] || samtools faidx "$G"
export LC_ALL=C
RC='function rc(s,   r, i, c) { r = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1); r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N") } return r }'

# 1) flanks, outward from the junction (5' flank reversed), soft-masked bases kept in lower case (seqs.raw.tsv); make_seqs SM writes
#    seqs.tsv in upper case with masked bases as N (SM=1) or as bases (SM=0). History: the first version always used SM=1; SM=0 alone
#    (2026-10-05) found the rsi MEG-RS twins (11 copies, 96-100 % identical flanks, all-N in both flanks under SM=1) but called
#    hot-spot copies twins; SOFTMASK=both (header) keeps both properties.
bedtools flank -i "$B" -g "$G.fai" -l $W -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f5.fa
bedtools flank -i "$B" -g "$G.fai" -l 0 -r $W -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null > f3.fa
for S in 5 3; do
    gawk -v S=$S 'function rev(s,   r, i) { r = ""; for (i = length(s); i >= 1; i--) r = r substr(s, i, 1); return r }
        /^>/ { n = substr($0, 2); sub(/\([+-]\)$/, "", n); next }
        { x = $0; gsub(/[^ACGTacgt]/, "N", x); print n "\t" S "\t" (S == 5 ? rev(x) : x) }' "f$S.fa"
done > seqs.raw.tsv
make_seqs() {   # SM: 1 = soft-masked bases become N; then upper case
    gawk -F'\t' -v OFS='\t' -v SM=$1 '{ x = $3; if (SM == 1) gsub(/[acgt]/, "N", x); print $1, $2, toupper(x) }' seqs.raw.tsv > seqs.tsv
}
make_seqs 0
# untestable: a copy with a missing or short flank (contig end) on either side
gawk -F'\t' -v M=$MINFL 'FILENAME == ARGV[1] { c[$4]; next } { if (length($3) >= M) ok[$1 SUBSEP $2] = 1 }
    END { for (n in c) if (!((n SUBSEP 5) in ok) || !((n SUBSEP 3) in ok)) print n }' "$B" seqs.tsv | sort > untestable.txt

# repeat filter: 20-mers seen more than GCAP times in the genome mask the flank bases they cover
GK=20
if [[ -n "${KCOUNT:-}" && -s "${KCOUNT:-}" ]]; then
    # jellyfish dump -c writes "KMER COUNT" separated by a SPACE; the first version split on tabs only, so $2 was empty,
    # rep20.tsv came out empty and the repeat filter was silently off in every run with a KCOUNT (rsi, 2026-10-01)
    # only the 20-mers that occur in a flank are kept: a genome's repeated 20-mers can be hundreds of millions (anoCar2)
    if [[ $FS9_AWK == 1 ]]; then
    gawk -F'\t' -v K=$GK "$RC"'{ s = $3; n = length(s); for (i = 1; i + K - 1 <= n; i++) { m = substr(s, i, K); if (m ~ /N/) continue
                                  r = rc(m); if (r < m) m = r; print m } }' seqs.tsv | sort -u > flank20.tsv
    gawk -v C=$GCAP 'FILENAME == ARGV[1] { want[$1]; next } NF >= 2 { any = 1; m = toupper($1); if ($2 + 0 > C && (m in want)) print m }
                     END { if (!any) exit 1 }' flank20.tsv "$KCOUNT" > rep20.tsv ||
        { echo "fs9: KCOUNT $KCOUNT has no 'KMER COUNT' lines: wrong file or format (jellyfish dump -c)" >&2; exit 1; }
    rm -f flank20.tsv
    else
    # the same rep20.tsv with the flank 20-mers held in memory as numbers (flank20.tsv was 1.3 GB on mau VES, left in the results)
    python3 -I "$HERE9/fs9_repmask.py" rep seqs.tsv "$KCOUNT" "$GCAP" > rep20.tsv ||
        { echo "fs9: KCOUNT $KCOUNT has no 'KMER COUNT' lines: wrong file or format (jellyfish dump -c)" >&2; exit 1; }
    fi
else
    SZ=$(gawk '{ s += $2 } END { print s + 0 }' "$G.fai")
    (( SZ <= 300000000 )) || { echo "fs9: genome > 300 Mb: give KCOUNT (jellyfish dump -c, k=20, canonical)" >&2; exit 1; }
    gawk -v K=$GK -v C=$GCAP "$RC"'
        function scan(   n, i, m, r) { n = length(seq); for (i = 1; i + K - 1 <= n; i++) { m = substr(seq, i, K); if (m ~ /N/) continue; r = rc(m); if (r < m) m = r; c[m]++ } }
        /^>/ { if (seq != "") scan(); seq = ""; next } { seq = seq toupper($0) }
        END { scan(); for (m in c) if (c[m] > C) print m }' "$G" > rep20.tsv
fi
rep_mask() {   # bases covered by a repeated 20-mer become N (seqs.tsv in place)
    if [[ $FS9_AWK != 1 ]]; then python3 -I "$HERE9/fs9_repmask.py" mask seqs.tsv rep20.tsv > seqs.masked.tsv; mv seqs.masked.tsv seqs.tsv; return; fi
    gawk -F'\t' -v K=$GK "$RC"'
        FILENAME == ARGV[1] { rep[$1]; next }
        { s = $3; n = length(s); for (i = 1; i <= n; i++) mk[i] = 0
          for (i = 1; i + K - 1 <= n; i++) { m = substr(s, i, K); if (m ~ /N/) continue; r = rc(m); if (r < m) m = r; if (m in rep) for (j = i; j < i + K; j++) mk[j] = 1 }
          o = ""; for (i = 1; i <= n; i++) o = o (mk[i] ? "N" : substr(s, i, 1)); print $1 "\t" $2 "\t" o }' rep20.tsv seqs.tsv > seqs.masked.tsv
    mv seqs.masked.tsv seqs.tsv
}

# 2-4) one pass: pass_run PATTERN MINKEYS IDMIN TAG
pass_run() {
    local PAT=$1 MK=$2 IDM=$3 TAG=$4
    if [[ $FS9_AWK != 1 ]]; then
        # steps 2-3 in memory, row for row the awk below (fs9_seeds.py): the awk writes one text row per seed shared by two copies
        # and sorts them (tbr VES, 620,807 copies: hits.p1.tsv 45 GB, hits.p2.tsv 92 GB, keys 3.6 GB each)
        python3 -I "$HERE9/fs9_seeds.py" seqs.tsv "$PAT" "$MK" --iw "$IW" --cap "$CAP" --band "$BAND" --start "$START" --slack 25 \
            --threads "$FS9_THREADS" | sort -u > "cand.$TAG.tsv"
    else
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
    rm -f "keys.$TAG.tsv" "hits.$TAG.tsv"
    fi
    # confirmation: score the diagonal and BAND diagonals around it (offset in copy a = offset in copy b + delta).
    # FS9_THREADS > 1: the candidate rows are cut into consecutive ranges, one awk each, and the outputs joined in order: the same
    # rows in the same order as one awk (one awk took 6-28 min per pass on the 340k-620k copy VES families)
    local NCAND P k lo hi rc=0; local -a pids=()
    NCAND=$(wc -l < "cand.$TAG.tsv"); P=$FS9_THREADS; (( NCAND < 20000 )) && P=1
    for (( k = 0; k < P; k++ )); do
    lo=$(( NCAND * k / P )); hi=$(( NCAND * (k + 1) / P ))
    gawk -F'\t' -v OFS='\t' -v IDM=$IDM -v ML=$MINLEN -v START=$START -v BAND=$BAND -v TAG=$TAG -v LO=$lo -v HI=$hi '
        FILENAME == ARGV[1] { sq[$1 SUBSEP $2] = $3; next }
        FNR <= LO || FNR > HI { next }
        { A =sq[$2 SUBSEP $1]; Bq = sq[$3 SUBSEP $1]; la = length(A); lb = length(Bq); bestid = 0; bestl = 0
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
          if (bestl >= ML) printf "%s\t%s\t%s\t%.1f\t%d\t%s\n", $2, $3, $1, bestid, bestl, TAG }' seqs.tsv "cand.$TAG.tsv" > "conf.$TAG.tsv.$k" &
    pids+=($!)
    done
    for k in "${pids[@]}"; do wait "$k" || rc=1; done
    (( rc == 0 )) || { echo "fs9: confirmation of $TAG failed" >&2; exit 1; }
    for (( k = 0; k < P; k++ )); do cat "conf.$TAG.tsv.$k"; done > "conf.$TAG.tsv"
    for (( k = 0; k < P; k++ )); do rm -f "conf.$TAG.tsv.$k"; done
}

find_pairs() {   # SM IDa IDb OUTFILE: pairs of one masking mode, pass 1 at IDa, pass 2 at IDb
    make_seqs $1; rep_mask
    pass_run 111111111111 $MK1 $2 p1
    pass_run 111010010100110111 $MK2 $3 p2
    # a pair found in both passes keeps its pass 1 row; the tier follows the identity
    gawk -F'\t' -v OFS='\t' -v ID1=$ID1 '{ key = $1 SUBSEP $2 SUBSEP $3; if (key in done) next; done[key] = 1
          print $1, $2, $3, $4, $5, ($4 >= ID1 ? "T1" : "T2"), "twin" }' conf.p1.tsv conf.p2.tsv | sort -k1,1 -k2,2 -k3,3 > "$4"
}
if [[ $SOFTMASK == both ]]; then
    find_pairs 0 $IDMASK $IDMASK pairs.masked95.body       # masked bases kept, near-identical flanks only
    find_pairs 1 $ID1 $ID2 pairs.strict.body               # masked bases as N, the usual tiers (its seqs.tsv stays)
    gawk -F'\t' '{ key = $1 SUBSEP $2 SUBSEP $3; if (!(key in done)) { done[key] = 1; print } }' pairs.strict.body pairs.masked95.body |
        sort -k1,1 -k2,2 -k3,3 > pairs.body
else
    find_pairs $SOFTMASK $ID1 $ID2 pairs.body
fi
# masked: a copy with fewer than MINFL testable (non-N) bases in a flank after masking has nothing to compare on that side; it is
# reported as "masked", never as "unique" (the first version called such copies unique). With SOFTMASK=both: the strict pass, so a
# masked copy was compared only at >= IDMASK.
gawk -F'\t' -v M=$MINFL '{ s = $3; nn = gsub(/N/, "", s); if (length(s) < M) print $1 }' seqs.tsv | sort -u > masked.txt
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
# element test (BODYCHECK=1, default): a pair stays a twin only if its two SINE bodies are more similar than random copies of the
# family (fs9_bodycheck.py; hot spots fail it, copies duplicated with their neighbourhood pass); flank-only results are kept
# as copy_status.flank.tsv and flank_groups.flank.tsv
if [[ ${BODYCHECK:-1} == 1 ]]; then
    python3 -I "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/fs9_bodycheck.py" "$G" "$B" . > bodycheck.log 2>&1 ||
        { cat bodycheck.log >&2; echo "fs9: body check failed" >&2; exit 1; }
fi
gawk -F'\t' '{ n[$2]++; t++ } END { printf "copies\t%d\n", t; for (k in n) printf "%s\t%d\t%.1f%%\n", k, n[k], 100 * n[k] / t }' copy_status.tsv | sort > summary.txt
cat summary.txt
