#!/usr/bin/env bash
# fs_lib.sh - functions shared by the flankscan stages that build candidate consensuses (stage 5, stage 6c).
# Sourced; expects OUT_DIR as the working directory and the variables RUN, T, FL, SHARE, SKIP (and NEAR,
# MINLEN for the end check) from the caller.

# fs_build_cons NAME - candidate NAME from cand/NAME.bed: the copies with FL bp flanks (lowercase), aligned
# (MAFFT L-INS-i), consensus over the element extended into the flanks while >= SHARE of the copies agree.
# Out: cand/NAME.aln.fa (consensus as row 1), cand/NAME.fa
fs_build_cons() {
    local NAME=$1
    # every copy with FL bp of genomic flank (clamped at contig ends), in the copy orientation:
    # flanks lowercase, element uppercase - the consensus is taken over the element and then extended
    # into the flanks as far as the copies keep agreeing (the element cut can stop short, e.g. inside a
    # simple-repeat tail that the A-tail rule does not follow: rsi r1_r3 lost 28 bp at its 3' end)
    gawk -F'\t' -v OFS='\t' -v F=$FL 'FNR == NR { len[$1] = $2; next }
        { ws = $2 - F; if (ws < 0) ws = 0; we = $3 + F; if (we > len[$1]) we = len[$1]
          l = $2 - ws; r = we - $3; if ($6 == "-") { t = l; l = r; r = t }   # 5 flank first in copy orientation
          print $1, ws, we, $1 ":" $2 "-" $3 "(" $6 ")|" l "|" r, 0, $6 }' \
        "$RUN/genome.clean.fa.fai" "cand/$NAME.bed" > cand/win.bed
    bedtools getfasta -fi "$RUN/genome.clean.fa" -bed cand/win.bed -s -nameOnly | seqkit seq -w 0 \
    | gawk '/^>/ { h = substr($0, 2); sub(/\([+-]\)$/, "", h); split(h, p, "|"); l = p[2]; r = p[3]
                   print ">" p[1]; next }
            { L = length($0); print tolower(substr($0, 1, l)) toupper(substr($0, l + 1, L - l - r)) tolower(substr($0, L - r + 1)) }' \
        > "cand/$NAME.copies.fa"
    mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet --thread "$T" \
        "cand/$NAME.copies.fa" > "cand/$NAME.mafft" 2> /dev/null
    # consensus: the element span = first..last column where >= half the copies have an UPPERCASE base;
    # inside it the majority base of every column where >= half the copies have a base (as before);
    # outside it, walk outward: a column with < SKIP of the copies occupied is passed over, a column
    # whose top base is carried by >= SHARE of ALL copies extends the consensus, the first other column
    # stops the walk. The alignment is written with the consensus (gapped) as row 1.
    gawk -v NAME="$NAME" -v ALN="cand/$NAME.aln.fa" -v SHARE=$SHARE -v SKIP=$SKIP '
        /^>/ { n++; h[n] = $0; next } { s[n] = s[n] $0 }
        function top(x,   i, b, occ) {                   # sets TB (top base), TC (its count), OC (occupied)
            delete k; occ = 0; TB = "-"; TC = 0
            for (i = 1; i <= n; i++) { b = toupper(substr(s[i], x, 1)); if (b != "-") { k[b]++; occ++ } }
            for (b in k) if (k[b] > TC || (k[b] == TC && b < TB)) { TC = k[b]; TB = b }
            OC = occ
        }
        END {
            L = length(s[1]); a = 0; z = 0
            for (x = 1; x <= L; x++) { up = 0
                for (i = 1; i <= n; i++) if (substr(s[i], x, 1) ~ /[ACGTN]/) up++
                if (up >= n / 2) { if (!a) a = x; z = x } }
            for (x = 1; x <= L; x++) g[x] = "-"
            for (x = a; x <= z; x++) { top(x); if (OC >= n / 2) g[x] = TB }
            for (x = a - 1; x >= 1; x--) { top(x); if (OC < SKIP * n) continue; if (TC >= SHARE * n) g[x] = TB; else break }
            for (x = z + 1; x <= L; x++) { top(x); if (OC < SKIP * n) continue; if (TC >= SHARE * n) g[x] = TB; else break }
            gcons = ""; cons = ""
            for (x = 1; x <= L; x++) { gcons = gcons g[x]; if (g[x] != "-") cons = cons g[x] }
            print ">" NAME > ALN; print gcons > ALN
            for (i = 1; i <= n; i++) { print h[i] > ALN; print s[i] > ALN }
            print ">" NAME; print cons
        }' "cand/$NAME.mafft" > "cand/$NAME.fa"
}

# fs_fold LIST STATUS - the same element built twice -> merged by best hit into groups; each group keeps the
# first of the list (the list is in rank order: largest peak first). LIST: name TAB cons_len, cand/NAME.fa must
# exist for each; STATUS (out): name TAB kept | same_as:ROOT
fs_fold() {
    local LIST=$1 STATUS=$2
cut -f1 "$LIST" | while read -r N; do cat "cand/$N.fa"; done > cand/all.fa
ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 cand/all.fa cand/all.fa 2> /dev/null > cand/self.m8 || true
gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { L[$1] = $2; ord[++n] = $1; next }             # list: name cons_len
    # m8: pid, q_s q_e ($7 $8), s_s s_e ($9 $10). Same element = the alignment covers >= 90 % of
    # BOTH: a shorter candidate contained in a longer one is a different element (r3 with its
    # internal repeat inside r1 + 39 bp + r3; r10 + r6 part inside r10 + 105 bp + group B)
    # Merge by BEST hit, not first hit: every candidate links to the one candidate it matches best
    # (highest identity among all pairs that pass), linked candidates form one group, and the group
    # keeps the candidate from the largest peak. First-hit in peak order merged rsi P42 (r8 + r8) into
    # P26 (r5h_r6, 90.6 %) although it is 99.0 % identical to P43 (r8 + r8, a smaller peak).
    $1 != $2 { if ($3 >= 90 && $8 - $7 + 1 >= 0.9 * L[$1] && $10 - $9 + 1 >= 0.9 * L[$2]) {
                   if ($3 > pid[$1, $2]) pid[$1, $2] = pid[$2, $1] = $3 } }
    function root(x) { while (up[x] != x) x = up[x]; return x }
    END {
        for (i = 1; i <= n; i++) { rank[ord[i]] = i; up[ord[i]] = ord[i] }
        for (i = 1; i <= n; i++) {                      # best partner of each candidate
            a = ord[i]; bb = ""; bv = 0
            for (j = 1; j <= n; j++) { b = ord[j]; if (b != a && ((a, b) in pid) && pid[a, b] > bv) { bv = pid[a, b]; bb = b } }
            if (bb != "") { ra = root(a); rb = root(bb); if (ra != rb) { if (rank[ra] < rank[rb]) up[rb] = ra; else up[ra] = rb } }
        }
        for (i = 1; i <= n; i++) { r = root(ord[i]); print ord[i], (r == ord[i] ? "kept" : "same_as:" r) }
    }' "$LIST" cand/self.m8 > "$STATUS"
    rm -f cand/all.fa cand/self.m8
}

# fs_endcheck NAME - are the ends of candidate NAME closed? The copies it was built from (cand/NAME.bed) with
# EFL bp of flank beyond each end (strand-aware), searched against the bank (cons.masked.fa, ssearch36
# E <= 1e-3); a copy's end is OPEN when its best hit lies within NEAR (100) bp of the junction (a linker of up to that length is allowed, rsi r1 + r3 has 39 bp) and is >= MINLEN bp.
# Out: cand/ends/NAME.h5, .h3 (one row per open copy: bed name, unit, cons start, cons end, flank lo, flank hi,
# bits; flank coordinates 1 = farthest from the element, EFL = next to it for 5', 1 = next to it for 3') and,
# on stdout, the ends.tsv row: name side5_share side5_unit side5_second side3_share side3_unit side3_second open
fs_endcheck() {
    local N=$1 EFL=${EFL:-250} NEAR=${NEAR:-100} MINLEN=${MINLEN:-40} OPENFRAC=${OPENFRAC:-0.5} S NC G="$RUN/genome.clean.fa"
    local B=cand/$N.bed; mkdir -p cand/ends
    gawk -F'\t' -v OFS='\t' '{ print $1, $2, $3, $4, 0, $6 }' "$B" > "cand/ends/$N.b6"
    bedtools flank -i "cand/ends/$N.b6" -g "$G.fai" -l $EFL -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null \
        | seqkit seq -w 0 | sed '/^>/s/([+-])$//' > "cand/ends/$N.f5.fa"
    bedtools flank -i "cand/ends/$N.b6" -g "$G.fai" -l 0 -r $EFL -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null \
        | seqkit seq -w 0 | sed '/^>/s/([+-])$//' > "cand/ends/$N.f3.fa"
    NC=$(wc -l < "$B")
    for S in 5 3; do
        ssearch36 -m 8 -E 1e-3 -Z 1000 -z 11 -T "$T" cons.masked.fa "cand/ends/$N.f$S.fa" 2> /dev/null \
        | gawk -F'\t' -v S=$S -v FL=$EFL -v NEAR=$NEAR -v ML=$MINLEN -v OFS='\t' '
            { lo = ($9 < $10 ? $9 : $10); hi = ($9 < $10 ? $10 : $9); if (hi - lo + 1 < ML) next
              near = (S == 5) ? (FL - hi <= NEAR) : (lo <= NEAR + 1)
              if (!near) next
              if (!($2 in be) || $12 > be[$2]) { be[$2] = $12; bu[$2] = $1; qs[$2] = ($7 < $8 ? $7 : $8); qe[$2] = ($7 < $8 ? $8 : $7); fl[$2] = lo; fh[$2] = hi } }
            END { for (k in bu) print k, bu[k], qs[k], qe[k], fl[k], fh[k], be[k] }' > "cand/ends/$N.h$S"
    done
    for S in 5 3; do
        gawk -F'\t' -v NC=$NC 'BEGIN { OFS = "\t" } { c[$2]++; n++ }
            END { t = ""; tn = 0; s2 = ""; sn = 0
                  for (u in c) { if (c[u] > tn) { s2 = t; sn = tn; t = u; tn = c[u] } else if (c[u] > sn) { s2 = u; sn = c[u] } }
                  printf "%.2f\t%s\t%s\n", (NC ? n / NC : 0), (t == "" ? "-" : t), (s2 == "" ? "-" : s2 ":" sn) }' "cand/ends/$N.h$S" > "cand/ends/$N.s$S"
    done
    paste <(printf "%s\n" "$N") "cand/ends/$N.s5" "cand/ends/$N.s3" \
    | gawk -F'\t' -v OFS='\t' -v F=$OPENFRAC '{ o = ($2 >= F ? "5" : "") ($5 >= F ? "3" : ""); print $1, $2, $3, $4, $5, $6, $7, (o == "" ? "-" : o) }'
}
