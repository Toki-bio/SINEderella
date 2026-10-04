#!/usr/bin/env bash
# fs_lib.sh - functions shared by the flankscan stages that build candidate consensuses (stage 5, stage 6c).
# Sourced; expects OUT_DIR as the working directory and the variables RUN, T, FL, SHARE, SKIP (and NEAR,
# MINLEN for the end check) from the caller.

# fs_build_cons NAME - candidate NAME from cand/NAME.bed. The copies are cut with FLW (300) bp of genomic flank
# (clamped at contig ends; flank lowercase, element uppercase) and aligned (MAFFT L-INS-i). Consensus:
#   1) over the element columns (first..last column where >= half the copies have an UPPERCASE base): the
#      majority base wherever >= half the copies have a base; then walked outward through the flank COLUMNS while
#      >= SHARE of all copies agree (a column occupied by < SKIP of the copies is passed over);
#   2) then each end is extended through the PROXIMAL FLANK: the TAILW (60) bases each copy has beyond that
#      column are taken degapped and aligned on their own (a short MAFFT job), and the consensus is walked along
#      those columns while >= TSHARE (0.4) of all copies carry the same base (looser than SHARE: the copies
#      differ in how many A they carry before the terminator, so no base reaches 60 % there). The columns of the big alignment are
#      unreliable in a flank (a poly-A of varying length, indels, 1-3 extra A before the terminator), a small
#      alignment of the proximal flank is not: it takes the consensus through the polyadenylation signal, the
#      terminator and the A tail (rsi chain C5: ...CCCCAATAAAATCTT + A).
# Out: cand/NAME.fa (the consensus), cand/NAME.aln.fa = the PLATE: consensus as row 1, then every copy as
#   [FLD (100) bp of flank, degapped, packed against the element, lowercase] [extension columns] [element
#   columns as aligned] [extension columns] [FLD bp of flank]; a flank shorter than FLD (contig end) is padded
#   with "-". No gap columns in the flanks.
fs_build_cons() {
    local NAME=$1 FLW=${FLW:-300} FLD=${FLD:-100} TAILW=${TAILW:-60} TSHARE=${TSHARE:-0.4} P="cand/$1"
    # every copy with FLW bp of genomic flank (clamped at contig ends), in the copy orientation
    gawk -F'\t' -v OFS='\t' -v F=$FLW 'FNR == NR { len[$1] = $2; next }
        { ws = $2 - F; if (ws < 0) ws = 0; we = $3 + F; if (we > len[$1]) we = len[$1]
          l = $2 - ws; r = we - $3; if ($6 == "-") { t = l; l = r; r = t }   # 5 flank first in copy orientation
          print $1, ws, we, $1 ":" $2 "-" $3 "(" $6 ")|" l "|" r, 0, $6 }' \
        "$RUN/genome.clean.fa.fai" "$P.bed" > cand/win.bed
    bedtools getfasta -fi "$RUN/genome.clean.fa" -bed cand/win.bed -s -nameOnly | seqkit seq -w 0 \
    | gawk '/^>/ { h = substr($0, 2); sub(/\([+-]\)$/, "", h); split(h, p, "|"); l = p[2]; r = p[3]
                   print ">" p[1]; next }
            { L = length($0); print tolower(substr($0, 1, l)) toupper(substr($0, l + 1, L - l - r)) tolower(substr($0, L - r + 1)) }' \
        > "$P.copies.fa"
    mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet --thread "$T" --threadit 0 \
        "$P.copies.fa" > "$P.mafft" 2> /dev/null
    # pass A: the element columns and the column walk; the proximal flank of every copy goes to its own alignment
    gawk -v SHARE=$SHARE -v SKIP=$SKIP -v TAILW=$TAILW -v ST="$P.state" -v RF="$P.right.fa" -v LF="$P.left.fa" '
        /^>/ { n++; h[n] = $0; next } { s[n] = s[n] $0 }
        function top(x,   i, b, occ) {
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
            a2 = 0; z2 = 0; for (x = 1; x <= L; x++) if (g[x] != "-") { if (!a2) a2 = x; z2 = x }
            gc = ""; for (x = a2; x <= z2; x++) gc = gc g[x]
            print n "\t" gc > ST                                                # state: copies, gapped core consensus
            for (i = 1; i <= n; i++) { lf = ""; rf = ""
                for (x = 1; x < a2; x++) { c = substr(s[i], x, 1); if (c != "-") lf = lf c }
                for (x = z2 + 1; x <= L; x++) { c = substr(s[i], x, 1); if (c != "-") rf = rf c }
                print i "\t" substr(s[i], a2, z2 - a2 + 1) "\t" lf "\t" rf > ST ".copies"
                if (length(rf) > 0) { print ">" i > RF; print substr(rf, 1, TAILW) > RF }
                if (length(lf) > 0) { print ">" i > LF; print substr(lf, (length(lf) > TAILW ? length(lf) - TAILW + 1 : 1)) > LF } }
            close(ST); close(ST ".copies"); close(RF); close(LF)
        }' "$P.mafft"
    : > "$P.right.aln"; : > "$P.left.aln"
    if [[ -s "$P.right.fa" ]]; then mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --quiet --thread "$T" --threadit 0 "$P.right.fa" > "$P.right.aln" 2> /dev/null; fi
    if [[ -s "$P.left.fa" ]]; then mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --quiet --thread "$T" --threadit 0 "$P.left.fa" > "$P.left.aln" 2> /dev/null; fi
    # pass B: extend through the proximal-flank alignments, pack the plate
    gawk -v NAME="$NAME" -v ALN="$P.aln.fa" -v SHARE=$TSHARE -v TAILA=${TAILA:-15} -v FLD=$FLD -v RAF="$P.right.aln" -v LAF="$P.left.aln" '
        function dashes(m,   d) { d = sprintf("%*s", m, ""); gsub(/ /, "-", d); return d }
        function readaln(file, arr,   l, id) { while ((getline l < file) > 0) { if (l ~ /^>/) id = substr(l, 2); else arr[id] = arr[id] toupper(l) }; close(file) }
        FILENAME == ARGV[1] { n = $1; gc = $2; next }
        FILENAME == ARGV[2] { elem[$1] = $2; Lf[$1] = $3; Rf[$1] = $4; next }
        END {
            readaln(RAF, RA); readaln(LAF, LA)
            e3 = ""; c3 = 0; run = 0; pend = ""; Lr = 0; for (id in RA) { Lr = length(RA[id]); break }
            for (x = 1; x <= Lr; x++) { delete kc; tc = 0; tb = ""
                oc = 0; for (id in RA) { b = substr(RA[id], x, 1); if (b != "-") { kc[b]++; oc++ } }
                if (oc < 0.5 * n) { pend = pend "-"; continue }          # sparse column (leading gaps of the local alignment): kept as a gap
                for (b in kc) if (kc[b] > tc || (kc[b] == tc && b < tb)) { tc = kc[b]; tb = b }
                if (tc >= SHARE * n) { e3 = e3 pend tb; pend = ""; c3 = length(e3); run = (tb == "A") ? run + 1 : 0; if (run >= TAILA) break } else break }
            e5 = ""; c5 = 0; pend = ""; Ll = 0; for (id in LA) { Ll = length(LA[id]); break }
            for (x = Ll; x >= 1; x--) { delete kc; tc = 0; tb = ""
                oc = 0; for (id in LA) { b = substr(LA[id], x, 1); if (b != "-") { kc[b]++; oc++ } }
                if (oc < 0.5 * n) { pend = "-" pend; continue }
                for (b in kc) if (kc[b] > tc || (kc[b] == tc && b < tb)) { tc = kc[b]; tb = b }
                if (tc >= SHARE * n) { e5 = tb pend e5; pend = ""; c5 = length(e5) } else break }
            cons = gc; gsub(/-/, "", cons); ec5 = e5; gsub(/-/, "", ec5); ec3 = e3; gsub(/-/, "", ec3)
            print ">" NAME > ALN; print dashes(FLD) e5 gc e3 dashes(FLD) > ALN
            for (i = 1; i <= n; i++) {
                x3 = (i in RA) ? substr(RA[i], 1, c3) : dashes(c3)
                u3 = x3; gsub(/-/, "", u3)
                rr = substr(Rf[i], length(u3) + 1)
                fl3 = (length(rr) >= FLD) ? substr(rr, 1, FLD) : rr dashes(FLD - length(rr))
                x5 = (c5 == 0) ? "" : ((i in LA) ? substr(LA[i], length(LA[i]) - c5 + 1) : dashes(c5))
                u5 = x5; gsub(/-/, "", u5)
                lr = substr(Lf[i], 1, length(Lf[i]) - length(u5))
                fl5 = (length(lr) >= FLD) ? substr(lr, length(lr) - FLD + 1) : dashes(FLD - length(lr)) lr
                print ">" i > ALN; print tolower(fl5) toupper(x5) elem[i] toupper(x3) tolower(fl3) > ALN }
            print ">" NAME; print ec5 cons ec3
        }' "$P.state" "$P.state.copies" > "$P.fa"
    # the plate rows carry the copy names of the alignment (pass B only knows indices)
    gawk 'FILENAME == ARGV[1] { if (/^>/) { n++; nm[n] = $0 }; next }
          /^>/ && FNR > 1 { k++; print nm[k]; next } { print }' "$P.copies.fa" "$P.aln.fa" > "$P.aln.tmp" && mv "$P.aln.tmp" "$P.aln.fa"
    rm -f "$P.mafft" "$P.state" "$P.state.copies" "$P.right.fa" "$P.left.fa" "$P.right.aln" "$P.left.aln"
}

# fs_fold LIST STATUS - the same element built twice -> merged by best hit into groups; each group keeps the
# first of the list (the list is in rank order: largest peak first). LIST: name TAB cons_len [TAB key], cand/NAME.fa must
# exist for each; STATUS (out): name TAB kept | same_as:ROOT
# Rules: (1) same element = the alignment of the two consensuses covers >= 90 % of BOTH, identity >= 90 %; a shorter
# candidate contained in a longer one is a different element (r3 with its internal repeat inside r1 + 39 bp + r3).
# (2) the part of one candidate that the other does not cover (>= XMIN = 40 bp) must not be a known unit: when it hits
# the bank (cons.masked.fa, >= 40 bp, E <= 1e-3) the candidates differ (the four-unit rsi chain, 632 bp, carries a
# 49 bp repeat of r3 that the three-unit chain, 583 bp, lacks; routes that differ only in flank-extended ends fold).
# (3) key (3rd list column, optional): homodimers carry "H:<family>" and merge only with homodimers of the same family -
# a tandem of one unit is its own element and is never absorbed by a composite of other units (rsi r8 + r8 into r5 + r6).
# Merge by BEST hit, not first hit: every candidate links to the one candidate it matches best (highest identity among
# the pairs that pass), linked candidates form one group, the group keeps the candidate from the largest peak.
fs_fold() {
    local LIST=$1 STATUS=$2 XMIN=${XMIN:-40}
    cut -f1 "$LIST" | while read -r N; do cat "cand/$N.fa"; done > cand/all.fa
    ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 cand/all.fa cand/all.fa 2> /dev/null > cand/self.m8 || true
    # pairs that pass rules 1 and 3, with the uncovered stretches (name start end) of both candidates
    gawk -F'\t' -v OFS='\t' -v XMIN=$XMIN '
        FILENAME == ARGV[1] { L[$1] = $2; K[$1] = $3; next }
        $1 != $2 && K[$1] == K[$2] && $3 >= 90 && $9 < $10 && $8 - $7 + 1 >= 0.9 * L[$1] && $10 - $9 + 1 >= 0.9 * L[$2] {
            print $1, $2, $3, ($7 - 1 >= XMIN ? $1 ":1-" $7 - 1 : "-"), (L[$1] - $8 >= XMIN ? $1 ":" $8 + 1 "-" L[$1] : "-"),
                  ($9 - 1 >= XMIN ? $2 ":1-" $9 - 1 : "-"), (L[$2] - $10 >= XMIN ? $2 ":" $10 + 1 "-" L[$2] : "-") }' "$LIST" cand/self.m8 > cand/fold.pairs
    # rule 2: an uncovered stretch that is a known unit separates the candidates
    : > cand/fold.block
    while IFS=$'\t' read -r A B PID X1 X2 X3 X4; do
        for X in "$X1" "$X2" "$X3" "$X4"; do
            [[ "$X" == - ]] && continue
            N=${X%%:*}; R=${X#*:}; S=${R%-*}; E=${R#*-}
            # (the string constant used to span a line break; gawk rejects that as an "unterminated string", so this
            #  gawk failed, fold.x.fa stayed empty and rule 2 never blocked a pair - found in the audit of 2026-10-05)
            gawk -v S=$S -v E=$E '!/^>/ { print ">x\n" substr($0, S, E - S + 1) }' "cand/$N.fa" > cand/fold.x.fa
            if [[ -n $(ssearch36 -m 8 -E 1e-3 -Z 1000 cand/fold.x.fa cons.masked.fa 2> /dev/null | gawk -F'\t' '$4 >= 40 { print "hit"; exit }') ]]; then
                printf "%s\t%s\n" "$A" "$B" >> cand/fold.block; break
            fi
        done
    done < <(gawk -F'\t' '$4 != "-" || $5 != "-" || $6 != "-" || $7 != "-"' cand/fold.pairs)
    gawk -F'\t' -v OFS='\t' '
        FILENAME == ARGV[1] { ord[++n] = $1; next }
        FILENAME == ARGV[2] { bl[$1, $2] = bl[$2, $1] = 1; next }
        { if (!(($1, $2) in bl) && $3 > pid[$1, $2]) pid[$1, $2] = pid[$2, $1] = $3 }
        function root(x) { while (up[x] != x) x = up[x]; return x }
        END {
            for (i = 1; i <= n; i++) { rank[ord[i]] = i; up[ord[i]] = ord[i] }
            for (i = 1; i <= n; i++) {                      # best partner of each candidate
                a = ord[i]; bb = ""; bv = 0
                for (j = 1; j <= n; j++) { b = ord[j]; if (b != a && ((a, b) in pid) && pid[a, b] > bv) { bv = pid[a, b]; bb = b } }
                if (bb != "") { ra = root(a); rb = root(bb); if (ra != rb) { if (rank[ra] < rank[rb]) up[rb] = ra; else up[ra] = rb } }
            }
            for (i = 1; i <= n; i++) { r = root(ord[i]); print ord[i], (r == ord[i] ? "kept" : "same_as:" r) }
        }' "$LIST" cand/fold.block cand/fold.pairs > "$STATUS"
    rm -f cand/all.fa cand/self.m8 cand/fold.pairs cand/fold.block cand/fold.x.fa
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
