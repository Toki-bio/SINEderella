#!/usr/bin/env bash
# fs5_build.sh OUT_DIR RUN_DIR [NBEST=60] [MINK=50] [THREADS=8]
#
# Stage 5 of flankscan: a candidate consensus for every structure stage 4 found - the whole element
# (first unit start -> last unit end, spacer included), built from its best copies. The way
# r10_groupB, r1_r3 and r5h_r6 were built by hand (Tal rsi/REFINEMENT.md §3, §12).
#
# 1) Peaks used: same strand, type composite / homodimer / piecewise, >= MINK copies. Inverted
#    pairs are not built (head-to-head is not one element read in one direction).
#    Every pair is seen twice, once from each unit (r3 side 5 r1 = r1 side 3 r3). In the ELEMENT's
#    orientation both peaks have the same upstream unit U, downstream unit D, U's consensus end,
#    D's consensus start and gap (each within TOL): the smaller one is a MIRROR of the larger and
#    is not built again (its copies count for the larger one).
# 2) Per copy of a peak, the element in window coordinates = from the outer end of one unit to the
#    outer end of the other, with the copy's and the partner's own A tails (stage 3) and the A run
#    after the element's 3' end (same >= 80 % A rule as stage 3). Ranked by main + partner
#    bitscore; copies whose genomic spans overlap are one element (kept once: the mirror peak sees
#    the same element from its other unit). The best NBEST are cut from the genome.
# 3) MAFFT L-INS-i; consensus = majority base of every column where >= half the copies have a base.
# 4) Candidates that are the same element twice (one ssearch36 alignment, >= 90 % identity,
#    covering >= 90 % of BOTH) are kept once: the one from the larger peak. Stage 6 counts the
#    copies of a folded candidate for the kept one.
#
# Out: OUT/candidates.fa       the kept candidate consensuses (name = U__D_Pn, Pn = the peak id)
#      OUT/candidates.tsv      name peak type upstream downstream n_peak mirrors elements n_used
#                              len_p10 len_med len_p90 cons_len status   (status kept / same_as:NAME)
#                              elements = distinct genomic elements among the peak's copies
#      OUT/cand/NAME.aln.fa    the alignment, consensus as row 1 (to look at); NAME.bed the copies used
set -euo pipefail
OUT=${1:?OUT_DIR}; RUN=$(readlink -f "${2:?RUN_DIR}"); NBEST=${3:-60}; MINK=${4:-50}; T=${5:-8}
TOL=10
HDR="name\tpeak\ttype\tupstream\tdownstream\tn_peak\tmirrors\telements\tn_used\tlen_p10\tlen_med\tlen_p90\tcons_len\tstatus\n"
cd "$OUT"; rm -rf cand; mkdir -p cand

# 1) peaks to build, mirrors folded in -> cand/peaks_used.tsv: peak type U D n mirrors
gawk -F'\t' -v OFS='\t' -v MINK=$MINK -v TOL=$TOL '
function abs(x) { return x < 0 ? -x : x }
NR > 1 && $4 == "same" && $5 ~ /^(composite|homodimer|piecewise)$/ && $6 >= MINK {
    n++; id[n] = $15; ty[n] = $5; cnt[n] = $6; gp[n] = $11
    if ($2 == 3) { U[n] = $1; D[n] = $3; ue[n] = $9;  ds[n] = $10 }    # the copy is upstream
    else         { U[n] = $3; D[n] = $1; ue[n] = $10; ds[n] = $9  }    # the partner is upstream
}
END {
    for (i = 1; i <= n; i++) o[i] = i                                   # larger peaks first
    for (i = 2; i <= n; i++) for (j = i; j > 1 && cnt[o[j]] > cnt[o[j-1]]; j--) { t = o[j]; o[j] = o[j-1]; o[j-1] = t }
    for (a = 1; a <= n; a++) {
        i = o[a]; mir = 0
        for (b = 1; b <= nk; b++) { k = K[b]
            if (U[k] == U[i] && D[k] == D[i] && abs(gp[k] - gp[i]) <= TOL && abs(ue[k] - ue[i]) <= TOL && abs(ds[k] - ds[i]) <= TOL) { mir = k; break } }
        if (mir) { M[mir] = M[mir] (M[mir] ? "," : "") id[i]; continue }
        K[++nk] = i
    }
    for (b = 1; b <= nk; b++) { k = K[b]; print id[k], ty[k], U[k], D[k], cnt[k], (M[k] ? M[k] : "-") }
}' peaks.tsv > cand/peaks_used.tsv
if [[ ! -s cand/peaks_used.tsv ]]; then
    : > candidates.fa; printf "$HDR" > candidates.tsv
    echo "fs5: no peak to build" >&2; exit 0
fi

# 2) the element span of every member copy -> cand/spans.tsv: peak score contig start end id strand
#    (genomic BED coordinates; strand = the element's genomic strand, both units have it)
gawk -F'\t' -v OFS='\t' '
function tailrun(s, x, step, base,   n, a, c, best) {          # as in fs3: the A run from x on
    a = 0; best = 0
    for (n = 1; x >= 1 && x <= length(s); n++) {
        c = toupper(substr(s, x, 1)); if (c == base) a++
        if (n >= 10 && a / n < 0.5) break
        if (a / n >= 0.8 && c == base) best = n
        x += step
    }
    return best
}
FILENAME == ARGV[1] { own[$1] = $1; n = split($6, m, ","); for (i = 1; i <= n; i++) if (m[i] != "-") own[m[i]] = $1; next }
FILENAME == ARGV[2] { if (FNR > 1 && ($1 in own)) { mem[$2, $3] = own[$1]; want[$2] = 1 }; next }   # members.tsv
FILENAME == ARGV[3] { if (FNR > 1 && ($1 in want)) { bits[$1, $3] = $12                          # units.tsv
                          if ($4 == "main") { ML[$1] = $10; MH[$1] = $11; mst[$1] = $6; mb[$1] = $12 } }; next }
FILENAME == ARGV[4] { if (FNR > 1 && ($1 in want)) { if ($3 == 3) mt[$1] = $24 + 0              # junctions.tsv
                          J[$1, $3] = $14 "\t" $15 "\t" $26 }; next }
FILENAME == ARGV[5] { if ($4 in want) { ctg[$4] = $1; ws[$4] = $2; we[$4] = $3; wst[$4] = $6 }; next }   # windows.bed
/^>/ { w = substr($1, 2); next }                                                                  # windows.fa
!(w in want) { next }
{
    for (sd = 3; sd <= 5; sd += 2) {
        if (!((w, sd) in mem)) continue
        split(J[w, sd], j, "\t")
        lo = ML[w]; hi = MH[w]
        if (mst[w] == "+") hi += mt[w]; else lo -= mt[w]            # the copy s own A tail
        if (j[1] < lo) lo = j[1]; if (j[2] > hi) hi = j[2]          # the partner (with its tail)
        if (mst[w] == "+") hi += tailrun($0, hi + 1, 1, "A")        # the A run after the element
        else               lo -= tailrun($0, lo - 1, -1, "T")
        # window -> genome: the window is the genome read on its own strand (wst); the element s
        # genomic strand = window strand x unit strand in the window
        if (wst[w] == "+") { gs = ws[w] + lo - 1; ge = ws[w] + hi } else { gs = we[w] - hi; ge = we[w] - lo + 1 }
        print mem[w, sd], mb[w] + bits[w, j[3]], ctg[w], gs, ge, w "_" sd, (wst[w] == mst[w] ? "+" : "-")
    }
}' cand/peaks_used.tsv members.tsv units.tsv junctions.tsv windows.bed <(seqkit seq -w 0 windows.fa) \
| sort -t$'\t' -k1,1 -k2,2gr > cand/spans.tsv

# 3) per peak: distinct elements, the best NBEST, alignment, majority consensus
: > cand/built.tsv
while IFS=$'\t' read -r P TY U D NP MIR; do
    NAME="${U}__${D}_${P}"
    # overlapping spans = one element (seen from both of its units); the best-scoring one is kept
    gawk -F'\t' -v OFS='\t' -v P="$P" -v NB=$NBEST '
        $1 != P { next }
        { for (i = 1; i <= nk; i++) if (c[i] == $3 && $4 < e[i] && $5 > s[i]) next
          nk++; c[nk] = $3; s[nk] = $4; e[nk] = $5
          if (nk <= NB) print $3, $4, $5, $6, $2, $7 }
        END { print nk + 0 > "cand/n.tmp" }' cand/spans.tsv > "cand/$NAME.bed"
    NEL=$(cat cand/n.tmp)
    bedtools getfasta -fi "$RUN/genome.clean.fa" -bed "cand/$NAME.bed" -s -nameOnly > "cand/$NAME.copies.fa"
    mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --quiet --thread "$T" "cand/$NAME.copies.fa" > "cand/$NAME.mafft" 2> /dev/null
    # majority consensus; the alignment written with the consensus (gapped) as row 1
    gawk -v NAME="$NAME" -v ALN="cand/$NAME.aln.fa" '
        /^>/ { n++; h[n] = $0; next } { s[n] = s[n] toupper($0) }
        END {
            L = length(s[1]); cons = ""; gcons = ""
            for (x = 1; x <= L; x++) {
                delete k; occ = 0; bb = "-"; bv = 0
                for (i = 1; i <= n; i++) { b = substr(s[i], x, 1); if (b != "-") { k[b]++; occ++ } }
                if (occ >= n / 2) { for (b in k) if (k[b] > bv || (k[b] == bv && b < bb)) { bv = k[b]; bb = b }
                                    cons = cons bb; gcons = gcons bb }
                else gcons = gcons "-"
            }
            print ">" NAME > ALN; print gcons > ALN
            for (i = 1; i <= n; i++) { print h[i] > ALN; print s[i] > ALN }
            print ">" NAME; print cons
        }' "cand/$NAME.mafft" > "cand/$NAME.fa"
    LENS=$(gawk -F'\t' '{ print $3 - $2 }' "cand/$NAME.bed" | sort -n \
           | gawk -v OFS='\t' '{ v[NR] = $1 } END { print v[int(NR * 0.1) + 1], v[int((NR + 1) / 2)], v[(int(NR * 0.9) > 0) ? int(NR * 0.9) : 1] }')
    CL=$(tail -1 "cand/$NAME.fa" | tr -d '\n' | wc -c)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NAME" "$P" "$TY" "$U" "$D" "$NP" "$MIR" "$NEL" \
        "$(wc -l < "cand/$NAME.bed")" "$LENS" "$CL" >> cand/built.tsv
    rm -f "cand/$NAME.mafft" cand/n.tmp
done < cand/peaks_used.tsv

# 4) the same element built twice -> keep the larger peak's (built.tsv is in peak-size order)
while read -r N; do cat "cand/$N.fa"; done < <(cut -f1 cand/built.tsv) > cand/all.fa
ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 cand/all.fa cand/all.fa 2> /dev/null > cand/self.m8 || true
gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { L[$1] = $13; ord[++n] = $1; next }             # built.tsv: cons_len
    # m8: pid, q_s q_e ($7 $8), s_s s_e ($9 $10). Same element = the alignment covers >= 90 % of
    # BOTH: a shorter candidate contained in a longer one is a different element (r3 with its
    # internal repeat inside r1 + 39 bp + r3; r10 + r6 part inside r10 + 105 bp + group B)
    $1 != $2 { if ($3 >= 90 && $8 - $7 + 1 >= 0.9 * L[$1] && $10 - $9 + 1 >= 0.9 * L[$2]) same[$1, $2] = same[$2, $1] = 1 }
    END { for (i = 1; i <= n; i++) { st[ord[i]] = "kept"
              for (j = 1; j < i; j++) if (st[ord[j]] == "kept" && same[ord[i], ord[j]]) { st[ord[i]] = "same_as:" ord[j]; break } }
          for (i = 1; i <= n; i++) print ord[i], st[ord[i]] }' cand/built.tsv cand/self.m8 > cand/status.tsv
{ printf "$HDR"; paste cand/built.tsv <(cut -f2 cand/status.tsv); } > candidates.tsv
gawk -F'\t' 'NR > 1 && $14 == "kept" { print $1 }' candidates.tsv | while read -r N; do cat "cand/$N.fa"; done > candidates.fa
rm -f cand/all.fa cand/self.m8 cand/built.tsv cand/status.tsv
echo "fs5: $(grep -c '^>' candidates.fa) candidates kept of $(($(wc -l < candidates.tsv) - 1)) built -> $OUT/candidates.fa" >&2
