#!/usr/bin/env bash
# fs4_junctions.sh OUT_DIR RUN_DIR [MINK=50]
#
# Stage 4 of flankscan: from the junctions of stage 3, which copies are parts of something larger.
#
# A) Junction peaks (per family). All partners within NEAR bp of a copy are grouped by
#    (family, side, partner family, same / opp strand). Inside a group, a PEAK is a set of copies
#    whose junction is the same within TOL bp on all three coordinates: the gap, the copy's
#    consensus position at the junction (main_j) and the partner's (p_j). Found one coordinate at a
#    time: the gap value with most copies within +-TOL, then among those the main_j value, then p_j.
#    A peak is kept when it has >= MINK copies and >= FOLD x the chance expectation
#        E = copies of the family x (copies of the partner / genome bp) x (2 TOL + 1) x 1/2 (strand)
#    = how many copies would have a random partner copy at that gap (it ignores the consensus
#    positions, so it overestimates chance - conservative). Then its copies are removed and the next
#    peak in the group is looked for, until none passes (one group can hold several peaks: TB after
#    a whole-TA linker dimer and TB after a piecewise TA head).
#    Peak type, with the DOWNSTREAM unit of the pair (same strand only):
#      composite   the downstream unit starts at its consensus 5' end (<= 15): a new copy -> the
#                  pair is one larger element the library holds as two pieces
#      homodimer   composite with partner = the family itself
#      piecewise   the downstream unit starts mid-consensus: one element that two consensuses each
#                  cover in part -> a consensus problem, not two copies
#      inverted    opposite strands (head-to-head / tail-to-tail)
#    Per peak: the linker consensus (column majority of the gaps with the commonest length, in the
#    copy's consensus orientation), mean identity of those gaps to it, and its A share - an A-rich
#    linker may be a real dimer linker (Alu-like) OR an insertion hotspot in a partner's A tail;
#    both are reported, the call is his.
# B) Per copy and side (the nearest partner from stage 3):
#      nested     partner on BOTH sides, same family and strand, both gaps <= ADJ, and the host
#                 continues across the copy: its facing positions differ by -30..+10 bp (a target
#                 site duplication makes the difference negative)
#      atail      side 5, same strand, the partner reaches its own A tail (p_j within 15 bp of it),
#                 owns >= 5 bp of A run right before the copy (p_tail, stage 3) and the gap left is
#                 <= ADJ: the copy sits in the partner's A tail
#      peak member   the peak's type (composite / homodimer / piecewise / inverted)
#      chance     any other partner <= ADJ bp away;  near: ADJ < gap <= NEAR
#    Copy class = the first that applies of: satellite (stage 2) > nested > piecewise > homodimer >
#    composite > inverted > atail > chance > near > single ("nomain" = no consensus covers the core).
#
# Out: OUT/peaks.tsv           family side partner rel type n n_group expected main_j p_j gap
#                              linker linker_id linker_A
#      OUT/copies.tsv          wid locus family class side5 side3   (side = partner:label:gap or -)
#      OUT/family_summary.tsv  per family: copies and % of each class
set -euo pipefail
OUT=${1:?OUT_DIR}; RUN=${2:?RUN_DIR}; MINK=${3:-50}
TOL=10; ADJ=30; NEAR=200; FOLD=10
GSIZE=$(gawk '{ s += $2 } END { print s }' "$RUN/genome.clean.fa.fai")
cd "$OUT"

gawk -F'\t' -v OFS='\t' -v MINK=$MINK -v TOL=$TOL -v ADJ=$ADJ -v NEAR=$NEAR -v FOLD=$FOLD -v GSIZE=$GSIZE '
function rc(s,   o, i, c) { o = ""; for (i = length(s); i >= 1; i--) { c = substr(s, i, 1)
    o = o (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c == "a" ? "t" : c == "c" ? "g" : c == "g" ? "c" : c == "t" ? "a" : c) }; return o }
# the value (integer) with the most of V[1..n] within +-TOL of it
function mode(n, V,   i, c, lo, hi, x, s, best, bv, bx) {
    delete c; lo = 1e9; hi = -1e9
    for (i = 1; i <= n; i++) { c[V[i]]++; if (V[i] < lo) lo = V[i]; if (V[i] > hi) hi = V[i] }
    s = 0; for (x = lo - TOL; x <= lo + TOL; x++) s += c[x]
    best = lo; bv = s
    for (x = lo + 1; x <= hi; x++) { s += c[x + TOL] - c[x - TOL - 1]; if (s > bv) { bv = s; best = x } }
    # the commonest single value inside that densest window (the window centre is not the answer:
    # with ties the first window covering a tight cluster is centred TOL bp to its left)
    bv = -1; for (x = best - TOL; x <= best + TOL; x++) if (c[x] > bv) { bv = c[x]; bx = x }
    return bx
}
FILENAME == ARGV[1] { if (FNR > 1) { N[$3]++; locus[$1] = $2; fam[$1] = $3; W[++nw] = $1 }; next }   # loci.tsv
FILENAME == ARGV[2] { ctail[$1] = $3; next }                                                          # cons.tsv
FILENAME == ARGV[3] { if (FNR > 1 && $3 == "satellite") sat[$1] = 1; next }                           # trf.tsv
FNR == 1 { next }                                                                                     # junctions.tsv
{
    w = $1; sd = $3
    if ($9 == "nomain") { nomain[w] = 1; next }
    if ($9 == "-" || $16 > NEAR) next
    r = ++nr
    Rw[r] = w; Rf[r] = $2; Rs[r] = sd; Rmst[r] = $5; Rmj[r] = $8; Rp[r] = $9; Rrel[r] = $10; Rpj[r] = $11
    Rgap[r] = $16; Rseq[r] = $23; Rpt[r] = $25
    at[w, sd] = r
    g = $2 SUBSEP sd SUBSEP $9 SUBSEP $10
    if (!(g in gn)) G[++ng] = g
    gm[g, ++gn[g]] = r
}
END {
    # ---- A) peaks
    print "family", "side", "partner", "rel", "type", "n", "n_group", "expected", "main_j", "p_j", "gap",
          "linker", "linker_id", "linker_A" > "peaks.tsv"
    for (gi = 1; gi <= ng; gi++) {
        g = G[gi]; split(g, K, SUBSEP); F = K[1]; sd = K[2]; P = K[3]; rel = K[4]
        delete alive; na = 0; for (j = 1; j <= gn[g]; j++) { alive[gm[g, j]] = 1; na++ }
        E = N[F] * (N[P] / GSIZE) * (2 * TOL + 1) * 0.5
        while (na >= MINK) {
            n = 0; delete V; delete R
            for (r in alive) { n++; R[n] = r + 0; V[n] = Rgap[r] }
            gmode = mode(n, V)
            m = 0; delete R2; for (i = 1; i <= n; i++) if (Rgap[R[i]] >= gmode - TOL && Rgap[R[i]] <= gmode + TOL) R2[++m] = R[i]
            delete V; for (i = 1; i <= m; i++) V[i] = Rmj[R2[i]]
            mmode = mode(m, V)
            k = 0; delete R3; for (i = 1; i <= m; i++) if (Rmj[R2[i]] >= mmode - TOL && Rmj[R2[i]] <= mmode + TOL) R3[++k] = R2[i]
            delete V; for (i = 1; i <= k; i++) V[i] = Rpj[R3[i]]
            pmode = mode(k, V)
            np = 0; delete PK; for (i = 1; i <= k; i++) if (Rpj[R3[i]] >= pmode - TOL && Rpj[R3[i]] <= pmode + TOL) PK[++np] = R3[i]
            if (np < MINK || np < FOLD * E) break
            # type from the downstream unit of the pair
            if (rel != "same") type = "inverted"
            else {
                dstart = (sd == 3) ? pmode : mmode          # side 3: partner is downstream; side 5: the copy is
                type = (dstart > 15) ? "piecewise" : (P == F ? "homodimer" : "composite")
            }
            npk++
            # linker: gaps with the commonest exact length, in the copy consensus orientation
            delete lc; bl = 0; bc = 0
            for (i = 1; i <= np; i++) { lc[Rgap[PK[i]]]++; if (lc[Rgap[PK[i]]] > bc) { bc = lc[Rgap[PK[i]]]; bl = Rgap[PK[i]] } }
            lk = "."; lid = "."; lA = "."
            if (bl > 0 && bl <= 300) {
                delete col; ns = 0
                for (i = 1; i <= np; i++) { r = PK[i]; if (Rgap[r] != bl || Rseq[r] == ".") continue
                    s = toupper(Rmst[r] == "-" ? rc(Rseq[r]) : Rseq[r]); S[++ns] = s
                    for (x = 1; x <= bl; x++) col[x, substr(s, x, 1)]++ }
                lk = ""; for (x = 1; x <= bl; x++) { bb = "N"; bv = 0
                    split("A C G T", B, " "); for (y = 1; y <= 4; y++) if (col[x, B[y]] > bv) { bv = col[x, B[y]]; bb = B[y] }
                    lk = lk bb }
                id = 0; for (i = 1; i <= ns; i++) { mt = 0; for (x = 1; x <= bl; x++) if (substr(S[i], x, 1) == substr(lk, x, 1)) mt++; id += mt / bl }
                lid = sprintf("%.2f", ns ? id / ns : 0); t = lk; lA = sprintf("%.2f", gsub(/A/, "", t) / bl)
            }
            print F, sd, P, rel, type, np, gn[g], sprintf("%.2f", E), mmode, pmode, gmode, lk, lid, lA > "peaks.tsv"
            for (i = 1; i <= np; i++) { lab[PK[i]] = type; delete alive[PK[i]]; na-- }
        }
    }
    # ---- B) per copy
    print "wid", "locus", "family", "class", "side5", "side3" > "copies.tsv"
    split("satellite nested piecewise homodimer composite inverted atail chance near single nomain", ORD, " ")
    for (i = 1; i <= 11; i++) rank[ORD[i]] = i
    for (i = 1; i <= nw; i++) {
        w = W[i]; F = fam[w]
        r5 = at[w, 5]; r3 = at[w, 3]
        for (sd = 5; sd >= 3; sd -= 2) {
            r = (sd == 5) ? r5 : r3; L = "single"
            if (r) {
                if (r in lab) L = lab[r]
                else if (sd == 5 && Rrel[r] == "same" && Rpj[r] >= ctail[Rp[r]] - 1 - 15 &&
                         Rpt[r] >= 5 && Rgap[r] <= ADJ) L = "atail"
                else L = (Rgap[r] <= ADJ) ? "chance" : "near"
            }
            lb[sd] = L; tx[sd] = r ? Rp[r] ":" L ":" Rgap[r] : "-"
        }
        # nested: the same host on both sides, continuing across the copy
        if (r5 && r3 && Rp[r5] == Rp[r3] && Rrel[r5] == "same" && Rrel[r3] == "same" && Rgap[r5] <= ADJ && Rgap[r3] <= ADJ) {
            d = Rpj[r3] - Rpj[r5]                           # sides are in the copy consensus orientation
            if (d >= -30 && d <= 10) { lb[5] = lb[3] = "nested"; tx[5] = Rp[r5] ":nested:" Rgap[r5]; tx[3] = Rp[r3] ":nested:" Rgap[r3] }
        }
        c = (lb[5] in rank && rank[lb[5]] < rank[lb[3]]) ? lb[5] : lb[3]
        if (sat[w]) c = "satellite"
        if (nomain[w]) c = "nomain"
        cnt[F, c]++
        print w, locus[w], F, c, tx[5], tx[3] > "copies.tsv"
    }
    # ---- family summary
    hdr = "family\tcopies"; for (i = 1; i <= 11; i++) hdr = hdr "\t" ORD[i]
    print hdr > "family_summary.tsv"; close("family_summary.tsv")
    for (F in N) { line = F "\t" N[F]; for (i = 1; i <= 11; i++) line = line "\t" sprintf("%.1f%%", 100 * cnt[F, ORD[i]] / N[F])
                   print line | "sort -t\"\t\" -k2,2nr >> family_summary.tsv" }
    close("sort -t\"\t\" -k2,2nr >> family_summary.tsv")
    printf "fs4: %d peaks; classes in copies.tsv, per family in family_summary.tsv\n", npk > "/dev/stderr"
}' loci.tsv cons.tsv trf.tsv junctions.tsv
