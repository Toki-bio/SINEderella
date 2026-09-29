#!/usr/bin/env bash
# fs3_partners.sh OUT_DIR CONSENSUSES.fa [THREADS=8]
#
# Stage 3 of flankscan: which consensus units lie in each window (stage 1 windows, flank repeats
# masked by stage 2), and what sits right next to the copy on each side.
#
# 1) Consensuses are masked before the search (length kept, so consensus coordinates stay valid):
#    - the 3' A tail: the longest 3' end that is >= 80 % A and starts with A (>= 6 bp); scanned from
#      the 3' end, stops once the A share falls below 50 %;
#    - low complexity (dustmasker).
#    Otherwise every copy would "match" every other copy's A tail.
# 2) ssearch36 -m 8 -E 1e-5 -z 11, every masked consensus vs every masked window, both strands
#    (-z 11 because the library is all homologs; default statistics return nothing). ssearch36
#    reports several alignments of one consensus in one window (a split or repeated partner).
# 3) Per window, hits are taken best bitscore first; the bases a better hit already covers belong to
#    it, and a weaker hit keeps only its parts outside every kept unit (>= MIN_ALN bp each; its
#    consensus coordinates follow linearly along the alignment). -> the window's UNITS, disjoint.
#    Why trim instead of drop: local alignments overrun a junction by 10-30 bp whenever the sequence
#    across it happens to score positive (toy: TC's alignment ran 31 bp into an inserted TA copy
#    and through its A tail). Dropping such a hit loses a real partner; trimming keeps it and the
#    stronger unit decides where the junction is. A hit spanning a whole inserted copy splits into
#    two pieces (a nested insertion).
#    The MAIN unit is the first kept unit overlapping the core by >= 30 bp. Every other unit is
#    upstream or downstream of it.
# 4) Per window and side, the NEAREST unit on that side is the partner. Sides are named in the MAIN
#    unit's consensus orientation: side 3 = beyond the main unit's consensus 3' end (downstream in
#    the window when the main unit is +, upstream when it is -). ov = bp trimmed off the partner's
#    facing end (how far its alignment overran into the main unit).
#
# Out: OUT/cons.masked.fa, OUT/cons.tsv (name length tail_start dust_bp; tail_start = length+1 if none)
#      OUT/hits.m8      raw ssearch36 hits
#      OUT/units.tsv    every kept unit: wid family k role unit strand cons_s cons_e tag win_s win_e bits pid
#                       role main / up / down (window orientation); tag h = starts <= 15 bp from its
#                       consensus 5' end, e = ends <= 15 bp before its A tail, mid = neither
#      OUT/junctions.tsv one row per window and side (5, 3); columns:
#        wid family side  main main_strand main_cons_s main_cons_e main_j   (main_j = main's consensus
#                                                                            coordinate at the junction)
#        partner rel      (partner family, same / opp strand as the main unit; "-" = none on that side,
#                          "nomain" = no unit covers the core)
#        p_j p_far p_tag  partner consensus coordinate at the junction (after trimming), at its far end
#        p_ws p_we        partner window coordinates (trimmed)
#        gap ov           bp between main and partner; ov = bp trimmed off the partner's facing end
#        free             bp of flank searched on that side without a partner (window end - main end);
#                         with a partner: the same up to the partner
#        clamp            bp of flank lost at the contig end on that side (from loci.tsv)
#        gapA gapT gapN   A / T share and masked-N count of the gap (window orientation)
#        gap_seq          the gap (window orientation) when 1..GAPSEQ bp, else "."
#        main_tail        bp of the copy's own A tail added to the main unit (side 3 row; 0 on side 5)
#        p_tail           bp of the partner's own A tail added to it (only when its 3' end faces the copy)
#      gap and free are counted between the units INCLUDING their own A tails: the gap is what
#      neither element explains (a linker, target site, or an extra A run).
set -euo pipefail
OUT=${1:?OUT_DIR}; CONS=$(readlink -f "${2:?CONSENSUSES.fa}"); T=${3:-8}
MIN_ALN=40; EVAL=1e-5; GAPSEQ=300
cd "$OUT"

# 1) mask consensus A tails and low complexity; record where each tail starts
dustmasker -in "$CONS" -outfmt fasta 2> /dev/null | seqkit seq -w 0 \
| gawk -v OFS='\t' '
    /^>/ { n = substr($1, 2); sub(/^lcl\|/, "", n); next }
    {
        u = toupper($0); L = length(u); a = 0; t = L + 1
        for (i = L; i >= 1; i--) {                       # walk in from the 3 end
            if (substr(u, i, 1) == "A") a++
            f = a / (L - i + 1)
            if (f < 0.5) break
            if (f >= 0.8 && substr(u, i, 1) == "A") t = i # most 5 start that keeps >= 80 % A
        }
        if (L - t + 1 < 6) t = L + 1                      # shorter than 6 bp: no tail
        s = $0; nd = gsub(/[a-z]/, "N", s)                # dustmasker lower-cases low complexity
        s = substr(s, 1, t - 1) sprintf("%*s", L - t + 1, ""); gsub(/ /, "N", s)
        print ">" n > "cons.masked.fa"; print s > "cons.masked.fa"
        print n, L, t, nd > "cons.tsv"
    }'

# 2) search
ssearch36 -m 8 -E $EVAL -z 11 -T "$T" cons.masked.fa windows.masked.fa > hits.m8 2> /dev/null

# 3) + 4) units and junctions
sort -t$'\t' -k2,2V -k12,12gr hits.m8 \
| gawk -F'\t' -v OFS='\t' -v MIN_ALN=$MIN_ALN -v GAPSEQ=$GAPSEQ '
function min(x, y) { return x < y ? x : y }
function max(x, y) { return x > y ? x : y }
# consensus coordinate at window position x, linear along the alignment lo..hi <-> qa..qb
function cpos(x) { return int(0.5 + (st == "+" ? qa + (x - lo) * (qb - qa) / max(1, hi - lo) \
                                             : qb - (x - lo) * (qb - qa) / max(1, hi - lo))) }
# length of the A run (base A, or T on the other strand) starting at x and walking by step: the
# longest stretch that is >= 80 % base and ends on base; gives up once it falls below 50 % after
# 10 bp; never reaches the position stop
function tailrun(s, x, step, base, stop,   n, a, c, best) {
    a = 0; best = 0
    for (n = 1; x >= 1 && x <= length(s) && x != stop; n++) {
        c = toupper(substr(s, x, 1)); if (c == base) a++
        if (n >= 10 && a / n < 0.5) break
        if (a / n >= 0.8 && c == base) best = n
        x += step
    }
    return best
}
function tag(f, a, b,   t) { t = (a <= 15 ? "h" : "") (b >= ctail[f] - 1 - 15 ? "e" : ""); return t == "" ? "mid" : t }
FILENAME == ARGV[1] { clen[$1] = $2; ctail[$1] = $3; next }                               # cons.tsv
FILENAME == ARGV[2] { if (FNR > 1) { W[++nw] = $1; fam[$1] = $3; cs[$1] = $8; ce[$1] = $9    # loci.tsv
                                     wl[$1] = $10; cl5[$1] = $11; cl3[$1] = $12 }; next }
FILENAME == ARGV[3] { if (/^>/) w = substr($1, 2); else seq[w] = $0; next }               # windows.fa
FILENAME == ARGV[4] { if (/^>/) w = substr($1, 2); else mseq[w] = $0; next }              # windows.masked.fa
{                                                                                          # hits, best first
    w = $2; q1 = $7; q2 = $8; s1 = $9; s2 = $10
    lo = min(s1, s2); hi = max(s1, s2); qa = min(q1, q2); qb = max(q1, q2)
    st = ((q1 > q2) != (s1 > s2)) ? "-" : "+"
    if (hi - lo + 1 < MIN_ALN) next
    # subtract every kept unit from lo..hi -> pieces PL[j]..PH[j]
    np = 1; PL[1] = lo; PH[1] = hi
    for (k = 1; k <= nk[w]; k++) {
        n2 = 0
        for (j = 1; j <= np; j++) {
            if (PH[j] < KL[w, k] || PL[j] > KH[w, k]) { n2++; L2[n2] = PL[j]; H2[n2] = PH[j]; continue }
            if (PL[j] < KL[w, k]) { n2++; L2[n2] = PL[j]; H2[n2] = KL[w, k] - 1 }
            if (PH[j] > KH[w, k]) { n2++; L2[n2] = KH[w, k] + 1; H2[n2] = PH[j] }
        }
        np = n2; for (j = 1; j <= np; j++) { PL[j] = L2[j]; PH[j] = H2[j] }
    }
    for (j = 1; j <= np; j++) {
        if (PH[j] - PL[j] + 1 < MIN_ALN) continue
        k = ++nk[w]
        KL[w, k] = PL[j]; KH[w, k] = PH[j]; KF[w, k] = $1; KS[w, k] = st; KBits[w, k] = $12; KP[w, k] = $3
        TL[w, k] = PL[j] - lo; TR[w, k] = hi - PH[j]            # bp trimmed off each window end
        a = cpos(PL[j]); b = cpos(PH[j]); KA[w, k] = min(a, b); KB[w, k] = max(a, b)
    }
}
END {
    print "wid", "family", "k", "role", "unit", "strand", "cons_s", "cons_e", "tag", "win_s", "win_e", "bits", "pid" > "units.tsv"
    print "wid", "family", "side", "main", "main_strand", "main_cons_s", "main_cons_e", "main_j", "partner", "rel",
          "p_j", "p_far", "p_tag", "p_ws", "p_we", "gap", "ov", "free", "clamp", "gapA", "gapT", "gapN", "gap_seq", "main_tail", "p_tail" > "junctions.tsv"
    for (i = 1; i <= nw; i++) {
        w = W[i]; m = 0
        for (k = 1; k <= nk[w]; k++)                      # kept in bitscore order: first core hit = main
            if (min(KH[w, k], ce[w]) - max(KL[w, k], cs[w]) + 1 >= 30) { m = k; break }
        if (!m) {
            for (k = 1; k <= nk[w]; k++) print w, fam[w], k, "other", KF[w, k], KS[w, k], KA[w, k], KB[w, k],
                tag(KF[w, k], KA[w, k], KB[w, k]), KL[w, k], KH[w, k], KBits[w, k], KP[w, k] > "units.tsv"
            for (sd = 5; sd >= 3; sd -= 2) print w, fam[w], sd, "-", ".", ".", ".", ".", "nomain", ".", ".", ".", ".",
                ".", ".", ".", ".", ".", (sd == 5 ? cl5[w] : cl3[w]), ".", ".", ".", ".", ".", "." > "junctions.tsv"
            continue
        }
        ML = KL[w, m]; MH = KH[w, m]; mst = KS[w, m]
        # the copy owns its own A tail: the consensus tail is masked, so the main unit stops where the
        # tail starts and a neighbouring alignment could overrun into it. When the main unit reaches
        # its consensus tail (tag e), extend it over the A run that follows (T run before it when -),
        # same >= 80 % rule as the consensus mask; neighbours are trimmed back from there.
        mt = 0
        if (tag(KF[w, m], KA[w, m], KB[w, m]) ~ /e/)
            mt = (mst == "+") ? tailrun(seq[w], MH + 1, 1, "A", -1) : tailrun(seq[w], ML - 1, -1, "T", -1)
        if (mst == "+") MH += mt; else ML -= mt
        up = dn = 0
        for (k = 1; k <= nk[w]; k++) {
            if (k != m && mt) {                           # trim what now lies inside the main unit + tail
                t = (mst == "+") ? MH - KL[w, k] + 1 : KH[w, k] - ML + 1
                if (t > 0 && KL[w, k] <= MH && KH[w, k] >= ML) {
                    sl = (KB[w, k] - KA[w, k]) / max(1, KH[w, k] - KL[w, k])
                    if (mst == "+") { KL[w, k] += t; TL[w, k] += t; if (KS[w, k] == "+") KA[w, k] = int(KA[w, k] + t * sl + 0.5); else KB[w, k] = int(KB[w, k] - t * sl + 0.5) }
                    else            { KH[w, k] -= t; TR[w, k] += t; if (KS[w, k] == "+") KB[w, k] = int(KB[w, k] - t * sl + 0.5); else KA[w, k] = int(KA[w, k] + t * sl + 0.5) }
                    if (KH[w, k] - KL[w, k] + 1 < MIN_ALN) { drop[w, k] = 1; continue }
                }
            }
            role = (k == m) ? "main" : (KL[w, k] < ML ? "up" : "down")
            if (role == "up"   && (!up || KH[w, k] > KH[w, up])) up = k
            if (role == "down" && (!dn || KL[w, k] < KL[w, dn])) dn = k
            print w, fam[w], k, role, KF[w, k], KS[w, k], KA[w, k], KB[w, k], tag(KF[w, k], KA[w, k], KB[w, k]),
                  KL[w, k], KH[w, k], KBits[w, k], KP[w, k] > "units.tsv"
        }
        # the two window sides; in the main units consensus orientation "up" is side 5 when main is +
        for (d = 0; d <= 1; d++) {                        # d = 0 upstream, 1 downstream (window)
            p = d ? dn : up
            sd = ((d == 1) == (mst == "+")) ? 3 : 5
            mj = ((d == 1) == (mst == "+")) ? KB[w, m] : KA[w, m]    # main end facing this side
            clamp = (d ? cl3[w] : cl5[w])                            # loci clamps are in window orientation
            if (!p) {
                print w, fam[w], sd, KF[w, m], mst, KA[w, m], KB[w, m], mj, "-", ".", ".", ".", ".", ".", ".", ".", ".",
                      (d ? wl[w] - MH : ML - 1), clamp, ".", ".", ".", ".", (sd == 3 ? mt : 0), "." > "junctions.tsv"
                continue
            }
            ps = KS[w, p]; pa = KA[w, p]; pb = KB[w, p]; pl = KL[w, p]; ph = KH[w, p]
            # the partner end facing the main unit: its window-left end when downstream, right end when upstream
            if (d) { ov = TL[w, p]; if (ps == "+") { pj = pa; pf = pb } else { pj = pb; pf = pa } }
            else   { ov = TR[w, p]; if (ps == "+") { pj = pb; pf = pa } else { pj = pa; pf = pb } }
            # the partner owns its own A tail too, when its 3 end (reached: tag e) faces the copy
            pt = 0
            if (tag(KF[w, p], pa, pb) ~ /e/) {
                if (!d && ps == "+") { pt = tailrun(seq[w], ph + 1, 1, "A", ML); ph += pt }
                if (d && ps == "-")  { pt = tailrun(seq[w], pl - 1, -1, "T", MH); pl -= pt }
            }
            gap = d ? pl - MH - 1 : ML - ph - 1
            gs = (gap > 0) ? substr(seq[w],  d ? MH + 1 : ph + 1, gap) : ""
            gm = (gap > 0) ? substr(mseq[w], d ? MH + 1 : ph + 1, gap) : ""
            nA = gsub(/[Aa]/, "&", gs); nT = gsub(/[Tt]/, "&", gs); nN = gsub(/N/, "&", gm)
            print w, fam[w], sd, KF[w, m], mst, KA[w, m], KB[w, m], mj, KF[w, p], (ps == mst ? "same" : "opp"),
                  pj, pf, tag(KF[w, p], pa, pb), pl, ph, gap, ov, gap, clamp,
                  (gap ? sprintf("%.2f", nA / gap) : "."), (gap ? sprintf("%.2f", nT / gap) : "."), nN,
                  (gap >= 1 && gap <= GAPSEQ ? gs : "."), (sd == 3 ? mt : 0), pt > "junctions.tsv"
        }
    }
}' cons.tsv loci.tsv <(seqkit seq -w 0 windows.fa) <(seqkit seq -w 0 windows.masked.fa) -
echo "fs3: $(grep -c . hits.m8) hits, $(($(wc -l < units.tsv) - 1)) units, $(gawk -F'\t' 'NR>1 && $9!="-" && $9!="nomain"' junctions.tsv | wc -l) partners -> $OUT/junctions.tsv" >&2
