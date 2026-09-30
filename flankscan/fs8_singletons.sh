#!/usr/bin/env bash
# fs8_singletons.sh OUT_DIR FAMILY [CONTROL_FAMILY] [THREADS=8]
#
# Stage 8 of flankscan: are the "single" copies of a family (stage 4: no bank consensus within 200 bp)
# real standalone insertions, or composites whose partner the strict search missed, or longer elements?
# For a family that is mostly part of composites (rsi r3: 98 %) its few singles decide whether it
# exists as a monomer.
#
# Groups (each measurement is run on all three, so every number has its reference):
#   S = single copies of FAMILY
#   B = FAMILY copies inside a composite (a partner was found: what a composite looks like)
#   A = single copies of CONTROL_FAMILY (a clean monomer: what standalone insertions look like)
# 1) Relaxed partner search: every bank consensus (masked, stage 3) against RELAX_BP of 5' and of 3'
#    flank (next to the copy and its own tail), ssearch36 E <= 1 per flank (-Z 1), >= 20 bp.
# 2) Completeness: where the copy starts / ends in its own consensus (full = within 15 bp of both
#    ends, the 3' end measured before the tail).
# 3) Target-site duplications as END MARKERS, detected exactly as ViewAlign does (MSA-viewer script.js
#    _findBestTsdInFlanks): 30 bases before the left boundary, 30 + 20 + 6 after the right boundary,
#    TSD 4-20 bp, <= 20 % mismatch (N = half), the 5' copy ending <= 4 bp from the left boundary and
#    the 3' copy starting <= 3 bp from the right boundary, score (1 - div) * sqrt(len)
#    - 0.04 * upstream offset - 0.13 * downstream start; best pair per copy.
#    One change (his note, 2026-09-30): a SINE's 5' end is sharp (the head starts at a fixed base) but
#    its 3' end is not (tails of simple motifs vary in length, so "unit + own tail" is only an estimate).
#    The 5' copy keeps ViewAlign's <= 4 bp; the 3' copy may start up to TSD_RSLACK (default 25) bp past
#    the estimated 3' end instead of 3, same per-bp penalty; the calibration below uses the same rule.
#    The distance found is kept with each TSD (unit_tsd = SEQ@bp) to tune the slack per family.
#    If a family makes TSDs, the TSD pair shows where the insertion really starts and ends:
#      unit: the TSD search at the unit's own ends (its start; its end after its own tail);
#      scan: the left boundary moved outward 1 bp at a time up to WIDE bp (right boundary fixed), and
#            the right boundary likewise; the boundary with the best score is where the insertion ends
#            (off5 / off3 = bp beyond the unit; 0 = the unit is the whole insertion).
#    ViewAlign's defaults are too relaxed (a 4 bp minimum calls short chance motifs): the minimum length
#    is calibrated on shuffled pairs (5' side of one copy with the 3' side of the next in its group) as
#    the smallest at which <= 5 % of them still give a TSD (tsd_calibration.tsv; TSD_MIN=n overrides).
#    Only when a family's TSDs are well above the shuffled rate are the offsets read as element ends.
#
# Out: OUT/singletons/FAMILY/copies.tsv   group wid start end full up_hit up_dist down_hit down_dist
#                                         near_len near_tsd wide_len wide_tsd off5 off3 bits
#      OUT/singletons/FAMILY/summary.tsv  per group: copies, full %, relaxed partner 5'/3' %, near TSD %
#                                         (and shuffled %), wide TSD % (and shuffled %), median off5/off3,
#                                         partners found
#      OUT/singletons/FAMILY/offsets.tsv  group off5 off3 count - where the wide TSDs put the element ends
#      OUT/singletons/FAMILY/plate.aln.fa the best single copies, flanks packed lowercase, FAMILY's
#                                         consensus + every stage-5 candidate containing it on top
set -euo pipefail
OUT=${1:?OUT_DIR}; FAM=${2:?FAMILY}; CTRL=${3:-}; T=${4:-8}
RELAX_BP=300; MINLEN=20; NEAR_MIN=8; WIDE=${WIDE:-150}; WIDE_MIN=8   # TSD_MIN (env) fixes the TSD minimum; default calibrated
 PLATE_N=100; PLATE_F=150; MAXG=1000    # copies per group (random, seed 42)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$OUT"; D=singletons/$FAM; mkdir -p "$D"

gawk -F'\t' -v OFS='\t' -v F="$FAM" -v C="$CTRL" 'NR > 1 {
        if ($3 == F && $4 == "single") print "S", $1
        else if ($3 == F && ($4 == "composite" || $4 == "homodimer")) print "B", $1
        else if (C != "" && $3 == C && $4 == "single") print "A", $1 }' copies.tsv   | gawk -F'	' -v M=$MAXG 'BEGIN { srand(42) } { print rand() "	" $0 }' | sort -t$'	' -k1,1g   | gawk -F'	' -v OFS='	' -v M=$MAXG '++n[$2] <= M { print $2, $3 }' > "$D/groups.tsv"   # <= MAXG per group, seed 42

seqkit seq -w 0 windows.fa > "$D/win.tmp"
# per copy: the unit and its own tail in window coordinates, flanks for the relaxed search, the
# 5' and 3' TSD windows (written once, reused for the real pairs and the shuffled pairs)
gawk -F'\t' -v OFS='\t' -v D="$D" -v R=$RELAX_BP -v WD=$WIDE '
    FILENAME == ARGV[1] { g[$2] = $1; next }
    FILENAME == ARGV[2] { if ($4 == "main" && ($1 in g)) { ms[$1] = $10; me[$1] = $11; mst[$1] = $6
                                 ca[$1] = $7; cb[$1] = $8; bits[$1] = $12; fam[$1] = $2 }; next }
    FILENAME == ARGV[3] { if (FNR > 1 && $3 == 3 && ($1 in g)) tl[$1] = ($24 == "." ? 0 : $24 + 0); next }
    FILENAME == ARGV[4] { tailst[$1] = $3; next }
    /^>/ { w = substr($1, 2); next }
    (w in g) && (w in ms) {
        s = toupper($0); a = ms[w]; b = me[w] + tl[w]              # windows are in copy orientation
        u0 = (a - R > 1 ? a - R : 1)
        print ">" w "|up" > (D "/flanks.fa");   print substr(s, u0, a - u0) > (D "/flanks.fa")
        print ">" w "|down" > (D "/flanks.fa"); print substr(s, b + 1, R) > (D "/flanks.fa")
        print g[w], w, a, b, s > (D "/tsdwin.tsv")          # unit start, end (with tail), window
        full = (ca[w] <= 15 && cb[w] >= tailst[fam[w]] - 1 - 15) ? 1 : 0
        print g[w], w, ca[w], cb[w], full, bits[w] > (D "/base.tsv")
    }' "$D/groups.tsv" units.tsv junctions.tsv cons.tsv "$D/win.tmp"

# composite hypotheses from stage 7 (hierarchy.tsv): if a single copy of FAMILY is really one part of a
# composite whose other part decayed, its TSD sits at a predictable distance beyond the unit:
#   FAMILY = part 2 -> 5' boundary moved out by part 1 length + gap; FAMILY = part 1 -> 3' out by gap + part 2
HYP=""
if [ -s hierarchy.tsv ]; then
    HYP=$(gawk -F'	' -v F="$FAM" 'BEGIN { sub(/_[0-9]+seqs$/, "", F) }
        NR > 1 && $12 == F { printf "%s5:%d:%s", (n++ ? "," : ""), $10 - $9 + 1 + $11, $2 }
        NR > 1 && $8 == F  { printf "%s3:%d:%s", (n++ ? "," : ""), $11 + $14 - $13 + 1, $2 }' hierarchy.tsv)
fi
echo "fs8: composite hypotheses for $FAM: ${HYP:-none}" >&2

# TSDs, ViewAlign's detector (MSA-viewer script.js _findBestTsdInFlanks), real and shuffled pairs
gawk -F'\t' -v OFS='\t' -v WD=$WIDE -v FIXED="${TSD_MIN:-0}" -v CAL="$D/tsd_calibration.tsv" -v HYP="$HYP" -v RSLK="${TSD_RSLACK:-25}" '
    function mm(x, y,   i, m, c1, c2) { m = 0
        for (i = 1; i <= length(x); i++) { c1 = substr(x, i, 1); c2 = substr(y, i, 1)
            if (c1 == "N" || c2 == "N") m += 0.5; else if (c1 != c2) m++ }
        return m }
    function best(up, down, mn,   L, cl, uo, us, ds, dq, dv, sc) {  # sets TL TS TD TSC; 0 if none
        TL = 0; TSC = -1e9
        for (L = mn; L <= 20; L++) {
            if (length(up) < L || length(down) < L) continue
            cl = length(up) - L
            for (uo = 0; uo <= (cl < 4 ? cl : 4); uo++) {
                us = substr(up, cl - uo + 1, L)
                for (ds = 0; ds <= ((length(down) - L) < RSLK ? length(down) - L : RSLK); ds++) {   # 3 side: RSLK bp slack (not RS: gawk record separator)
                    dq = substr(down, ds + 1, L); dv = mm(us, dq) / L
                    if (dv > 0.20) continue
                    sc = (1 - dv) * sqrt(L) - uo * 0.04 - ds * 0.08 - ds * 0.05
                    if (sc > TSC || (sc == TSC && L > TL)) { TSC = sc; TL = L; TS = us; TD = dv; TDS = ds }
                } } }
        return TL }
    function upw(s, l) { return substr(s, (l - 30 > 1 ? l - 30 : 1), (l - 30 > 1 ? 30 : l - 1)) }   # 30 before l
    function dnw(s, r) { return substr(s, r + 1, 56) }                                             # from r + 1
    { gr[NR] = $1; id[NR] = $2; A[NR] = $3 + 0; B[NR] = $4 + 0; S[NR] = $5 }
    function partner(r,   q) { q = r + 1; while (q <= NR && gr[q] != gr[r]) q++
                                if (q > NR) { q = 1; while (q < r && gr[q] != gr[r]) q++ }
                                return q }
    END {
        # the ViewAlign defaults (4-20 bp, <= 20 % mismatch) call short motifs by chance; the minimum
        # length is calibrated here: the smallest at which <= 5 % of shuffled pairs still give a TSD
        MINL = (FIXED > 0 ? FIXED : 0)
        for (mn = 4; mn <= 14; mn++) {
            ns = 0; hs = 0
            for (r = 1; r <= NR; r++) { q = partner(r); if (q == r) continue
                ns++; if (best(upw(S[r], A[r]), dnw(S[q], B[q]), mn)) hs++ }
            print mn, ns, hs, (ns ? sprintf("%.1f", 100 * hs / ns) : "-") > CAL
            if (!MINL && ns && hs / ns <= 0.05) MINL = mn
            if (mn == 4 && ns > 100 && hs == 0) {        # random pairs always give some 4 bp matches
                print "fs8: TSD detector found nothing on " ns " shuffled pairs at 4 bp - broken input or parsing, stop" > "/dev/stderr"
                exit 3 }
        }
        if (!MINL) MINL = 14
        print "chosen_min_len", MINL > CAL
        # the scan keeps the best of up to 2 x WD boundaries, so chance finds something almost every
        # time at the single-position minimum: its minimum is calibrated the same way, on the scan
        # itself run on shuffled pairs (up to 200), as the smallest length found in <= 5 % of them
        SCANL = 0; nsc = 0
        for (r = 1; r <= NR && nsc < 200; r++) { q = partner(r); if (q == r) continue; nsc++
            s = S[r] substr(S[q], B[q] + 1); a = A[r]; b = length(S[r]) - (length(S[r]) - B[r])   # 3 side from copy q
            s = substr(S[r], 1, B[r]) substr(S[q], B[q] + 1); b = B[r]; ml = 0
            for (k = 0; k <= WD && a - k > 31; k++) { if (best(upw(s, a - k), dnw(s, b), MINL) > ml) ml = TL }
            for (k = 1; k <= WD && b + k + 56 <= length(s); k++) { if (best(upw(s, a), dnw(s, b + k), MINL) > ml) ml = TL }
            smax[nsc] = ml }
        for (L = MINL; L <= 20 && !SCANL; L++) { h = 0; for (i = 1; i <= nsc; i++) if (smax[i] >= L) h++
            print "scan_min", L, nsc, h, (nsc ? sprintf("%.1f", 100 * h / nsc) : "-") > CAL
            if (nsc && h / nsc <= 0.05) SCANL = L }
        if (!SCANL) SCANL = 21
        print "chosen_scan_min_len", SCANL > CAL
        # hypothesis windows: +-10 bp around a predicted boundary = 21 positions; their own minimum,
        # calibrated on shuffled pairs with the window at 100 bp out on the 5 side
        nh = (HYP == "" ? 0 : split(HYP, hy, ","))
        HYPL = 0; nsc = 0
        for (r = 1; r <= NR && nsc < 500; r++) { q = partner(r); if (q == r || A[r] - 110 < 32) continue; nsc++
            s = substr(S[r], 1, B[r]) substr(S[q], B[q] + 1); ml = 0
            for (k = 90; k <= 110; k++) if (best(upw(s, A[r] - k), dnw(s, B[r]), MINL) > ml) ml = TL
            hmax[nsc] = ml }
        for (L = MINL; L <= 20 && !HYPL; L++) { h = 0; for (i = 1; i <= nsc; i++) if (hmax[i] >= L) h++
            print "hyp_min", L, nsc, h, (nsc ? sprintf("%.1f", 100 * h / nsc) : "-") > CAL
            if (nsc && h / nsc <= 0.05) HYPL = L }
        if (!HYPL) HYPL = 21
        print "chosen_hyp_min_len", HYPL > CAL
        for (r = 1; r <= NR; r++) {
            s = S[r]; a = A[r]; b = B[r]
            best(upw(s, a), dnw(s, b), MINL); ul = TL; uts = (TL ? TS "@" TDS : "-"); usc = TSC   # @ = bp past the 3 end
            bl = 0; bsc = usc; b5 = 0; b3 = 0; bts = uts                  # scan: left outward, then right outward
            for (k = 1; k <= WD && a - k > 31; k++) { best(upw(s, a - k), dnw(s, b), SCANL); if (TL && TSC > bsc) { bsc = TSC; b5 = k; b3 = 0; bl = TL; bts = TS } }
            for (k = 1; k <= WD && b + k + 56 <= length(s); k++) { best(upw(s, a), dnw(s, b + k), SCANL); if (TL && TSC > bsc) { bsc = TSC; b5 = 0; b3 = k; bl = TL; bts = TS } }
            if (!bl) { bl = ul }
            hs = ""
            for (h = 1; h <= nh; h++) { split(hy[h], hp, ":"); ml = 0
                for (k = hp[2] - 10; k <= hp[2] + 10; k++) {
                    if (hp[1] == "5") { if (a - k < 32) continue; x = best(upw(s, a - k), dnw(s, b), HYPL) }
                    else { if (b + k + 56 > length(s)) continue; x = best(upw(s, a), dnw(s, b + k), HYPL) }
                    if (x > ml) ml = x }
                if (ml) hs = hs (hs == "" ? "" : ";") hp[3] ":" ml }
            print "real", gr[r], id[r], ul, uts, bl, bts, b5, b3, (hs == "" ? "-" : hs)
            q = partner(r)
            if (q != r) { best(upw(s, a), dnw(S[q], B[q]), MINL); print "shuf", gr[r], id[r], TL, "-", "-", "-", "-", "-" }
        } }' "$D/tsdwin.tsv" > "$D/tsd.tsv"

# relaxed partner search, E per flank
ssearch36 -m 8 -E 1 -Z 1 -z 11 -T "$T" cons.masked.fa "$D/flanks.fa" 2> /dev/null \
  | gawk -F'\t' -v OFS='\t' -v M=$MINLEN -v R=$RELAX_BP '
      { lo = ($9 < $10 ? $9 : $10); hi = ($9 < $10 ? $10 : $9); if (hi - lo + 1 < M) next
        split($2, p, "|"); d = (p[2] == "up") ? R - hi : lo - 1
        if (!($2 in be) || $11 < be[$2]) { be[$2] = $11; bd[$2] = d; bq[$2] = $1 } }
      END { for (k in be) { split(k, p, "|"); print p[1], p[2], bq[k], bd[k] } }' > "$D/relaxed.tsv"

gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { h[$1, $2] = $3; dd[$1, $2] = $4; next }
    FILENAME == ARGV[2] { if ($1 == "real") { t[$3] = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 }; next }
    { u = (($2, "up") in h); v = (($2, "down") in h)
      print $1, $2, $3, $4, $5, (u ? h[$2, "up"] : "-"), (u ? dd[$2, "up"] : "-"),
            (v ? h[$2, "down"] : "-"), (v ? dd[$2, "down"] : "-"), t[$2], $6 }' \
    "$D/relaxed.tsv" "$D/tsd.tsv" "$D/base.tsv" > "$D/copies.body"
{ printf "group\twid\tstart\tend\tfull\tup_hit\tup_dist\tdown_hit\tdown_dist\tunit_tsd_len\tunit_tsd\tbest_tsd_len\tbest_tsd\toff5\toff3\thypotheses\tbits\n"
  cat "$D/copies.body"; } > "$D/copies.tsv"

# summary per group; offsets histogram of the wide TSDs (10 bp bins)
gawk -F'\t' -v OFS='\t' -v NM=$NEAR_MIN -v WM=$WIDE_MIN -v D="$D" '
    function med(arr, n,   i, j, t) { for (i = 2; i <= n; i++) { t = arr[i]; for (j = i - 1; j >= 1 && arr[j] > t; j--) arr[j + 1] = arr[j]; arr[j + 1] = t }
                                    return n ? arr[int((n + 1) / 2)] : "-" }
    FILENAME == ARGV[1] { if ($1 == "shuf") { sn[$2]++; if ($4 > 0) sne[$2]++ }; next }
    FNR == 1 { next }
    { g = $1; n[g]++; fu[g] += $5
      if ($6 != "-") { up[g]++; tu[g, $6]++ }; if ($8 != "-") { dn[g]++; td[g, $8]++ }
      if ($10 > 0) ne[g]++
      if ($12 > 0 && ($14 > 0 || $15 > 0)) wi[g]++
      o5[g, ++no[g]] = $14; o3[g, no[g]] = $15
      hb[g, int($14 / 10) * 10, int($15 / 10) * 10]++
      if ($16 != "-") { nx = split($16, hh, ";"); for (x = 1; x <= nx; x++) { split(hh[x], hq, ":"); hyc[g, hq[1]]++; hyn[hq[1]] } } }
    END {
        nm["S"] = "singles of the family"; nm["B"] = "family copies in composites"; nm["A"] = "singles of the control"
        print "group", "what", "copies", "full_pct", "partner5_pct", "partner3_pct", "tsd_at_unit_pct", "tsd_shuffled_pct",
              "tsd_moved_out_pct", "-", "median_off5", "median_off3", "partners5", "partners3"
        for (g in n) {
            u = ""; for (k in tu) { split(k, q, SUBSEP); if (q[1] == g) u = u q[2] ":" tu[k] " " }
            d = ""; for (k in td) { split(k, q, SUBSEP); if (q[1] == g) d = d q[2] ":" td[k] " " }
            delete A5; delete A3; for (i = 1; i <= no[g]; i++) { A5[i] = o5[g, i]; A3[i] = o3[g, i] }
            printf "%s\t%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%s\t%s\t%s\n", g, nm[g], n[g],
                100 * fu[g] / n[g], 100 * up[g] / n[g], 100 * dn[g] / n[g],
                100 * ne[g] / n[g], (sn[g] ? 100 * sne[g] / sn[g] : 0), 100 * wi[g] / n[g], "-",
                med(A5, no[g]), med(A3, no[g]), u, d }
        print "group", "hypothesis", "copies_with_tsd_there", "pct" > (D "/hypotheses.tsv")
        for (g in n) for (hname in hyn) printf "%s\t%s\t%d\t%.1f\n", g, hname, hyc[g, hname], 100 * hyc[g, hname] / n[g] > (D "/hypotheses.tsv")
        print "group", "off5_bin", "off3_bin", "copies" > (D "/offsets.tsv")
        for (k in hb) { split(k, q, SUBSEP); print q[1], q[2], q[3], hb[k] > (D "/offsets.tsv") }
    }' "$D/tsd.tsv" "$D/copies.tsv" > "$D/summary.tsv"

# the plate
gawk -F'\t' '$1 == "S"' "$D/copies.body" | sort -t$'\t' -k17,17gr | head -$PLATE_N | cut -f2 > "$D/plate.ids"
gawk -v F=$PLATE_F 'FILENAME == ARGV[1] { want[$1]; next }
     FILENAME == ARGV[2] { if ($4 == "main" && ($1 in want)) { a[$1] = $10; b[$1] = $11 }; next }
     /^>/ { w = substr($1, 2); next }
     (w in a) { s = $0; lo = (a[w] - F > 1 ? a[w] - F : 1)
                print ">" w; print tolower(substr(s, lo, a[w] - lo)) toupper(substr(s, a[w], b[w] - a[w] + 1)) tolower(substr(s, b[w] + 1, F)) }' \
     "$D/plate.ids" units.tsv "$D/win.tmp" > "$D/plate.copies.fa"
{ seqkit grep -p "$FAM" cons.masked.fa 2> /dev/null
  if [ -s candidates.fa ]; then seqkit grep -r -p "(^|__)${FAM}(__|_P)" candidates.fa 2> /dev/null || true; fi; } > "$D/plate.cons.fa"
if [ -f "$HERE/../tools/compare_cons.py" ]; then
    python3 "$HERE/../tools/compare_cons.py" <(printf ">none\nN\n"; cat "$D/plate.copies.fa") "$D/plate.aln.fa" "$D/plate.cons.fa" > /dev/null 2>&1 \
        || echo "fs8: plate alignment failed" >&2
fi
rm -f "$D/win.tmp" "$D/copies.body"
echo "fs8: $FAM singletons -> $OUT/$D/summary.tsv" >&2
