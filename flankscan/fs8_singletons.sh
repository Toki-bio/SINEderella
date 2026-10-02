#!/usr/bin/env bash
# fs8_singletons.sh OUT_DIR FAMILY [CONTROL_FAMILY] [THREADS=8]
#
# Stage 8 of flankscan: are the "single" copies of a family (stage 4: no bank consensus within 200 bp)
# real standalone insertions, or parts of composites whose other part decayed, or longer elements?
# For a family that is mostly part of composites (rsi r3: 98 %) its few singles decide whether it
# exists as a monomer.
#
# The same three groups are measured, so every number has its reference:
#   S = single copies of FAMILY
#   B = FAMILY copies inside a composite (what a composite looks like)
#   A = single copies of CONTROL_FAMILY, a clean monomer (what standalone insertions look like)
# Each group: up to NALN copies (random, seed 42).
#
# Per group, the way ViewAlign does it (MSA-viewer script.js, TSD finder, auto mode):
#   1) the copies with FL bp of flank (lowercase; the unit and its own tail uppercase) are aligned, the
#      family consensus as the first row (shown, not scored);
#   2) the element's ends are found ONCE, from the alignment: column conservation (top-base share x
#      coverage, columns with < 35 % of copies ignored), the first / last window of W columns whose mean
#      reaches the threshold, trimmed inward to a conserved column (_findSineBoundaryColumns, auto).
#      If the singles are really a composite with a decayed head, the copies stay conserved upstream of
#      the unit and the 5' boundary moves there by itself - no per-copy boundary search;
#   3) one TSD search per copy at those ends (_findBestTsdInFlanks, ported exactly: 30 bases before the
#      5' end, 56 after the 3' end, TSD 4-20 bp, <= 20 % mismatch, the 5' copy <= 4 bp from the end).
#      His two corrections: the 3' end of a SINE is fuzzy (simple-motif tails), so the 3' copy may start
#      up to TSD_RSLACK (45) bp past the 3' end instead of 3; and the default minimum (4 bp) counts chance
#      motifs, so the minimum is calibrated: the smallest length found in <= 5 % of shuffled pairs
#      (5' side of one copy with the 3' side of the next, same ends).
# Plus, per copy: a relaxed partner search (every bank consensus vs RELAX_BP of each flank, E <= 1 per
# flank, >= 20 bp) and completeness (starts / ends within 15 bp of its consensus ends).
#
# Out: OUT/singletons/FAMILY/<group>.aln.fa  the group alignment (consensus row 1), to look at
#      (default flank in the alignment 100 bp, auto-extended; the TSD flanks and plate flanks are cut from the raw copy, 400 bp)
#      OUT/singletons/FAMILY/ends.tsv        group copies columns left right, median bp from the 5' end to
#                                            the unit start (positive = the element starts that far
#                                            upstream of the unit) and from the unit end to the 3' end
#      OUT/singletons/FAMILY/copies.tsv      group wid full up_hit down_hit tsd_len tsd slack3 off5 off3
#      OUT/singletons/FAMILY/summary.tsv     per group: copies, full %, partner 5' / 3' %, TSD %, shuffled %,
#                                            5' / 3' end offsets
#      OUT/singletons/FAMILY/tsd_calibration.tsv   group min_len pairs shuffled_hits shuffled_pct real_pct excess_points
#                                            (real_pct = share of the copies with a TSD of at least min_len; excess = real - shuffled,
#                                            the enrichment curve: a minimum chosen only to silence the shuffled pairs hides short real TSDs)
set -euo pipefail
OUT=${1:?OUT_DIR}; FAM=${2:?FAMILY}; CTRL=${3:-}; T=${4:-8}
FL=${FL:-100}; FULLFL=${FULLFL:-400}; FLD=${FLD:-100}; NALN=${NALN:-200}; RSLK=${TSD_RSLACK:-45}; ETRIM=${ENDTRIM:-0}; EXT3=${EXT3:-0}; EXT3_THR=${EXT3_THR:-0.7}; RELAX_BP=300; MINLEN=20
cd "$OUT"; D=singletons/$FAM; rm -rf "$D"; mkdir -p "$D"

# groups, NALN each at most (seed 42)
gawk -F'\t' -v OFS='\t' -v F="$FAM" -v C="$CTRL" 'NR > 1 {
        if ($3 == F && $4 == "single") print "S", $1
        else if ($3 == F && ($4 == "composite" || $4 == "homodimer")) print "B", $1
        else if (C != "" && $3 == C && $4 == "single") print "A", $1 }' copies.tsv \
  | gawk 'BEGIN { srand(42) } { print rand() "\t" $0 }' | sort -t$'\t' -k1,1g \
  | gawk -F'\t' -v OFS='\t' -v M=$NALN '++n[$2] <= M { print $2, $3 }' > "$D/groups.tsv"

# per copy: the unit (+ own tail) with FL bp flanks, the flanks for the partner search, completeness
seqkit seq -w 0 windows.fa > "$D/win.tmp"
gawk -F'\t' -v OFS='\t' -v D="$D" -v FL=$FL -v R=$RELAX_BP '
    FILENAME == ARGV[1] { g[$2] = $1; next }
    FILENAME == ARGV[2] { if ($4 == "main" && ($1 in g)) { ms[$1] = $10; me[$1] = $11; ca[$1] = $7; cb[$1] = $8; fam[$1] = $2 }; next }
    FILENAME == ARGV[3] { if (FNR > 1 && $3 == 3 && ($1 in g)) tl[$1] = ($24 == "." ? 0 : $24 + 0); next }
    FILENAME == ARGV[4] { tailst[$1] = $3; next }
    /^>/ { w = substr($1, 2); next }
    (w in g) && (w in ms) {
        s = toupper($0); a = ms[w]; b = me[w] + tl[w]              # windows are in copy orientation
        print g[w], w, a, b > (D "/unit.tsv")                     # the unit (+ own tail) in the window
        u0 = (a - R > 1 ? a - R : 1)
        print ">" w "|up" > (D "/flanks.fa");   print substr(s, u0, a - u0) > (D "/flanks.fa")
        print ">" w "|down" > (D "/flanks.fa"); print substr(s, b + 1, R) > (D "/flanks.fa")
        print g[w], w, (ca[w] <= 15 && cb[w] >= tailst[fam[w]] - 1 - 15) ? 1 : 0 > (D "/base.tsv")
    }' "$D/groups.tsv" units.tsv junctions.tsv cons.tsv "$D/win.tmp"

# one group's copies with fl bp flanks (lowercase), the unit uppercase
extract() {   # extract GROUP FL
    gawk -v G="$1" -v FL="$2" 'FILENAME == ARGV[1] { if ($1 == G) { a[$2] = $3; b[$2] = $4 }; next }
        /^>/ { w = substr($1, 2); next }
        (w in a) { s = toupper($0); l0 = (a[w] - FL > 1 ? a[w] - FL : 1)
                   print ">" w; print tolower(substr(s, l0, a[w] - l0)) substr(s, a[w], b[w] - a[w] + 1) tolower(substr(s, b[w] + 1, FL)) }' \
        "$D/unit.tsv" "$D/win.tmp"
}

# relaxed partner search, E per flank
ssearch36 -m 8 -E 1 -Z 1 -z 11 -T "$T" cons.masked.fa "$D/flanks.fa" 2> /dev/null \
  | gawk -F'\t' -v OFS='\t' -v M=$MINLEN '
      { lo = ($9 < $10 ? $9 : $10); hi = ($9 < $10 ? $10 : $9); if (hi - lo + 1 < M) next
        if (!($2 in be) || $11 < be[$2]) { be[$2] = $11; bq[$2] = $1 } }
      END { for (k in be) { split(k, p, "|"); print p[1], p[2], bq[k] } }' > "$D/relaxed.tsv"

# per group: align, ends from conservation, one TSD search per copy
REF=$(seqkit grep -p "$FAM" cons.masked.fa 2> /dev/null | seqkit seq -s -w 0 | tr -d 'N')
: > "$D/ends.tsv"; : > "$D/tsd.tsv"; : > "$D/tsd_calibration.tsv"
for G in S B A; do
    grep -q "^$G	" "$D/unit.tsv" || continue
    FLC=$FL
    while :; do                              # extend the flank while an end sits close to the flank end
        rm -f "$D/$G".*.tmp "$D/$G.plate.aln.fa" "$D/$G.proposed.fa"
        { printf ">REF_%s\n%s\n" "$FAM" "$REF"; extract "$G" "$FLC"; } > "$D/$G.in.fa"
        extract "$G" "$FULLFL" | seqkit seq -w 0 > "$D/$G.full.fa"       # the same copies with long raw flanks (not aligned)
        mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --preservecase --quiet --thread "$T" "$D/$G.in.fa" 2> /dev/null \
            | seqkit seq -w 0 > "$D/$G.aln.fa"
        gawk -v OFS='\t' -v G=$G -v RSLK=$RSLK -v ETRIM=$ETRIM -v EXT3=$EXT3 -v EXT3_THR=$EXT3_THR -v CAL="$D/$G.cal.tmp" -v ENDS="$D/$G.ends.tmp" -v PLATE="$D/$G.plate.aln.fa" -v PROP="$D/$G.proposed.fa" -v FLC=$FLC -v FLD=$FLD -v FULL="$D/$G.full.fa" '
            BEGIN { while ((getline ln < FULL) > 0) { if (ln ~ /^>/) fk = substr(ln, 2); else FS_[fk] = ln } }
            /^>/ { n++; h[n] = substr($1, 2); next } { s[n] = $0 }
            function isb(c) { c = toupper(c); return c ~ /^[ACGTN]$/ }
            # ---- ViewAlign _findBestTsdInFlanks (ported; 3 side slack RSLK) ----
            function mm(x, y,   i, m, c1, c2) { m = 0
                for (i = 1; i <= length(x); i++) { c1 = substr(x, i, 1); c2 = substr(y, i, 1)
                    if (c1 == "N" || c2 == "N") m += 0.5; else if (c1 != c2) m++ }
                return m }
            function best(up, down, mn,   L, cl, uo, us, ds, dq, dv, sc) {      # sets TL TS TDS
                TL = 0; TSC = -1e9
                for (L = mn; L <= 20; L++) {
                    if (length(up) < L || length(down) < L) continue
                    cl = length(up) - L
                    for (uo = 0; uo <= (cl < 4 ? cl : 4); uo++) {
                        us = substr(up, cl - uo + 1, L)
                        for (ds = 0; ds <= ((length(down) - L) < RSLK ? length(down) - L : RSLK); ds++) {
                            dq = substr(down, ds + 1, L); dv = mm(us, dq) / L
                            if (dv > 0.20) continue
                            sc = (1 - dv) * sqrt(L) - uo * 0.04 - ds * 0.08 - ds * 0.05
                            if (sc > TSC || (sc == TSC && L > TL)) { TSC = sc; TL = L; TS = us; TDS = ds } } } }
                return TL }
            END {
                L = length(s[1]); nc = n - 1                               # row 1 = consensus, not scored
                # ---- ViewAlign _columnConservationScores + auto-mode ends ----
                minc = int(nc * 0.35 + 0.999); if (minc < 3) minc = 3
                for (x = 1; x <= L; x++) { delete k; v = 0; top = 0
                    for (i = 2; i <= n; i++) { c = toupper(substr(s[i], x, 1)); if (c ~ /^[ACGT]$/) { k[c]++; v++ } }
                    for (c in k) if (k[c] > top) { top = k[c]; mj[x] = c }
                    sc[x] = (v >= minc) ? (top / v) * (v / nc) : 0 }
                W = int(L / 24); if (W < 8) W = 8; if (W > 16) W = 16
                thr = (nc < 8) ? 0.68 : 0.58; left = 0; right = 0
                for (st = 1; st + W - 1 <= L && !left; st++) { m = 0; for (x = st; x < st + W; x++) m += sc[x]; if (m / W >= thr) left = st }
                for (st = L - W + 1; st >= 1 && !right; st--) { m = 0; for (x = st; x < st + W; x++) m += sc[x]; if (m / W >= thr) right = st + W - 1 }
                if (!left || right < left) { print G, nc, L, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", FLC >> ENDS; exit }
                for (kk = 0; kk < W - 1 && left < right && sc[left] < thr; kk++) left++
                for (kk = 0; kk < W - 1 && right > left && sc[right] < thr; kk++) right--
                # per copy: 30 bases before the 5 end, 56 from the 3 end + 1; bp between each end and the unit
                for (i = 2; i <= n; i++) {
                    # the flanks are cut from the raw copy, not from the alignment: n5 / m3 = copy bases before the left column / up to the right column
                    fr = FS_[h[i]]; f5 = 0; while (f5 < length(fr) && substr(fr, f5 + 1, 1) ~ /[a-z]/) f5++
                    a5 = 0; for (x = 1; x <= L; x++) { c = substr(s[i], x, 1); if (c ~ /[A-Z]/) break; if (c ~ /[a-z]/) a5++ }
                    off[i] = f5 - a5; nb5[i] = 0; for (x = 1; x < left; x++) if (substr(s[i], x, 1) != "-") nb5[i]++
                    mb3[i] = 0; for (x = 1; x <= right; x++) if (substr(s[i], x, 1) != "-") mb3[i]++
                    if (ETRIM) {       # per-copy ends: the last / first base that agrees with the column majority and sits in a stretch of >= 6 agreeing of 8 bases
                        for (x = right; x > left; x--) { if (toupper(substr(s[i], x, 1)) != mj[x]) continue
                            mt = 0; tt = 0; for (y = x; y > left && tt < 8; y--) { d = toupper(substr(s[i], y, 1)); if (d ~ /^[ACGT]$/) { tt++; if (d == mj[y]) mt++ } }
                            if (tt >= 6 && mt >= 6) break }
                        mb3[i] = 0; for (y = 1; y <= x; y++) if (substr(s[i], y, 1) != "-") mb3[i]++
                        for (x = left; x < right; x++) { if (toupper(substr(s[i], x, 1)) != mj[x]) continue
                            mt = 0; tt = 0; for (y = x; y < right && tt < 8; y++) { d = toupper(substr(s[i], y, 1)); if (d ~ /^[ACGT]$/) { tt++; if (d == mj[y]) mt++ } }
                            if (tt >= 6 && mt >= 6) break }
                        nb5[i] = 0; for (y = 1; y < x; y++) if (substr(s[i], y, 1) != "-") nb5[i]++ }
                    p5 = off[i] + nb5[i]; up[i] = toupper(substr(fr, (p5 > 30 ? p5 - 29 : 1), (p5 > 30 ? 30 : p5)))
                    dn[i] = toupper(substr(fr, off[i] + mb3[i] + 1, RSLK + 31))
                    o5 = 0; seen = 0; for (x = left; x <= L && !seen; x++) { c = substr(s[i], x, 1); if (c ~ /[ACGTN]/) seen = 1; else if (c ~ /[acgtn]/) o5++ }
                    o3 = 0; seen = 0; for (x = right; x >= 1 && !seen; x--) { c = substr(s[i], x, 1); if (c ~ /[ACGTN]/) seen = 1; else if (c ~ /[acgtn]/) o3++ }
                    O5[i - 1] = o5; O3[i - 1] = o3 }
                # chance: shuffled pairs (5 side of copy i, 3 side of copy i + 1)
                MINL = 0
                for (mn = 4; mn <= 14; mn++) { hs = 0; ns = 0
                    for (i = 2; i <= n; i++) { j = (i < n) ? i + 1 : 2; if (j == i) continue; ns++; if (best(up[i], dn[j], mn)) hs++ }
                    hr = 0; for (i = 2; i <= n; i++) if (best(up[i], dn[i], mn)) hr++        # the same search on the real pairs
                    print G, mn, ns, hs, (ns ? sprintf("%.1f", 100 * hs / ns) : "-"), (nc ? sprintf("%.1f", 100 * hr / nc) : "-"),
                          (ns && nc ? sprintf("%.1f", 100 * hr / nc - 100 * hs / ns) : "-") >> CAL
                    if (mn == 4 && ns > 50 && hs == 0) { print "fs8: TSD detector found nothing on shuffled pairs - broken, stop" > "/dev/stderr"; exit 3 }
                    if (!MINL && ns && hs / ns <= 0.05) { MINL = mn; SHUF = 100 * hs / ns } }
                if (!MINL) { MINL = 14; SHUF = 0 }
                print G, "chosen_min_len", MINL >> CAL
                for (i = 2; i <= n; i++) { best(up[i], dn[i], MINL)
                    print G, h[i], TL, (TL ? TS : "-"), (TL ? TDS : "-"), O5[i - 1], O3[i - 1] }
                # column majority (all copies with a base): for the carrying test and the proposed consensus
                for (x = 1; x <= L; x++) { delete k; v = 0; top = 0; tb = "-"
                    for (i = 2; i <= n; i++) { c = toupper(substr(s[i], x, 1)); if (c ~ /^[ACGT]$/) { k[c]++; v++ } }
                    for (c in k) if (k[c] > top) { top = k[c]; tb = c }
                    maj[x] = tb; cov[x] = v / nc }
                # per copy: flank bases beyond each end (is the flank long enough?), and whether the copy carries
                # the stretch between its unit and a moved end (identity to the column majority >= 0.6)
                c5 = 0; c3 = 0; n5 = 0; n3 = 0
                for (i = 2; i <= n; i++) {
                    r5 = 0; for (x = 1; x < left; x++) if (isb(substr(s[i], x, 1))) r5++
                    r3 = 0; for (x = right + 1; x <= L; x++) if (isb(substr(s[i], x, 1))) r3++
                    R5[i - 1] = r5; R3[i - 1] = r3
                    m = 0; t = 0; for (x = left; x <= L; x++) { c = substr(s[i], x, 1); if (c ~ /[ACGTN]/) break; if (c ~ /[acgt]/) { t++; if (toupper(c) == maj[x]) m++ } }
                    if (t >= 10) { n5++; if (m / t >= 0.6) c5++ }
                    m = 0; t = 0; for (x = right; x >= 1; x--) { c = substr(s[i], x, 1); if (c ~ /[ACGTN]/) break; if (c ~ /[acgt]/) { t++; if (toupper(c) == maj[x]) m++ } }
                    if (t >= 10) { n3++; if (m / t >= 0.6) c3++ } }
                asort(R5); asort(R3)
                # the plate: element columns as aligned, flanks packed (FLD bases, no gap columns)
                pc = ""; for (x = left; x <= right; x++) pc = pc (cov[x] >= 0.5 ? maj[x] : "-")
                # EXT3: the proposed consensus is continued past the right end while the bases that the copies carry
                # there (the packed flank, the same bases the plate shows) agree: majority share >= EXT3_THR among the
                # copies that have a base, those are >= half of the copies; at most FLD - 1 bases. 0 = off (default).
                ext = ""
                if (EXT3) for (xk = 1; xk < FLD; xk++) { delete ek; ev = 0; et = 0; eb = ""
                    for (i = 2; i <= n; i++) { c = toupper(substr(FS_[h[i]], off[i] + mb3[i] + xk, 1)); if (c ~ /^[ACGT]$/) { ek[c]++; ev++ } }
                    for (c in ek) if (ek[c] > et) { et = ek[c]; eb = c }
                    if (ev < 0.5 * nc || et / ev < EXT3_THR) break
                    ext = ext eb }
                pc = pc ext
                pad = sprintf("%*s", FLD, ""); gsub(/ /, "-", pad)
                print ">" h[1] > PLATE; print pad substr(s[1], left, right - left + 1) pad > PLATE
                print ">proposed_" G > PLATE; print pad pc substr(pad, 1, FLD - length(ext)) > PLATE
                for (i = 2; i <= n; i++) {
                    fr = FS_[h[i]]; p5 = off[i] + nb5[i]
                    u = toupper(substr(fr, (p5 > FLD ? p5 - FLD + 1 : 1), (p5 > FLD ? FLD : p5)))
                    d = toupper(substr(fr, off[i] + mb3[i] + 1, FLD))
                    print ">" h[i] > PLATE; print substr(pad, 1, FLD - length(u)) u substr(s[i], left, right - left + 1) d substr(pad, 1, FLD - length(d)) > PLATE }
                gsub(/-/, "", pc); print ">proposed_" G > PROP; print pc > PROP
                asort(O5); asort(O3)
                print G, nc, L, left, right, O5[int((nc + 1) / 2)], O3[int((nc + 1) / 2)], MINL, sprintf("%.1f", SHUF),
                      R5[int((nc + 1) / 2)], R3[int((nc + 1) / 2)], (n5 ? sprintf("%.0f", 100 * c5 / n5) : "-"),
                      (n3 ? sprintf("%.0f", 100 * c3 / n3) : "-"), FLC >> ENDS
            }' "$D/$G.aln.fa" > "$D/$G.tsd.tmp"
        R5=$(cut -f10 "$D/$G.ends.tmp"); R3=$(cut -f11 "$D/$G.ends.tmp")
        if [[ "$R5" =~ ^[0-9]+$ && "$R3" =~ ^[0-9]+$ ]] && (( (R5 < 60 || R3 < 60) && FLC < 1000 )); then
            FLC=$(( FLC * 2 > 1000 ? 1000 : FLC * 2 )); echo "fs8: $FAM group $G: an end within 60 bp of the flank end - flank $FLC" >&2; continue
        fi
        break
    done
    cat "$D/$G.ends.tmp" >> "$D/ends.tsv"; cat "$D/$G.cal.tmp" >> "$D/tsd_calibration.tsv"; cat "$D/$G.tsd.tmp" >> "$D/tsd.tsv"
    rm -f "$D/$G".*.tmp
done

# per copy table and per group summary
gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { h[$1, $2] = $3; next }
    FILENAME == ARGV[2] { t[$2] = $3 OFS $4 OFS $5 OFS $6 OFS $7; next }
    ($2 in t) { print $1, $2, $3, (($2, "up") in h ? h[$2, "up"] : "-"), (($2, "down") in h ? h[$2, "down"] : "-"), t[$2] }' \
    "$D/relaxed.tsv" "$D/tsd.tsv" "$D/base.tsv" > "$D/copies.body"
{ printf "group\twid\tfull\tup_hit\tdown_hit\ttsd_len\ttsd\tslack3\toff5\toff3\n"; cat "$D/copies.body"; } > "$D/copies.tsv"
gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { e5[$1] = $6; e3[$1] = $7; ml[$1] = $8; sh[$1] = $9; k5[$1] = $12; k3[$1] = $13; fl[$1] = $14; next }
    FNR == 1 { next }
    { g = $1; n[g]++; fu[g] += $3; if ($4 != "-") u[g]++; if ($5 != "-") d[g]++; if ($6 > 0) t[g]++ }
    END { nm["S"] = "singles of the family"; nm["B"] = "family copies in composites"; nm["A"] = "singles of the control"
          print "group", "what", "copies", "full_pct", "partner5_pct", "partner3_pct", "tsd_pct", "tsd_shuffled_pct", "tsd_min_len",
                "median_bp_5end_before_unit", "median_bp_3end_after_unit", "carry5_pct", "carry3_pct", "flank_bp"
          for (g in n) printf "%s\t%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", g, nm[g], n[g], 100 * fu[g] / n[g],
                100 * u[g] / n[g], 100 * d[g] / n[g], 100 * t[g] / n[g], sh[g], ml[g], e5[g], e3[g], k5[g], k3[g], fl[g] }' \
    "$D/ends.tsv" "$D/copies.tsv" > "$D/summary.tsv"
rm -f "$D/copies.body" "$D"/*.in.fa "$D"/*.full.fa "$D/win.tmp"
echo "fs8: $FAM -> $OUT/$D/summary.tsv" >&2
