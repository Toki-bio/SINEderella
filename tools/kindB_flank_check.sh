#!/usr/bin/env bash
# kindB_flank_check.sh - are the regularly spaced runs of satellite_stage (kind B) real tandem arrays? Short pieces only.
#
#   kindB_flank_check.sh HITS.tsv GENOME.fa OUTDIR
#   env: FL=150 (flank length)  MIN_ID=85  BORDER_LO=55  MAX_UNITS=12  THREADS=8  SSEARCH=ssearch36  MINFL=60  PIECE=150  TARGET_MAX=40000
#
# HITS.tsv (tab separated, 0-based BED coordinates, hits of one run together and sorted by start): run_id contig start end
# Output OUTDIR/verdicts.tsv: run_id  units_tested  median_identity  frac_ge_min  verdict(ARRAY|COPIES|-)  method(flank|piece|-)
#
# Why: the earlier check cut the DNA between hit k and hit k+1 (2-24 kb) and aligned whole units against each other. No other step of
# SINEderella compares long sequences. In an array every hit sits at the same place of a repeating unit, so the ~150 bp right after
# hit k (and right before it) are copies of the 150 bp after hit k+1 / k+2 (lag 2 covers dimeric arrays). Between chance copies flanks
# are unrelated. Per hit: best score over the two flanks and the two lags, score = identity x min(aligned / shorter piece, 1); run
# verdict = median over hits >= MIN_ID (same rule as the unit check). Runs whose median falls between BORDER_LO and MIN_ID get a second
# opinion: a PIECE (150 bp) from the middle of stretch k is searched in the stretches k+1 and k+2 (150 bp against a few kb); same rule.
# Never more than ~150 bp against a stretch. Needs gawk, samtools, ssearch36.
set -euo pipefail
HITS=$(readlink -f "$1"); GENOME=$(readlink -f "$2"); mkdir -p "$3"; OUT=$(readlink -f "$3")
FL=${FL:-150}; MIN_ID=${MIN_ID:-85}; BORDER_LO=${BORDER_LO:-55}; MAX_UNITS=${MAX_UNITS:-12}; THREADS=${THREADS:-8}
SSEARCH=${SSEARCH:-ssearch36}; MINFL=${MINFL:-60}; PIECE=${PIECE:-150}; TARGET_MAX=${TARGET_MAX:-40000}; UMAX=${UMAX:-400}
mkdir -p "$OUT/m" "$OUT/p"
cd "$OUT" || exit 1
rm -f m/* p/*
[[ -s "$GENOME.fai" ]] || samtools faidx "$GENOME"

# 1. flank regions: for hit k of a run, F = the FL bp after it, B = the FL bp before it, cut where the neighbouring hit starts/ends
gawk -F'\t' -v OFS='\t' -v FL="$FL" -v MINFL="$MINFL" -v UMAX="$UMAX" -v MAXH=$((MAX_UNITS + 1)) '
    FNR == NR { len[$1] = $2; next }
    function flush(   k, s, e, ps, pe, ns) {
        for (k = 0; k < m; k++) {
            # after
            s = E[k] + 1; e = E[k] + FL
            if (k + 1 < m && S[k + 1] < e) e = S[k + 1]
            if (e > len[C]) e = len[C]
            if (e - s + 1 >= MINFL) { print C ":" s "-" e > "regions.txt"; print id "_" k "_F" > "names.txt" }
            # before
            s = S[k] - FL + 1; e = S[k]
            if (k > 0 && E[k - 1] + 1 > s) s = E[k - 1] + 1
            if (s < 1) s = 1
            if (e - s + 1 >= MINFL) { print C ":" s "-" e > "regions.txt"; print id "_" k "_B" > "names.txt" }
        }
        # a short unit (hit start to the next hit start, <= UMAX bp: the hits nearly touch, there is no flank to test) is itself a short piece
        for (k = 0; k + 1 < m; k++)
            if (S[k + 1] - S[k] >= 50 && S[k + 1] - S[k] <= UMAX) { print C ":" (S[k] + 1) "-" S[k + 1] > "regions.txt"; print id "_" k "_U" > "names.txt" }
        if (started) print id, m > "nhits.tsv"
    }
    { if (!started || ($1 "") != id) { flush(); started = 1; id = $1 ""; C = $2; m = 0 }
      if (m < MAXH) { S[m] = $3; E[m] = $4; m++ } }
    END { flush() }' "$GENOME.fai" "$HITS"
: > /dev/null
[[ -s regions.txt ]] || { echo "kindB_flank_check: no flank regions" >&2; : > verdicts.tsv; exit 0; }

# 2. cut them, name them, one FASTA per run
samtools faidx "$GENOME" -r regions.txt \
| gawk 'FNR == NR { n[FNR] = $0; next } /^>/ { i++; print ">" n[i]; next } { print }' names.txt - > flanks.fa
gawk '/^>/ { if (nm != "") print nm "\t" L; nm = substr($0, 2); L = 0; next } { L += length($0) } END { if (nm != "") print nm "\t" L }' flanks.fa > lens.tsv
gawk -v RS='>' 'NR > 1 { split($0, a, "\n"); split(a[1], p, "_"); f = "m/" p[1] ".fa"; printf ">%s", $0 >> f; close(f) }' flanks.fa

# 3. all-against-all inside each run (26 pieces of 150 bp), the lag-1/lag-2 pairs of one flank kind are read from it
# ssearch36 exits 0 when it finds nothing, but also when it cannot open a file (only "***" on stderr): both are checked. Until
# 2026-10-19 a failed search was "|| true": its run scored 0 and was verdicted COPIES, so a real array stayed in (GLM audit)
ss_check='rc=$?; if [ "$rc" -ne 0 ] || grep -q "\*\*\*" "$b.err"; then echo "kindB_flank_check: ssearch36 failed on $1 (exit $rc): $(head -c 300 "$b.err")" >&2; exit 255; fi'
( cd m && ls *.fa | xargs -P "$THREADS" -I{} sh -c 'b=$(basename "$1" .fa); "$2" -m 8 -E 10 -z -1 -T 1 "$1" "$1" > "$b.m8" 2> "$b.err"; '"$ss_check" _ {} "$SSEARCH" )

# 4. per-hit best score, per-run median
cat m/*.m8 2>/dev/null | gawk -F'\t' -v OFS='\t' -v MIN_ID="$MIN_ID" '
    FILENAME == "lens.tsv" { L[$1] = $2; next }
    { split($1, a, "_"); split($2, b, "_")
      if (a[1] != b[1] || a[3] != b[3]) next
      d = b[2] - a[2]; if (d < 0) { d = -d; k = b[2] } else k = a[2]
      if (d != 1 && d != 2) next
      mn = L[$1] < L[$2] ? L[$1] : L[$2]; c = $4 / mn; if (c > 1) c = 1
      sc = $3 * c; key = a[1] SUBSEP k
      if (!(key in best) || sc > best[key]) best[key] = sc }
    END { for (key in best) { split(key, z, SUBSEP); print z[1], z[2], best[key] > "flank_scores.tsv" } }' lens.tsv -
: >> flank_scores.tsv
gawk -F'\t' -v OFS='\t' -v MIN_ID="$MIN_ID" -v BORDER_LO="$BORDER_LO" '
    FILENAME == "names.txt" { split($0, a, "_"); has[a[1], a[2], a[3]] = 1; next }
    FILENAME == "nhits.tsv" { m[$1] = $2; ord[++nr] = $1; next }
    { sc[$1, $2] = $3 }
    function elig(r, k,   t) {      # could hit/unit k be compared with k+1 or k+2 at all (pieces long enough to cut)?
        for (t in KINDS) if (((r, k, t) in has) && (((r, k + 1, t) in has) || ((r, k + 2, t) in has))) return 1
        return 0 }
    END {
        KINDS["F"]; KINDS["B"]; KINDS["U"]
        for (i = 1; i <= nr; i++) { r = ord[i]; n = m[r] - 1; c0 = 0
            for (k = 0; k < n; k++) if (elig(r, k)) v[++c0] = ((r, k) in sc) ? sc[r, k] : 0
            if (c0 < 2) { print r, c0, "-", "-", "-", "-"; delete v; continue }
            c = asort(v, w); med = (c % 2) ? w[(c + 1) / 2] : (w[c / 2] + w[c / 2 + 1]) / 2
            ge = 0; for (k = 1; k <= c; k++) if (v[k] >= MIN_ID) ge++
            verdict = med >= MIN_ID ? "ARRAY" : (med >= BORDER_LO ? "BORDER" : "COPIES")
            printf "%s\t%d\t%.1f\t%.2f\t%s\tflank\n", r, c0, med, ge / c0, verdict
            delete v; delete w } }' names.txt nhits.tsv flank_scores.tsv > verdicts.flank.tsv

# 5. second opinion for the borderline runs: a 150 bp piece from the middle of stretch k against stretches k+1 and k+2
awk -F'\t' '$5 == "BORDER" || $5 == "ARRAY" { print $1 }' verdicts.flank.tsv > border.ids
if [[ -s border.ids ]]; then
    gawk -F'\t' -v OFS='\t' -v PIECE="$PIECE" -v TMAX="$TARGET_MAX" -v MAXH=$((MAX_UNITS + 1)) '
        FILENAME == "border.ids" { want[$1] = 1; next }
        FILENAME ~ /fai$/ { len[$1] = $2; next }
        function flush(   k, mid, s, e, st, c, w) {
            if (!started || !(id in want)) return
            for (k = 0; k + 1 < m; k++) st[k + 1] = S[k + 1] - S[k]
            c = asort(st, w)
            if (c < 1 || w[int((c + 1) / 2)] < 2 * PIECE) return      # short units: the units themselves were compared, no interior piece to cut
            for (k = 0; k + 1 < m; k++) {
                e = S[k + 1]; s = S[k]
                if (e - s >= 2 * PIECE) { mid = int((s + e) / 2); print C ":" (mid - PIECE / 2 + 1) "-" (mid + PIECE / 2) > "pregions.txt"; print id "_" k "_P" > "pnames.txt" }
                if (e - s >= 50) { if (e - s > TMAX) e = s + TMAX; print C ":" (s + 1) "-" e > "pregions.txt"; print id "_" k "_T" > "pnames.txt" }
            }
        }
        { if (!started || ($1 "") != id) { flush(); started = 1; id = $1 ""; C = $2; m = 0 }
          if (m < MAXH) { S[m] = $3; m++ } }
        END { flush() }' border.ids "$GENOME.fai" "$HITS"
    if [[ -s pregions.txt ]]; then
        samtools faidx "$GENOME" -r pregions.txt \
        | gawk 'FNR == NR { n[FNR] = $0; next } /^>/ { i++; print ">" n[i]; next } { print }' pnames.txt - > pieces.fa
        gawk -v RS='>' 'NR > 1 { split($0, a, "\n"); split(a[1], p, "_"); f = "p/" p[1] "." (p[3] == "P" ? "q" : "t") ".fa"; printf ">%s", $0 >> f; close(f) }' pieces.fa
        ( cd p && ls *.q.fa | xargs -P "$THREADS" -I{} sh -c 'b=$(basename "$1" .q.fa); [ -s "$b.t.fa" ] || exit 0; "$2" -m 8 -E 10 -z -1 -T 1 "$1" "$b.t.fa" > "$b.m8" 2> "$b.err"; '"$ss_check" _ {} "$SSEARCH" )
        cat p/*.m8 2>/dev/null | gawk -F'\t' -v OFS='\t' -v PIECE="$PIECE" '
            { split($1, a, "_"); split($2, b, "_"); d = b[2] - a[2]; if (a[1] != b[1] || (d != 1 && d != 2)) next
              c = $4 / PIECE; if (c > 1) c = 1; sc = $3 * c; key = a[1] SUBSEP a[2]
              if (!(key in best) || sc > best[key]) best[key] = sc }
            END { for (key in best) { split(key, z, SUBSEP); print z[1], z[2], best[key] } }' > piece_scores.tsv || true
    fi
    : >> piece_scores.tsv; : >> pnames.txt
    gawk -F'\t' -v OFS='\t' '
        FILENAME == "pnames.txt" { split($0, a, "_"); if (a[3] == "P") has[a[1], a[2]] = 1; next }
        FILENAME == "verdicts.flank.tsv" { nn[$1] = $2; next }
        { sc[$1, $2] = $3 }
        END { while ((getline r < "border.ids") > 0) { n = 0
                for (k = 0; k <= nn[r] - 2; k++) if ((r, k) in has) v[++n] = ((r, k) in sc) ? sc[r, k] : 0
                if (n == 0) { print r, 0, 0; continue }
                c = asort(v, w); med = (c % 2) ? w[(c + 1) / 2] : (w[c / 2] + w[c / 2 + 1]) / 2
                print r, med, n; delete v; delete w } }' pnames.txt verdicts.flank.tsv piece_scores.tsv > piece_median.tsv
fi
: >> piece_median.tsv
gawk -F'\t' -v OFS='\t' -v MIN_ID="$MIN_ID" '
    FILENAME == "piece_median.tsv" { pm[$1] = $2; pn[$1] = $3; next }
    { if ($5 == "ARRAY" || $5 == "BORDER") {
          if (($1 in pn) && pn[$1] > 0) {                 # the interior of the unit must repeat as well, not only the flanks
              if (pm[$1] >= MIN_ID) { $6 = ($5 == "ARRAY") ? "flank+piece" : "piece"; $5 = "ARRAY"; if ($3 + 0 < pm[$1]) $3 = sprintf("%.1f", pm[$1]) }
              else $5 = "COPIES"
          } else if ($5 == "BORDER") $5 = "COPIES"         # no piece could be cut (short units): the flank/unit score alone decides
      }
      print }' piece_median.tsv verdicts.flank.tsv > verdicts.tsv
echo "kindB_flank_check: $(wc -l < verdicts.tsv) runs; $(awk -F'\t' '$5 == "ARRAY"' verdicts.tsv | wc -l) ARRAY ($(awk -F'\t' '$6 == "piece"' verdicts.tsv | wc -l) only by the piece test)" >&2
