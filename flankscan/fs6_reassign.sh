#!/usr/bin/env bash
# fs6_reassign.sh OUT_DIR RUN_DIR [THREADS=8] [ACCEPT=70]
#
# Stage 6 of flankscan: re-assign ONCE with the candidates of stage 5 added to the library, and
# measure what each candidate explains. The check done by hand for r10_groupB / r1_r3 / r5h_r6
# (Tal rsi/REFINEMENT.md §10, §12: accept when >= 70 % of its copies read as one full unit).
#
# Which copies: every copy of every family that is a part (upstream or downstream unit) of a kept
# candidate - not only the peak members: a copy called "single r10" may still be the left part of
# r10 + group B with a partner too weak to be seen (a MISASSIGNED copy).
# How: stage 3 again (fs3_partners.sh, unchanged) on those windows with library = the run's
# consensuses + candidates.fa; the new MAIN unit of a copy (best hit covering its core) is its
# re-assignment. A copy "reads as one full unit" when its main unit is the candidate and spans it
# from its 5' end to its tail (tag he).
# Copies are counted as ELEMENTS: a two-part element is two copies of the run (one per part), but
# both parts now have the same main unit - counted once (same contig, strand, main unit start in
# the genome within 10 bp).
#
# Out: OUT/reassign/            stage 3 files for these copies (units.tsv, junctions.tsv ...)
#      OUT/reassign.tsv         per candidate: name peak_elements now_main full pct_full verdict
#                               from_full from_part
#                               peak_elements = elements of its peak (+ mirrors); now_main = of those,
#                               main unit is now the candidate; full = ... and it reads as one full
#                               unit; verdict accept (pct_full >= ACCEPT) / check;
#                               from_full = ALL elements (peak members or not) that now read as one
#                               full candidate unit, by old family - the misassigned ones show up
#                               here; from_part = main is the candidate but only in part (a single
#                               TB copy scores about as well against TB__TB as against TB: not a
#                               finding, listed so it is visible)
#      OUT/reassign_moves.tsv   old_family new_main elements full in_peak   (every pair, incl. unchanged)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT=$(readlink -f "${1:?OUT_DIR}"); RUN=$(readlink -f "${2:?RUN_DIR}"); T=${3:-8}; ACCEPT=${4:-70}
RE="$OUT/reassign"; rm -rf "$RE"; mkdir -p "$RE"
cd "$OUT"
if [[ ! -s candidates.fa ]]; then
    printf "name\tpeak_elements\tnow_main\tfull\tpct_full\tverdict\tfrom\n" > reassign.tsv
    printf "old_family\tnew_main\telements\tfull\tin_peak\n" > reassign_moves.tsv
    echo "fs6: no candidates" >&2; exit 0
fi

# the copies: all copies of the parts' families
gawk -F'\t' 'FILENAME == ARGV[1] { if (FNR > 1 && $14 == "kept") { fam[$4]; fam[$5] }; next }
             FNR == 1 || ($3 in fam)' candidates.tsv loci.tsv > "$RE/loci.tsv"
tail -n +2 "$RE/loci.tsv" | cut -f1 > "$RE/wids.txt"
seqkit grep -f "$RE/wids.txt" windows.fa        > "$RE/windows.fa"        2> /dev/null
seqkit grep -f "$RE/wids.txt" windows.masked.fa > "$RE/windows.masked.fa" 2> /dev/null
gawk 'FNR == NR { w[$1]; next } ($4 in w)' "$RE/wids.txt" windows.bed > "$RE/windows.bed"
cat "$RUN/consensuses.clean.fa" candidates.fa > "$RE/library.fa"
echo "fs6: $(wc -l < "$RE/wids.txt") copies of $(tail -n +2 candidates.tsv | gawk -F'\t' '$14=="kept"{f[$4];f[$5]} END{print length(f)}') families vs $(grep -c '^>' "$RE/library.fa") consensuses" >&2
bash "$HERE/fs3_partners.sh" "$RE" "$RE/library.fa" "$T"

gawk -F'\t' -v OFS='\t' -v ACCEPT=$ACCEPT -v RE="$RE" '
FILENAME == ARGV[1] { if (FNR == 1) next                                          # candidates.tsv
                      to = ($14 == "kept") ? $1 : (($14 ~ /^same_as:/) ? substr($14, 9) : "")   # same_as:NAME -> NAME
                      if ($14 == "kept") { C[++nc] = $1; cand[$1] = 1 }
                      if (to == "") next                                            # extended:CHAIN - its copies belong to the chain
                      pk[$2] = to; n = split($7, m, ","); for (i = 1; i <= n; i++) if (m[i] != "-") pk[m[i]] = to
                      next }
FILENAME == ARGV[2] { if (FNR > 1 && ($1 in pk)) inpk[$2] = pk[$1]; next }       # members.tsv: wid -> candidate
FILENAME == ARGV[3] { ctg[$4] = $1; ws[$4] = $2; we[$4] = $3; wst[$4] = $6; next }  # windows.bed
FILENAME == ARGV[4] { if (FNR > 1) { W[++nw] = $1; old[$1] = $3 }; next }          # reassign/loci.tsv
FNR > 1 && $4 == "main" {                                                          # reassign/units.tsv
    w = $1
    g = (wst[w] == "+") ? ws[w] + $10 - 1 : we[w] - $11                            # main unit start in the genome
    new[w] = $5; full[w] = ($9 == "he"); key[w] = ctg[w] SUBSEP (wst[w] == $6 ? "+" : "-") SUBSEP $5 SUBSEP int(g / 10)
}
END {
    for (i = 1; i <= nw; i++) {
        w = W[i]; nm = (w in new) ? new[w] : "nomain"
        # one element = one key; a key seen from a neighbouring 10-bp bin is the same element too
        k = key[w]; split(k, q, SUBSEP)
        if (w in new) { if ((q[1] SUBSEP q[2] SUBSEP q[3] SUBSEP q[4] - 1) in seen || (q[1] SUBSEP q[2] SUBSEP q[3] SUBSEP q[4] + 1) in seen || k in seen) continue
                        seen[k] = 1 }
        mv[old[w], nm]++; if (full[w]) mf[old[w], nm]++
        if (w in inpk) { mp[old[w], nm]++; pe[inpk[w]]++; if (nm == inpk[w]) { pm[nm]++; if (full[w]) pf[nm]++ } }
        if (nm in cand) { if (full[w]) fr[nm, old[w]]++; else fp[nm, old[w]]++ }
    }
    print "old_family", "new_main", "elements", "full", "in_peak" > RE "/../reassign_moves.tsv"; close(RE "/../reassign_moves.tsv")
    for (k in mv) { split(k, q, SUBSEP); print q[1], q[2], mv[k], mf[k] + 0, mp[k] + 0 | "sort -t\"\t\" -k1,1 -k3,3nr >> " RE "/../reassign_moves.tsv" }
    print "name", "peak_elements", "now_main", "full", "pct_full", "verdict", "from_full", "from_part"
    for (i = 1; i <= nc; i++) {
        c = C[i]; f = ""; h = ""
        for (k in fr) { split(k, q, SUBSEP); if (q[1] == c) f = f (f ? ", " : "") q[2] " " fr[k] }
        for (k in fp) { split(k, q, SUBSEP); if (q[1] == c) h = h (h ? ", " : "") q[2] " " fp[k] }
        p = pe[c] ? 100 * pf[c] / pe[c] : 0
        print c, pe[c] + 0, pm[c] + 0, pf[c] + 0, sprintf("%.1f", p), (p >= ACCEPT ? "accept" : "check"), (f ? f : "-"), (h ? h : "-")
    }
}' candidates.tsv members.tsv windows.bed "$RE/loci.tsv" "$RE/units.tsv" > reassign.tsv
echo "fs6: per candidate in $OUT/reassign.tsv, family moves in $OUT/reassign_moves.tsv" >&2
