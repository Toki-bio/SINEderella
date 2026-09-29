#!/usr/bin/env bash
# fs1_extract.sh RUN_DIR OUT_DIR [FLANK=1000]
#
# Stage 1 of flankscan: every assigned copy with FLANK bp on each side, extracted ONCE; all later
# stages (TRF, partner search, junctions) read these windows and never re-derive coordinates.
#
# In : RUN_DIR/results/assignment_full.tsv (Status "assigned"; the Sequence column is the locus as
#      "contig:start-end(strand)", BED coordinates: 0-based start, end exclusive), RUN_DIR/genome.clean.fa
# Out: OUT_DIR/loci.tsv    one row per copy (columns below)
#      OUT_DIR/windows.fa  one window per copy, in the COPY's orientation (5' flank, core, 3' flank);
#                          header ">wN" = the loci.tsv key
#
# loci.tsv columns
#   wid        window key (w1, w2, ...)
#   locus      the Sequence id as in assignment_full.tsv
#   family     assigned subfamily
#   contig start end strand       the core, BED coordinates
#   core_s core_e                 the core inside the window, 1-based, copy orientation
#   win_len                       window length
#   clamp5 clamp3                 bp of flank MISSING because the contig ends (0 = full FLANK);
#                                 a short flank is not "no partner found"
#   gap5 gap3                     distance (bp) to the nearest other assigned copy on the 5' / 3'
#                                 side (-1 = none on this contig); negative-overlap shows as 0
#   nb5 nb3                       family of that neighbouring copy (- = none)
# Neighbouring copies overlap each other's windows; gap5/gap3 let later stages count an adjacency
# once (from the copy with the lower wid) instead of twice.
set -euo pipefail
RUN=${1:?RUN_DIR}; OUT=${2:?OUT_DIR}; F=${3:-1000}
mkdir -p "$OUT"
G="$RUN/genome.clean.fa"
[[ -s "$G.fai" ]] || samtools faidx "$G"

# 1) assigned loci, sorted by contig and start: contig start end strand family locus
gawk -F'\t' 'NR > 1 && $5 == "assigned" {
        if (match($1, /^(.+):([0-9]+)-([0-9]+)\(([^)]*)\)$/, m)) {
            st = (m[4] == "-") ? "-" : "+"          # "+,-" loci are extracted as + by SINEderella
            print m[1] "\t" m[2] "\t" m[3] "\t" st "\t" $2 "\t" $1
        }
    }' "$RUN/results/assignment_full.tsv" | sort -k1,1 -k2,2n > "$OUT/loci.sorted.tmp"

# 2) neighbours, windows, clamps, core position - one pass per contig in sorted order
gawk -F'\t' -v F="$F" -v OUT="$OUT" '
    FNR == NR { clen[$1] = $2; next }                      # genome.clean.fa.fai: contig length
    { n++; c[n]=$1; s[n]=$2; e[n]=$3; st[n]=$4; fam[n]=$5; id[n]=$6 }   # loci, sorted
    END {
        print "wid\tlocus\tfamily\tcontig\tstart\tend\tstrand\tcore_s\tcore_e\twin_len\tclamp5\tclamp3\tgap5\tgap3\tnb5\tnb3" > (OUT "/loci.tsv")
        k = 0
        for (i = 1; i <= n; i++) {
            # left / right neighbour on the same contig (genomic)
            gl = -1; fl = "-"; gr = -1; fr = "-"
            if (i-1 in c && c[i-1] == c[i]) { gl = s[i] - e[i-1]; if (gl < 0) gl = 0; fl = fam[i-1] }
            if (i+1 in c && c[i+1] == c[i]) { gr = s[i+1] - e[i]; if (gr < 0) gr = 0; fr = fam[i+1] }
            if (!(c[i] in clen)) { print "fs1: contig " c[i] " is not in the genome index - stop" > "/dev/stderr"; exit 1 }
            L = clen[c[i]]
            ws = s[i] - F; if (ws < 0) ws = 0
            we = e[i] + F; if (we > L) we = L
            cl = F - (s[i] - ws); cr = F - (we - e[i])           # genomic left / right clamp
            if (st[i] == "+") { cs = s[i] - ws + 1; ce = e[i] - ws; c5 = cl; c3 = cr; g5 = gl; g3 = gr; n5 = fl; n3 = fr }
            else               { cs = we - e[i] + 1; ce = we - s[i]; c5 = cr; c3 = cl; g5 = gr; g3 = gl; n5 = fr; n3 = fl }
            k++
            print c[i] "\t" ws "\t" we "\tw" k "\t0\t" st[i] > (OUT "/windows.bed")
            print "w" k "\t" id[i] "\t" fam[i] "\t" c[i] "\t" s[i] "\t" e[i] "\t" st[i] "\t" cs "\t" ce "\t" we-ws "\t" c5 "\t" c3 "\t" g5 "\t" g3 "\t" n5 "\t" n3 > (OUT "/loci.tsv")
        }
    }
' "$G.fai" "$OUT/loci.sorted.tmp"
rm -f "$OUT/loci.sorted.tmp"

# 3) windows in copy orientation; header ">wN" (bedtools adds "(+)" / "(-)" - removed)
bedtools getfasta -fi "$G" -bed "$OUT/windows.bed" -s -nameOnly | sed '/^>/s/([+-])$//' > "$OUT/windows.fa"
echo "fs1: $(($(wc -l < "$OUT/loci.tsv") - 1)) copies, flank $F bp -> $OUT/windows.fa" >&2
