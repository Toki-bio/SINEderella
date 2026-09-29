#!/usr/bin/env bash
# fs2_trf.sh OUT_DIR
#
# Stage 2 of flankscan: tandem repeats in every window (stage 1), classified by WHERE they sit
# relative to the core - "is there a repeat" alone means nothing (every SINE has an A tail).
#
# TRF settings: trf windows.fa 2 5 7 80 10 20 2000 -h -ngs
#   match 2, mismatch 5, indel 7, PM 80, PI 10, min score 20, periods up to 2000 bp - from 2 bp
#   dimers to SINE-length units (a satellite made of SINE copies has a period near the SINE's length).
# Kept: copies >= 4 for periods < 50 bp (short periods: random (AT)n noise below that), copies >= 2
#   for periods >= 50 bp (two SINE-length units already make a dimer array).
#
# Classes (window coordinates, core = core_s..core_e from loci.tsv, TOL = 15 bp):
#   tail       starts within 30 bp before / TOL after the core's 3' end and ends no earlier than TOL
#              before it: the SINE's own tail ((A)n, (TA)n, (CA)n ...). The assigned core usually
#              already contains the A tail, so a tail that ENDS at the core end is still a tail, not
#              "core"; one that runs past the end is a candidate boundary extension, not contamination
#   head       the mirror image at the core's 5' end ((A)n head = the copy sits in a host's A tail)
#   satellite  starts >= 50 bp before the core and ends >= 50 bp after it: the copy is a unit of
#              (or sits inside) a tandem array - flag the copy
#   core       inside the core (internal repeat of the element)
#   partial    overlaps the core in any other way
#   flank5 / flank3   entirely in a flank: background, masked with N for the partner search
#
# Out: OUT/trf.tsv            wid family class start end period copies pct_match score motif
#      OUT/trf_summary.tsv    per family: copies, % with each class, top tail motifs
#      OUT/windows.masked.fa  windows with flank5 / flank3 repeats replaced by N
set -euo pipefail
OUT=${1:?OUT_DIR}
cd "$OUT"
trf windows.fa 2 5 7 80 10 20 2000 -h -ngs > trf.raw 2> /dev/null || true   # trf's exit code is not an error code

gawk -F'\t' -v OFS='\t' '
    FNR == NR { if (FNR > 1) { fam[$1]=$3; cs[$1]=$8; ce[$1]=$9 }; next }      # loci.tsv
    /^@/ { w = substr($1, 2); next }
    {
        split($0, f, " ")                   # start end period copies consize pmatch pindel score A C G T entropy motif seq ...
        a = f[1]+0; b = f[2]+0; per = f[3]+0; cop = f[4]+0
        if (cop < (per < 50 ? 4 : 2)) next
        s = cs[w]; e = ce[w]
        if      (a >= e - 30 && a <= e + 15 && b >= e - 15)            cls = "tail"
        else if (b <= s + 30 && b >= s - 15 && a <= s + 15)            cls = "head"
        else if (a <= s - 50 && b >= e + 50)                           cls = "satellite"
        else if (a >= s - 15 && b <= e + 15)                           cls = "core"
        else if (b < s - 15)                                           cls = "flank5"
        else if (a > e + 15)                                           cls = "flank3"
        else                                                           cls = "partial"
        print w, fam[w], cls, a, b, per, cop, f[6], f[8], (length(f[14]) > 40 ? substr(f[14],1,40) "..." : f[14])
    }' loci.tsv trf.raw | sort -k1,1V -k4,4n > trf.body
{ printf "wid\tfamily\tclass\tstart\tend\tperiod\tcopies\tpct_match\tscore\tmotif\n"; cat trf.body; } > trf.tsv
rm -f trf.body

# per family: share of copies with each class; the commonest tail motifs
gawk -F'\t' '
    FNR == NR { if (FNR > 1) { n[$3]++ }; next }                            # loci.tsv: copies per family
    FNR > 1 { k = $1 SUBSEP $3; if (!(k in seen)) { seen[k]; c[$2, $3]++ }
              if ($3 == "tail") tm[$2, (length($10) <= 12 ? $10 : "long:" $6)]++ }
    END {
        for (f in n) {
            line = f "\t" n[f]
            split("tail head satellite core partial flank5 flank3", C, " ")
            for (i = 1; i <= 7; i++) line = line "\t" sprintf("%.1f%%", 100 * c[f, C[i]] / n[f])
            # top 3 tail motifs
            delete best
            for (k in tm) { split(k, q, SUBSEP); if (q[1] == f) best[q[2]] = tm[k] }
            m = ""; for (r = 1; r <= 3; r++) { bk = ""; bv = 0; for (x in best) if (best[x] > bv) { bv = best[x]; bk = x }
                                               if (bk == "") break; m = m (m ? ", " : "") bk " " bv; delete best[bk] }
            print line "\t" m
        }
    }' loci.tsv trf.tsv | sort -t$'\t' -k2,2nr > trf_summary.body
{ printf "family\tcopies\ttail\thead\tsatellite\tcore\tpartial\tflank5\tflank3\ttop_tail_motifs\n"; cat trf_summary.body; } > trf_summary.tsv
rm -f trf_summary.body

# mask flank background repeats (N) for the partner search; everything else unchanged
gawk -F'\t' '
    FNR == NR { if (FNR > 1 && ($3 == "flank5" || $3 == "flank3")) { m[$1] = m[$1] " " $4 "-" $5 }; next }
    /^>/ { w = substr($0, 2); print; next }
    {
        s = $0
        n = split(m[w], iv, " ")
        for (i = 1; i <= n; i++) { split(iv[i], ab, "-"); L = ab[2] - ab[1] + 1
            s = substr(s, 1, ab[1] - 1) sprintf("%*s", L, "") substr(s, ab[2] + 1); }
        gsub(/ /, "N", s); print s
    }' trf.tsv <(seqkit seq -w 0 windows.fa) > windows.masked.fa
echo "fs2: $(($(wc -l < trf.tsv) - 1)) repeats kept; summary in $OUT/trf_summary.tsv" >&2
