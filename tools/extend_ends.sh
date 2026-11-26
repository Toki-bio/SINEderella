#!/usr/bin/env bash
# extend_ends.sh - extend merged loci whose element goes on past the end of the hits (step 1, before extraction).
#
#   extend_ends.sh MERGED_BED LABELED_BED GENOME OUT_BED OUT_DIR [--win 100] [--k 10] [--step 10] [--min-n 20]
#                  [--min-frac 0.02] [--share 0.2] [--min-shared 10] [--bridge 2] [--flag-frac 0.2] [--threads N]
#
# Bash/awk rewrite (2026-11-26) of tools/extend_ends.py, same outputs (docs/END_EXTENSION.md). --threads defaults to $THREADS,
# else nproc. Needs bedtools, sort, split, xargs and an awk (mawk when present, else awk; EXTEND_AWK overrides; gawk and
# mawk give the same result).
#
# Why. A locus is the union of the hits of all bank consensuses, so it ends where the longest matching consensus stops. When
# the genome's element is longer than every bank consensus that hit it, the copy is cut there (hla shorties, 2026-10-29).
#
# How. Loci are grouped by the label of their highest-scoring hit (first such hit in bedtools intersect order on ties). For
# every locus the WIN bp past each end are read in the element's orientation (never past the neighbouring locus or the contig
# end); 5' windows are reversed (not complemented) so that index 0 is next to the locus end. Per group and side, every K-mer
# is counted once per window it occurs in; K-mers with a non-ACGT base, made only of A/T, or with >= K-3 copies of one base
# are not counted. A K-mer is shared when it occurs in at least max(MIN_N, ceil(MIN_FRAC x windows of the group and side))
# windows; fewer than MIN_SHARED shared K-mers: the group/side extends nothing. Each end is walked outward in STEP-bp steps
# (a step = the K-mers starting in it); a step is element when it has >= 3 counted K-mers and >= SHARE of them are shared;
# the walk stops at the (BRIDGE+1)-th miss in a row; the new end is the far end of the last shared K-mer (at most WIN). A
# walk ending at >= WIN-STEP is left as it is (at_cap_not_extended); an extended end gets the A run that follows (>= 8 A per
# 10 bp, then base by base, never past WIN). Two neighbours that both grew into the gap between them meet in its middle.
#
# Parallel design. bedtools intersect runs per bin of contigs (THREADS bins). The windows are read once with bedtools
# getfasta (THREADS chunks of the window list, sorted by group, side and position) into pieces per group, side and chunk,
# kept in OUT_DIR/extend_*/pieces until the end (about 100 bytes per window). Stage 1: the group-wide statistic (windows per
# K-mer) needs every window of a group and side, so each job reads all windows of one group and side but counts only the
# K-mers with one two-base prefix (16 jobs per group and side); the shared K-mers of the 16 prefixes are pooled and the
# MIN_SHARED rule is applied to the pooled count. Stage 2 (the walk, per window) runs per piece. No group-wide number is
# computed from a part of the group. The new coordinates are one awk pass in (contig, start) order.
#
# Out: OUT_BED (same columns as MERGED_BED, ends moved), OUT_DIR/ends_extension.tsv, OUT_DIR/loci_extended.tsv, one line
# on stdout: EXTEND loci=.. extended5=.. extended3=.. flagged=.. -> OUT_DIR/ends_extension.tsv
set -euo pipefail
export LC_ALL=C

usage() { awk 'NR>=2 && /^#/{sub(/^# ?/, ""); print; next} NR>=2{exit}' "$0" >&2; exit 2; }

WIN=100; K=10; STEP=10; MIN_N=20; MIN_FRAC=0.02; SHARE=0.2; MIN_SHARED=10; BRIDGE=2; FLAG_FRAC=0.2
THR="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
POS=()
while [[ $# -gt 0 ]]; do
    a="$1"; v=""
    case "$a" in
        --*=*) v="${a#*=}"; a="${a%%=*}" ;;
        --*) [[ $# -ge 2 ]] || { echo "extend_ends.sh: $a needs a value" >&2; exit 2; }; v="$2"; shift ;;
        -h) usage ;;
        *) POS+=("$a"); shift; continue ;;
    esac
    shift
    case "$a" in
        --win) WIN="$v" ;; --k) K="$v" ;; --step) STEP="$v" ;; --min-n) MIN_N="$v" ;; --min-frac) MIN_FRAC="$v" ;;
        --share) SHARE="$v" ;; --min-shared) MIN_SHARED="$v" ;; --bridge) BRIDGE="$v" ;; --flag-frac) FLAG_FRAC="$v" ;;
        --threads) THR="$v" ;; --help) usage ;;
        *) echo "extend_ends.sh: unknown option $a" >&2; exit 2 ;;
    esac
done
[[ ${#POS[@]} -eq 5 ]] || usage
MERGED="${POS[0]}"; LABELED="${POS[1]}"; GENOME="${POS[2]}"; OUT_BED="${POS[3]}"; OUT_DIR="${POS[4]}"
for x in "$WIN" "$K" "$STEP" "$MIN_N" "$MIN_SHARED" "$BRIDGE" "$THR"; do
    [[ "$x" =~ ^[0-9]+$ ]] || { echo "extend_ends.sh: not a whole number: $x" >&2; exit 2; }
done
[[ "$THR" -ge 1 ]] || THR=1
AWK="${EXTEND_AWK:-$(command -v mawk || command -v awk)}"
[[ -f "$GENOME.fai" ]] || { echo "extend_ends.sh: no index $GENOME.fai" >&2; exit 1; }
mkdir -p "$OUT_DIR"

NLOCI=$("$AWK" '/[^[:space:]]/{n++} END{print n+0}' "$MERGED")
if [[ "$NLOCI" -eq 0 ]]; then
    : > "$OUT_BED"; echo "EXTEND loci=0 extended5=0 extended3=0 flagged=0"; exit 0
fi
TMP=$(mktemp -d "$OUT_DIR/extend_XXXXXXXX")
tick() { [[ -z "${EXTEND_VERBOSE:-}" ]] || echo "[extend_ends $(date +%T.%N | cut -c1-12)] $*" >&2; }
tick start
export AWK TMP GENOME K STEP WIN SHARE BRIDGE MIN_N MIN_FRAC MIN_SHARED

# ---- awk programs ------------------------------------------------------------------------------------------------------
# ok(): K-mer counted? (memo per process; the K-mer is upper case)
cat > "$TMP/ok.awk" <<'EOF'
function ok(km,   t, a, c, g, x, m) {
    if (km in OK) return OK[km]
    if (km ~ /[^ACGT]/) return OK[km] = 0
    t = km; a = gsub(/A/, "", t); c = gsub(/C/, "", t); g = gsub(/G/, "", t); x = K - a - c - g
    m = a; if (c > m) m = c; if (g > m) m = g; if (x > m) m = x
    return OK[km] = !((a + x == K) || (m >= K - 3))
}
EOF
# stage 1: windows of one group and side (lines i<TAB>seq), only the K-mers that start with the prefix CL (the K-mers of a
# window fall into disjoint prefix classes, so the per-window counts of the classes are independent) -> K-mers in >= thr
# windows (file SH) and one line on stdout: KEY n thr count
cat > "$TMP/count.awk" <<'EOF'
{
    u = toupper($2); L = length(u) - K + 1; n++
    j = 0
    while ((p = index(substr(u, j + 1), CL)) > 0) {
        j += p
        if (j > L) break
        km = substr(u, j, K)
        if (!(km in seen) && ok(km)) { seen[km]; cnt[km]++ }
    }
    delete seen
}
END {
    x = MIN_FRAC * n; thr = int(x); if (thr < x) thr++
    if (thr < MIN_N) thr = MIN_N
    ns = 0
    for (km in cnt) if (cnt[km] >= thr) { print km > SH; ns++ }
    close(SH)
    print KEY "\t" n "\t" thr "\t" ns
}
EOF
# stage 2: walk every window of a piece (lines i<TAB>seq) over the shared K-mers in SH; out: KEY<TAB>i<TAB>bp (-1 = at the cap)
cat > "$TMP/walk.awk" <<'EOF'
BEGIN { while ((getline km < SH) > 0) sh[km] = 1; close(SH) }
{
    u = toupper($2); len = length(u); L = len - K + 1
    last = 0; miss = 0
    for (s0 = 0; s0 < L; s0 += STEP) {
        hi = s0 + STEP; if (hi > L) hi = L
        nc = 0; ns = 0; li = -1
        for (j = s0; j < hi; j++) {
            km = substr(u, j + 1, K)
            if (ok(km)) { nc++; if (km in sh) { ns++; li = j } }
        }
        if (nc < 3 || ns / nc < SHARE) { if (++miss > BRIDGE) break; continue }
        miss = 0; last = li + K
    }
    b = (last < WIN) ? last : WIN
    if (b >= WIN - STEP) { print KEY "\t" $1 "\t-1"; next }
    if (b > 0) {
        p = b
        while (p + 10 <= len && p + 10 <= WIN) { t = substr(u, p + 1, 10); if (gsub(/A/, "", t) < 8) break; p += 10 }
        while (p < len && p < WIN && substr(u, p + 1, 1) == "A") p++
        b = p; if (b > len) b = len
        print KEY "\t" $1 "\t" b
    }
}
EOF

# ---- loci, groups, windows -------------------------------------------------------------------------------------------
# loci.bed: chrom start end i 0 strand (i = 0-based number of the non-blank line)
"$AWK" -F'\t' -v OFS='\t' '/[^[:space:]]/{print $1, $2, $3, n++, 0, $6}' "$MERGED" > "$TMP/loci.bed"
# label of the highest-scoring labeled hit on the locus strand (first in bedtools order when equal). The contigs are dealt
# into NB bins (largest first, each to the least loaded bin) and bedtools intersect runs per bin: a bin's hit file keeps the
# hits of its contigs in file order, so every locus meets its hits in the same order as with the whole file (bedtools keeps
# one interval tree per contig).
mkdir -p "$TMP/bins"
NB=$THR; [[ "$NB" -le 256 ]] || NB=256          # awk keeps 2 x NB files open (mawk: about 1000 at most)
"$AWK" -F'\t' '{c[$1]++} END {for (k in c) print c[k] "\t" k}' "$TMP/loci.bed" | sort -t$'\t' -k1,1nr -k2,2 \
  | "$AWK" -F'\t' -v N="$NB" '{m = 1; for (j = 2; j <= N; j++) if (load[j] < load[m]) m = j; load[m] += $1; print $2 "\t" m}' \
  > "$TMP/bins/contig_bin.tsv"
"$AWK" -F'\t' -v D="$TMP/bins" 'FNR == 1 {fi++} fi == 1 {bin[$1] = $2; next}
    fi == 2 {print > (D "/a." bin[$1]); next} ($1 in bin) {print > (D "/b." bin[$1])}' \
    "$TMP/bins/contig_bin.tsv" "$TMP/loci.bed" "$LABELED"
best1() {   # $1 = bin number
    [[ -f "$TMP/bins/b.$1" ]] || : > "$TMP/bins/b.$1"
    bedtools intersect -s -wa -wb -a "$TMP/bins/a.$1" -b "$TMP/bins/b.$1" \
      | "$AWK" -F'\t' '{ sc = ($11 == "" || $11 == ".") ? 0 : $11 + 0
                         if (!($4 in b) || sc > b[$4]) { b[$4] = sc; g[$4] = $10 } }
                       END { for (i in g) print i "\t" g[i] }' > "$TMP/bins/lab.$1"
}
export -f best1
"$AWK" -F'\t' '{print $2}' "$TMP/bins/contig_bin.tsv" | sort -u -n | tr '\n' '\0' | xargs -0 -n 1 -P "$THR" bash -c 'best1 "$1"' _
find "$TMP/bins" -type f -name 'lab.*' -exec cat {} + > "$TMP/lab.tsv"
rm -rf "$TMP/bins"
tick intersect done
# group ids in sorted name order ("-" = no labeled hit)
{ "$AWK" -F'\t' '{print $2}' "$TMP/lab.tsv"
  [[ "$(wc -l < "$TMP/lab.tsv")" -eq "$NLOCI" ]] || echo "-"; } | sort -u | "$AWK" '{print NR "\t" $0}' > "$TMP/groups.tsv"
# sorted.bed: chrom start end i gid strand, in (chrom, start, i) order
"$AWK" -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] {gid[$2] = $1; next} FILENAME == ARGV[2] {lab[$1] = $2; next}
    {print $1, $2, $3, $4, gid[($4 in lab) ? lab[$4] : "-"], $6}' "$TMP/groups.tsv" "$TMP/lab.tsv" "$TMP/loci.bed" \
  | sort -t$'\t' -k1,1 -k2,2n -k4,4n > "$TMP/sorted.bed"
tick loci sorted, groups
# windows: room up to the neighbouring locus (any strand) and the contig end; name gid_side_i; sorted by gid, side, position
"$AWK" -F'\t' -v OFS='\t' -v W="$WIN" '
FILENAME == ARGV[1] { size[$1] = $2; next }
function emit(   lr, rr, r3, r5) {
    if (!(c in size)) { print "extend_ends.sh: contig " c " is not in the genome index" > "/dev/stderr"; bad = 1; exit 1 }
    lr = s - (pc == c ? pe : 0); if (lr > W) lr = W
    rr = (nc == c ? nx : size[c]) - e; if (rr > W) rr = W
    if (st == "+") { r3 = rr; r5 = lr } else { r3 = lr; r5 = rr }
    if (r3 > 0) print c, (st == "+" ? e : s - r3), (st == "+" ? e + r3 : s), g "_3_" i, 0, st, g, 3
    if (r5 > 0) print c, (st == "+" ? s - r5 : e), (st == "+" ? s : e + r5), g "_5_" i, 0, st, g, 5
}
{   nc = $1; nx = $2 + 0
    if (FNR > 1) { emit(); pc = c; pe = e }
    c = $1; s = $2 + 0; e = $3 + 0; i = $4; g = $5; st = $6 }
END { if (bad) exit 1; nc = "\001"; emit() }' "$GENOME.fai" "$TMP/sorted.bed" \
  | sort -t$'\t' -k7,7n -k8,8n -k1,1 -k2,2n > "$TMP/win.bed"

tick windows listed
# ---- read the windows (THR chunks in parallel) into pieces p.GID.SIDE.CHUNK (lines i<TAB>seq, 5' windows reversed) ----
mkdir -p "$TMP/chunks" "$TMP/pieces" "$TMP/res"
split -n l/"$THR" -d -a 4 "$TMP/win.bed" "$TMP/chunks/c"
rm -f "$TMP/win.bed"
cat > "$TMP/split.awk" <<'EOF'
function rev(s,   r, j) { r = ""; for (j = length(s); j > 0; j--) r = r substr(s, j, 1); return r }
{   n = $1; sub(/\([^(]*$/, "", n); split(n, f, "_")
    key = f[1] "." f[2]
    if (key != cur) { if (cur != "") close(out); cur = key; out = DIR "/p." key "." CH }
    print f[3] "\t" (f[2] == "5" ? rev($2) : $2) > out }
EOF
getfa() {   # $1 = chunk bed
    local ch; ch=$(basename "$1")
    bedtools getfasta -s -tab -nameOnly -fi "$GENOME" -bed "$1" | "$AWK" -F'\t' -v DIR="$TMP/pieces" -v CH="$ch" -f "$TMP/split.awk"
}
export -f getfa
find "$TMP/chunks" -type f -name 'c*' -print0 | xargs -0 -n 1 -P "$THR" bash -c 'getfa "$1"' _
rm -rf "$TMP/chunks"

tick windows read
# ---- stage 1: count per group and side (all windows of the key) --------------------------------------------------------
#      one job per group, side and K-mer prefix (16 two-base prefixes; 4 one-base prefixes when K = 1), so that the largest
#      group is spread over 16 jobs; the classes are summed per group and side afterwards.
count1() {  # $1 = key gid.side, $2 = prefix class
    cat "$TMP/pieces/p.$1".* | "$AWK" -F'\t' -v KEY="$1" -v CL="$2" -v SH="$TMP/res/sh.$1.$2" -v K="$K" -v MIN_N="$MIN_N" \
        -v MIN_FRAC="$MIN_FRAC" -f "$TMP/ok.awk" -f "$TMP/count.awk" > "$TMP/res/st.$1.$2"
}
export -f count1
if [[ "$K" -ge 2 ]]; then CLASSES="AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT"; else CLASSES="A C G T"; fi
# keys gid.side, largest first (so that its jobs start first); keys.txt in (gid, side) order
find "$TMP/pieces" -type f -name 'p.*' -printf '%s %f\n' \
  | "$AWK" '{n = $2; sub(/^p\./, "", n); sub(/\.[^.]*$/, "", n); s[n] += $1} END {for (n in s) print s[n] "\t" n}' > "$TMP/keysize.tsv"
"$AWK" -F'\t' '{print $2}' "$TMP/keysize.tsv" | sort -t. -k1,1n -k2,2n > "$TMP/keys.txt"
sort -t$'\t' -k1,1nr -k2,2 "$TMP/keysize.tsv" | "$AWK" -F'\t' -v C="$CLASSES" 'BEGIN {nc = split(C, cl, " ")}
    {for (j = 1; j <= nc; j++) printf "%s\t%s\n", $2, cl[j]}' \
  | tr '\t\n' ' \0' | xargs -0 -n 1 -P "$THR" bash -c 'count1 ${1% *} ${1#* }' _
# per key: windows, threshold, shared K-mers (all classes; 0 when fewer than MIN_SHARED); the shared K-mers to res/sh.KEY
while IFS= read -r key; do
    read -r n thr ns < <("$AWK" -F'\t' 'NR == 1 {n = $2; t = $3} {s += $4} END {print n, t, s}' "$TMP/res/st.$key".*)
    if [[ "$ns" -lt "$MIN_SHARED" || "$ns" -eq 0 ]]; then ns=0; else cat "$TMP/res/sh.$key".* > "$TMP/res/sh.$key"; fi
    printf '%s\t%s\t%s\t%s\n' "$key" "$n" "$thr" "$ns"
done < "$TMP/keys.txt" > "$TMP/stats.tsv"

tick stage 1 counted
# ---- stage 2: walk, per piece, for the keys with shared K-mers ---------------------------------------------------------
walk1() {   # $1 = piece file
    local b key; b=$(basename "$1"); key=${b#p.}; key=${key%.*}
    "$AWK" -F'\t' -v KEY="$key" -v SH="$TMP/res/sh.$key" -v K="$K" -v STEP="$STEP" -v WIN="$WIN" -v SHARE="$SHARE" \
        -v BRIDGE="$BRIDGE" -f "$TMP/ok.awk" -f "$TMP/walk.awk" "$1" > "$TMP/res/w.$b"
}
export -f walk1
"$AWK" -F'\t' '$4 > 0 {print $1}' "$TMP/stats.tsv" | while IFS= read -r key; do
    find "$TMP/pieces" -type f -name "p.$key.*" -print0
done | xargs -0 -r -n 1 -P "$THR" bash -c 'walk1 "$1"' _
find "$TMP/res" -type f -name 'w.*' -exec cat {} + > "$TMP/walk.tsv"

tick stage 2 walked
# ---- per group and side table ------------------------------------------------------------------------------------------
"$AWK" -F'\t' -v MIN_N="$MIN_N" -v FF="$FLAG_FRAC" -v OUT="$OUT_DIR/ends_extension.tsv" -v FL="$TMP/nflag" '
FILENAME == ARGV[1] { name[$1] = $2; next }
FILENAME == ARGV[2] { if ($3 < 0) cap[$1]++; else { ne[$1]++; h[$1, $3]++; if ($3 >= 10) ng[$1]++; if ($3 > mx[$1]) mx[$1] = $3 } next }
function median(k, m,   b, c, lo, hi) {   # statistics.median of the values in h[k, 1..]: int for odd m, a.0 / a.5 for even m
    c = 0; lo = hi = ""
    for (b = 1; c < m; b++) {
        if (!((k, b) in h)) continue
        c += h[k, b]
        if (lo == "" && c >= int((m + 1) / 2)) lo = b
        if (hi == "" && c >= int(m / 2) + 1) hi = b
    }
    if (m % 2) return lo
    return ((lo + hi) % 2) ? sprintf("%d.5", int((lo + hi) / 2)) : sprintf("%d.0", (lo + hi) / 2)
}
{   split($1, f, "."); k = $1; n = $2; e = ne[k] + 0; g = ng[k] + 0
    fl = (n >= MIN_N && g >= FF * n) ? "flag" : ""; nf += (fl != "")
    printf "%s\t%s\047\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%s\n", name[f[1]], f[2], n, $3, $4, e, g, (e ? median(k, e) : 0), mx[k] + 0, cap[k] + 0, fl > OUT
}
BEGIN { print "group\tside\twindows\tshared_threshold\tshared_kmers\textended\textended_ge10\tmedian_bp\tmax_bp\tat_cap_not_extended\tflag" > OUT }
END { close(OUT); print nf + 0 > FL }' "$TMP/groups.tsv" "$TMP/walk.tsv" "$TMP/stats.tsv"

tick table written
# ---- new coordinates (neighbours that both grew into the gap meet in its middle), outputs -------------------------------
"$AWK" -F'\t' -v OFS='\t' -v MOVED="$OUT_DIR/loci_extended.tsv" -v OUT="$OUT_BED" -v CNT="$TMP/counts" '
FILENAME == ARGV[1] { name[$1] = $2; next }
FILENAME == ARGV[2] { if ($3 > 0) { if ($1 ~ /[.]5$/) { a5[$2 + 0] = $3; n5++ } else { a3[$2 + 0] = $3; n3++ } } next }   # numeric keys: much faster in mawk
FILENAME == ARGV[3] {      # sorted loci: chrom start end i gid strand
    s = $2 + 0; e = $3 + 0; i = $4 + 0
    b5 = (i in a5) ? a5[i] : 0; b3 = (i in a3) ? a3[i] : 0
    if ($6 == "+") { ns = s - b5; ne = e + b3 } else { ns = s - b3; ne = e + b5 }
    if (FNR > 1 && pc == $1 && pne > ns) {
        mid = int((pe + s) / 2)
        x = (pne < mid) ? pne : mid; pne = (pe > x) ? pe : x
        x = (ns > mid) ? ns : mid; ns = (s < x) ? s : x
    }
    if (FNR > 1 && (pns != ps || pne != pe)) { NS[pi] = pns; NE[pi] = pne; GG[pi] = pg }
    pc = $1; ps = s; pe = e; pi = i; pns = ns; pne = ne; pg = $5
    next
}
FNR == 1 && !done { if (pns != ps || pne != pe) { NS[pi] = pns; NE[pi] = pne; GG[pi] = pg }; done = 1
                    print "chrom\tstart\tend\tstrand\tnew_start\tnew_end\tgroup\tadded_5\tadded_3" > MOVED }
/[^[:space:]]/ {
    i = n++
    if (i in NS) {
        s = $2 + 0; e = $3 + 0; ns = NS[i]; ne = NE[i]
        if ($6 == "+") { b5 = s - ns; b3 = ne - e } else { b5 = ne - e; b3 = s - ns }
        print $1, s, e, $6, ns, ne, name[GG[i]], b5, b3 > MOVED
        $2 = ns; $3 = ne
    } else if ($2 !~ /^(0|[1-9][0-9]*)$/ || $3 !~ /^(0|[1-9][0-9]*)$/) { $2 = $2 + 0; $3 = $3 + 0 }   # written as integers, as before
    print > OUT
}
END { close(OUT); close(MOVED); print n5 + 0, n3 + 0 > CNT }' "$TMP/groups.tsv" "$TMP/walk.tsv" "$TMP/sorted.bed" "$MERGED"

tick coordinates written
read -r E5 E3 < "$TMP/counts"
NF_=$(cat "$TMP/nflag")
rm -rf "$TMP"
echo "EXTEND loci=$NLOCI extended5=$E5 extended3=$E3 flagged=$NF_ -> $OUT_DIR/ends_extension.tsv"
