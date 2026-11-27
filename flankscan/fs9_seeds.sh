#!/usr/bin/env bash
# fs9_seeds.sh SEQS.tsv PATTERN MK [--iw 80] [--cap 20] [--band 3] [--start 20] [--slack 25] [--threads 1]
#
# Candidate pairs of one fs9_twins.sh pass in awk + sort, the same rows as fs9_seeds.py (rows unsorted on stdout; the caller
# pipes them through sort -u). Written to test whether awk can do this step at the speed of the numpy version (2026-11-27).
# The rules (fs9_seeds.py header): seeds = the PATTERN's 1-positions of windows 1..IW (no N in the window, > 2 distinct bases);
# a bucket = all rows of one seed in the order of the whole line "seed offset copy side" (offset compared as text); rows whose
# (seed, offset) occurs more than CAP times are not paired; two rows pair when they are of one side and two copies; per pair
# (side, a, b): best = most seeds within delta +- BAND of one of its deltas, bd = the delta of the first seed (bucket order)
# reaching best; kept when best >= MK and the smallest offset of a and of b are <= START + SLACK.
#
# How (three stages, each in THREADS processes; temporary files under TMPDIR, about 50 bytes per seed row in all):
#  A  seeds: the rows in THREADS chunks, one awk each; every seed row is written as one fixed-width line
#     SEED OFFSETRANK COPYRANK SIDE (copy rank = place of the name in byte order) into one of 64 files by the seed's first 3 bases
#  B  buckets: per 3-base prefix, sort (byte order of the fixed-width line = the bucket order) and one awk: drops the (seed, offset)
#     runs longer than CAP and the rows with no partner, numbers the buckets (serial) and the rows inside them (position), and
#     writes each row to one of Q residue files by copy rank mod Q
#  C  pairs: the copy pairs fall into Q(Q+1)/2 classes by the residues of their two copies; one awk per class reads only the rows
#     of its two residues (sort -m of the two residue streams: bucket order kept), pairs them, and counts per (pair, delta) in
#     a hash with NUMERIC keys (mawk hashes numbers about 7x faster than strings). The count's value also holds the order in which
#     the (pair, delta) was first seen = the bucket order, which decides bd on ties. Memory per class: about Q^-2 of all seed pairs.
# Env: FS9_SEEDS_AWK (default mawk if on PATH, else gawk; the same rows), FS9_SORT (default gnusort if on PATH, else sort),
#      FS9_Q (default 16), FS9_SORT_MEM (default 512M per sort).
set -euo pipefail
export LC_ALL=C
SQ=${1:?SEQS.tsv}; PAT=${2:?PATTERN}; MK=${3:?MK}; shift 3
IW=80; CAP=20; BAND=3; START=20; SLACK=25; TH=1
while (( $# )); do
    case $1 in
        --iw) IW=$2 ;; --cap) CAP=$2 ;; --band) BAND=$2 ;; --start) START=$2 ;; --slack) SLACK=$2 ;; --threads) TH=$2 ;;
        *) echo "fs9_seeds.sh: unknown option $1" >&2; exit 2 ;;
    esac; shift 2
done
[[ $PAT =~ ^1[01]*$ ]] || { echo "fs9_seeds.sh: PATTERN must be 0/1 and start with 1" >&2; exit 2; }
[[ $TH =~ ^[1-9][0-9]*$ ]] || TH=1
AWK=${FS9_SEEDS_AWK:-$(command -v mawk || command -v gawk)}
SORT=${FS9_SORT:-$(command -v gnusort || command -v sort)}
Q=${FS9_Q:-16}; SMEM=${FS9_SORT_MEM:-512M}
(( IW + BAND < 128 && IW < 128 )) || { echo "fs9_seeds.sh: IW + BAND must be < 128 (delta and count fields)" >&2; exit 2; }
TD=$(mktemp -d "${TMPDIR:-/tmp}/fs9_seeds.XXXXXX"); trap 'rm -rf "$TD"' EXIT

cut -f1 "$SQ" | "$SORT" -u > "$TD/names"                 # copy rank = line number - 1 (byte order, as sorted(set(names)))
NC=$(wc -l < "$TD/names"); (( NC > 0 )) || exit 0
SIDES=$(cut -f2 "$SQ" | "$SORT" -u | "$AWK" '{ s = s (NR > 1 ? "," : "") $0; n++ } END { if (n > 2) exit 1; print s }') ||
    { echo "fs9_seeds.sh: more than two flank sides in $SQ" >&2; exit 1; }
"$AWK" -v n="$NC" 'BEGIN { if (n * n * 512 >= 2 ^ 53) exit 1 }' || { echo "fs9_seeds.sh: $NC copies: too many for the 53-bit pair keys" >&2; exit 1; }
WC=${#NC}; WO=${#IW}; WEIGHT=$(( $(tr -cd 1 <<< "$PAT" | wc -c) )); PL=$(( WEIGHT < 3 ? WEIGHT : 3 ))
mkdir -p "$TD/A" "$TD/B" "$TD/C"

# ---- A: seed rows, by row chunk
NR_=$(wc -l < "$SQ"); P=$TH; (( NR_ >= P )) || P=1
"$AWK" -v P="$P" -v N="$NR_" -v D="$TD" '{ print > sprintf("%s/q%03d", D, int((FNR - 1) * P / N)) }' "$SQ"
seedgen() {   # CHUNK INDEX
    "$AWK" -F'\t' -v PAT="$PAT" -v IW="$IW" -v WO="$WO" -v WC="$WC" -v PL="$PL" -v NAMES="$TD/names" -v OUT="$TD/A/$2" '
    BEGIN { L = length(PAT); nr = 0
            for (j = 1; j <= L; j++) if (substr(PAT, j, 1) == "1") { if (j > 1 && substr(PAT, j - 1, 1) == "1") rl[nr]++; else { rs[++nr] = j - 1; rl[nr] = 1 } }
            for (o = 1; o <= IW; o++) { r = 0; for (o2 = 1; o2 <= IW; o2++) if ((o2 "") < (o "")) r++; ORK[o] = sprintf("%0" WO "d", r) }   # offsets ranked as text
            while ((getline l < NAMES) > 0) RK[l] = sprintf("%0" WC "d", nn++); close(NAMES)
            np = split("A C G T", B, " "); F[""] = ""; for (k = 1; k <= PL; k++) { split("", G); for (p in F) for (b = 1; b <= 4; b++) G[p B[b]] = 1; split("", F); for (p in G) F[p] = OUT "." p } }
    { if (length($2) != 1) { print "fs9_seeds.sh: flank side must be one character: " $2 > "/dev/stderr"; exit 2 }
      s = substr($3, 1, IW + L - 1); n = length(s); hasn = index(s, "N"); c = RK[$1] $2
      for (i = 1; i <= IW && i + L - 1 <= n; i++) {
          if (hasn && index(substr(s, i, L), "N")) continue
          if (nr == 1) m = substr(s, i, rl[1]); else { m = ""; for (k = 1; k <= nr; k++) m = m substr(s, i + rs[k], rl[k]) }
          if ((index(m, "A") > 0) + (index(m, "C") > 0) + (index(m, "G") > 0) + (index(m, "T") > 0) <= 2) continue
          print m ORK[i] c > F[substr(m, 1, PL)] } }' "$1"
}
pids=(); k=0; for c in "$TD"/q[0-9][0-9][0-9]; do seedgen "$c" "$k" & pids+=($!); k=$((k + 1)); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done; (( rc == 0 )) || { echo "fs9_seeds.sh: seed step failed" >&2; exit 1; }
rm -f "$TD"/q[0-9][0-9][0-9]

# ---- B: buckets, by seed prefix (prefix index in byte order = the order of the buckets)
"$AWK" -v PL="$PL" 'BEGIN { F[""] = 1; for (k = 1; k <= PL; k++) { split("", G); for (p in F) { G[p "A"]; G[p "C"]; G[p "G"]; G[p "T"] }; split("", F); for (p in G) F[p] = 1 }
                             for (p in F) print p }' | "$SORT" | "$AWK" '{ print NR - 1, $1 }' > "$TD/prefixes"
buckets() {   # PREFIX_INDEX PREFIX
    local f; local -a fs=(); for f in "$TD"/A/*."$2"; do [[ -e $f ]] && fs+=("$f"); done
    (( ${#fs[@]} )) || return 0
    cat "${fs[@]}" | "$SORT" -S "$SMEM" -T "$TD" | "$AWK" -v LS="$WEIGHT" -v WO="$WO" -v WC="$WC" -v IW="$IW" -v CAP="$CAP" -v Q="$Q" \
        -v BASE=$(( $1 * 16777216 )) -v OUT="$TD/B/$(printf %02d "$1")" '
        BEGIN { for (o = 1; o <= IW; o++) { r = 0; for (o2 = 1; o2 <= IW; o2++) if ((o2 "") < (o "")) r++; OFF[sprintf("%0" WO "d", r)] = o }
                WF = length(IW ""); LK = LS + WO; for (r = 0; r < Q; r++) RF[r] = OUT "." r }
        function flush(   i, j, t, m, keep, k) {
            if (n < 2) { n = 0; return }
            i = 1; m = 0; split("", seen); split("", nd)
            while (i <= n) { k = substr(R[i], 1, LK); j = i; while (j < n && substr(R[j + 1], 1, LK) == k) j++
                keep = (j - i + 1 <= CAP)                            # a (seed, offset) seen more than CAP times is not paired
                if (keep) for (t = i; t <= j; t++) { m++; cp[m] = substr(R[t], LK + 1, WC); sd[m] = substr(R[t], LK + WC + 1, 1); of[m] = OFF[substr(R[t], LS + 1, WO)]
                                                     if (!((sd[m] SUBSEP cp[m]) in seen)) { seen[sd[m] SUBSEP cp[m]]; nd[sd[m]]++ } }
                i = j + 1 }
            n = 0
            if (m < 2) return
            if (m >= 2048) { print "fs9_seeds.sh: more than 2047 rows in one seed bucket" > "/dev/stderr"; exit 3 }
            for (t = 1; t <= m; t++) if (nd[sd[t]] >= 2) break
            if (t > m) return                                        # no side with two copies: nothing pairs
            if (++ser >= 16777216) { print "fs9_seeds.sh: more than 2^24 buckets with one prefix" > "/dev/stderr"; exit 3 }
            S = sprintf("%010d", BASE + ser)
            for (t = 1; t <= m; t++) if (nd[sd[t]] >= 2) printf "%s%04d%0" WF "d%s%s\n", S, t, of[t], cp[t], sd[t] > RF[cp[t] % Q] }
        { k = substr($0, 1, LS); if (k != cur) { flush(); cur = k } R[++n] = $0 }
        END { flush() }'
    rm -f "${fs[@]}"
}
export -f buckets; export TD SORT SMEM AWK WEIGHT WO WC IW CAP Q
xargs -P "$TH" -n 2 bash -c 'buckets "$@"' _ < "$TD/prefixes"
# a failing buckets() makes xargs return 123: set -e stops the script

# ---- C: pairs, by residue class (r1 <= r2)
"$AWK" -v Q="$Q" 'BEGIN { for (a = 0; a < Q; a++) for (b = a; b < Q; b++) print a, b }' > "$TD/classes"
pairs() {   # R1 R2
    local -a f1=( "$TD"/B/*."$1" ) f2=( "$TD"/B/*."$2" )
    [[ -e ${f1[0]} ]] || f1=(); [[ -e ${f2[0]} ]] || f2=()
    if [[ $1 == "$2" ]]; then (( ${#f1[@]} )) || return 0; cat "${f1[@]}"
    else (( ${#f1[@]} && ${#f2[@]} )) || return 0; "$SORT" -m <(cat "${f1[@]}") <(cat "${f2[@]}"); fi |
    "$AWK" -v WF="${#IW}" -v WC="$WC" -v Q="$Q" -v SAME=$([[ $1 == "$2" ]] && echo 1 || echo 0) -v NC="$NC" -v MK="$MK" -v BAND="$BAND" \
        -v LIM=$(( START + SLACK )) -v SIDES="$SIDES" -v NAMES="$TD/names" '
        BEGIN { ns = split(SIDES, SV, ","); for (i = 1; i <= ns; i++) SB[SV[i]] = i - 1; PO = 14 + WF; PS = PO + WC }
        function flush(   b, i, j, n, x0, ci, cj, a, bb, oa, ob, pk, k, v) {
            for (b = 0; b < 2; b++) { n = N[b]; x0 = b * 4096
                for (i = x0 + 1; i < x0 + n; i++) { ci = C[i]
                    for (j = i + 1; j <= x0 + n; j++) { cj = C[j]
                        if (ci == cj || (!SAME && RS[i] == RS[j])) continue
                        if (ci < cj) { a = ci; bb = cj; oa = O[i]; ob = O[j] } else { a = cj; bb = ci; oa = O[j]; ob = O[i] }
                        pk = (a * NC + bb) * 2 + b; k = pk * 256 + oa - ob + 128
                        if ((v = c[k]) == 0) c[k] = (++idx) * 128 + 1; else c[k] = v + 1   # first-seen order * 128 + count (count <= IW)
                        if (oa <= LIM) la[pk]; if (ob <= LIM) lb[pk] } }
                N[b] = 0 } }
        { s = substr($0, 1, 10); if (s != cur) { flush(); cur = s }
          b = SB[substr($0, PS + 1, 1)]; n = b * 4096 + (++N[b]); C[n] = substr($0, PO + 1, WC) + 0; O[n] = substr($0, 15, WF) + 0; RS[n] = C[n] % Q }   # numeric keys: side * 4096 + row
        END { flush()
              for (k in c) { k += 0; pk = int(k / 256); if (!(pk in la) || !(pk in lb)) continue
                  s = 0; for (dd = -BAND; dd <= BAND; dd++) if ((k + dd) in c) s += c[k + dd] % 128
                  if (s < MK) continue
                  id = int(c[k] / 128)
                  if (!(pk in BE) || s > BE[pk] || (s == BE[pk] && id < BI[pk])) { BE[pk] = s; BI[pk] = id; BD[pk] = k % 256 - 128 } }
              if (length(BE) == 0) exit 0
              while ((getline l < NAMES) > 0) NM[nn++] = l; close(NAMES)
              for (pk in BE) { pk += 0; ab = int(pk / 2); a = int(ab / NC); bb = ab - a * NC
                  printf "%s\t%s\t%s\t%d\t%d\n", SV[pk % 2 + 1], NM[a], NM[bb], BD[pk], BE[pk] } }' > "$TD/C/$1.$2"
}
export -f pairs; export NC MK BAND START SLACK SIDES
xargs -P "$TH" -n 2 bash -c 'pairs "$@"' _ < "$TD/classes"
cat "$TD"/C/*
