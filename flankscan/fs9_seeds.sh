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
# Env: SD_AWK (the awk; default mawk >= 1.3.4, else gawk; the same rows), FS9_SORT (default gnusort if on PATH, else sort),
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
# the awk (pipeline rule 2026-11-27): $SD_AWK when set, else mawk if its version is >= 1.3.4 (1.3.3 crashed on KIT), else gawk
AWK=${SD_AWK:-}
if [[ -z $AWK ]]; then
    m=$(command -v mawk || true); v=""; [[ -n $m ]] && v=$("$m" -W version 2>&1 | head -n 1)
    if [[ $v =~ mawk\ ([0-9]+)\.([0-9]+)\.([0-9]+) ]] && (( BASH_REMATCH[1] * 1000000 + BASH_REMATCH[2] * 1000 + BASH_REMATCH[3] >= 1003004 )); then AWK=$m
    else AWK=$(command -v gawk) || { echo "fs9_seeds.sh: needs mawk >= 1.3.4 or gawk" >&2; exit 1; }; fi
fi
echo "fs9_seeds.sh: awk $AWK ($("$AWK" -W version 2>&1 | head -n 1))" >&2
SORT=${FS9_SORT:-$(command -v gnusort || command -v sort)}
Q=${FS9_Q:-16}; SMEM=${FS9_SORT_MEM:-512M}
(( IW + BAND < 128 && IW < 128 )) || { echo "fs9_seeds.sh: IW + BAND must be < 128 (delta and count fields)" >&2; exit 2; }
TD=$(mktemp -d "${TMPDIR:-/tmp}/fs9_seeds.XXXXXX"); trap 'rm -rf "$TD"' EXIT
stamp() { [[ -z ${FS9_SEEDS_LOG:-} ]] || echo "$1 $(date +%s.%N) $(du -sk "$TD" | cut -f1)" >> "$FS9_SEEDS_LOG"; }   # FS9_SEEDS_LOG=file: stage times (diagnostics)
stamp start

# copy rank of every row = place of its name among the distinct names in byte order (as sorted(set(names))); made with sort, not
# with an awk hash of the names (mawk hashes 10 characters of a long string: names that differ in the middle share hash slots)
TAB=$'\t'
cut -f1 "$SQ" | "$AWK" '{ print NR "\t" $0 }' | "$SORT" -t "$TAB" -k2,2 -S "$SMEM" -T "$TD" |
    "$AWK" -F'\t' -v NM="$TD/names" '{ if (NR == 1 || $2 != p) { r++; p = $2; print $2 > NM } print $1 "\t" r - 1 }' |
    "$SORT" -t "$TAB" -k1,1n -S "$SMEM" -T "$TD" | cut -f2 > "$TD/rank"
[[ -s $TD/names ]] || exit 0
NC=$(wc -l < "$TD/names")
SIDES=$(cut -f2 "$SQ" | "$SORT" -u | "$AWK" '{ s = s (NR > 1 ? "," : "") $0; n++ } END { if (n > 2) exit 1; print s }') ||
    { echo "fs9_seeds.sh: more than two flank sides in $SQ" >&2; exit 1; }
"$AWK" -v n="$NC" 'BEGIN { if (n * n * 512 >= 2 ^ 53) exit 1 }' || { echo "fs9_seeds.sh: $NC copies: too many for the 53-bit pair keys" >&2; exit 1; }
WC=${#NC}; WO=${#IW}; WEIGHT=$(( $(tr -cd 1 <<< "$PAT" | wc -c) )); PL=$(( WEIGHT < 3 ? WEIGHT : 3 ))
mkdir -p "$TD/A" "$TD/B" "$TD/C"

stamp ranks
# ---- A: seed rows, by row chunk
NR_=$(wc -l < "$SQ"); P=$TH; (( NR_ >= P )) || P=1
paste "$TD/rank" <(cut -f2,3 "$SQ") | "$AWK" -v P="$P" -v N="$NR_" -v D="$TD" '{ print > sprintf("%s/q%03d", D, int((NR - 1) * P / N)) }'
rm -f "$TD/rank"
seedgen() {   # CHUNK INDEX
    "$AWK" -F'\t' -v PAT="$PAT" -v IW="$IW" -v WO="$WO" -v WC="$WC" -v PL="$PL" -v OUT="$TD/A/$2" '
    BEGIN { L = length(PAT); nr = 0
            for (j = 1; j <= L; j++) if (substr(PAT, j, 1) == "1") { if (j > 1 && substr(PAT, j - 1, 1) == "1") rl[nr]++; else { rs[++nr] = j - 1; rl[nr] = 1 } }
            for (o = 1; o <= IW; o++) { r = 0; for (o2 = 1; o2 <= IW; o2++) if ((o2 "") < (o "")) r++; ORK[o] = sprintf("%0" WO "d", r) }   # offsets ranked as text
            np = split("A C G T", B, " "); F[""] = ""; for (k = 1; k <= PL; k++) { split("", G); for (p in F) for (b = 1; b <= 4; b++) G[p B[b]] = 1; split("", F); for (p in G) F[p] = OUT "." p } }
    { if (length($2) != 1) { print "fs9_seeds.sh: flank side must be one character: " $2 > "/dev/stderr"; exit 2 }
      s = substr($3, 1, IW + L - 1); n = length(s); hasn = index(s, "N"); c = sprintf("%0" WC "d", $1) $2   # $1: copy rank
      for (i = 1; i <= IW && i + L - 1 <= n; i++) {
          if (hasn && index(substr(s, i, L), "N")) continue
          if (nr == 1) m = substr(s, i, rl[1]); else { m = ""; for (k = 1; k <= nr; k++) m = m substr(s, i + rs[k], rl[k]) }
          if ((index(m, "A") > 0) + (index(m, "C") > 0) + (index(m, "G") > 0) + (index(m, "T") > 0) <= 2) continue
          print m ORK[i] c > F[substr(m, 1, PL)] } }' "$1"
}
pids=(); k=0; for c in "$TD"/q[0-9][0-9][0-9]; do seedgen "$c" "$k" & pids+=($!); k=$((k + 1)); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done; (( rc == 0 )) || { echo "fs9_seeds.sh: seed step failed" >&2; exit 1; }
rm -f "$TD"/q[0-9][0-9][0-9]

stamp A
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

stamp B
# ---- C: pairs, by residue class (r1 <= r2)
"$AWK" -v Q="$Q" 'BEGIN { for (a = 0; a < Q; a++) for (b = a; b < Q; b++) print a, b }' > "$TD/classes"
pairs() {   # R1 R2
    local -a f1=( "$TD"/B/*."$1" ) f2=( "$TD"/B/*."$2" )
    [[ -e ${f1[0]} ]] || f1=(); [[ -e ${f2[0]} ]] || f2=()
    if [[ $1 == "$2" ]]; then (( ${#f1[@]} )) || return 0; cat "${f1[@]}"
    else (( ${#f1[@]} && ${#f2[@]} )) || return 0; "$SORT" -m <(cat "${f1[@]}") <(cat "${f2[@]}"); fi |
    "$AWK" -v WF="${#IW}" -v WC="$WC" -v Q="$Q" -v R1="$1" -v SAME=$([[ $1 == "$2" ]] && echo 1 || echo 0) -v NC="$NC" -v MK="$MK" -v BAND="$BAND" \
        -v LIM=$(( START + SLACK )) -v SIDES="$SIDES" -v NAMES="$TD/names" '
        BEGIN { ns = split(SIDES, SV, ","); for (i = 1; i <= ns; i++) SB[SV[i]] = i - 1; PO = 14 + WF; PS = PO + WC; HP = 2147483647 }
        # The pair keys k = ((a * NC + b) * 2 + side) * 256 + delta + 128 are up to 2^48. mawk hashes a number that is not a 31-bit
        # integer as its decimal string, and its string hash looks at 10 characters only: 10^7 such keys took 40 us each (e1-2,
        # 2026-11-27). So the arrays here are indexed by k mod (2^31 - 1), a 31-bit integer (hashed well by mawk and gawk), with
        # the full key in K[] and linear probing when two keys share a slot (open addressing; few slots are shared at 10^7 keys).
        function pslot(x,   h) { h = x % HP; while (h in PK) { if (PK[h] == x) return h; if (++h >= HP) h = 0 } PK[h] = x; return h }   # pair slot (made if new)
        function kfind(x,   h) { h = x % HP; while (h in K) { if (K[h] == x) return h; if (++h >= HP) h = 0 } return -1 }
        function flush(   b, i, j, n, x0, ci, cj, a, bb, oa, ob, k, h, f, v, t, n1, n2, c1, c2, q, qs, qe) {
            for (b = 0; b < 2; b++) { n = N[b]; x0 = b * 4096; N[b] = 0
                if (n < 2) continue
                # partners of row i: the later rows of the bucket (SAME) or the later rows of the other residue, so that the pairs
                # come in the order (i, j), i < j, of the bucket rows without visiting the pairs of one residue (other classes)
                if (SAME) for (t = 1; t <= n; t++) IX[t] = x0 + t
                else { n1 = 0; n2 = 0; for (i = x0 + 1; i <= x0 + n; i++) if (RR[i] == R1) IX[++n1] = i; else IX[10000 + (++n2)] = i
                       if (!n1 || !n2) continue }
                c1 = 0; c2 = 0
                for (i = x0 + 1; i <= x0 + n; i++) { ci = C[i]
                    if (SAME) { qs = i - x0 + 1; qe = n } else if (RR[i] == R1) { c1++; qs = 10001 + c2; qe = 10000 + n2 } else { c2++; qs = 1 + c1; qe = n1 }
                    for (q = qs; q <= qe; q++) { j = IX[q]; cj = C[j]
                        if (ci == cj) continue
                        if (ci < cj) { a = ci; bb = cj; oa = O[i]; ob = O[j] } else { a = cj; bb = ci; oa = O[j]; ob = O[i] }
                        k = ((a * NC + bb) * 2 + b) * 256 + oa - ob + 128
                        h = k % HP; f = 0
                        while (h in K) { if (K[h] == k) { f = 1; break } if (++h >= HP) h = 0 }
                        # V = (a hit with offset_a <= LIM) * 256 + (offset_b <= LIM) * 128 + count (<= IW < 128); SL[t] = slot of the
                        # t-th new (pair, delta): t is the first-seen order (= the bucket order of the reference), which decides bd on ties
                        if (f) v = V[h] + 1; else { K[h] = k; SL[++nk] = h; v = 1 }
                        if (oa <= LIM && int(v / 256) % 2 == 0) v += 256
                        if (ob <= LIM && int(v / 128) % 2 == 0) v += 128
                        V[h] = v } } } }
        { s = substr($0, 1, 10); if (s != cur) { flush(); cur = s }
          b = SB[substr($0, PS + 1, 1)]; n = b * 4096 + (++N[b]); C[n] = substr($0, PO + 1, WC) + 0; O[n] = substr($0, 15, WF) + 0; RR[n] = C[n] % Q }   # numeric keys: side * 4096 + row
        END { flush()
              # (no "for (x in A)" here: mawk turns every key of A into a string for it, with the poor string hash)
              # per pair: seeds in all, and whether some seed has offset_a <= LIM and some offset_b <= LIM (PF 2 + 1)
              for (t = 1; t <= nk; t++) { h = SL[t]; v = V[h]; hp = pslot(int(K[h] / 256)); HS[t] = hp; PT[hp] += v % 128; x = PF[hp] + 0
                  if (int(v / 256) % 2 && int(x / 2) % 2 == 0) x += 2; if (int(v / 128) % 2 && x % 2 == 0) x += 1; PF[hp] = x }
              # best = most seeds within delta +- BAND; ties: the (pair, delta) seen first (t ascending, strict >)
              for (t = 1; t <= nk; t++) { h = SL[t]; k = K[h]; hp = HS[t]; if (PF[hp] != 3 || PT[hp] < MK) continue
                  s = 0; for (dd = -BAND; dd <= BAND; dd++) { y = k + dd; g = y % HP        # kfind(y), written out
                      while (g in K) { if (K[g] == y) { s += V[g] % 128; break } if (++g >= HP) g = 0 } }
                  if (s < MK) continue
                  if (!(hp in BE)) { PL[++np] = hp; BE[hp] = s; BD[hp] = k % 256 - 128 } else if (s > BE[hp]) { BE[hp] = s; BD[hp] = k % 256 - 128 } }
              if (!np) exit 0
              while ((getline l < NAMES) > 0) NM[nn++] = l; close(NAMES)
              for (t = 1; t <= np; t++) { hp = PL[t]; pk = PK[hp]; ab = int(pk / 2); a = int(ab / NC); bb = ab - a * NC
                  printf "%s\t%s\t%s\t%d\t%d\n", SV[pk % 2 + 1], NM[a], NM[bb], BD[hp], BE[hp] } }' > "$TD/C/$1.$2"
}
export -f pairs; export NC MK BAND START SLACK SIDES
xargs -P "$TH" -n 2 bash -c 'pairs "$@"' _ < "$TD/classes"
stamp C
cat "$TD"/C/*
