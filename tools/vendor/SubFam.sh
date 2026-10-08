#!/usr/bin/env bash
# SubFam: reduce a large set of related sequences (TE copies, satellite monomers,
# multigene family members) to an alignment of consensus sequences.
#
#   1. order all sequences along a k-mer guide tree (UPGMA; neighbours end up adjacent)
#   2. cut the ordered list into consecutive chunks of N neighbours
#   3. align each chunk (MAFFT) and call its plurality consensus
#   4. align the consensus sequences (MAFFT L-INS-i)
#
# Each consensus averages N closely related copies: private mutations of single
# copies cancel out, variants shared by a group survive. No identity threshold
# is needed, so it works equally for young and very old (highly diverged) families.
#
# Original idea by Toki, improved by Nikita (parallel chunk alignment).
# The k-mer guide tree is a port of kmer-tree.js from ViewAlign (Toki-bio/MSA-viewer);
# the consensus and MSF writers reproduce EMBOSS cons and seqret output exactly.

set -euo pipefail

VERSION=1.2.0
N=50            # sequences per chunk
PLURALITY=0.36  # fraction of chunk sequences that must agree for a consensus base
K=6            # k-mer size for the ordering tree
THREADS=$(nproc)
OUTDIR=.
PREFIX=
ORDER=kmer
BOTH=
PARTTREE=
KEEP_N=
KEEP_TMP=
RELCOV=
MINCOV=3       # with -c: fewest covering sequences for a consensus base

usage() {
    cat <<EOF
SubFam $VERSION - reduce related sequences to an alignment of chunk consensuses

Usage: $(basename "$0") [options] <sequences.fasta>

Options:
  -n INT    sequences per chunk (default: $N)
  -p FLOAT  plurality: fraction of a chunk that must agree on a base (default: $PLURALITY)
  -k INT    k-mer size of the ordering tree, 3-12 (default: $K, best at all divergences tested)
  -t INT    threads (default: all, $THREADS)
  -o DIR    output directory (default: current directory)
  -x STR    output prefix (default: input file name without extension)
  -r        sequences may be on both strands: strand-independent k-mers, then each copy
            is oriented like its neighbour
  -m        order with the MAFFT guide tree (mafft --retree 0 --reorder) instead of k-mers
  -P        order with MAFFT PartTree (implies -m; for very large inputs, >20,000 sequences)
  -c        coverage-relative plurality, for truncated copies (fragments, 5'-truncated LINEs):
            a base needs -p of the sequences that span that column (internal gaps count,
            end gaps do not) instead of -p of the whole chunk, with at least 3 spanning
            sequences. Default: -p of the whole chunk, which cuts a consensus down to the
            region most copies cover.
  -a        keep no-consensus positions as N (default: drop them from the consensus)
  -K        keep intermediate files (chunks, chunk alignments)
  -v        print version
  -h        show this help

Outputs (in DIR):
  PREFIX.cons.fasta   unaligned consensus sequences, one per chunk
  PREFIX.aln.fasta    aligned consensus sequences (FASTA)
  PREFIX.msf          the same alignment in MSF format
  PREFIX.chunks.tsv   input id, consensus it went into, strand (- = reverse-complemented
                      by -r, relative to the first sequence of the order)

Requires: mafft, awk, and for the k-mer ordering a C compiler (kmer_order.c) or Python >= 3.6 with numpy (not needed with -m/-P).
EOF
}

die() { echo "SubFam: $*" >&2; exit 1; }
elapsed() { local s=$(( $(date +%s) - $1 )); echo "  done in $((s / 3600))h $(((s / 60) % 60))m $((s % 60))s"; }

while getopts ":n:p:k:t:o:x:rmPcaKvh" opt; do
    case $opt in
        n) N=$OPTARG ;;
        p) PLURALITY=$OPTARG ;;
        k) K=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        x) PREFIX=$OPTARG ;;
        r) BOTH=1 ;;
        m) ORDER=mafft ;;
        P) ORDER=mafft; PARTTREE=--parttree ;;
        c) RELCOV=1 ;;
        a) KEEP_N=1 ;;
        K) KEEP_TMP=1 ;;
        v) echo "SubFam $VERSION"; exit 0 ;;
        h) usage; exit 0 ;;
        :) die "option -$OPTARG needs a value" ;;
        *) die "unknown option -$OPTARG (see -h)" ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 1 ] || { usage >&2; exit 1; }
IN=$1

[ -s "$IN" ] || die "input file '$IN' not found or empty"
for tool in mafft awk; do
    command -v "$tool" >/dev/null 2>&1 || die "'$tool' not found in PATH"
done
[[ $N =~ ^[0-9]+$ ]] && [ "$N" -ge 2 ] || die "-n must be an integer >= 2"
[[ $K =~ ^[0-9]+$ ]] && [ "$K" -ge 3 ] && [ "$K" -le 12 ] || die "-k must be an integer from 3 to 12"
[[ $THREADS =~ ^[0-9]+$ ]] && [ "$THREADS" -ge 1 ] || die "-t must be a positive integer"
awk -v p="$PLURALITY" 'BEGIN { exit !(p > 0 && p <= 1) }' || die "-p must be in (0, 1]"

PREFIX=${PREFIX:-$(basename "${IN%.*}")}
NAME=${PREFIX//[[:space:]]/_}   # FASTA ids cannot contain spaces
mkdir -p "$OUTDIR"
WORK=$(mktemp -d "$OUTDIR/.subfam_${PREFIX}_XXXXXX")
[ -n "$KEEP_TMP" ] || trap 'rm -rf "$WORK"' EXIT

# k-mer ordering program: the C version (kmer_order.c, compiled into $WORK for this run, so there is
# no stale binary and nothing is written outside the output directory), else the embedded Python
# version (needs Python >= 3.6 with numpy; KIT's default python3 is 3.5). KMER_IMPL=c|py forces one.
KMER_BIN=
PY=
if [ "$ORDER" = kmer ]; then
    HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)
    if [ "${KMER_IMPL:-}" != py ] && [ -f "$HERE/kmer_order.c" ]; then
        for cc in "${CC:-}" cc gcc clang; do
            [ -n "$cc" ] && command -v "$cc" >/dev/null 2>&1 || continue
            for flags in "-fopenmp" ""; do
                if "$cc" -O2 -ffp-contract=off $flags -o "$WORK/kmer_order" "$HERE/kmer_order.c" -lm 2>/dev/null; then
                    KMER_BIN="$WORK/kmer_order"; break 2
                fi
            done
        done
    fi
    if [ -z "$KMER_BIN" ] && [ "${KMER_IMPL:-}" != c ]; then
        for cand in "${PYTHON:-}" python3 python3.12 python3.11 python3.10 python3.9 python3.8 python3.7 python3.6; do
            [ -n "$cand" ] || continue
            if command -v "$cand" >/dev/null 2>&1 &&
               "$cand" -c 'import sys, numpy; sys.exit(sys.version_info < (3, 6))' 2>/dev/null; then PY=$cand; break; fi
        done
    fi
    [ -n "$KMER_BIN$PY" ] || die "k-mer ordering needs a C compiler (kmer_order.c) or Python >= 3.6 with numpy; or use -m"
fi

# ---------------------------------------------------------------------------------------
# Embedded programs
# ---------------------------------------------------------------------------------------
KMER_ORDER_PY=$(cat <<'PY'
import sys
import numpy as np

# k-mer guide-tree order: a port of kmer-tree.js (ViewAlign, Toki-bio/MSA-viewer, MIT).
# Weighted Jaccard distance on k-mer counts, UPGMA with ViewAlign's tie-breaking, and at
# every merge the two clusters are flipped so that their closest ends meet.
#   args: in.fasta out.fasta K CANONICAL(0/1) ORIENT(0/1)
# Reverse-complemented sequences (ORIENT=1) get MAFFT's "_R_" prefix on their id.

CODE = np.full(256, -1, dtype=np.int64)
for ch, v in zip("ACGTUacgtu", (0, 1, 2, 3, 3, 0, 1, 2, 3, 3)):
    CODE[ord(ch)] = v


def read_fasta(path):
    names, seqs, cur = [], [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip())
                cur = []
                seqs.append(cur)
            elif cur is not None:
                cur.append(line.strip())
    return names, ["".join(s) for s in seqs]


def kmer_codes(seq, k, canonical):
    """k-mer codes over the A/C/G/T/U letters of seq (others skipped, as in kmer-tree.js)."""
    v = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    v = v[v >= 0]
    if len(v) < k:
        return np.zeros(0, dtype=np.int64)
    w = np.lib.stride_tricks.sliding_window_view(v, k)
    h = (w << (2 * np.arange(k - 1, -1, -1))).sum(1)
    if canonical:
        r = ((3 - w) << (2 * np.arange(k))).sum(1)
        h = np.where(r < h, r, h)
    return h


def distances(seqs, k, canonical):
    """float32 weighted Jaccard distance matrix: 1 - sum(min) / sum(max) of k-mer counts."""
    n = len(seqs)
    codes = [kmer_codes(s, k, canonical) for s in seqs]
    total = np.array([len(c) for c in codes], dtype=np.float64)
    # sparse count table over k-mers seen in >= 2 sequences (others never contribute to sum(min))
    rows = np.concatenate([np.full(len(c), i) for i, c in enumerate(codes)]) if n else np.zeros(0, int)
    allc = np.concatenate(codes) if n else np.zeros(0, int)
    pair = np.unique(rows * (1 << 26) + np.searchsorted(np.unique(allc), allc), return_counts=True)
    r_idx, c_idx, cnt = pair[0] >> 26, pair[0] & ((1 << 26) - 1), pair[1]
    seen = np.bincount(c_idx)
    keep = seen[c_idx] >= 2
    r_idx, c_idx, cnt = r_idx[keep], c_idx[keep], cnt[keep]
    _, c_idx = np.unique(c_idx, return_inverse=True)
    m = int(c_idx.max()) + 1 if len(c_idx) else 0
    shared = np.zeros((n, n), dtype=np.float32)          # integer counts, exact below 2^24
    # sum(min(a, b)) = sum over t >= 1 of [a >= t][b >= t]: a matrix product per count level
    block = max(1, min(m, (1 << 24) // max(n, 1)))
    for t in range(1, int(cnt.max()) + 1 if len(cnt) else 1):
        sel = cnt >= t
        rr, cc = r_idx[sel], c_idx[sel]
        for lo in range(0, m, block):
            s = (cc >= lo) & (cc < lo + block)
            if not s.any():
                continue
            B = np.zeros((n, block), dtype=np.float32)
            B[rr[s], cc[s] - lo] = 1
            for r in range(0, n, 1024):                  # row blocks: no n x n temporary
                shared[r:r + 1024] += B[r:r + 1024] @ B.T
    d = shared                                            # converted in place, a row block at a time
    for lo in range(0, n, 256):
        sh = shared[lo:lo + 256].astype(np.float64)
        union = total[lo:lo + 256, None] + total[None, :] - sh
        with np.errstate(invalid="ignore", divide="ignore"):
            d[lo:lo + 256] = np.where(union > 0, 1 - sh / union, 1.0)
    np.fill_diagonal(d, 0)
    profiles = [np.unique(c, return_counts=True) for c in codes]
    return d, profiles, total


def pair_distance(profiles, total, a, b):
    """One entry of the distance matrix, recomputed exactly (float32, as stored)."""
    (ka, ca), (kb, cb) = profiles[a], profiles[b]
    _, ia, ib = np.intersect1d(ka, kb, assume_unique=True, return_indices=True)
    shared = float(np.minimum(ca[ia], cb[ib]).sum())
    union = total[a] + total[b] - shared
    return np.float32(1 - shared / union if union > 0 else 1.0) if a != b else np.float32(0)


def upgma_order(cd, dist):
    """Leaf order of kmer-tree.js upgma() (average linkage), same tie-breaking.
    cd: distance matrix, overwritten; dist(a, b): original distance of two sequences."""
    n = len(cd)
    if n < 2:
        return list(range(n))
    clusters = [[i] for i in range(n)]
    active = np.ones(n, dtype=bool)
    nnD = np.full(n, np.inf)
    nnJ = np.full(n, -1)

    def scan(i):
        row = np.where(active[i + 1:], cd[i, i + 1:], np.float32(np.inf))
        if len(row) and np.isfinite(row.min()):
            j = int(row.argmin())
            nnD[i], nnJ[i] = row[j], i + 1 + j
        else:
            nnD[i], nnJ[i] = np.inf, -1

    for i in range(n):
        scan(i)
    idx = np.arange(n)
    for _ in range(n - 1):
        cand = np.where(active & (nnJ >= 0), nnD, np.inf)
        ci = int(cand.argmin())
        if not np.isfinite(cand[ci]):
            break
        cj = int(nnJ[ci])
        A, B = clusters[ci], clusters[cj]
        opts = [dist(A[-1], B[0]), dist(A[-1], B[-1]), dist(A[0], B[0]), dist(A[0], B[-1])]
        best = int(np.argmin(opts))                      # first minimum, like a stable sort
        A = A[::-1] if best >= 2 else A
        B = B[::-1] if best in (1, 3) else B
        sI, sJ = len(A), len(B)
        clusters[ci], clusters[cj] = A + B, None
        active[cj] = False
        upd = active.copy()
        upd[ci] = False
        new = ((cd[ci].astype(np.float64) * sI + cd[cj].astype(np.float64) * sJ) / (sI + sJ)).astype(np.float32)
        cd[ci, upd] = new[upd]
        cd[upd, ci] = new[upd]
        scan(ci)
        others = active.copy()
        others[ci] = False
        rescan = others & ((nnJ == ci) | (nnJ == cj))
        for i in np.flatnonzero(rescan):
            scan(int(i))
        low = others & ~rescan & (idx < ci)
        v = cd[low, ci].astype(np.float64)
        win = (v < nnD[low]) | ((v == nnD[low]) & (ci < nnJ[low]))
        li = np.flatnonzero(low)[win]
        nnD[li], nnJ[li] = v[win], ci
    return next(c for i, c in enumerate(clusters) if active[i])


COMP = str.maketrans("ACGTUacgtuRYKMBVDHrykmbvdh", "TGCAAtgcaaYRMKVBHDyrmkvbhd")


def main():
    src, out, k, canonical, orient = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] == "1", sys.argv[5] == "1"
    names, seqs = read_fasta(src)
    d, profiles, total = distances(seqs, k, canonical)
    order = upgma_order(d, lambda a, b: pair_distance(profiles, total, a, b))
    flip = [False] * len(seqs)
    if orient:
        # orient each sequence like its predecessor in the order: its nearest neighbour
        prev = None
        for i in order:
            if prev is not None:
                ref = set(kmer_codes(prev, k, False).tolist())
                fwd = sum(c in ref for c in kmer_codes(seqs[i], k, False).tolist())
                rev = sum(c in ref for c in kmer_codes(seqs[i].translate(COMP)[::-1], k, False).tolist())
                flip[i] = rev > fwd
            prev = seqs[i].translate(COMP)[::-1] if flip[i] else seqs[i]
    with open(out, "w") as fh:
        for i in order:
            if flip[i]:
                fh.write(f">_R_{names[i]}\n{seqs[i].translate(COMP)[::-1]}\n")
            else:
                fh.write(f">{names[i]}\n{seqs[i]}\n")


if __name__ == "__main__":
    main()
PY
)

CONS_AWK=$(cat <<'AWK'
# Plurality consensus of an aligned nucleotide FASTA: EMBOSS cons (embConsCalc) with its
# defaults (EDNAFULL matrix, unit weights), reimplemented. Variables: plur (plurality,
# absolute count), setcase (default: half the number of sequences), name.
# rel=1 (SubFam -c): the threshold is frac * (sequences spanning the column) instead of plur,
# needs mincov spanning sequences, and end gaps do not take part in the base choice.
BEGIN {
    sym = "A T G C S W R Y K M B V H D N U"
    nsym = split(sym, S, " ")
    split("5 -4 -4 -4 -4 1 1 -4 -4 1 -4 -1 -1 -1 -2 -4 " \
          "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5 " \
          "-4 -4 5 -4 1 -4 1 -4 1 -4 -1 -1 -4 -1 -2 -4 " \
          "-4 -4 -4 5 1 -4 -4 1 -4 1 -1 -1 -1 -4 -2 -4 " \
          "-4 -4 1 1 -1 -4 -2 -2 -2 -2 -1 -1 -3 -3 -1 -4 " \
          "1 1 -4 -4 -4 -1 -2 -2 -2 -2 -3 -3 -1 -1 -1 1 " \
          "1 -4 1 -4 -2 -2 -1 -4 -2 -2 -3 -1 -3 -1 -1 -4 " \
          "-4 1 -4 1 -2 -2 -4 -1 -2 -2 -1 -3 -1 -3 -1 1 " \
          "-4 1 1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1 1 " \
          "1 -4 -4 1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1 -4 " \
          "-4 -1 -1 -1 -1 -3 -3 -1 -1 -3 -1 -2 -2 -2 -1 -1 " \
          "-1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1 -4 " \
          "-1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1 -1 " \
          "-1 -1 -1 -4 -3 -1 -1 -3 -1 -3 -2 -2 -2 -1 -1 -1 " \
          "-2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 " \
          "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5", V, " ")       # EDNAFULL
    for (a = 1; a <= nsym; a++) {
        isym[S[a]] = 1
        for (b = 1; b <= nsym; b++) M[S[a], S[b]] = V[(a - 1) * nsym + b]
    }
}
/^>/ { n++; next }
{ gsub(/[ \t\r]/, ""); gsub(/[.~]/, "-"); gsub(/[Xx]/, "N"); seq[n] = seq[n] toupper($0) }   # as EMBOSS reads DNA
END {
    if (setcase == "") setcase = n / 2
    for (i = 1; i <= n; i++) if (length(seq[i]) > L) L = length(seq[i])
    if (rel) for (i = 1; i <= n; i++) {          # span of each sequence: first..last residue
        t = seq[i]; fi[i] = 0; la[i] = -1
        if (match(t, /[^-]/)) { fi[i] = RSTART; sub(/-+$/, "", t); la[i] = length(t) }
    }
    out = ""
    for (k = 1; k <= L; k++) {
        delete cnt; ncov = 0
        for (i = 1; i <= n; i++) {
            c = substr(seq[i], k, 1); col[i] = c
            if (c in isym) cnt[c]++
            if (rel && k >= fi[i] && k <= la[i]) ncov++
        }
        if (rel) {
            need = frac * ncov; need = (need == int(need)) ? need : int(need) + 1
            mc = (mincov < n) ? mincov : n
            if (ncov < mc) { out = out "N"; continue }
        }
        delete sc; delete mt
        for (s in cnt) {                       # score of one sequence with symbol s, and its +ve matches
            x = -M[s, s]; y = 0
            for (c in cnt) { x += cnt[c] * M[s, c]; if (M[s, c] > 0) y += cnt[c] }
            sc[s] = x; mt[s] = y
        }
        hi = 0; max = "unset"
        for (i = 1; i <= n; i++) {             # first sequence with the top score; a gap yields ties
            if (rel && !(col[i] in isym) && (k < fi[i] || k > la[i])) continue   # end gap: not a candidate
            v = (col[i] in isym) ? sc[col[i]] : 0
            if (max == "unset" || v > max || (v == max && col[hi] == "-")) { hi = i; max = v }
        }
        r = col[hi]; m = (r in isym) ? mt[r] : 0
        res = (m >= (rel ? need : plur)) ? r : "N"
        if (m <= setcase) res = tolower(res)
        out = out res
    }
    print ">" name
    for (i = 1; i <= length(out); i += 60) print substr(out, i, 60)
}
AWK
)

MSF_AWK=$(cat <<'AWK'
# Aligned FASTA -> GCG MSF, as EMBOSS seqret -osformat2 msf writes it (nucleotides).
# Variables: file (name printed in the header), date (dd/mm/yy).
function check(s,   i, c, t) {                 # GCG checksum
    t = 0; c = 0
    for (i = 1; i <= length(s); i++) { if (++c > 57) c = 1; t += c * ORD[toupper(substr(s, i, 1))] }
    return t % 10000
}
BEGIN { for (i = 32; i < 127; i++) ORD[sprintf("%c", i)] = i }
/^>/ { n++; split(substr($0, 2), w, /[ \t]/); name[n] = w[1]; next }
{ gsub(/[ \t\r]/, ""); seq[n] = seq[n] $0 }
END {
    L = 0; W = 10
    for (i = 1; i <= n; i++) { if (length(seq[i]) > L) L = length(seq[i]); if (length(name[i]) > W) W = length(name[i]) }
    tot = 0
    for (i = 1; i <= n; i++) {
        s = seq[i]
        gsub(/[ .~-]/, ".", s)                   # internal gaps '.', end gaps and padding '~'
        if (match(s, /^\.+/)) s = sprintf("%" RLENGTH "s", "") substr(s, RLENGTH + 1)
        if (match(s, /\.+$/)) s = substr(s, 1, RSTART - 1) sprintf("%" RLENGTH "s", "")
        while (length(s) < L) s = s " "
        gsub(/ /, "~", s)
        seq[i] = s; ck[i] = check(s); tot = (tot + ck[i]) % 10000
    }
    printf "!!NA_MULTIPLE_ALIGNMENT 1.0\n\n"
    printf "  %s MSF: %d Type: N %s CompCheck: %4d ..\n\n", file, L, date, tot
    for (i = 1; i <= n; i++) printf "  Name: %-" W "s Len: %d  Check: %4d Weight: 1.00\n", name[i], length(seq[i]), ck[i]
    printf "\n//\n\n"
    for (p = 1; p <= L; p += 50) {
        e = p + 49; if (e > L) e = L
        if (e == L) {
            g = e - p - length(p "")
            if (g < 0) printf "%" W "s %d %d\n", " ", p, e      # EMBOSS: unsigned width wraps
            else if (g >= length(e "")) printf "%" W "s %d %" g "d\n", " ", p, e
            else printf "           %d\n", p
        } else printf "           %-25d%25d\n", p, e
        for (i = 1; i <= n; i++) printf "%-" W "s %s\n", name[i], substr(seq[i], p, e - p + 1)
        printf "\n"
    }
}
AWK
)

# ---------------------------------------------------------------------------------------
T0=$(date +%s)
TOTAL=$(grep -c '^>' "$IN" || true)
[ "$TOTAL" -ge 1 ] || die "no FASTA records in '$IN'"
NCHUNK=$(( TOTAL / N )); [ "$NCHUNK" -ge 1 ] || NCHUNK=1
echo "SubFam $VERSION: $TOTAL sequences -> $NCHUNK chunks of ~$N"

T1=$(date +%s)
if [ "$ORDER" = kmer ]; then
    echo "Ordering sequences along the $K-mer guide tree"
    if [ -n "$KMER_BIN" ]; then
        OMP_NUM_THREADS=$THREADS "$KMER_BIN" "$IN" "$WORK/ordered.fasta" "$K" "${BOTH:-0}" "${BOTH:-0}"
    else
        "$PY" -c "$KMER_ORDER_PY" "$IN" "$WORK/ordered.fasta" "$K" "${BOTH:-0}" "${BOTH:-0}"
    fi
else
    echo "Ordering sequences along the MAFFT guide tree"
    mafft --thread "$THREADS" --nuc --quiet --retree 0 --reorder $PARTTREE ${BOTH:+--adjustdirection} \
        "$IN" > "$WORK/ordered.fasta"
fi
elapsed "$T1"

# Split into NCHUNK chunks of N consecutive sequences; the remainder (< N)
# joins the last chunk so that no sequence is lost.
awk -v n="$N" -v nc="$NCHUNK" -v dir="$WORK" -v pre="$NAME" -v map="$WORK/chunks.tsv" '
    BEGIN { fmt = "%s_%0" length(nc "") "d" }
    /^>/ {
        c = int(i++ / n); if (c >= nc) c = nc - 1
        name = sprintf(fmt, pre, c + 1)
        out = dir "/" name ".fasta"
        id = substr($1, 2); sub(/^_R_/, "", id)   # strip the reverse-complement tag
        print id "\t" name "\t" (substr($1, 2, 3) == "_R_" ? "-" : "+") > map
    }
    { print > out }
' "$WORK/ordered.fasta"
rm -f "$WORK/ordered.fasta"

echo "Aligning chunks and calling consensus sequences"
T1=$(date +%s)
export PLURALITY KEEP_N CONS_AWK RELCOV MINCOV
find "$WORK" -maxdepth 1 -name "${NAME}_*.fasta" -print0 | sort -z |
    xargs -0 -P "$THREADS" -I {} sh -c '
        set -e; f=$1; name=$(basename "$f" .fasta)
        n=$(grep -c "^>" "$f")
        plur=$(awk -v n="$n" -v p="$PLURALITY" "BEGIN { x = n * p; print (x == int(x)) ? x : int(x) + 1 }")
        if [ "$n" -ge 2 ]; then mafft --thread 1 --nuc --quiet "$f" > "$f.aln"; else cp "$f" "$f.aln"; fi
        awk -v plur="$plur" -v name="$name" -v rel="$RELCOV" -v frac="$PLURALITY" -v mincov="$MINCOV" "$CONS_AWK" "$f.aln" |
            awk -v keep="$KEEP_N" "!/^>/ && !keep { gsub(/[Nn]/, \"\") } 1" > "$f.cons"
    ' _ {}
find "$WORK" -maxdepth 1 -name "${NAME}_*.fasta.cons" -print0 | sort -z | xargs -0 cat > "$OUTDIR/$PREFIX.cons.fasta"
cp "$WORK/chunks.tsv" "$OUTDIR/$PREFIX.chunks.tsv"
elapsed "$T1"

echo "Aligning consensus sequences"
T1=$(date +%s)
if [ "$NCHUNK" -ge 2 ]; then
    mafft --thread "$THREADS" --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet \
        "$OUTDIR/$PREFIX.cons.fasta" > "$OUTDIR/$PREFIX.aln.fasta"
else
    cp "$OUTDIR/$PREFIX.cons.fasta" "$OUTDIR/$PREFIX.aln.fasta"
fi
awk -v file="$OUTDIR/$PREFIX.msf" -v date="$(date +%d/%m/%y)" "$MSF_AWK" "$OUTDIR/$PREFIX.aln.fasta" > "$OUTDIR/$PREFIX.msf"
elapsed "$T1"

[ -z "$KEEP_TMP" ] || echo "Intermediate files kept in $WORK"
echo "Job completed: $OUTDIR/$PREFIX.{cons.fasta,aln.fasta,msf,chunks.tsv}"
elapsed "$T0"
