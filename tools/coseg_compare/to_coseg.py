#!/usr/bin/env python3
"""Align labelled copies to ONE reference and write COSEG input (NAME.seqs, NAME.ins, NAME.names).

usage: to_coseg.py REF.fa COPIES.fa OUTPREFIX [drop|keep] [JOBS]

Pure Python (no Biopython). Affine-gap global alignment with free end gaps on both sequences
(match 2, mismatch -3, gap open -5, gap extend -2: the settings used for the Konkel Alu run).
  drop (default): copies missing more than 5 reference bases at either end are left out, as in
                  Price et al. 2004 (the list of dropped copies goes to OUTPREFIX.dropped).
  keep:           all copies are written; missing ends are '-' like internal deletions.
JOBS > 1 aligns copies in parallel (order of the output does not change).
COSEG row format: one character per reference position (lower-case base, '-' = deletion or
missing end); an insertion after reference position p is marked '+' in the row after that
position and recorded in the .ins line as "p:bases".
"""
import sys

MATCH, MISMATCH, GAP_OPEN, GAP_EXT = 2, -3, -5, -2
NEG = -10 ** 9
_REF = ""


def read_fa(path):
    recs, cur = [], None
    for line in open(path):
        if line.startswith(">"):
            cur = [line[1:].strip(), []]
            recs.append(cur)
        elif cur is not None:
            cur[1].append(line.strip())
    return [(n, "".join(s).upper()) for n, s in recs]


def align(ref, q):
    """Gotoh with free end gaps. Returns the alignment columns as (ref_idx or None, q_idx or None)."""
    n, m = len(ref), len(q)
    # M: ends in a (mis)match; X: ends in a ref base against a gap in q; Y: ends in a q base against a gap in ref
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):          # free leading gap in q (reference bases before the copy starts)
        X[i][0] = 0
    for j in range(1, m + 1):          # free leading overhang of q
        Y[0][j] = 0
    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            s = MATCH if ri == q[j - 1] else MISMATCH
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + s
            X[i][j] = max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXT, Y[i - 1][j] + GAP_OPEN)
            Y[i][j] = max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXT, X[i][j - 1] + GAP_OPEN)
    # free trailing gaps: the best score anywhere on the last row or column, then free end
    best, bi, bj, bs = NEG, n, m, "M"
    for j in range(m + 1):
        for s, T in (("M", M), ("X", X), ("Y", Y)):
            if T[n][j] > best:
                best, bi, bj, bs = T[n][j], n, j, s
    for i in range(n + 1):
        for s, T in (("M", M), ("X", X), ("Y", Y)):
            if T[i][m] > best:
                best, bi, bj, bs = T[i][m], i, m, s
    cols = [(None, j - 1) for j in range(m, bj, -1)]       # trailing q overhang
    cols += [(i - 1, None) for i in range(n, bi, -1)]      # trailing reference bases
    i, j, s = bi, bj, bs
    while i > 0 or j > 0:
        if i == 0:
            cols.append((None, j - 1)); j -= 1; continue
        if j == 0:
            cols.append((i - 1, None)); i -= 1; continue
        if s == "M":
            cols.append((i - 1, j - 1))
            sc = M[i][j] - (MATCH if ref[i - 1] == q[j - 1] else MISMATCH)
            s = "M" if M[i - 1][j - 1] == sc else ("X" if X[i - 1][j - 1] == sc else "Y")
            i -= 1; j -= 1
        elif s == "X":
            cols.append((i - 1, None))
            v = X[i][j]
            s = "M" if M[i - 1][j] + GAP_OPEN == v else ("X" if X[i - 1][j] + GAP_EXT == v else "Y")
            i -= 1
        else:
            cols.append((None, j - 1))
            v = Y[i][j]
            s = "M" if M[i][j - 1] + GAP_OPEN == v else ("Y" if Y[i][j - 1] + GAP_EXT == v else "X")
            j -= 1
    cols.reverse()
    return cols


def to_row(cols, ref, q):
    """Row over the reference, insertions, and how many reference bases are missing at each end."""
    L = len(ref)
    row = ["-"] * L
    ins = {}
    # internal span = from the first to the last column where the copy has a base against a reference base
    matched = [k for k, (a, b) in enumerate(cols) if a is not None and b is not None]
    if not matched:
        return row, ins, L, L
    first, last = matched[0], matched[-1]
    pend, p = [], None
    for k in range(first, last + 1):
        a, b = cols[k]
        if a is not None and b is not None:
            if pend:
                ins[p] = ins.get(p, "") + "".join(pend); pend = []
            row[a] = q[b].lower(); p = a + 1
        elif a is None:
            pend.append(q[b].lower())
        # a reference base against a gap stays '-'
    lead = cols[first][0] if cols[first][0] is not None else 0
    trail = L - 1 - cols[last][0]
    return row, ins, lead, trail


def init(ref):
    global _REF
    _REF = ref


def one(rec):
    """Worker: align one copy; returns (name, row string with '+' marks, ins string, lead, trail)."""
    name, q = rec
    cols = align(_REF, q)
    row, ins, lead, trail = to_row(cols, _REF, q)
    out_row, insl = [], []
    for i, c in enumerate(row):
        out_row.append(c)
        if (i + 1) in ins:
            out_row.append("+")
            insl.append("%d:%s" % (i + 1, ins[i + 1]))
    return name, "".join(out_row), " ".join(insl), lead, trail


def main():
    global _REF
    ref_fa, copies_fa, out = sys.argv[1:4]
    mode = sys.argv[4] if len(sys.argv) > 4 else "drop"
    jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    _REF = "".join(s for _, s in read_fa(ref_fa)[:1])
    recs = read_fa(copies_fa)
    if jobs > 1:
        import multiprocessing
        with multiprocessing.Pool(jobs, initializer=init, initargs=(_REF,)) as pool:
            results = pool.map(one, recs, chunksize=max(1, len(recs) // (jobs * 8)))
    else:
        results = [one(r) for r in recs]
    kept = dropped = 0
    fs = open(out + ".seqs", "w")
    fi = open(out + ".ins", "w")
    fn = open(out + ".names", "w")
    fd = open(out + ".dropped", "w")
    for name, row, insl, lead, trail in results:
        if mode == "drop" and (lead > 5 or trail > 5):
            dropped += 1
            fd.write("%s\t%d\t%d\n" % (name, lead, trail))
            continue
        kept += 1
        fs.write(row + "\n")
        fi.write(insl + "\n")
        fn.write(name + "\n")
    for f in (fs, fi, fn, fd):
        f.close()
    print("reference %d bp; %s: kept %d, dropped %d" % (len(_REF), mode, kept, dropped))


if __name__ == "__main__":
    main()
