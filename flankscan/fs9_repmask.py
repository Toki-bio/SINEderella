#!/usr/bin/env python3
"""fs9_repmask.py rep SEQS.tsv KCOUNT GCAP > rep20.tsv      repeated 20-mers that occur in the flanks
   fs9_repmask.py mask SEQS.tsv REP20.tsv > seqs.masked.tsv  flank bases covered by a repeated 20-mer become N

The repeat filter of fs9_twins.sh, the same rows as its awk version, in numbers instead of strings.
  rep : the canonical 20-mers (the smaller of the k-mer and its reverse complement, A < C < G < T) of the flanks (windows with
        an N skipped) that KCOUNT ("KMER COUNT", whitespace separated, e.g. jellyfish dump -c) lists with COUNT > GCAP; the
        k-mer of KCOUNT is upper-cased and must equal a flank k-mer as written. Exit 1 if KCOUNT has no line with two fields.
        Written one per line in the order of KCOUNT, as the awk did.
        The awk version first wrote every distinct flank 20-mer to flank20.tsv (1.3 GB on mau VES, 425,267 copies, kept in
        the results) and sorted it; here the set is held in memory as 40-bit numbers.
  mask: seqs.tsv rows (name TAB side TAB sequence) in input order, every base covered by a 20-mer of rep20.tsv as N.
"""
import sys
import numpy as np

K = 20
LUT = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate(b"ACGT"): LUT[c] = i

def read_seqs(path):
    names, rows = [], []
    with open(path, "rb") as fh:
        for line in fh:
            f = line.rstrip(b"\n").split(b"\t")
            names.append(f[:2]); rows.append(f[2] if len(f) > 2 else b"")
    return names, rows

def kmers(rows):
    """per chunk of rows: canonical 20-mer codes of every window without N: (row index, window start, code), int64 (2 bits/base)"""
    if not rows: return
    W = max(len(s) for s in rows)
    if W < K: return
    CH = 200_000
    for r0 in range(0, len(rows), CH):
        part = rows[r0:r0 + CH]
        X = LUT[np.frombuffer(b"".join(s.ljust(W, b"N") for s in part), dtype=np.uint8)].reshape(len(part), W)
        nw = W - K + 1
        isn = np.concatenate([np.zeros((len(X), 1), np.int32), np.cumsum(X == 4, axis=1, dtype=np.int32)], axis=1)
        valid = (isn[:, K:K + nw] - isn[:, 0:nw]) == 0
        f = np.zeros((len(X), nw), np.int64); r = np.zeros((len(X), nw), np.int64)
        for j in range(K):
            c = X[:, j:j + nw].astype(np.int64); c[c > 3] = 0
            f = f * 4 + c
            r = r + ((3 - c) << (2 * j))
        code = np.minimum(f, r)
        rr, ii = np.nonzero(valid)
        yield rr + r0, ii, code[rr, ii]

def sorted_unique(x):
    """np.unique by sorting (numpy 2.4's hash-based np.unique took 19 of 33 s on 92k copies)"""
    x = np.sort(x)
    return x[np.concatenate([[True], x[1:] != x[:-1]])] if len(x) else x

def encode(s):
    v = 0
    for ch in s:
        x = int(LUT[ch])
        if x > 3: return -1
        v = v * 4 + x
    return v

mode = sys.argv[1] if len(sys.argv) > 1 else ""
w = sys.stdout.buffer
if mode == "rep" and len(sys.argv) == 5:
    _, rows = read_seqs(sys.argv[2])
    parts = [sorted_unique(c) for _, _, c in kmers(rows)]; del rows
    flank = sorted_unique(np.concatenate(parts)) if parts else np.zeros(0, np.int64); del parts
    import re
    C = float(sys.argv[4]); anyline = False
    NUM = re.compile(rb"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
    def awknum(x):                     # awk's $2 + 0: the leading number of the field, 0 if none
        try: return float(x)
        except ValueError:
            m_ = NUM.match(x); return float(m_.group(0)) if m_ else 0.0
    P4 = 4 ** np.arange(K - 1, -1, -1, dtype=np.int64)
    def flush(ms):
        if not ms: return
        X = LUT[np.frombuffer(b"".join(ms), dtype=np.uint8)].reshape(len(ms), K)
        ok = (X < 4).all(axis=1)
        v = (X.astype(np.int64) * P4).sum(axis=1)
        j = np.searchsorted(flank, v); j[j >= len(flank)] = 0
        hit = ok & (flank[j] == v) if len(flank) else np.zeros(len(ms), bool)
        w.write(b"".join(ms[i] + b"\n" for i in np.flatnonzero(hit).tolist()))
    with open(sys.argv[3], "rb") as fh:
        ms = []
        for line in fh:
            f = line.split()
            if len(f) < 2: continue
            anyline = True
            m = f[0].upper()
            if len(m) != K or not awknum(f[1]) > C: continue
            ms.append(m)
            if len(ms) >= 1_000_000: flush(ms); ms = []
        flush(ms)
    sys.exit(0 if anyline else 1)
elif mode == "mask" and len(sys.argv) == 4:
    names, rows = read_seqs(sys.argv[2])
    rep = []
    with open(sys.argv[3], "rb") as fh:
        for line in fh:
            m = line.rstrip(b"\n")
            if len(m) == K:
                v = encode(m)
                if v >= 0: rep.append(v)
    rep = sorted_unique(np.array(rep, dtype=np.int64))
    cover = {}
    if len(rep):
        for rr, ii, cc in kmers(rows):
            j = np.searchsorted(rep, cc); j[j >= len(rep)] = 0
            hit = rep[j] == cc
            for r_, i_ in zip(rr[hit].tolist(), ii[hit].tolist()): cover.setdefault(r_, []).append(i_)
    for n, (nm, s) in enumerate(zip(names, rows)):
        if n in cover:
            b = bytearray(s)
            for i_ in cover[n]: b[i_:i_ + K] = b"N" * K
            s = bytes(b)
        w.write(nm[0] + b"\t" + (nm[1] if len(nm) > 1 else b"") + b"\t" + s + b"\n")
else:
    sys.exit(__doc__)
