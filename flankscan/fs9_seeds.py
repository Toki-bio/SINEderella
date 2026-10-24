#!/usr/bin/env python3
"""fs9_seeds.py SEQS.tsv PATTERN MK [--iw 80] [--cap 20] [--band 3] [--start 20] [--slack 25] [--threads 1]

Candidate pairs of one fs9_twins.sh pass (seeds -> pairs of copies sharing seeds -> >= MK seeds on one diagonal band), the
same rows as the awk pipeline it replaces (seed rows keys.TAG.tsv -> every pair of a seed bucket hits.TAG.tsv -> stable sort
by pair -> one pair at a time). Rows go to stdout unsorted; fs9_twins.sh pipes them through sort -u as before.

Why: on tbr VES (620,807 copies) the awk version wrote hits.p1.tsv 45 GB and hits.p2.tsv 92 GB (one text row per seed shared
by two copies) plus a sort of the same size, and took hours. Here the seed rows are held as numbers (about 20 bytes each) and
the pairs of one block of copies are counted in memory and dropped; nothing but the candidate rows is written.

Rules reproduced exactly (fs9_twins.sh pass_run, checked row for row against its output, tests/test_fs9_seeds.py):
  seeds   windows i = 1..IW of the first IW + len(PATTERN) - 1 bases of a flank; a window with an N is skipped; the seed is
          the window's bases at the 1-positions of PATTERN; a seed with at most 2 distinct bases is skipped
  bucket  all rows of one seed; a row whose (seed, offset) occurs more than CAP times in the bucket (both sides counted) is not
          paired; two rows pair when they are of one flank side and of two different copies; copy a = the smaller name
          (byte order), delta = offset in a - offset in b
  pair    one candidate per (side, a, b): best = the largest number of the pair's seeds within delta +- BAND of one of its deltas;
          bd = the delta of the first seed (in the awk's row order) that reaches best; kept when best >= MK and the smallest
          offset of a and of b are both <= START + SLACK
  order   the awk's row order, which decides bd on ties, is: seeds in byte order; inside a seed the rows sorted by the whole
          line "seed TAB offset TAB copy TAB side" (sort's last-resort comparison: offsets compare as strings); the pairs of a
          bucket in the order i < j of those rows. Here every pair carries (position of row i, position of row j) as its key.
Output rows: side TAB copy_a TAB copy_b TAB bd TAB best.
"""
import argparse, os, sys
import numpy as np

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("seqs"); ap.add_argument("pattern"); ap.add_argument("mk", type=int)
ap.add_argument("--iw", type=int, default=80); ap.add_argument("--cap", type=int, default=20)
ap.add_argument("--band", type=int, default=3); ap.add_argument("--start", type=int, default=20)
ap.add_argument("--slack", type=int, default=25)
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--block-hits", type=int, default=8_000_000, help="seed pairs counted at once by one worker (memory ~ 60 bytes each)")
a = ap.parse_args()

PAT = a.pattern
if not PAT or set(PAT) - {"0", "1"} or PAT[0] != "1":
    sys.exit("fs9_seeds: PATTERN must be 0/1 and start with 1")
L = len(PAT); POS = [j for j, c in enumerate(PAT) if c == "1"]; W2 = a.iw + L - 1

# ---- flanks: name, side, sequence (seqs.tsv of fs9_twins.sh: upper case ACGTN)
names, sides, rows = [], [], []
with open(a.seqs, "rb") as fh:
    for line in fh:
        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 2: continue
        names.append(f[0]); sides.append(f[1]); rows.append(f[2] if len(f) > 2 else b"")
R = len(rows)
if R == 0: sys.exit(0)
uniq = sorted(set(names))                      # bytes sort = the awk's string comparison in the C locale
rank = {n: i for i, n in enumerate(uniq)}
NC = len(uniq)
copy_of_row = np.fromiter((rank[n] for n in names), dtype=np.int64, count=R)
side_vals = sorted(set(sides))                 # "3" < "5": the side is the last field of the sorted line
if len(side_vals) > 2: sys.exit("fs9_seeds: more than two flank sides in %s" % a.seqs)
side_of_row = np.fromiter((side_vals.index(s) for s in sides), dtype=np.int64, count=R)

LUT = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate(b"ACGT"): LUT[c] = i
ok_chars = set(b"ACGTN")
buf = bytearray()
for s in rows:
    if set(s[:W2]) - ok_chars:
        sys.exit("fs9_seeds: unexpected character in a flank (only ACGTN): %r" % s[:W2])
    buf += s[:W2].ljust(W2, b"N")              # padding as N: a window past the end of the flank is skipped, as i + L - 1 > n is
del rows
S = LUT[np.frombuffer(bytes(buf), dtype=np.uint8)].reshape(R, W2); del buf

# offsets 1..IW ranked as strings ("1" < "10" < "11" < ... < "2"): the order of sort's whole-line comparison inside a seed
offrank = np.zeros(a.iw + 1, dtype=np.int64)
for r_, o in enumerate(sorted(range(1, a.iw + 1), key=str)): offrank[o] = r_
POP = np.array([bin(x).count("1") for x in range(16)], dtype=np.int8)

seed_l, off_l, row_l = [], [], []
CH = 100_000
for r0 in range(0, R, CH):
    X = S[r0:r0 + CH]
    isn = np.concatenate([np.zeros((len(X), 1), np.int32), np.cumsum(X == 4, axis=1, dtype=np.int32)], axis=1)
    i0 = np.arange(a.iw)
    valid = (isn[:, i0 + L] - isn[:, i0]) == 0
    seed = np.zeros((len(X), a.iw), dtype=np.int64); msk = np.zeros((len(X), a.iw), dtype=np.int64)
    for p in POS:
        c = X[:, i0 + p].astype(np.int64)
        seed = seed * 4 + np.where(c < 4, c, 0); msk |= np.where(c < 4, 1 << c, 0)
    valid &= POP[msk] > 2
    rr, ii = np.nonzero(valid)
    seed_l.append(seed[rr, ii]); off_l.append(ii + 1); row_l.append(rr + r0)
del S
seed = np.concatenate(seed_l); off = np.concatenate(off_l).astype(np.int64); row = np.concatenate(row_l).astype(np.int64)
del seed_l, off_l, row_l
cop = copy_of_row[row]; sid = side_of_row[row]; del row

# ---- bucket order of the awk (seed, then the rest of the line), CAP per (seed, offset)
o = np.lexsort((sid, cop, offrank[off], seed))
seed, off, cop, sid = seed[o], off[o], cop[o], sid[o]; del o
N = len(seed)
if N == 0: sys.exit(0)
brk = np.flatnonzero((seed[1:] != seed[:-1]) | (off[1:] != off[:-1])) + 1
st = np.concatenate([[0], brk]); ln = np.diff(np.concatenate([st, [N]]))
keep = np.repeat(ln <= a.cap, ln)
seed, off, cop, sid = seed[keep], off[keep], cop[keep], sid[keep]; del keep, brk, st, ln
K = len(seed)
if K == 0: sys.exit(0)
pos = np.arange(K, dtype=np.int64)            # position in the awk's row order (rows not paired removed; order kept)

# ---- partner index: rows grouped by (seed, side), copies ascending; partners of a row = rows of its group with a larger copy
o2 = np.lexsort((cop, sid, seed))
seed2, sid2, cop2, off2, pos2 = seed[o2], sid[o2], cop[o2], off[o2], pos[o2]; del seed, off, cop, sid, pos, o2
g_brk = np.flatnonzero((seed2[1:] != seed2[:-1]) | (sid2[1:] != sid2[:-1])) + 1
c_brk = np.flatnonzero((seed2[1:] != seed2[:-1]) | (sid2[1:] != sid2[:-1]) | (cop2[1:] != cop2[:-1])) + 1
def run_end(brks):
    ends = np.concatenate([brks, [K]]); lens = np.diff(np.concatenate([[0], ends])); return np.repeat(ends, lens)
ge = run_end(g_brk); ub = run_end(c_brk); del g_brk, c_brk
del seed2
cnt = ge - ub                                   # partners of each row
hits_per_copy = np.bincount(cop2, weights=cnt, minlength=NC)
order_by_copy = np.argsort(cop2, kind="stable")  # rows of one copy together
copy_start = np.searchsorted(cop2[order_by_copy], np.arange(NC + 1))

# blocks of whole copies of about --block-hits seed pairs
blocks, acc, b0 = [], 0, 0
for c in range(NC):
    acc += hits_per_copy[c]
    if acc >= a.block_hits:
        blocks.append((b0, c + 1)); b0, acc = c + 1, 0
if b0 < NC: blocks.append((b0, NC))
LIM = a.start + a.slack

def do_block(blk):
    c0, c1 = blk
    T = order_by_copy[copy_start[c0]:copy_start[c1]]
    T = T[cnt[T] > 0]
    if len(T) == 0: return b""
    n = cnt[T]; tot = int(n.sum())
    t_rep = np.repeat(T, n)
    first = np.cumsum(n) - n
    p = np.repeat(ub[T] - first, n) + np.arange(tot, dtype=np.int64)
    ca, cb = cop2[t_rep], cop2[p]
    ofa, ofb = off2[t_rep], off2[p]
    pa, pb = pos2[t_rep], pos2[p]
    ordk = np.minimum(pa, pb) * K + np.maximum(pa, pb)
    pk = (ca * NC + cb) * 2 + sid2[t_rep]
    key = pk * 256 + (ofa - ofb + 128)
    del t_rep, p, pa, pb, ca, cb
    o = np.argsort(key, kind="stable")
    key, ordk, ofa, ofb = key[o], ordk[o], ofa[o], ofb[o]; del o
    gb = np.flatnonzero(key[1:] != key[:-1]) + 1
    gs = np.concatenate([[0], gb])
    gkey = key[gs]; gcnt = np.diff(np.concatenate([gs, [tot]]))
    gfirst = np.minimum.reduceat(ordk, gs); gma = np.minimum.reduceat(ofa, gs); gmb = np.minimum.reduceat(ofb, gs)
    del key, ordk, ofa, ofb
    cum = np.concatenate([[0], np.cumsum(gcnt)])
    s = cum[np.searchsorted(gkey, gkey + a.band, "right")] - cum[np.searchsorted(gkey, gkey - a.band, "left")]
    gpk = gkey >> 8
    pb_ = np.flatnonzero(gpk[1:] != gpk[:-1]) + 1
    ps = np.concatenate([[0], pb_]); plen = np.diff(np.concatenate([ps, [len(gkey)]]))
    best = np.maximum.reduceat(s, ps); ma = np.minimum.reduceat(gma, ps); mb = np.minimum.reduceat(gmb, ps)
    BIG = np.iinfo(np.int64).max
    cand = np.where(s == np.repeat(best, plen), gfirst, BIG)
    firstbest = np.minimum.reduceat(cand, ps)
    sel = cand == np.repeat(firstbest, plen)                       # one group per pair (keys of seed pairs are unique)
    good = (best >= a.mk) & (ma <= LIM) & (mb <= LIM)
    g_ok = sel & np.repeat(good, plen)
    out = []
    for k_, d_, b_ in zip(gkey[g_ok].tolist(), ((gkey[g_ok] & 255) - 128).tolist(), np.repeat(best, plen)[g_ok].tolist()):
        pk_ = k_ >> 8; sd = pk_ & 1; ab = pk_ >> 1
        out.append(b"%s\t%s\t%s\t%d\t%d\n" % (side_vals[sd], uniq[ab // NC], uniq[ab % NC], d_, b_))
    return b"".join(out)

w = sys.stdout.buffer
if a.threads > 1 and len(blocks) > 1:
    import multiprocessing as mp
    with mp.get_context("fork").Pool(min(a.threads, len(blocks))) as pool:
        for chunk in pool.imap_unordered(do_block, blocks):
            w.write(chunk)
else:
    for blk in blocks: w.write(do_block(blk))
w.flush()
