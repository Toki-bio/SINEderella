#!/usr/bin/env python3
"""Slide the element edge out along ungapped flanks until flanks are unique.

Does not use alignment columns. Each copy's flank is the ungapped sequence
outside the current consensus, with position 0 touching the consensus edge.
Two flanks "share" when the first WIN bases agree at >= THR (the same cut
verdict.py uses: 0.55 over up to 70 bp, at least 35 bp compared).

Groups are single-linkage, same as flank_uniqueness.py. A depth d moves the
edge d bases out; the flank used for the test is whatever remains beyond d.

Stops are reported for three readings of "no 2-5 copies share":
  largest cluster <= 1   no two copies share
  largest cluster <= 2   no three copies share
  largest cluster <= 4   no five copies share
"""
import sys

WIN = 70
THR = 0.55
MINLEN = 35
OCC = 0.50
AGREE = 0.45
MISS = 8
# Largest shared group still allowed. 4 means no five copies share a flank.
CLUSTER_CAP = 4
# Extension starts only if the group sharing flank at the current edge holds at
# least MIN_SHARED of the measurable copies, and at least MIN_MEASURED copies can be
# measured. Calibrated 2026-09-27 on rsi peel + 13 bat alignments: real missing element
# 0.30-0.97 (rsi r3 top100 0.74, r6 rand100 0.30, rmi VES rand100 0.97); noise <= 0.24
# (hla Rhin-1 rand100 0.24 -> 50 bp of scattered majority bases, read as AAAAAAC).
MIN_SHARED = 0.30
MIN_MEASURED = 20
# A run of this many A's is polyA, not element. The 3' A-tail already inside
# the seed consensus is left alone; this only stops an extension, and a run
# that is already the 5' end of the consensus.
POLYA = 5


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:].split()[0]
            names.append(cur)
            buf = []
        else:
            buf.append(line)
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def consensus_index(names):
    for i, n in enumerate(names):
        if ":" not in n:
            return i
    return 0


def element_bounds(seq):
    upper = [i for i, c in enumerate(seq) if c.isupper()]
    if upper:
        return upper[0], upper[-1]
    nz = [i for i, c in enumerate(seq) if c not in "-."]
    return (nz[0], nz[-1]) if nz else (0, len(seq) - 1)


def flanks_of(seqs, ci, skip=()):
    """skip: row indices that are not copies (the seed-as-searched row)."""
    lo, hi = element_bounds(seqs[ci])
    lefts, rights = [], []
    for i, s in enumerate(seqs):
        if i == ci or i in skip:
            continue
        left = "".join(c for c in s[:lo] if c not in "-.").upper()
        right = "".join(c for c in s[hi + 1:] if c not in "-.").upper()
        lefts.append(left[::-1])
        rights.append(right)
    return lefts, rights, lo, hi


def pair_identity(a, b):
    m = min(len(a), len(b), WIN)
    if m < MINLEN:
        return 0.0
    same = 0
    for i in range(m):
        if a[i] == b[i]:
            same += 1
    return same / float(m)


def largest_cluster(flanks):
    """Single-linkage at THR. Copies shorter than MINLEN are not measured."""
    idx = [i for i, f in enumerate(flanks) if len(f) >= MINLEN]
    n = len(idx)
    if n < 2:
        return 1, n, []
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a in range(n):
        fa = flanks[idx[a]]
        for b in range(a + 1, n):
            if pair_identity(fa, flanks[idx[b]]) >= THR:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
    sizes = {}
    for i in range(n):
        r = find(i)
        sizes[r] = sizes.get(r, 0) + 1
    ordered = sorted(sizes.values(), reverse=True)
    return ordered[0], n, ordered[:6]


def slide(flanks):
    if not flanks:
        return []
    max_d = max(len(f) for f in flanks)
    half = len(flanks) / 2.0
    rows = []
    for d in range(0, max_d + 1):
        rem = [f[d:] for f in flanks]
        largest, n_ok, top = largest_cluster(rem)
        rows.append((d, largest, n_ok, top))
        if n_ok < half:
            break
    return rows


def first_depth(rows, cap, n_copies):
    """Smallest d where the largest shared group is <= cap and most copies remain."""
    half = n_copies / 2.0
    for d, largest, n_ok, _top in rows:
        if n_ok >= half and largest <= cap:
            return d
    return None


def conservation_tail(flanks, d):
    """The long window goes unique while its first bases are still shared.

    From d, keep bases that half the copies carry and that agree, and swallow
    fewer than MISS misses. This is one ungapped position at a time, anchored
    at the element edge, not an alignment column.
    """
    n = len(flanks)
    if n == 0:
        return d
    limit = max(len(f) for f in flanks)
    misses = 0
    edge = d
    k = d
    while k < limit and misses < MISS:
        bases = [f[k] for f in flanks if len(f) > k and f[k] in "ACGT"]
        occ = len(bases) / float(n)
        agree = 0.0
        if bases:
            counts = {}
            for b in bases:
                counts[b] = counts.get(b, 0) + 1
            agree = max(counts.values()) / float(len(bases))
        if occ >= OCC and agree >= AGREE:
            edge = k + 1
            misses = 0
        else:
            misses += 1
        k += 1
    return edge


def trim_outer_polya(flanks, edge, min_run=POLYA):
    """Drop a terminal A-run at the far end of the extension.

    Index 0 touches the old element, so the far base is at edge-1. Gaps do not
    break the run. A shorter run is a real base and stays.
    """
    if edge <= 0:
        return 0
    ext = majority(flanks, edge)
    k = edge
    while k > 0 and ext[k - 1] == "-":
        k -= 1
    run = 0
    while k > 0 and ext[k - 1] in ("A", "-"):
        if ext[k - 1] == "A":
            run += 1
        k -= 1
        if run >= min_run and (k == 0 or ext[k - 1] not in ("A", "-")):
            return k
    return edge


def stop_before_polya(flanks, edge, min_run=POLYA):
    """Do not extend through a polyA run. The run itself stays flank."""
    if edge <= 0:
        return 0
    ext = majority(flanks, edge)
    run = 0
    start = 0
    for k in range(edge):
        if ext[k] == "A":
            if run == 0:
                start = k
            run += 1
            if run >= min_run:
                # Gaps immediately outside the old edge are not sequence.
                # Extending across them would pull the polyA back in.
                if all(ext[i] == "-" for i in range(start)):
                    return 0
                return start
        elif ext[k] != "-":
            run = 0
    return edge


def extension_edge(flanks, cap=CLUSTER_CAP, polya="outer"):
    """How far the element can grow. 0 means the current edge stands.

    polya='stop' refuses to cross a run of A's.
    polya='outer' only cuts an A-run that is the far end of the extension.

    If the extracted flank never becomes unique, the shared sequence is still
    the element. Leaving it as flank is how r10's 3' tail and the r3 SubFam
    head stayed outside the consensus.
    """
    rows = slide(flanks)
    # Missing element is carried by a large share of copies (rsi r3 head: 74 of 100
    # top copies). A small group sharing flank (copies sitting in the same older
    # repeat, say) is not element, and an edge that is already unique must not creep
    # outward on per-position agreement alone (it did by 1-16 bp on most files).
    if not rows or rows[0][2] < MIN_MEASURED or rows[0][1] < MIN_SHARED * rows[0][2]:
        return 0
    d = first_depth(rows, cap, len(flanks))
    if d is None:
        d = 0
    edge = conservation_tail(flanks, d)
    if polya == "stop":
        edge = stop_before_polya(flanks, edge)
    else:
        edge = trim_outer_polya(flanks, edge)
    return edge


def majority(flanks, d):
    """Bases at distance 0..d-1 from the old edge. Distance 0 touches the edge."""
    out = []
    n = len(flanks)
    for k in range(d):
        bases = [f[k] for f in flanks if len(f) > k and f[k] in "ACGT"]
        if not bases:
            out.append("-")
            continue
        counts = {}
        for b in bases:
            counts[b] = counts.get(b, 0) + 1
        top = max(counts, key=counts.get)
        occ = len(bases) / float(n)
        agree = counts[top] / float(len(bases))
        out.append(top if occ >= OCC and agree >= AGREE else "-")
    return "".join(out)


def summarize(name, flanks):
    rows = slide(flanks)
    n = len(flanks)
    print("== %s  copies=%d  median_flank=%d ==" % (
        name, n, sorted(len(f) for f in flanks)[n // 2]))
    if not rows:
        print("  no flanks")
        return
    d0, L0, n0, top0 = rows[0]
    print("  at the current edge: largest shared group %d of %d measured  sizes %s" % (
        L0, n0, top0))
    for cap, label in ((1, "no two share"), (2, "no three share"), (4, "no five share")):
        d = first_depth(rows, cap, n)
        if d is None:
            print("  %s: not reached inside the extracted flank" % label)
            continue
        edge = conservation_tail(flanks, d)
        ext = majority(flanks, edge)
        called = ext.replace("-", "")
        print("  %s: window unique at %d bp, edge at %d bp after the shared tail" % (
            label, d, edge))
        print("    consensus letters %d/%d" % (len(called), edge))
        seq5 = called[::-1]
        print("    5' end: %s" % seq5[:90])
        print("    junction into the old edge: %s" % seq5[-40:])
    # a short trace so the drop is visible
    print("  depth  largest  measured")
    marked = set()
    for d, largest, n_ok, _top in rows:
        cross = None
        for cap in (4, 2, 1):
            if largest <= cap and cap not in marked:
                marked.add(cap)
                cross = cap
        if d % 25 == 0 or cross is not None:
            print("  %5d  %8d  %9d" % (d, largest, n_ok))
        if 1 in marked:
            break


def _seq(seed, n):
    x = seed & 0x7FFFFFFF
    out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        out.append("ACGT"[(x >> 16) % 4])
    return "".join(out)


def _edge(flanks, cap):
    rows = slide(flanks)
    d = first_depth(rows, cap, len(flanks))
    if d is None:
        return None
    return conservation_tail(flanks, d)


def synthetic():
    shared = _seq(11, 80)
    lefts = [shared + _seq(1000 + i * 997, 90) for i in range(40)]
    print("SYNTHETIC shared 80 bp then unique")
    summarize("synthetic", lefts)
    got = _edge(lefts, 4)
    if got != 80:
        raise SystemExit("synthetic shared-80 failed: edge %s" % got)

    # shared element, then a 5' polyA. The A-run must not enter the consensus.
    shared = ("CGTGAC" * 20)[:70]
    lefts = [shared + ("A" * 12) + _seq(8000 + i * 997, 40) for i in range(40)]
    edge = extension_edge(lefts, polya="stop")
    if edge != 70:
        raise SystemExit("synthetic polya failed: edge %s" % edge)

    a, b = _seq(21, 60), _seq(22, 60)
    lefts = []
    for i in range(6):
        lefts.append(a + _seq(3000 + i * 997, 80))
    for i in range(6):
        lefts.append(b + _seq(4000 + i * 997, 80))
    for i in range(28):
        lefts.append(_seq(5000 + i * 997, 140))
    print("SYNTHETIC two groups of 6 share 60 bp, 28 already unique")
    summarize("two-groups", lefts)
    got = _edge(lefts, 4)
    # two different shared sequences do not make one majority, so the edge
    # must stop inside the shared block and must not run into the unique tails
    if got is None or not (34 <= got <= 60):
        raise SystemExit("synthetic two-groups failed: edge %s" % got)

    # The extract ends while every copy still shares the flank. That sequence
    # is element. A terminal polyA is not.
    shared = ("CGTGAC" * 20)[:100]
    lefts = [shared + ("A" * 15) for _ in range(40)]
    edge = extension_edge(lefts, polya="stop")
    if edge != 100:
        raise SystemExit("synthetic still-shared failed: edge %s" % edge)


def main(argv):
    synthetic()
    for path in argv[1:]:
        names, seqs = read_fa(path)
        ci = consensus_index(names)
        lefts, rights, lo, hi = flanks_of(seqs, ci)
        print("FILE %s  consensus=%s  element cols %d-%d" % (
            path, names[ci], lo, hi))
        summarize("left", lefts)
        summarize("right", rights)


if __name__ == "__main__":
    main(sys.argv)
