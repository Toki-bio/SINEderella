#!/usr/bin/env python3
"""merge_hits.py - merge the hits of all consensuses into loci, strand-aware (replaces `bedtools merge -c 4,5,6 -o max,max,distinct`).

    sorted BED (chrom start end identity score strand) on stdin  ->  merged BED (same columns) on stdout

Why. The plain merge ignores the strand: a DIP hit on "-" that touches a B1 hit on "+" became ONE interval with strand "+,-". `bedtools
getfasta -s` flips only a pure "-", so that interval was written forward: the DIP part came out reverse-complemented, in a piece that
holds two elements (found in Sicista, 2026-10-13: 2.0 % of the DIP copies, 98 % of the reverse-complemented outliers of the SubFam
alignment). Here a cluster of overlapping/book-ended hits is resolved like this:

    * all hits on one strand: the old result (max identity, max score, the strand);
    * hits on both strands: the hits of each strand are merged first. Two opposite-strand intervals that overlap by at least
      --min-overlap (default 0.5) of the shorter one are ONE locus seen from both strands (a consensus matching its own reverse
      complement, a low-scoring spurious reverse hit): they become one interval (union), strand of the highest-scoring hit.
      Otherwise they are different elements next to each other (or nested insertions): each keeps its own interval and strand,
      so every extracted piece is a single element in its own orientation.

The strand column is therefore always "+" or "-". --old gives the previous behaviour (strands joined, "+,-").
A summary goes to stderr: clusters, mixed-strand clusters, merged as one locus, split into pieces.
"""
from __future__ import annotations

import argparse
import sys


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float("-inf")


def merge_same_strand(hits):
    """hits: list of (start, end, ident_str, score_str, strand) of ONE strand, sorted; overlapping or book-ended hits are merged"""
    out = []
    for s, e, idn, sc, st in hits:
        if out and s <= out[-1][1]:
            o = out[-1]
            out[-1] = [o[0], max(o[1], e), idn if fnum(idn) > fnum(o[2]) else o[2], sc if fnum(sc) > fnum(o[3]) else o[3], st]
        else:
            out.append([s, e, idn, sc, st])
    return out


def resolve_cluster(hits, min_overlap, stats):
    """hits: all hits of one cluster (any strand). Returns [(start, end, ident, score, strand)] sorted by start."""
    strands = {h[4] for h in hits}
    if len(strands) == 1:
        iv = merge_same_strand(sorted(hits))
        return [tuple(x) for x in iv]
    stats["mixed"] += 1
    ivs = []
    for st in sorted(strands):
        ivs += merge_same_strand(sorted(h for h in hits if h[4] == st))
    n = len(ivs)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            a, b = ivs[i], ivs[j]
            if a[4] == b[4]:
                continue
            ov = min(a[1], b[1]) - max(a[0], b[0])
            if ov > 0 and ov >= min_overlap * min(a[1] - a[0], b[1] - b[0]):
                parent[find(i)] = find(j)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), []).append(ivs[i])
    out = []
    for members in comp.values():
        if len(members) == 1:
            out.append(tuple(members[0]))
            continue
        best = max(members, key=lambda m: fnum(m[3]))
        idn = max((m[2] for m in members), key=fnum)
        sc = max((m[3] for m in members), key=fnum)
        out.append((min(m[0] for m in members), max(m[1] for m in members), idn, sc, best[4]))
    out.sort(key=lambda x: (x[0], x[1]))
    if len(out) == 1:
        stats["one_locus"] += 1
    else:
        stats["split"] += 1
        stats["pieces"] += len(out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--min-overlap", type=float, default=0.5, help="opposite-strand intervals overlapping by this share of the shorter one are one locus (default 0.5)")
    ap.add_argument("--old", action="store_true", help="the previous behaviour: strands joined as '+,-' (bedtools merge -c 4,5,6 -o max,max,distinct)")
    a = ap.parse_args(argv)
    stats = {"clusters": 0, "mixed": 0, "one_locus": 0, "split": 0, "pieces": 0}
    w = sys.stdout.write
    chrom = None
    cl = []
    cend = -1

    def flush():
        if not cl:
            return
        stats["clusters"] += 1
        if a.old:
            s = min(h[0] for h in cl)
            e = max(h[1] for h in cl)
            idn = max((h[2] for h in cl), key=fnum)
            sc = max((h[3] for h in cl), key=fnum)
            st = ",".join(sorted({h[4] for h in cl}))
            w("%s\t%d\t%d\t%s\t%s\t%s\n" % (chrom, s, e, idn, sc, st))
            return
        for s, e, idn, sc, st in resolve_cluster(cl, a.min_overlap, stats):
            w("%s\t%d\t%d\t%s\t%s\t%s\n" % (chrom, s, e, idn, sc, st))

    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        h = (s, e, f[3], f[4], f[5])
        if c != chrom or s > cend:
            flush()
            cl = []
            chrom = c
            cend = -1
        cl.append(h)
        cend = max(cend, e)
    flush()
    print("merge_hits: %d clusters; %d with hits on both strands: %d merged as one locus, %d split into %d pieces"
          % (stats["clusters"], stats["mixed"], stats["one_locus"], stats["split"], stats["pieces"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
