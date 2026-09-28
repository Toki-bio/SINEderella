#!/usr/bin/env python3
"""colq.py ALN... - self-quality of each alignment, reference-free: fraction of residues that equal their
column's majority (columns with >=3 residues), over the element span (CONSENSUS row first..last letter).
Higher = more consistent columns. Lets two modes be compared when Q between them is low: equal scores =
equally good alternative alignments; lower = the mode made worse columns."""
import sys
from collections import Counter
def read(p):
    d, n, cur = [], None, []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n: d.append((n, "".join(cur)))
            n = l[1:]; cur = []
        else: cur.append(l)
    d.append((n, "".join(cur))); return d
for p in sys.argv[1:]:
    rows = read(p)
    cr = next(s for n, s in rows if "CONSENSUS" in n)
    let = [j for j, c in enumerate(cr) if c not in "-."]
    lo, hi = let[0], let[-1]
    seqs = [s.upper() for n, s in rows if "CONSENSUS" not in n]
    tot = ok = 0
    for j in range(lo, hi + 1):
        col = [s[j] for s in seqs if s[j] not in "-."]
        if len(col) < 3: continue
        m = Counter(col).most_common(1)[0][1]
        tot += len(col); ok += m
    print("%-50s width %5d  elem cols %5d  majority-agreement %.3f" % (p.split("sets")[-1], len(cr), hi - lo + 1, ok / tot))
