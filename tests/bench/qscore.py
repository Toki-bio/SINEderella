#!/usr/bin/env python3
"""qscore.py REF.aln TEST.aln [--elem] - fraction of the reference's aligned residue pairs (between different
rows, same column) that the test alignment reproduces (Q score, sampled over row pairs). Rows are
matched by name; MAFFT's _R_ prefix is stripped.
--elem: count only reference pairs in the element columns (first..last letter of the CONSENSUS row of
REF). Unrelated flank has no true alignment, so whole-row Q on +600 bp plates mostly measures noise."""
import random, sys


def read(p):
    d, name, cur = {}, None, []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if name:
                d[name] = "".join(cur)
            name = l[1:].split()[0].replace("_R_", "", 1); cur = []
        else:
            cur.append(l)
    d[name] = "".join(cur)
    return d


def colmap(s):
    """residue index -> column"""
    m, k = [], 0
    for j, c in enumerate(s):
        if c not in "-.":
            m.append(j)
    return m


ref, test = read(sys.argv[1]), read(sys.argv[2])
lo, hi = 0, 10 ** 9
if "--elem" in sys.argv:
    cr = next(v for k, v in ref.items() if "CONSENSUS" in k)
    let = [j for j, c in enumerate(cr) if c not in "-."]
    lo, hi = let[0], let[-1]
names = [n for n in ref if n in test]
random.seed(1)
pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
pairs = random.sample(pairs, min(400, len(pairs)))
tot = ok = 0
for a, b in pairs:
    ra, rb, ta, tb = colmap(ref[a]), colmap(ref[b]), colmap(test[a]), colmap(test[b])
    if len(ra) != len(ta) or len(rb) != len(tb):
        continue   # reversed differently - skip pair
    col_b_ref = {c: i for i, c in enumerate(rb)}
    col_b_test = {c: i for i, c in enumerate(tb)}
    for i, c in enumerate(ra):
        j = col_b_ref.get(c)
        if j is None or not lo <= c <= hi:
            continue
        tot += 1
        ok += col_b_test.get(ta[i]) == j
print("%.3f" % (ok / float(tot)) if tot else "nan")
