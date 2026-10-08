#!/usr/bin/env python3
"""Score partitions of labelled copies.

usage: score.py LABELS.tsv  SPEC:LABEL [SPEC:LABEL ...]
SPEC is  ASSIGN@NAMES  for COSEG (its .assign file and the .names written by to_coseg.py),
or a SubFam .chunks.tsv (copy id, chunk, strand), or any two-column TSV (copy id, group).

Three views of every method (a method may leave copies unplaced: Price's rule drops truncated
copies, the peel leaves a residue):
  placed    metrics on the copies the method put in a group (coverage = placed / all copies);
  common    the same metrics on the copies placed by every method listed (a like-for-like set,
            small when one method places few copies);
  all       ARI over ALL copies, each unplaced copy counting as a group of its own, so a method
            cannot gain by discarding hard copies.
Columns: method, placed, groups | purity, homogeneity, completeness, V, ARI (placed) |
         ARIcom, Vcom, groups (common) | ARIall.
"""
import collections
import sys
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure

lab = dict(l.rstrip("\n").split("\t")[:2] for l in open(sys.argv[1]))


def load(spec):
    path, label = spec.rsplit(":", 1)
    asg = {}
    if "@" in path:                                   # COSEG: ASSIGN@NAMES (one assignment per name, same order)
        path, names_path = path.split("@", 1)
        names = [l.strip() for l in open(names_path)]
        rows = [l.split() for l in open(path)]
        assert len(rows) == len(names), "%s: %d assignments for %d names" % (path, len(rows), len(names))
        for n, r in zip(names, rows):
            asg[n] = "sf" + r[-1]
    else:                                             # SubFam chunks.tsv or two-column TSV
        for l in open(path):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] in lab:
                asg[f[0]] = f[1]
    return label, asg


methods = [load(s) for s in sys.argv[2:]]
allids = list(lab)
common = [i for i in allids if all(i in a for _, a in methods)]


def metrics(ids, asg):
    if not ids:
        return 0, 0.0, 0.0, 0.0, 0.0, 0.0
    tr = [lab[i] for i in ids]
    cl = [asg[i] for i in ids]
    h, c, v = homogeneity_completeness_v_measure(tr, cl)
    pur = sum(max(collections.Counter(t for t, k in zip(tr, cl) if k == g).values()) for g in set(cl)) / len(ids)
    return len(set(cl)), pur, h, c, v, adjusted_rand_score(tr, cl)


def ari_all(asg):
    tr = [lab[i] for i in allids]
    cl = [asg.get(i, "u%d" % k) for k, i in enumerate(allids)]
    return adjusted_rand_score(tr, cl)


print("copies %d, truth groups %d, placed by all methods %d" % (len(allids), len(set(lab.values())), len(common)))
print("%-20s %6s %6s | %6s %6s %6s %6s %6s | %6s %6s %6s | %6s" % (
    "method", "placed", "groups", "purity", "homog", "compl", "V", "ARI", "ARIcom", "Vcom", "grpcom", "ARIall"))
for label, asg in methods:
    ids = [i for i in allids if i in asg]
    g, p, h, c, v, a = metrics(ids, asg)
    gc, pc, hc, cc, vc, ac = metrics(common, asg)
    print("%-20s %6d %6d | %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %6.3f %6d | %6.3f" % (
        label, len(ids), g, p, h, c, v, a, ac, vc, gc, ari_all(asg)))
