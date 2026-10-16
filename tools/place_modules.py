#!/usr/bin/env python3
"""Place the declared modules of each composite on the composite consensus, in order and without overlap.

Usage: place_modules.py BANK.fa LAYOUT.tsv > composites_explicit.tsv
LAYOUT.tsv: composite <TAB> unit:from-to,unit:from-to,...   (coordinates on the UNIT; 'unit' alone = whole unit)
Each module part is aligned into the composite (ssearch36, all HSPs >= 30 bp, E <= 1e-3); one HSP per module is
chosen so that the modules follow the declared order without overlapping by more than 10 bp, maximising the summed
bitscore. Output: composite <TAB> unit:start-end,unit:start-end,... on COMPOSITE coordinates, plus a report on stderr.
"""
import itertools
import os
import subprocess
import sys
import tempfile

bank, layout = sys.argv[1], sys.argv[2]
seqs, name = {}, None
for l in open(bank):
    l = l.strip()
    if l.startswith(">"):
        name = l[1:].split()[0]
        seqs[name] = []
    elif name:
        seqs[name].append(l)
seqs = {k: "".join(v) for k, v in seqs.items()}
tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "."))


def hsps(part, comp):
    q, l = os.path.join(tmp, "q.fa"), os.path.join(tmp, "l.fa")
    open(q, "w").write(">q\n%s\n" % part)
    open(l, "w").write(">c\n%s\n" % seqs[comp])
    out = subprocess.run(["ssearch36", "-m", "8", "-E", "1e-3", "-z", "-1", "-T", "1", "q.fa", "l.fa"], cwd=tmp,
                         capture_output=True, text=True).stdout
    res = []
    for line in out.splitlines():
        f = line.split("\t")
        a, b = sorted((int(f[8]), int(f[9])))
        if b - a + 1 >= 30 and int(f[8]) <= int(f[9]):          # forward orientation only
            res.append((a, b, float(f[11]), float(f[2])))
    return res


for line in open(layout):
    if not line.strip() or line.startswith("#"):
        continue
    comp, mods = line.rstrip("\n").split("\t")[:2]
    cand = []
    for m in mods.split(","):
        u, _, rng = m.partition(":")
        s = seqs[u]
        if rng:
            f, t = (int(x) for x in rng.split("-"))
            s = s[f - 1:t]
        cand.append((u, rng or "all", hsps(s, comp)))
    best, bestsc = None, -1
    for combo in itertools.product(*[c[2] for c in cand]):
        if all(combo[k][0] > combo[k - 1][1] - 10 for k in range(1, len(combo))):
            sc = sum(h[2] for h in combo)
            if sc > bestsc:
                best, bestsc = combo, sc
    if best is None:
        print("%s: no ordered placement of %s" % (comp, mods), file=sys.stderr)
        continue
    print("%s\t%s" % (comp, ",".join("%s:%d-%d" % (cand[k][0], best[k][0], best[k][1]) for k in range(len(cand)))))
    print("%s (%d bp): %s" % (comp, len(seqs[comp]), "  ".join("%s[%s] -> %d-%d (%.1f %%)" % (cand[k][0], cand[k][1], best[k][0], best[k][1], best[k][3])
                                                           for k in range(len(cand)))), file=sys.stderr)
