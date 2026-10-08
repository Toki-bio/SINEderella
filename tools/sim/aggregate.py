#!/usr/bin/env python3
"""Collect scores.txt of every OUTDIR/<scenario>_s<seed>/ into one table (mean over seeds).

usage: aggregate.py OUTDIR     -> TSV on stdout, one row per scenario and method.
Columns follow score.py: placed groups | purity homog compl V ARI | ARIcom Vcom grpcom | ARIall.
Runs whose scores.txt is missing a method are averaged over the runs that have it ("runs" says how many).
"""
import collections
import glob
import os
import re
import sys

rows = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*_s*", "scores.txt"))):
    scen = re.sub(r"_s\d+$", "", os.path.basename(os.path.dirname(f)))
    for line in open(f):
        parts = line.split("|")
        if len(parts) != 4 or line.startswith("method"):
            continue
        a = parts[0].split()
        b = parts[1].split()
        c = parts[2].split()
        d = parts[3].split()
        if len(a) < 3 or len(b) < 5 or len(c) < 3 or len(d) < 1:
            continue
        # placed, groups, purity, homog, compl, V, ARI, ARIcom, Vcom, grpcom, ARIall
        rows[(scen, a[0])].append([float(a[1]), float(a[2])] + [float(x) for x in b] + [float(x) for x in c] + [float(d[0])])

print("scenario\tmethod\truns\tplaced\tgroups\tpurity\thomog\tcompl\tV\tARI\tARIcom\tVcom\tARIall")
for (scen, meth), v in sorted(rows.items()):
    n = len(v)
    m = [sum(x[i] for x in v) / n for i in range(len(v[0]))]
    print("%s\t%s\t%d\t%.0f\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (
        scen, meth, n, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[10]))
