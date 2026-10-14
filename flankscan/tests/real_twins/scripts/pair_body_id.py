#!/usr/bin/env python3
"""Body identity of fs9 twin pairs: independent of the flank test.
usage: pair_body_id.py GENOME.fa COPIES.bed flank_groups.tsv OUT.tsv
For each pair (copy_a, copy_b) in flank_groups.tsv: the two SINE bodies (copy orientation) aligned with ssearch36
(both strands, best HSP); identity and aligned length. Out: copy_a copy_b side flank_identity tier class body_identity body_alnlen"""
import subprocess, sys, os, tempfile
G, BED, FG, OUT = sys.argv[1:5]
loc = {}
for l in open(BED):
    f = l.split("\t"); loc[f[3]] = (f[0], int(f[1]), int(f[2]), f[5].strip())
pairs = []
for i, l in enumerate(open(FG)):
    if i == 0: continue
    f = l.rstrip("\n").split("\t"); pairs.append(f)
need = sorted({p[0] for p in pairs} | {p[1] for p in pairs})
tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(OUT)))
reg = os.path.join(tmp, "r.bed")
with open(reg, "w") as fh:
    for c in need:
        ch, s, e, st = loc[c]; fh.write(f"{ch}\t{s}\t{e}\t{c}\t0\t{st}\n")
fa = subprocess.run(["bedtools", "getfasta", "-fi", G, "-bed", reg, "-s", "-nameOnly"], capture_output=True, text=True).stdout
seq = {}; name = None
for l in fa.splitlines():
    if l.startswith(">"): name = l[1:].split("(")[0]; seq[name] = ""
    else: seq[name] += l.strip().upper()
done = {}
with open(OUT, "w") as out:
    out.write("copy_a\tcopy_b\tside\tflank_identity\ttier\tclass\tbody_identity\tbody_alnlen\n")
    for p in pairs:
        key = (p[0], p[1])
        if key not in done:
            a = os.path.join(tmp, "a.fa"); b = os.path.join(tmp, "b.fa")
            open(a, "w").write(">a\n" + seq[p[0]] + "\n"); open(b, "w").write(">b\n" + seq[p[1]] + "\n")
            r = subprocess.run(["ssearch36", "-m", "8", "-E", "10", "-b", "1", "-d", "0", a, b], capture_output=True, text=True).stdout
            best = None
            for l in r.splitlines():
                if l.startswith("#") or not l.strip(): continue
                f = l.split("\t"); best = (float(f[2]), int(f[3])); break
            done[key] = best or (0.0, 0)
        bi, bl = done[key]
        out.write("\t".join([p[0], p[1], p[2], p[3], p[5], p[6], "%.1f" % bi, str(bl)]) + "\n")
