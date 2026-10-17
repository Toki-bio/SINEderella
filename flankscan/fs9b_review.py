#!/usr/bin/env python3
"""fs9b_review.py COMPOSITE.fa CUT.fa OUT.aln.fa [THREADS=4]

Review alignment of one fs9b group, anchored on its candidate composite (partner + family consensus, from fs9b):
every member is laid on the composite's positions (mafft --addfragments --keeplength), so the columns are exactly the
composite (a free MAFFT of members with long flanks spreads them over many gap columns: rsi r6 vs P26, 376 positions over
1,066 columns). What lies outside the composite in a member (the genomic flank beyond the partner, and after the element)
is shown lower case, degapped and packed against that member's own first / last aligned base. Bases a member has in
addition to the composite are not shown in the columns and are counted in the row name (ins=).

CUT.fa: one record per member, the member's partner-side flank + element + other-side flank, in the element's orientation
(fs9b writes it). The member's span on the composite is taken from ssearch36 (forward strand).
"""
import os
import subprocess
import sys
import tempfile

comp_fa, cut_fa, out = sys.argv[1:4]
T = sys.argv[4] if len(sys.argv) > 4 else "4"
FL = int(os.environ.get("REVIEW_FLANK", 100))


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]
            d[n] = []
        elif n:
            d[n].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}


comp = read_fa(comp_fa)
cname, cseq = next(iter(comp.items()))
cut = read_fa(cut_fa)
tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "."))
key = {"m%d" % i: n for i, n in enumerate(cut)}
with open(os.path.join(tmp, "c.fa"), "w") as fh:
    fh.write(">c\n%s\n" % cseq)
with open(os.path.join(tmp, "q.fa"), "w") as fh:
    for k, n in key.items():
        fh.write(">%s\n%s\n" % (k, cut[n]))
o = subprocess.run(["ssearch36", "-3", "-z", "11", "-m", "8", "-E", "10", "-T", T, "c.fa", "q.fa"], cwd=tmp,
                   capture_output=True, text=True).stdout
span = {}
for line in o.splitlines():
    f = line.split("\t")
    if f[1] not in span:
        span[f[1]] = (int(f[8]), int(f[9]))                 # span on the member (1-based)
core, left, right = {}, {}, {}
for k, n in key.items():
    if k not in span:
        continue
    a, b = span[k]
    core[k] = cut[n][a - 1:b]
    left[k] = cut[n][:a - 1][-FL:].lower()
    right[k] = cut[n][b:][:FL].lower()
with open(os.path.join(tmp, "f.fa"), "w") as fh:
    for k, s in core.items():
        fh.write(">%s\n%s\n" % (k, s))
aln = subprocess.run(["mafft", "--quiet", "--thread", T, "--threadit", "0", "--addfragments", "f.fa", "--keeplength", "c.fa"],
                     cwd=tmp, capture_output=True, text=True, check=True).stdout
recs, cur = {}, None
for line in aln.splitlines():
    if line.startswith(">"):
        cur = line[1:]
        recs[cur] = []
    else:
        recs[cur].append(line.upper())
W = len(cseq)
rows = []
for k in core:
    s = "".join(recs[k])
    fc = next((j for j, ch in enumerate(s) if ch != "-"), None)
    if fc is None:
        continue
    lc = W - 1 - next(j for j, ch in enumerate(s[::-1]) if ch != "-")
    ins = len(core[k]) - sum(1 for ch in s if ch != "-")
    rows.append((fc, k, s, lc, ins))
rows.sort()
L = max([len(left[k]) - fc for fc, k, _, _, _ in rows] + [0])
R = max([len(right[k]) - (W - 1 - lc) for _, k, _, lc, _ in rows] + [0])
L, R = max(L, 0), max(R, 0)
with open(out, "w") as fo:
    fo.write(">CONS_%s\n%s\n" % (cname, "-" * L + cseq + "-" * R))
    for fc, k, s, lc, ins in rows:
        row = left[k].rjust(L + fc, "-") + s[fc:lc + 1] + right[k].ljust(W - 1 - lc + R, "-")
        fo.write(">%s|ins=%d\n%s\n" % (key[k], ins, row))
print("%s: %d members on %d composite columns (+%d / +%d flank columns)" % (out, len(rows), W, L, R))
