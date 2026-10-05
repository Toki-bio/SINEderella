#!/usr/bin/env python3
"""retally.py VOTEDIR CONS.fa FAMILIES.tsv [--label NAME]

Re-tally the ten cycles written by vote_full.sh three ways, on the same scores:
  flat     today's step 2: the consensus with the highest bitscore wins the cycle; firm = the same consensus wins all 10
  family   layer 1: a family's score in a cycle = the best score of any of its members; firm = the same family wins all 10
  subfam   layer 2, only for family-firm copies: the best member of that family wins the cycle; firm = the same member all 10
FAMILIES.tsv: consensus <tab> family (a consensus missing from the file is its own family).
Composite names (two or more unit names joined, e.g. r1_r3_P1, r1_r2_r3_r3_C11) are recognised by the _P<n> / _C<n> suffix.
For copies whose vote is split between two members, the pair is classified:
  composite-unit  one member is a composite and the copy's alignment covers < 60 % of the composite (a lone unit tying with the
                  composite that contains it: local alignment gives both the same score)
  composite-both  one member is a composite and the copy covers >= 60 % of it
  monomers        both members are single units (subfamily-like or length versions)
Also: near-tie cycles (top two within 2 % of each other) and the threshold question (0.45 x 10th-best summed score per family
versus per member, as step 2 does now).
"""
import collections
import os
import re
import sys

vd, cons_fa, fam_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
label = sys.argv[sys.argv.index("--label") + 1] if "--label" in sys.argv else os.path.basename(vd.rstrip("/"))

clen, name = {}, None
for l in open(cons_fa):
    l = l.strip()
    if l.startswith(">"):
        name = l[1:].split()[0]; clen[name] = 0
    elif name:
        clen[name] += len(l)
fam = {c: c for c in clen}
if os.path.exists(fam_tsv):
    for l in open(fam_tsv):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] in fam:
            fam[f[0]] = f[1]
COMP = re.compile(r"_[PC]\d+$")
is_comp = {c: bool(COMP.search(c)) for c in clen}

# cycles[c][copy] = {cons: (bits, cs, ce)}
cycles = []
for n in range(1, 11):
    p = os.path.join(vd, "cyc_%d.tsv" % n)
    d = collections.defaultdict(dict)
    if os.path.exists(p):
        for l in open(p):
            f = l.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            d[f[0]][f[1]] = (float(f[2]), int(f[3]), int(f[4]))
    cycles.append(d)
copies = set()
for d in cycles:
    copies |= set(d)


def winner(scores):
    if not scores:
        return None, 0.0, None, 0.0
    s = sorted(scores.items(), key=lambda x: -x[1])
    return s[0][0], s[0][1], (s[1][0] if len(s) > 1 else None), (s[1][1] if len(s) > 1 else 0.0)


res = {}
near_tie = collections.Counter()
for cp in copies:
    fw, famw, subw, bits_member, bits_fam, cov = [], [], [], collections.defaultdict(float), collections.defaultdict(float), {}
    for d in cycles:
        sc = d.get(cp, {})
        flat = {c: v[0] for c, v in sc.items()}
        w, b1, w2, b2 = winner(flat)
        fw.append(w)
        if w and w2 and b1 > 0 and (b1 - b2) / b1 < 0.02:
            near_tie[(fam[w] == fam[w2])] += 1
        fs = {}
        for c, v in sc.items():
            fs[fam[c]] = max(fs.get(fam[c], 0.0), v[0])
        fw_, _, _, _ = winner(fs)
        famw.append(fw_)
        for c, v in sc.items():
            bits_member[c] += v[0]
            cov[c] = max(cov.get(c, 0.0), (v[2] - v[1] + 1) / float(max(clen.get(c, 1), 1)))
        for f_, v in fs.items():
            bits_fam[f_] += v
    flat_firm = len(set(fw)) == 1 and fw[0] is not None
    fam_firm = len(set(famw)) == 1 and famw[0] is not None
    sub_firm, sub_pair = None, None
    if fam_firm:
        F = famw[0]
        sw = []
        for d in cycles:
            sc = {c: v[0] for c, v in d.get(cp, {}).items() if fam[c] == F}
            sw.append(winner(sc)[0])
        sub_firm = len(set(sw)) == 1
        if not sub_firm:
            cnt = collections.Counter(x for x in sw if x)
            top = [c for c, _ in cnt.most_common(2)]
            sub_pair = tuple(sorted(top)) if len(top) == 2 else None
    flat_pair = None
    if not flat_firm:
        cnt = collections.Counter(x for x in fw if x)
        top = [c for c, _ in cnt.most_common(2)]
        flat_pair = tuple(sorted(top)) if len(top) == 2 else None
    res[cp] = (flat_firm, fw[0] if flat_firm else None, fam_firm, famw[0] if fam_firm else None, sub_firm, sub_pair, flat_pair,
               dict(bits_member), dict(bits_fam), cov)

n = len(res)
flat_firm = sum(1 for r in res.values() if r[0])
fam_firm = sum(1 for r in res.values() if r[2])
sub_firm = sum(1 for r in res.values() if r[2] and r[4])
print("=== %s: %d copies with a hit; families: %s" % (label, n, ", ".join(sorted(set(fam.values())))))
print("flat firm (today's vote)                %7d  %5.1f %%" % (flat_firm, 100.0 * flat_firm / n))
print("family firm (layer 1)                   %7d  %5.1f %%" % (fam_firm, 100.0 * fam_firm / n))
print("  of these, subfamily firm (layer 2)    %7d  %5.1f %% of all" % (sub_firm, 100.0 * sub_firm / n))
print("  family firm, subfamily split          %7d" % (fam_firm - sub_firm))
print("flat firm but NOT family firm           %7d   (must be 0 unless a family's best member changes)" % sum(1 for r in res.values() if r[0] and not r[2]))
print("near-tie cycles (top two within 2 %%): same family %d, different families %d" % (near_tie[True], near_tie[False]))


def pair_class(pair, cov):
    a, b = pair
    if is_comp[a] or is_comp[b]:
        cm = a if is_comp[a] else b
        return "composite-unit" if cov.get(cm, 0.0) < 0.6 else "composite-both"
    return "monomers"


for title, idx in (("flat split pairs (today)", 6), ("layer-2 split pairs (inside a firm family)", 5)):
    pc, cls = collections.Counter(), collections.Counter()
    for r in res.values():
        p = r[idx]
        if p:
            pc[p] += 1
            cls[(fam[p[0]] == fam[p[1]], pair_class(p, r[9]))] += 1
    tot = sum(pc.values())
    print("--- %s: %d copies" % (title, tot))
    for (same, c), v in sorted(cls.items(), key=lambda x: -x[1]):
        print("    %-15s %-22s %6d  %5.1f %%" % ("same family" if same else "across families", c, v, 100.0 * v / max(tot, 1)))
    for p, v in pc.most_common(12):
        print("    %6d  %s | %s" % (v, p[0], p[1]))

# threshold: today 0.45 x Nth-best (N = min(10, count)) of the summed score among the member's own firm copies; the same rule per family
def thr(vals):
    v = sorted(vals, reverse=True)
    return 0.45 * v[min(10, len(v)) - 1] if v else 0.0

by_m, by_f = collections.defaultdict(list), collections.defaultdict(list)
for cp, r in res.items():
    if r[0]:
        by_m[r[1]].append(r[7][r[1]])
    if r[2]:
        by_f[r[3]].append(r[8][r[3]])
tm = {m: thr(v) for m, v in by_m.items()}
tf = {f: thr(v) for f, v in by_f.items()}
rej_m = sum(1 for r in res.values() if r[0] and r[7][r[1]] < tm[r[1]])
rej_f = sum(1 for r in res.values() if r[2] and r[8][r[3]] < tf[r[3]])
print("--- threshold 0.45 x 10th best: per member (today) rejects %d of the flat-firm copies; per family it would reject %d of the family-firm"
      % (rej_m, rej_f))
rf = collections.Counter()
for r in res.values():
    if r[2] and r[8][r[3]] < tf[r[3]] and r[4]:
        # which member did the family threshold reject (layer-2 firm member)?
        best = max((c for c in r[7] if fam[c] == r[3]), key=lambda c: r[7][c])
        rf[best] += 1
for m, v in rf.most_common(8):
    print("    rejected by the family threshold although subfamily-firm: %6d  %s (%d bp)" % (v, m, clen.get(m, 0)))
