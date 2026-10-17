#!/usr/bin/env python3
"""Toy genome for flankscan stage 9b (unknown partners from one-sided shared flanks).

Usage: make_toy_partners.py OUTDIR [seed=5]
Writes OUTDIR/genome.fa (8 x 1 Mb), OUTDIR/bank.fa (units U and V), OUTDIR/partners.fa (partner Y, as a partner
library), OUTDIR/truth.tsv (contig start0 end strand class), where the span is the U element of each copy.

Classes (U = 200 bp unit, every copy 6-12 % diverged, random strand, >= 5 kb apart):
  single   300  U alone
  X+U       25  an invented partner X (250 bp, in no bank or library) joined to U's 5' end, X 4-8 % diverged per copy
                -> stage 9b must report a 5' group of ~25 with an UNKNOWN partner of ~250 bp
  Y+U       15  partner Y (200 bp) joined to U's 5' end; Y is in partners.fa -> reported as known partner Y
  hotspot   20  U inserted at position 200 of a copy of an older repeat R (400 bp, 160 copies at 15 %): both flanks
                shared with other hot-spot copies -> two-sided, must NOT be reported
  segdup     8  4 pairs: U with 1 kb on each side duplicated (2 % between the two) -> twins, must NOT be reported
  X+U dup    6  3 X + U loci duplicated with 1 kb on each side: they share the X side AND the other side, so 9b must
                remove them from the X group (removed_duplicates = 6); the owner found such a block in rsi (NC_142512)
  V        600  a second family, alone
  Bk+U/V 15+25  partner Bk (180 bp, unknown) next to copies of BOTH families in proportion to their size -> may be grouped,
                but the summary must call it "all families: background", while X must be "family-specific"
"""
import os
import random
import sys

out = sys.argv[1]
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
B = "ACGT"


def rnd(n, gc=0.42):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=w, k=n))


def mut(s, p):
    return "".join(rng.choice(B.replace(c, "")) if rng.random() < p else c for c in s)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


U = rnd(185, 0.55) + "A" * 15
V = rnd(190, 0.50) + "T" * 4 + "A" * 12          # a second, unrelated family: makes family specificity testable
X = rnd(250, 0.35)
Y = rnd(200, 0.50)
R = rnd(400, 0.45)
BK = rnd(180, 0.40)                              # background: next to copies of BOTH families, in proportion to their size
os.makedirs(out, exist_ok=True)
open(os.path.join(out, "bank.fa"), "w").write(">U\n%s\n>V\n%s\n" % (U, V))
open(os.path.join(out, "partners.fa"), "w").write(">x.Y\n%s\n" % Y)
open(os.path.join(out, "partner_X_truth.fa"), "w").write(">X\n%s\n" % X)

ctg = {"chr%d" % i: list(rnd(1_000_000)) for i in range(1, 9)}
slots = [(c, p) for c in ctg for p in range(5000, 995000, 6000)]
rng.shuffle(slots)
slots = iter(slots)
truth = []
inserts = {c: [] for c in ctg}                 # contig -> [(pos, seq, el_offset, el_len, strand, class)]


def place(seq, el_off, el_len, cls, strand=None):
    c, p = next(slots)
    st = strand or rng.choice("+-")
    s = seq if st == "+" else rc(seq)
    off = el_off if st == "+" else len(seq) - el_off - el_len
    inserts[c].append((p, s, off, el_len, st, cls))


for _ in range(300):
    e = mut(U, rng.uniform(0.06, 0.12)); place(e, 0, len(e), "single")
for _ in range(600):
    e = mut(V, rng.uniform(0.06, 0.12)); place(e, 0, len(e), "V")
for fam, n in ((U, 15), (V, 25)):
    for _ in range(n):
        k = mut(BK, rng.uniform(0.04, 0.08)); e = mut(fam, rng.uniform(0.06, 0.12))
        place(k + e, len(k), len(e), "Bk+U" if fam is U else "Bk+V")
for _ in range(25):
    x = mut(X, rng.uniform(0.04, 0.08)); e = mut(U, rng.uniform(0.06, 0.12)); place(x + e, len(x), len(e), "X+U")
for _ in range(15):
    y = mut(Y, rng.uniform(0.04, 0.08)); e = mut(U, rng.uniform(0.06, 0.12)); place(y + e, len(y), len(e), "Y+U")
for k in range(160):                           # old repeat R; the first 20 carry a U at position 200 (hot spot)
    r = mut(R, 0.15)
    if k < 20:
        e = mut(U, rng.uniform(0.06, 0.12)); place(r[:200] + e + r[200:], 200, len(e), "hotspot", "+")
    else:
        place(r, 0, 0, "R")
for k in range(3):                             # X + U loci duplicated with 1 kb each side (twice): must leave the X group
    left, right = rnd(1000), rnd(1000); e = mut(X, 0.05) + mut(U, rng.uniform(0.06, 0.12))
    for _ in range(2):
        s = mut(left + e + right, 0.01); place(s, 1000 + len(X), len(s) - 2000 - len(X), "X+U dup")
for k in range(4):                             # segmental duplications: one locus with 1 kb each side, twice
    left, right = rnd(1000), rnd(1000); e = mut(U, rng.uniform(0.06, 0.12))
    for _ in range(2):
        s = mut(left + e + right, 0.01); place(s, 1000, len(s) - 2000, "segdup")

with open(os.path.join(out, "genome.fa"), "w") as fg, open(os.path.join(out, "truth.tsv"), "w") as ft:
    ft.write("contig\tstart0\tend\tstrand\tclass\n")
    for c, base in ctg.items():
        parts, prev, shift = [], 0, 0
        for p, s, off, el, st, cls in sorted(inserts[c]):
            parts.append("".join(base[prev:p]))
            parts.append(s)
            if cls != "R":
                ft.write("%s\t%d\t%d\t%s\t%s\n" % (c, p + shift + off, p + shift + off + el, st, cls))
            shift += len(s)
            prev = p
        parts.append("".join(base[prev:]))
        g = "".join(parts)
        fg.write(">%s\n" % c)
        for i in range(0, len(g), 60):
            fg.write(g[i:i + 60] + "\n")
print("toy written to", out)
