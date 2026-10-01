#!/usr/bin/env python3
"""Toy genome for the flank-twin stage (fs9): every case planted at known positions.

usage: make_toy_flank.py OUTDIR [SEED]
writes OUTDIR/genome.fa (LINE-like repeat copies in lower case = soft mask), copies.bed (the SINE copies, bed6, strand
= element orientation), truth.tsv (copy TAB class TAB group TAB identity):

  class  meaning                                             expected from fs9
  U      unique insertion into random sequence               no twin
  HOT    insertion at one offset of a LINE-like repeat       no twin (masked repeat; hot spot, not duplication)
  S      copy inside a duplicated 3 kb region (twins)        twin; tier 1 at >= 85 % identity, tier 2 at 70-85 %
  ARR    unit of a tandem array (units 600 bp apart)         array (not twin)
  R      two independent insertions into a duplicated region dup-region only (not twin)
  END    20 bp from a contig end                              untestable
  LC     (TA)n next to the junction, shared by two copies    no twin
"""
import os, random, sys

out = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
os.makedirs(out, exist_ok=True); rnd = random.Random(seed)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "t", "c": "g", "g": "c", "t": "a"}
def rc(s): return "".join(COMP[c] for c in reversed(s))
def rseq(n): return "".join(rnd.choice("ACGT") for _ in range(n))
def mutate(s, sub, indel=0.0, want_map=False):
    r = []; mp = []
    for c in s:
        x = rnd.random(); mp.append(len(r))
        if x < indel / 2: continue                       # deletion
        if x < indel: r.append(c)                        # insertion before c
        r.append(rnd.choice([b for b in "ACGT" if b != c.upper()]) if rnd.random() < sub else c)
    r = "".join(r)
    return (r, mp) if want_map else r

NC, CL = 6, 1000000
contig = [list(rseq(CL)) for _ in range(NC)]
# slots of 4 kb per contig, handed out in random order
slots = [(c, i * 4000) for c in range(NC) for i in range(1, CL // 4000 - 1)]
rnd.shuffle(slots)
def slot(): return slots.pop()
def slot_on(c):
    for i, (cc, p) in enumerate(slots):
        if cc == c: return slots.pop(i)
copies = []      # (contig, start, end, strand, name, class, group, identity)
def put(c, pos, seq):
    contig[c][pos:pos + len(seq)] = list(seq)

SINE = rseq(190) + "A" * 12
def sine_copy(): return mutate(SINE[:-12], 0.05) + "A" * rnd.randint(8, 16)
def insert_sine(c, pos, strand="+"):                      # element + 10 bp target-site duplication; returns (start, end) of the element
    s = sine_copy(); tsd = "".join(contig[c][pos:pos + 10])
    el = s if strand == "+" else rc(s)
    full = tsd + el + tsd
    put(c, pos, full); return pos + 10, pos + 10 + len(el)
n = [0]
def name(cls): n[0] += 1; return "%s%03d" % (cls, n[0])

# 1) LINE-like repeat: 300 copies of 600 bp, 10 % diverged, lower case (the soft mask a real genome carries)
LINE = rseq(600)
for _ in range(300):
    c, p = slot(); put(c, p, mutate(LINE, 0.10).lower())
# 2) hot spot: 40 SINEs inserted at offset 300 of further LINE copies
for _ in range(40):
    c, p = slot(); put(c, p, mutate(LINE, 0.10).lower())
    s, e = insert_sine(c, p + 300); copies.append((c, s, e, "+", name("HOT"), "HOT", "-", ""))
# 3) unique insertions
for _ in range(300):
    c, p = slot(); pos = p + rnd.randint(500, 3000)
    st = rnd.choice("+-"); s, e = insert_sine(c, pos, st); copies.append((c, s, e, st, name("U"), "U", "-", ""))
# 4) segmental duplications carrying a SINE (3 kb, SINE at 1500): same contig / cross contig / inverted, identities 98 .. 70 %
specs = [("same", "+", 0.02), ("cross", "+", 0.05), ("cross", "-", 0.05), ("cross", "+", 0.10), ("same", "+", 0.15),
         ("cross", "-", 0.15), ("cross", "+", 0.22), ("cross", "-", 0.27)]
for gi, (where, ori, div) in enumerate(specs, 1):
    c1, p1 = slot(); c2, p2 = slot()
    if where == "same": c2, p2 = slot_on(c1)
    s1, e1 = insert_sine(c1, p1 + 1500)
    seg = "".join(contig[c1][p1:p1 + 3000]); el = "".join(contig[c1][s1:e1])
    dup, mp = mutate(seg, div, 0.01, True)
    ds, de = mp[s1 - p1], mp[e1 - p1 - 1] + 1               # element span in the duplicate (before orientation)
    if ori == "-": dup = rc(dup); ds, de = len(dup) - de, len(dup) - ds
    put(c2, p2, dup); s2, e2 = p2 + ds, p2 + de
    ident = "%.0f" % (100 * (1 - div))
    copies.append((c1, s1, e1, "+", name("S"), "S", "S%d" % gi, ident)); copies.append((c2, s2, e2, "+" if ori == "+" else "-", name("S"), "S", "S%d" % gi, ident))
# 5) tandem array: 10 units of 600 bp (98 % identical), each carries a SINE
c, p = slot(); unit = rseq(300) + SINE[:-12] + "A" * 10 + rseq(90)
for u in range(10):
    uu = mutate(unit, 0.02); put(c, p + u * 400, uu[:400])
    k = uu.find(SINE[:-12][:30]); s = p + u * 400 + k; copies.append((c, s, s + 190, "+", name("ARR"), "ARR", "A1", "98"))
# 6) class R: a 3 kb region duplicated, a SINE inserted independently into each at different offsets
c1, p1 = slot(); c2, p2 = slot(); reg = "".join(contig[c1][p1:p1 + 3000]); put(c2, p2, mutate(reg, 0.03))
s, e = insert_sine(c1, p1 + 800); copies.append((c1, s, e, "+", name("R"), "R", "R1", ""))
s, e = insert_sine(c2, p2 + 1900); copies.append((c2, s, e, "+", name("R"), "R", "R1", ""))
# 7) contig ends
for c in range(4):
    s, e = insert_sine(c, 20); copies.append((c, s, e, "+", name("END"), "END", "-", ""))
# 8) low complexity: (TA)n in the proximal flank of ten copies (two of them share a longer one)
for i in range(10):
    c, p = slot(); put(c, p + 400, "TA" * 40); s, e = insert_sine(c, p + 480)
    copies.append((c, s, e, "+", name("LC"), "LC", "-", ""))

with open(os.path.join(out, "genome.fa"), "w") as g:
    for c in range(NC):
        g.write(">ctg%d\n" % (c + 1)); s = "".join(contig[c])
        for i in range(0, len(s), 80): g.write(s[i:i + 80] + "\n")
with open(os.path.join(out, "copies.bed"), "w") as b, open(os.path.join(out, "truth.tsv"), "w") as t:
    t.write("copy\tclass\tgroup\tidentity\n")
    for c, s, e, st, nm, cl, grp, ident in copies:
        b.write("ctg%d\t%d\t%d\t%s\t0\t%s\n" % (c + 1, s, e, nm, st)); t.write("%s\t%s\t%s\t%s\n" % (nm, cl, grp, ident))
print("toy written: %d copies" % len(copies), {k: sum(1 for x in copies if x[5] == k) for k in sorted({x[5] for x in copies})})
