"""Toy SINEderella run dir for step8a / align_for_publish: seconds, and every new branch exercised.

TOYS  6 dispersed firm copies (ctgA) + a 5-unit tandem array on ctgB (5 kb spacing, each unit =
      element + the same 300 bp downstream, 2 % mutated) + 4 soft copies  -> array marking/ordering,
      soft top-up, continuation 'unresolved' -> re-extraction loop (rand100 and/or top100)
TOYB  8 dispersed copies with a shared 40 bp tail after the element     -> continuation 'ends'
TOYC  3 copies, soft only                                               -> soft-only plates
TOYA  8 array units sharing 300 bp + 4 independent copies -> NO re-extraction (arrays are skipped)
TOYL  8 copies sharing 200 bp downstream (past the 70 bp flank)       -> continuation re-extraction loop
Dispersed copies sit 60 kb apart (a 12 kb spacing made array_order call everything an array).
"""
import os, random, sys

random.seed(7)
root = sys.argv[1]
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))


def mut(s, p):
    return "".join(random.choice(B.replace(c, "")) if random.random() < p else c for c in s)


cons = {"TOYS": rnd(140) + "A" * 12, "TOYB": rnd(150) + "A" * 10, "TOYC": rnd(130) + "A" * 10,
        "TOYL": rnd(145) + "A" * 10, "TOYA": rnd(138) + "A" * 10}
downA = rnd(300)
longL = rnd(200)
tailB = rnd(40)
arr_down = rnd(300)
ctg = {"ctgA": list(rnd(3000000)), "ctgB": list(rnd(100000)), "ctgC": list(rnd(100000))}
firm, soft = [], []   # (ctg, start1, end1, strand, sf, score)


def put(c, pos, seq):
    ctg[c][pos:pos + len(seq)] = list(seq)


pos = 5000
for i in range(6):                                     # TOYS dispersed, firm
    e = mut(cons["TOYS"], 0.08); put("ctgA", pos, e); firm.append(("ctgA", pos + 1, pos + len(e), "+", "TOYS", 250 - i)); pos += 60000
for i in range(5):                                     # TOYS tandem array, firm, top scores
    p = 20000 + i * 5000
    e = mut(cons["TOYS"], 0.02); d = mut(arr_down, 0.02)
    put("ctgB", p, e + d); firm.append(("ctgB", p + 1, p + len(e), "+", "TOYS", 300 - i))
for i in range(4):                                     # TOYS soft
    e = mut(cons["TOYS"], 0.12); put("ctgA", pos, e); soft.append(("ctgA", pos + 1, pos + len(e), "+", "TOYS", 150 - i)); pos += 60000
for i in range(8):                                     # TOYB with shared tail
    e = mut(cons["TOYB"], 0.08); t = mut(tailB, 0.05); put("ctgA", pos, e + t)
    firm.append(("ctgA", pos + 1, pos + len(e), "+", "TOYB", 240 - i)); pos += 60000
for i in range(3):                                     # TOYC soft only
    e = mut(cons["TOYC"], 0.08); put("ctgA", pos, e); soft.append(("ctgA", pos + 1, pos + len(e), "+", "TOYC", 140 - i)); pos += 60000

for i in range(8):                                     # TOYL: shared 200 bp downstream
    e = mut(cons["TOYL"], 0.06); d = mut(longL, 0.03); put("ctgA", pos, e + d)
    firm.append(("ctgA", pos + 1, pos + len(e), "+", "TOYL", 230 - i)); pos += 60000
for i in range(8):                                     # TOYA: array-dominated (8 units share 300 bp)
    p = 10000 + i * 6000
    e = mut(cons["TOYA"], 0.02); d = mut(downA, 0.02); put("ctgC", p, e + d)
    firm.append(("ctgC", p + 1, p + len(e), "+", "TOYA", 320 - i))
for i in range(4):                                     # TOYA: 4 independent copies
    e = mut(cons["TOYA"], 0.08); put("ctgA", pos, e); firm.append(("ctgA", pos + 1, pos + len(e), "+", "TOYA", 200 - i)); pos += 60000
assert pos < 3000000 - 1000
s2 = os.path.join(root, "step2", "step2_output")
os.makedirs(s2, exist_ok=True)
with open(os.path.join(root, "genome.clean.fa"), "w") as fh:
    for c, s in ctg.items():
        fh.write(">%s\n%s\n" % (c, "".join(s)))
for fn in ("consensuses.clean.fa", "consensuses.clean.fa.pre_publish.bak"):
    with open(os.path.join(root, fn), "w") as fh:
        for k, v in cons.items():
            fh.write(">%s\n%s\n" % (k, v))
G = {c: "".join(s) for c, s in ctg.items()}
with open(os.path.join(s2, "assigned.fasta"), "w") as fh:
    for c, a, b, st, sf, sc in firm:
        fh.write(">%s:%d-%d(%s)|%s|%d\n%s\n" % (c, a, b, st, sf, sc, G[c][a - 1:b]))
with open(os.path.join(s2, "unassigned.tsv"), "w") as fh:
    fh.write("SeqID\tSoft_Subfamily\tSearch_Score\tAll_Queries\tReason\n")
    for c, a, b, st, sf, sc in soft:
        fh.write("%s:%d-%d(%s)\t%s\t%d\t%s\tno_unanimous\n" % (c, a, b, st, sf, sc, sf))
print("toy run dir:", root, "firm", len(firm), "soft", len(soft))
