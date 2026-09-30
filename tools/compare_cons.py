#!/usr/bin/env python3
"""compare_cons.py COPIES_ALN OUT_ALN CONS.fa [CONS.fa ...]

Plate-style comparison: the copies' uppercase element parts are realigned together with several
candidate consensuses (MAFFT L-INS-i, plate settings); each copy's lowercase flanks are packed against
the element (5' flank right-justified, 3' flank left-justified), never scattered by the aligner.
Every consensus row is named with the mean identity of the copies to it over the columns where both
have a base, so the best-matching consensus can be read off the row names.
"""
import subprocess
import sys
import tempfile

src, out, cons_files = sys.argv[1], sys.argv[2], sys.argv[3:]


def read(p):
    recs = []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            recs.append([line[1:].split()[0], ""])
        elif recs:
            recs[-1][1] += line
    return recs


copies = read(src)[1:]                       # row 1 of the input is its own consensus - dropped
parts = []
for name, s in copies:
    u = s.replace("-", "")
    i = 0
    while i < len(u) and u[i].islower():
        i += 1
    j = len(u)
    while j > i and u[j - 1].islower():
        j -= 1
    parts.append((name, u[:i], u[i:j], u[j:]))
cons = []
for f in cons_files:
    for n, s in read(f):
        cons.append((n, s.replace("-", "").upper()))

with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
    for n, s in cons:
        fh.write(">C_%s\n%s\n" % (n, s))
    for k, (n, f5, el, f3) in enumerate(parts):
        fh.write(">E%d\n%s\n" % (k, el))
    tmp = fh.name
aln = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--nuc",
                      "--preservecase", "--quiet", "--thread", "16", tmp], capture_output=True, text=True).stdout
A = {}
for block in aln.split(">")[1:]:
    head, _, body = block.partition("\n")
    A[head.split()[0]] = "".join(body.split())

w5 = max(len(p[1]) for p in parts)
w3 = max(len(p[3]) for p in parts)
rows = []
for n, s in cons:
    a = A["C_" + n]
    idents = []
    for k in range(len(parts)):
        e = A["E%d" % k]
        both = [(x, y) for x, y in zip(a, e) if x != "-" and y != "-"]
        if both:
            idents.append(sum(x.upper() == y.upper() for x, y in both) / len(both))
    mean = 100 * sum(idents) / len(idents)
    rows.append((mean, ">%s mean_copy_identity=%.1f%%" % (n, mean), "-" * w5 + a.upper() + "-" * w3))
rows.sort(key=lambda r: -r[0])
with open(out, "w") as fh:
    for _, h, s in rows:
        fh.write(h + "\n" + s + "\n")
    for k, (n, f5, el, f3) in enumerate(parts):
        fh.write(">%s\n%s%s%s\n" % (n, "-" * (w5 - len(f5)) + f5, A["E%d" % k], f3 + "-" * (w3 - len(f3))))
print(out)
for m, h, _ in rows:
    print("  %s" % h[1:])
