#!/usr/bin/env python3
"""composite_pieces.py UNITS.tsv BANK.fa OUTDIR [MIN_LEN=50]

The stretches of the bank's composites as search queries of their own, for step 1 (2026-10-26, docs/PARTIAL_COPIES.md).

sear keeps a hit only when it covers >= 80 % of its query, so a copy of one unit of a composite (rsi: the head that
r1_r3_P1, r1_r3_P2 and C11 share) is found only if it is also 80 % of some other bank entry. Searching every composite
stretch of >= MIN_LEN bp (rows of kind "composite" in UNITS.tsv, from composite_units.py --bank) finds them: rsi without
r1_9seqs, all 21 r1 loci, the 3 the bank missed included.

Writes, per stretch, OUTDIR/pieces.part_<composite>__piece_<start>-<end>.fa whose header's first word is the COMPOSITE's
name, so every place that labels hits by the query header (step 1, the orchestrator's rebuild, --exclude) gives piece hits
the composite's label; step 2 votes with the bank only. OUTDIR/pieces.named.fa holds the same stretches under their file
names (<composite>__piece_<start>-<end>) for the satellite screen, which names queries by file.
Prints the number of pieces."""
import os
import re
import sys

units, bank, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
min_len = int(sys.argv[4]) if len(sys.argv) > 4 else 50

seq, name = {}, None
for line in open(bank):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seq[name] = []
    elif name:
        seq[name].append(line)
seq = {k: "".join(v) for k, v in seq.items()}

os.makedirs(outdir, exist_ok=True)
n = 0
with open(os.path.join(outdir, "pieces.named.fa"), "w") as named:
    for i, line in enumerate(open(units)):
        f = line.rstrip("\r\n").split("\t")
        if i == 0 or len(f) < 5 or f[1] != "composite" or f[0] not in seq:
            continue
        s, e = int(f[2]), int(f[3])
        if e - s + 1 < min_len:
            continue
        piece = seq[f[0]][s - 1:e]
        tag = "%s__piece_%d-%d" % (re.sub(r"[^\w.#+-]", "_", f[0]), s, e)
        with open(os.path.join(outdir, "pieces.part_%s.fa" % tag), "w") as fh:
            fh.write(">%s piece_%d-%d unit=%s\n%s\n" % (f[0], s, e, f[4], piece))
        named.write(">%s\n%s\n" % (tag, piece))
        n += 1
print(n)
