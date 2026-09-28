#!/usr/bin/env python3
"""array_order.py LOCI.tsv [--mark-only] [--limit N] > OUT.tsv

step8a loci rows: subfamily, score, ctg, start, end, strand[, soft]. Loci that sit in a TANDEM
cluster - >= MIN_COPIES on one contig with neighbours <= GAP bp apart - get column 8 = "array"
(and column 7 "-" where it was empty); step8a then marks their plate rows " [array]".

Without --mark-only (top100) the order is also changed: every independent locus keeps its rank,
the best-ranked member of each cluster keeps its place, and the other members move to the end, in
their original order. step8a takes the first 100, so the top100 plate SELECTS as many independent
insertions as the family has before any second copy of the same repeated unit. (It does not set
the rows' display order: MAFFT --reorder puts plate rows in guide-tree order.)

Clusters are looked for only among the first --limit rows (default LIMIT; step8a passes 100 for
rand100): the candidates that can reach the plate. Over ALL loci of an abundant family the rule is
meaningless - 100,000 Rhin-1 copies in a 2 Gb genome sit ~20 kb apart on average, so nearly every
copy has neighbours within GAP and was marked "array" (found on a toy run, 2026-09-28, before the
bat republish). Among a few hundred top-ranked copies spread over the genome, chance neighbours
within GAP are rare; a real array still shows, because its near-identical units rank together.

Why (bat corpus, 2026-09-28): the MEG-RS top100 plates of vmu, tbr, fho, tni and cse held 99, 98,
95, 87 and 82 of 100 copies from tandem clusters (spacing ~1-5 kb, near-identical flanks). Being
near-identical they top the bitscore ranking, and the plate then showed one repeated unit ~100
times: its shared flanks read as a "continuation" of the element and hid what the dispersed copies
look like.
"""
import sys

GAP = 50000     # hla MEG-RL: an array with a ~27 kb period
MIN_COPIES = 3
LIMIT = 300


def main(argv):
    path = argv[1]
    mark_only = "--mark-only" in argv
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else LIMIT
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    for r in rows:
        while len(r) < 8:
            r.append("-")
    cand = range(min(limit, len(rows)))
    order = sorted(cand, key=lambda i: (rows[i][2], int(rows[i][3])))
    cluster = {}
    k = 0
    while k < len(order):
        j = k
        while (j + 1 < len(order) and rows[order[j + 1]][2] == rows[order[j]][2]
               and int(rows[order[j + 1]][3]) - int(rows[order[j]][3]) <= GAP):
            j += 1
        if j - k + 1 >= MIN_COPIES:
            for m in order[k:j + 1]:
                cluster[m] = k
                rows[m][7] = "array"
        k = j + 1
    if mark_only:
        out = range(len(rows))
    else:
        seen, head, tail = set(), [], []
        for i in range(len(rows)):
            c = cluster.get(i)
            if c is None or c not in seen:
                head.append(i)
                if c is not None:
                    seen.add(c)
            else:
                tail.append(i)
        out = head + tail
    for i in out:
        sys.stdout.write("\t".join(rows[i]) + "\n")


if __name__ == "__main__":
    main(sys.argv)
