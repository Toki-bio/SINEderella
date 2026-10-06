#!/usr/bin/env python3
"""array_order.py LOCI.tsv [--mark-only] [--limit N] [--twins copy_status.tsv] > OUT.tsv

--twins: flankscan stage 9 output for the family (results/flank_twins/<family>/copy_status.tsv): twin copies get column 8 =
"twin" (an array mark wins) and are clustered by their twin group, so the top100 takes one copy per group first.

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

Second rule (2026-10-02), over ALL loci of the family: REGULAR spacing. Near-identical array units rank first, so the
first LIMIT rows can consist of array copies only and the rule above then finds no independent copy to put first (rsi
MEG-RS: 89 % of 1 715 copies in 30 arrays of ~2 kb spacing, 77 of the top 100 from two arrays). A tandem run is >= REG_MIN
consecutive copies on one contig, each gap <= REG_GAP (6 kb) and within a factor REG_RATIO (5) of the run's median gap; chance
neighbours in a dispersed family have gaps of very different size and almost never form such a run. The union of both
rules is used. `regular_runs` is also used by tools/array_flag.py.
"""
import os
import re
import sys

GAP = 50000     # hla MEG-RL: an array with a ~27 kb period
MIN_COPIES = 3
LIMIT = 300
REG_GAP = 6000     # regular-spacing rule (all loci)
REG_MIN = 5
REG_RATIO = 5.0
WIDE_GAP = 30000   # arrays of long period (rsi MEG-RS, 7 kb unit, was missed at 6 kb): only runs of >= WIDE_MIN copies count
WIDE_MIN = 10


def regular_runs(loci, reg_gap=None, reg_min=None):
    """loci: list of (contig, start). Returns {index: run id} for the copies in tandem runs of regular spacing."""
    reg_gap = REG_GAP if reg_gap is None else reg_gap
    reg_min = REG_MIN if reg_min is None else reg_min
    order = sorted(range(len(loci)), key=lambda i: (loci[i][0], loci[i][1]))
    runs, out = 0, {}
    k = 0
    while k < len(order):
        j, gaps = k, []
        while j + 1 < len(order) and loci[order[j + 1]][0] == loci[order[j]][0]:
            g = loci[order[j + 1]][1] - loci[order[j]][1]
            if g > reg_gap or g <= 0:
                break
            if gaps:
                med = sorted(gaps)[len(gaps) // 2]
                if g > REG_RATIO * med or g * REG_RATIO < med:
                    break
            gaps.append(g)
            j += 1
        if j - k + 1 >= reg_min:
            for m in order[k:j + 1]:
                out[m] = runs
            runs += 1
        k = j + 1 if j - k + 1 >= reg_min else (j if j > k else k + 1)
    return out


def regular_runs_wide(loci):
    """regular_runs plus long-period runs: >= WIDE_MIN copies with gaps up to WIDE_GAP replace the narrow runs they contain.
    Returns {index: run id}."""
    out = regular_runs(loci)
    nxt = max(out.values(), default=-1) + 1
    by_run = {}                       # narrow run id -> its members; one pass, instead of scanning all of out for every wide run
    for i, r in out.items():          # (that scan was quadratic: 6.7 h estimated for the 1.13 M DIP hits of Sicista, 2026-10-06)
        by_run.setdefault(r, []).append(i)
    wide = {}
    for i, r in regular_runs(loci, WIDE_GAP, WIDE_MIN).items():
        wide.setdefault(r, []).append(i)
    for members in wide.values():
        old = {out[i] for i in members if i in out}
        for r in old:
            for i in by_run.pop(r, ()):
                if out.get(i) == r:
                    del out[i]
        for i in members:
            out[i] = nxt
        nxt += 1
    return out


def read_twins(path):
    """flankscan stage 9 copy_status.tsv (copy status n_partners group): {(contig, start, end): group} for twin copies.
    Copy names are contig:start-end(strand), the names of assignment_full.tsv / assigned.fasta."""
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[1] not in ("twin1", "twin2"):
            continue
        m = re.match(r"^(.+):(\d+)-(\d+)\([+\-,]+\)$", f[0])
        if m:
            out[(m.group(1), m.group(2), m.group(3))] = f[3] if f[3] != "-" else f[0]
    return out


def main(argv):
    path = argv[1]
    mark_only = "--mark-only" in argv
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else LIMIT
    twins = read_twins(argv[argv.index("--twins") + 1]) if "--twins" in argv and os.path.exists(argv[argv.index("--twins") + 1]) else {}
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
    reg = regular_runs_wide([(r[2], int(r[3])) for r in rows])      # all loci, regular spacing (6 kb tier + 30 kb tier)
    for i, rid in reg.items():
        rows[i][7] = "array"
        cluster[i] = ("r", rid)
    # copies whose flanks are shared with other copies (flankscan stage 9, --twins): not independent insertions either.
    # Marked "twin" (an array mark wins) and clustered by their twin group, so the plate takes one copy per group first.
    for i, r in enumerate(rows):
        g = twins.get((r[2], r[3], r[4]))
        if g is None:
            continue
        if r[7] != "array":
            r[7] = "twin"
        if i not in cluster:
            cluster[i] = ("t", g)
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
