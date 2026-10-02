#!/usr/bin/env python3
"""array_flag.py RUN_DIR: how much of each family sits in tandem arrays?

For every family of results/assignment_full.tsv (firmly assigned copies) the copies are checked with the regular-spacing rule of
tools/array_order.py (>= 5 consecutive copies on one contig, gaps <= 10 kb and within a factor 2 of the run's median). Writes
results/array_flag.tsv: family, copies, copies_in_arrays, pct_in_arrays, arrays, largest_array, median_spacing_bp, flag.
flag ARRAY = at least 20 % of the copies are in arrays: the copies of that family are not independent insertions, their flanks
align, and the family cannot be judged as a dispersed SINE from them (the report replaces "Strong SINE" by "Tandem array").
"""
import collections
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import array_order as ao  # noqa: E402

ARRAY_PCT = 20.0
LOC = re.compile(r"(.+):(\d+)-(\d+)\(([+\-,]+)\)$")


def main():
    run = os.path.abspath(sys.argv[1])
    src = os.path.join(run, "results", "assignment_full.tsv")
    if not os.path.exists(src):
        print("array_flag: missing %s, skipped" % src)
        return 0
    fam = collections.defaultdict(list)
    for line in list(open(src))[1:]:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[4] != "assigned":
            continue
        m = LOC.match(f[0])
        if m:
            fam[f[1]].append((m.group(1), int(m.group(2))))
    out = os.path.join(run, "results", "array_flag.tsv")
    rows = []
    for name, loci in fam.items():
        runs = ao.regular_runs(loci)
        by = collections.defaultdict(list)
        for i, rid in runs.items():
            by[rid].append(loci[i])
        gaps = []
        for v in by.values():
            v.sort(key=lambda x: x[1])
            gaps += [v[i + 1][1] - v[i][1] for i in range(len(v) - 1)]
        pct = 100.0 * len(runs) / len(loci)
        rows.append((name, len(loci), len(runs), pct, len(by), max((len(v) for v in by.values()), default=0),
                     int(statistics.median(gaps)) if gaps else 0, "ARRAY" if pct >= ARRAY_PCT else "-"))
    rows.sort(key=lambda r: -r[3])
    with open(out, "w") as o:
        o.write("family\tcopies\tcopies_in_arrays\tpct_in_arrays\tarrays\tlargest_array\tmedian_spacing_bp\tflag\n")
        for r in rows:
            o.write("%s\t%d\t%d\t%.1f\t%d\t%d\t%d\t%s\n" % r)
    print("array_flag: %d families, %d flagged ARRAY -> results/array_flag.tsv" % (len(rows), sum(1 for r in rows if r[7] == "ARRAY")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
