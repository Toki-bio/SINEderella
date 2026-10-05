#!/usr/bin/env python3
"""array_flag.py RUN_DIR: how much of each family sits in tandem arrays?

For every family of results/assignment_full.tsv (firmly assigned copies) the copies are checked with the regular-spacing rule of
tools/array_order.py (>= 5 consecutive copies on one contig, gaps <= 6 kb and within a factor 5 of the run's median). Writes
results/array_flag.tsv: family, copies, copies_in_arrays, pct_in_arrays, null_pct, excess_pct, arrays, largest_array,
median_spacing_bp, flag.
The share in regular runs is compared with the chance null of tools/satellite_screen.null_b (the same copies spread at random over
the contigs): a hit-dense family forms regular runs by chance (tbr VES, 620 000 copies, one per 3 kb: 56 % of its copies in runs,
almost all of them chance; found on the release check of 2026-10-05, where VES was flagged "Tandem array" in the report). flag ARRAY =
the EXCESS over the null is at least 20 points (rsi MEG-RS: 92.8 % against a null near 0). An ARRAY family's copies are not independent
insertions, their flanks align, and the family cannot be judged as a dispersed SINE from them (the report replaces "Strong SINE" by
"Tandem array"). The satellite stage inside step 1 is the finer instrument (it checks the unit sequence of every run).
"""
import collections
import os
import re
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import array_order as ao  # noqa: E402
import satellite_screen as ss  # noqa: E402

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
            fam[f[1]].append((m.group(1), int(m.group(2)), int(m.group(3))))
    out = os.path.join(run, "results", "array_flag.tsv")
    rows = []
    for name, loci in fam.items():
        runs = ao.regular_runs([(c, s) for c, s, e in loci])
        by = collections.defaultdict(list)
        for i, rid in runs.items():
            by[rid].append(loci[i])
        gaps = []
        for v in by.values():
            v.sort(key=lambda x: x[1])
            gaps += [v[i + 1][1] - v[i][1] for i in range(len(v) - 1)]
        pct = 100.0 * len(runs) / len(loci)
        null = ss.null_b([(c, s, e, "+") for c, s, e in loci])[0] if runs else 0.0
        excess = pct - null
        rows.append((name, len(loci), len(runs), pct, null, excess, len(by), max((len(v) for v in by.values()), default=0),
                     int(statistics.median(gaps)) if gaps else 0, "ARRAY" if excess >= ARRAY_PCT else "-"))
    rows.sort(key=lambda r: -r[5])
    with open(out, "w") as o:
        o.write("family\tcopies\tcopies_in_arrays\tpct_in_arrays\tnull_pct\texcess_pct\tarrays\tlargest_array\tmedian_spacing_bp\tflag\n")
        for r in rows:
            o.write("%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%d\t%d\t%d\t%s\n" % r)
    print("array_flag: %d families, %d flagged ARRAY -> results/array_flag.tsv" % (len(rows), sum(1 for r in rows if r[7] == "ARRAY")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
