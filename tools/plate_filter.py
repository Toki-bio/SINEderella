#!/usr/bin/env python3
"""plate_filter.py LOCI.tsv LOCI_VERDICT.tsv [KEEP=SINE] [GREY_OUT.tsv] > FILTERED.tsv

Keeps only the plate candidates whose copy has verdict KEEP (default SINE) in loci_verdict.tsv (tools/flank_verdict.py).
With GREY_OUT, copies with verdict grey (flanks masked / untestable; owner 2026-10-18: "grey ones should still be aligned
and shown") are written there with column 7 = "grey": step8a uses them like the old soft top-up - after the SINE loci,
only when a subfamily has fewer than 100 - and marks their rows " [grey]". Copies with verdict not_SINE never pass.

LOCI.tsv: step8a's loci table (subfamily, bitscore, contig, start, end, strand[, ...]). Its coordinates can be wider than
the copy's own (the publish border loop extends element coordinates), so a row matches a verdict copy on the same contig
and strand whose interval lies inside the row's interval (5 bp slack); a row matching no copy is dropped.
Prints "PLATE_FILTER kept=.. grey=.. dropped=.." on stderr.
"""
import bisect
import re
import sys

loci, verdict = sys.argv[1], sys.argv[2]
keep = sys.argv[3] if len(sys.argv) > 3 else "SINE"
grey_out = sys.argv[4] if len(sys.argv) > 4 else None
NAME = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")

cop = {}                                     # (contig, strand) -> sorted [(start, end, verdict)]
with open(verdict) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        m = NAME.match(p[0])
        if m:
            cop.setdefault((m.group(1), m.group(4)), []).append((int(m.group(2)), int(m.group(3)), p[3]))
for v in cop.values():
    v.sort()
starts = {k: [a for a, _, _ in v] for k, v in cop.items()}


def verdict_of(ctg, strand, a, b):
    v = cop.get((ctg, strand))
    if not v:
        return None
    i = bisect.bisect_left(starts[(ctg, strand)], a - 5)
    best = None
    while i < len(v) and v[i][0] <= b + 5:
        if v[i][1] <= b + 5:
            if v[i][2] == keep:
                return keep
            best = best or v[i][2]
        i += 1
    return best


n_keep = n_grey = n_drop = 0
go = open(grey_out, "w") if grey_out else None
with open(loci) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        vd = verdict_of(p[2], p[5], int(p[3]), int(p[4])) if len(p) >= 6 else None
        if vd == keep:
            sys.stdout.write(line)
            n_keep += 1
        elif vd == "grey" and go:
            go.write("\t".join(p[:6] + ["grey"]) + "\n")
            n_grey += 1
        else:
            n_drop += 1
if go:
    go.close()
sys.stderr.write("PLATE_FILTER kept=%d grey=%d dropped=%d (verdict %s)\n" % (n_keep, n_grey, n_drop, keep))
