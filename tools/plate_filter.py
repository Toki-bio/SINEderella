#!/usr/bin/env python3
"""plate_filter.py LOCI.tsv LOCI_VERDICT.tsv [KEEP=SINE] > FILTERED.tsv

Keeps only the plate candidates whose copy has verdict KEEP (default SINE) in loci_verdict.tsv (tools/flank_verdict.py):
the report's alignments show interspersed copies only, never copies with shared, masked or untested flanks (docs/AIM.md).

LOCI.tsv: step8a's loci table (subfamily, bitscore, contig, start, end, strand[, ...]). Its coordinates can be wider than
the copy's own (the publish border loop extends element coordinates), so a row matches a verdict copy on the same contig
and strand whose interval lies inside the row's interval (5 bp slack); a row matching no kept copy is dropped.
Prints "PLATE_FILTER kept=.. dropped=.." on stderr.
"""
import bisect
import re
import sys

loci, verdict = sys.argv[1], sys.argv[2]
keep = sys.argv[3] if len(sys.argv) > 3 else "SINE"
NAME = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")

ok = {}                                      # (contig, strand) -> sorted [(start, end)]
with open(verdict) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or p[3] != keep:
            continue
        m = NAME.match(p[0])
        if m:
            ok.setdefault((m.group(1), m.group(4)), []).append((int(m.group(2)), int(m.group(3))))
for v in ok.values():
    v.sort()
starts = {k: [a for a, _ in v] for k, v in ok.items()}


def kept(ctg, strand, a, b):
    v = ok.get((ctg, strand))
    if not v:
        return False
    i = bisect.bisect_left(starts[(ctg, strand)], a - 5)
    while i < len(v) and v[i][0] <= b + 5:
        if v[i][1] <= b + 5:
            return True
        i += 1
    return False


n_keep = n_drop = 0
with open(loci) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 6 and kept(p[2], p[5], int(p[3]), int(p[4])):
            sys.stdout.write(line)
            n_keep += 1
        else:
            n_drop += 1
sys.stderr.write("PLATE_FILTER kept=%d dropped=%d (verdict %s)\n" % (n_keep, n_drop, keep))
