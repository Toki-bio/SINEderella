#!/usr/bin/env python3
"""BED of the loci an earlier run's consensuses hit, for a family-level mask.

Usage: mask_from_run.py RUN_DIR [NAME ...] -o mask.bed

Reads RUN_DIR/genome.clean_step1/searches/all_hits.labeled.bed (column 4 = the consensus that found
the hit; "a,b" is accepted when several did). With NAMEs, a locus is kept when any of its labels is
one of them; pass every subfamily consensus of a family to mask the family. Without NAMEs every
hit is kept. Hits, not assignments: a copy is masked even if it later failed the vote. Coordinates
are those of genome.clean.fa, so the mask applies to a run on the same sanitised genome.
"""
import argparse
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("names", nargs="*")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    cands = [os.path.join(a.run_dir, "genome.clean_step1", "searches", "all_hits.labeled.bed"),
             os.path.join(a.run_dir, "genome.clean_step1", "all_hits.labeled.bed")]
    src = next((p for p in cands if os.path.exists(p)), None)
    if src is None:
        sys.exit("mask_from_run: missing %s" % cands[0])
    want = set(a.names)
    seen, n_in, rows = set(), 0, []
    with open(src) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            n_in += 1
            labels = set(f[3].split(","))
            seen |= labels
            if not want or labels & want:
                rows.append((f[0], int(f[1]), int(f[2])))
    missing = want - seen
    if missing:
        sys.exit("mask_from_run: no hits labelled %s in %s (labels there: %s)"
                 % (", ".join(sorted(missing)), src, ", ".join(sorted(seen))))
    rows.sort()
    with open(a.out, "w") as fh:
        for c, s, e in rows:
            fh.write("%s\t%d\t%d\n" % (c, s, e))
    sys.stderr.write("mask_from_run: %d of %d hit intervals (%s)\n"
                     % (len(rows), n_in, ", ".join(sorted(want)) if want else "all consensuses"))


if __name__ == "__main__":
    main()
