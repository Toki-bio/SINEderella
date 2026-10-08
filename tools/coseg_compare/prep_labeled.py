#!/usr/bin/env python3
"""Labelled copy set from the owner's hand-curated alignments (KIT /data/W/toki/SINE_disc/aln_c).

usage: prep_labeled.py ALN_DIR SPECIES OUTPREFIX
Reads ALN_DIR/POS__SPECIES__<group>[_<N>seqs].aln.fa (a CONSENSUS_ row plus 100 copies per group),
degaps and upper-cases the copies, drops the consensus row, and writes
  OUTPREFIX.fa      copies, ids  <group>|<k>
  OUTPREFIX.labels  id <TAB> group (the owner's call)
"""
import glob
import os
import re
import sys

aln_dir, sp, out = sys.argv[1:4]
files = sorted(glob.glob(os.path.join(aln_dir, "POS__%s__*.aln.fa" % sp)))
if not files:
    sys.exit("no POS__%s__* files in %s" % (sp, aln_dir))
n = 0
with open(out + ".fa", "w") as fa, open(out + ".labels", "w") as lb:
    for f in files:
        group = re.sub(r"(_\d+seqs)?\.aln\.fa$", "", os.path.basename(f).split("__", 2)[2])
        recs, cur = [], None
        for line in open(f):
            if line.startswith(">"):
                cur = [line[1:].strip(), []]
                recs.append(cur)
            elif cur is not None:
                cur[1].append(line.strip())
        k = 0
        for name, parts in recs:
            if name.upper().startswith("CONSENSUS"):
                continue
            k += 1
            seq = "".join(parts).replace("-", "").replace(".", "").upper()
            cid = "%s|%d" % (group, k)
            fa.write(">%s\n%s\n" % (cid, seq)); lb.write("%s\t%s\n" % (cid, group)); n += 1
        print(group, k, file=sys.stderr)
print("%d copies in %d groups" % (n, len(files)), file=sys.stderr)
