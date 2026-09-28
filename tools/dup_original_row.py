#!/usr/bin/env python3
"""dup_original_row.py ALIGNMENT.aln.fa SUBFAMILY ORIGINAL.fa

Give a freshly built step8a plate its original-consensus row by COPYING row 1, not by realigning.

step8a aligns the copies together with the consensus (row 1), so the consensus sits exactly over the
copies. The publish chain then rebuilds row 1 from the copies. The original used to be put back at the
very end with `mafft --add` (tools/add_seed_row.py) - into a plate whose flanks had been packed and
whose consensus row had been rewritten, and MAFFT then misplaced it: on rsi r1_9seqs the original's 3'
part landed ~100 columns away, in columns no copy occupies, although in the step8a alignment it sat on
the copies at support 0.93.

So, right after step8a's MAFFT: row 1 is renamed <subfamily>_extended (the chain rebuilds it) and an
identical row named <subfamily> is inserted as row 2; every chain step skips it
(SINE-discriminator fix_alignments.is_seed) and moves it along with the rest. Done only when row 1 is
the original exactly (case and gaps ignored); if the publish consensus differs from the one searched
(border-loop widening), nothing is done and add_seed_row.py adds the original with mafft --add as before.
Prints one status line. Idempotent.
"""
import re
import sys


def read_fa(path):
    names, seqs = [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
        elif seqs:
            seqs[-1].append(line.strip())
    return names, ["".join(s) for s in seqs]


def main(argv):
    aln, sf, orig_fa = argv[1], argv[2], argv[3]
    names, seqs = read_fa(aln)
    if len(names) > 1 and names[1] == sf:
        print("%s: already has the original row" % aln)
        return
    on, os_ = read_fa(orig_fa)
    orig = {n.split()[0]: re.sub(r"[^A-Za-z]", "", s).upper() for n, s in zip(on, os_)}.get(sf)
    row1 = re.sub(r"[^A-Za-z]", "", seqs[0]).upper() if seqs else ""
    n1 = names[0].split()[0] if names else ""
    if orig is None or row1 != orig or n1.replace("_R_", "", 1) != sf:
        print("%s: row 1 is not the original %s - left for add_seed_row.py" % (aln, sf))
        return
    names = [sf + "_extended", sf] + names[1:]
    seqs = [seqs[0], seqs[0]] + seqs[1:]
    with open(aln, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            fh.write(">%s\n%s\n" % (n, s))
    print("%s: original copied to row 2" % aln)


if __name__ == "__main__":
    main(sys.argv)
