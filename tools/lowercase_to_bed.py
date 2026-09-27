#!/usr/bin/env python3
"""BED of the lowercase runs of a soft-masked FASTA, so a soft mask can feed mask_genome.py.

Usage: lowercase_to_bed.py GENOME.fa [--min-len 1] > lowercase.bed

Note that an NCBI or RepeatMasker soft mask covers every repeat class, not only SINE families;
use it when that is the intended mask.
"""
import argparse
import re
import sys

LOW = re.compile(r"[a-z]+")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("--min-len", type=int, default=1)
    a = ap.parse_args()
    out = sys.stdout
    name, pos = None, 0
    run = None  # [start, end) of a lowercase run that may continue on the next line

    def flush():
        if run is not None and run[1] - run[0] >= a.min_len:
            out.write("%s\t%d\t%d\n" % (name, run[0], run[1]))

    with open(a.genome) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, pos, run = line[1:].split()[0], 0, None
                continue
            seq = line.rstrip("\n\r")
            for m in LOW.finditer(seq):
                s, e = pos + m.start(), pos + m.end()
                if run is not None and run[1] == s:
                    run[1] = e
                else:
                    flush()
                    run = [s, e]
            pos += len(seq)
    flush()


if __name__ == "__main__":
    main()
