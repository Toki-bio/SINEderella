#!/usr/bin/env python3
"""Overlay border-loop rebuilt consensuses onto consensuses.clean.fa for publish."""
import glob
import os
import sys


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cur is not None:
                    seqs.append("".join(buf))
                cur = line[1:].split()[0]
                names.append(cur)
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return dict(zip(names, seqs))


def main():
    run_root = sys.argv[1]
    base = os.path.join(run_root, "consensuses.clean.fa")
    out = os.path.join(run_root, "consensuses.publish.fa")
    cons = read_fa(base)
    n = 0
    for p in glob.glob(os.path.join(run_root, "rebuilt_consensus", "*.fa")):
        sf = os.path.basename(p).replace(".fa", "")
        one = read_fa(p)
        if sf in one:
            cons[sf] = one[sf]
            n += 1
    with open(out, "w") as fh:
        for sf in sorted(cons):
            fh.write(">%s\n%s\n" % (sf, cons[sf]))
    print("wrote %s (%d rebuilt overlays)" % (out, n))


if __name__ == "__main__":
    main()
