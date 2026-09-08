#!/usr/bin/env python3
"""Merge border-loop adjusted element coordinates into assigned.fasta headers."""
import glob
import os
import sys


def read_bed_map(border_root):
    m = {}
    for bed in glob.glob(os.path.join(border_root, "*", "loci.adjusted.bed")):
        with open(bed) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                ctg, s, e, name, _sc, strand = parts[:6]
                m[name.split()[0]] = (ctg, int(s) + 1, int(e), strand)
    return m


def rewrite_header(hdr, bed_map):
    key = hdr.split()[0]
    if key not in bed_map:
        return hdr
    ctg, start, end, strand = bed_map[key]
    start = max(1, start)
    end = max(start, end)
    rest = key.split("|", 1)
    if len(rest) != 2:
        return hdr
    sf_bits = rest[1]
    return "%s:%d-%d(%s)|%s" % (ctg, start, end, strand, sf_bits)


def main():
    run_root = sys.argv[1]
    assigned = os.path.join(run_root, "step2/step2_output/assigned.fasta")
    out = os.path.join(run_root, "step2/step2_output/assigned.publish.fasta")
    bed_map = read_bed_map(os.path.join(run_root, "border_loop"))
    n = 0
    with open(assigned) as fin, open(out, "w") as fout:
        name, buf = None, []
        for line in fin:
            if line.startswith(">"):
                if name is not None:
                    fout.write(">%s\n%s\n" % (name, "".join(buf)))
                raw = line[1:].strip()
                new = rewrite_header(raw, bed_map)
                if new != raw:
                    n += 1
                name = new
                buf = []
            else:
                buf.append(line.rstrip("\n\r"))
        if name is not None:
            fout.write(">%s\n%s\n" % (name, "".join(buf)))
    print("wrote %s (%d headers updated from border loop)" % (out, n))


if __name__ == "__main__":
    main()
