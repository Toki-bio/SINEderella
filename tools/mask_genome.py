#!/usr/bin/env python3
"""Write a copy of a genome with BED intervals turned to N, for a masked search.

Usage: mask_genome.py GENOME.fa MASK.bed [MASK.bed ...] -o genome.masked.fa [--stats mask.stats.tsv]

The copy keeps every header, every sequence length and the case of every base outside the mask,
so a coordinate reported on the masked copy is the same coordinate on the original. The search
runs on the copy; extraction and reporting stay on the original. See docs/MASKING.md for why this
is N rather than lowercase: ssearch36 folds a DNA library to upper case, and a rebuilt -S masks
whole library sequences, not hits.

BED chromosome names are matched to the FASTA as given, or after the `_` -> `@U@` rewrite that
SINEderella's sanitize_fasta applies to genome headers, so a BED made against the original
assembly still matches genome.clean.fa. The write is refused (and the output removed) if any
header or length differs, if a BED name matches no sequence, if an interval runs past the end of
its sequence, or if the number of bases turned to N differs from the merged BED length.
"""
import argparse
import os
import sys


def read_bed(paths, names):
    """Merged intervals per sequence; exits on an unknown name or an interval past the end."""
    by = {}
    for path in paths:
        with open(path) as fh:
            for ln, line in enumerate(fh, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    sys.exit("%s:%d: fewer than 3 columns" % (path, ln))
                chrom, s, e = f[0], int(f[1]), int(f[2])
                if chrom not in names:
                    alt = chrom.replace("_", "@U@")
                    if alt not in names:
                        sys.exit("%s:%d: sequence %r is not in the genome" % (path, ln, chrom))
                    chrom = alt
                if s < 0 or e > names[chrom] or s >= e:
                    sys.exit("%s:%d: interval %s:%d-%d outside 0..%d" % (path, ln, chrom, s, e, names[chrom]))
                by.setdefault(chrom, []).append((s, e))
    merged = {}
    for chrom, iv in by.items():
        iv.sort()
        out = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = out
    return merged


def fasta_index(path):
    """Sequence name (first word) -> length, in file order, reading the FASTA once."""
    names, order, cur, n = {}, [], None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    names[cur] = n
                cur, n = line[1:].split()[0], 0
                order.append(cur)
            else:
                n += len(line.rstrip("\n\r"))
    if cur is not None:
        names[cur] = n
    return names, order


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("beds", nargs="+")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--stats", default=None)
    a = ap.parse_args()

    names, order = fasta_index(a.genome)
    mask = read_bed(a.beds, names)
    want = sum(e - s for iv in mask.values() for s, e in iv)

    tmp = a.out + ".tmp"
    masked = 0
    lengths = {}
    with open(a.genome) as fin, open(tmp, "w") as fout:
        cur, pos, iv, k = None, 0, [], 0
        for line in fin:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = pos
                cur, pos = line[1:].split()[0], 0
                iv, k = mask.get(cur, []), 0
                fout.write(line)
                continue
            seq = line.rstrip("\n\r")
            L = len(seq)
            if iv and k < len(iv):
                chars = None
                end = pos + L
                j = k
                while j < len(iv) and iv[j][0] < end:
                    s, e = max(iv[j][0], pos), min(iv[j][1], end)
                    if s < e:
                        if chars is None:
                            chars = list(seq)
                        chars[s - pos:e - pos] = "N" * (e - s)
                        masked += e - s
                    if iv[j][1] <= end:
                        j += 1
                    else:
                        break
                k = j
                if chars is not None:
                    seq = "".join(chars)
            fout.write(seq + "\n")
            pos += L
        if cur is not None:
            lengths[cur] = pos

    problems = []
    if list(lengths) != order:
        problems.append("sequence order or headers differ")
    for n in order:
        if lengths.get(n) != names[n]:
            problems.append("length differs for %s" % n)
    if masked != want:
        problems.append("masked %d bases, BED covers %d" % (masked, want))
    if problems:
        os.remove(tmp)
        sys.exit("mask_genome: refusing to write %s: %s" % (a.out, "; ".join(problems[:5])))
    os.replace(tmp, a.out)

    total = sum(names.values())
    msg = "masked\t%d\tof\t%d\tbp\t(%.3f%%)\tin\t%d\tintervals\n" % (
        masked, total, 100.0 * masked / max(total, 1), sum(len(v) for v in mask.values()))
    sys.stderr.write("mask_genome: " + msg)
    if a.stats:
        with open(a.stats, "w") as fh:
            fh.write("masked_bp\tgenome_bp\tfraction\tintervals\tbed_files\n")
            fh.write("%d\t%d\t%.6f\t%d\t%s\n" % (masked, total, masked / max(total, 1),
                                                sum(len(v) for v in mask.values()), ",".join(a.beds)))


if __name__ == "__main__":
    main()
