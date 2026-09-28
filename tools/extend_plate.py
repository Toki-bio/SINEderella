#!/usr/bin/env python3
"""extend_plate.py PLATE.aln EXTRACTED.fa ADD5 ADD3 OUT.aln [--threads N]

One continuation round of step8a WITHOUT re-aligning the plate: the copies were re-extracted with
ADD5 / ADD3 more bp per side (EXTRACTED.fa, same names as the plate rows, element in genomic/strand
orientation as bedtools getfasta -s gives it). Only the new segments are aligned - among themselves,
MAFFT L-INS-i on ~100 x ADD bp - and the block is appended to the plate on that side. Row 0 (the
consensus) and any row without a new segment (contig end) get gaps there.

MAFFT --adjustdirection may have reversed a copy in the plate (name "_R_..."): its genomic 3' segment
then belongs on the plate's 5' side, reverse-complemented, and vice versa.

Why (his review, 2026-09-28): each continuation round re-extracted and re-aligned the WHOLE plate
(L-INS-i, cost ~N^2 L^2), ~50x one plate; only the added part needs aligning, and one full alignment
at the final size. The decision (continuation.py --need) reads each copy's own bases past the edge and
the column majority, which the appended block provides.
"""
import os, re, subprocess, sys, tempfile

GAPS = "-."
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fa(p):
    n, s, cur = [], [], []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n:
                s.append("".join(cur))
            n.append(l[1:]); cur = []
        elif l:
            cur.append(l)
    if n:
        s.append("".join(cur))
    return n, s


LOC = re.compile(r"^(?:_R_)?(\S+):(\d+)-(\d+)\(([+-])\)")


def locus(name):
    """plate/getfasta row -> (contig, start, end, strand), or None. The coordinates are the EXTRACTED
    region, so they change when the flanks grow: rows are matched by containment, not by name
    (matching by name found nothing, and the loop ran to the cap on the toy)."""
    m = LOC.match(name.split()[0].replace("@U@", "_"))   # plates restore "_", raw extractions keep "@U@"
    return (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) if m else None


def align_block(segs, threads):
    """segs: {row_index: seq}; returns {row_index: aligned seq} (same width)"""
    if not segs:
        return {}, 0
    if len(segs) == 1:
        (i, s), = segs.items()
        return {i: s}, len(s)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        f = os.path.join(d, "seg.fa")
        with open(f, "w") as fh:
            for i, s in segs.items():
                fh.write(">r%d\n%s\n" % (i, s))
        out = subprocess.run(["mafft", "--localpair", "--maxiterate", "2", "--ep", "0.123", "--nuc",
                              "--preservecase", "--quiet", "--thread", str(threads), f],
                             capture_output=True, text=True, check=True).stdout
    n, s = [], []
    for l in out.splitlines():
        if l.startswith(">"):
            n.append(int(l[2:])); s.append("")
        else:
            s[-1] += l.strip()
    return dict(zip(n, s)), (len(s[0]) if s else 0)


def main(argv):
    plate, extracted, add5, add3, out = argv[1], argv[2], int(argv[3]), int(argv[4]), argv[5]
    threads = int(argv[argv.index("--threads") + 1]) if "--threads" in argv else 4
    pn, ps = read_fa(plate)
    en, es = read_fa(extracted)
    by_ctg = {}
    for n, s in zip(en, es):
        loc = locus(n)
        if loc:
            by_ctg.setdefault((loc[0], loc[3]), []).append((loc[1], loc[2], s))
    left, right = {}, {}       # plate side -> {row: new segment in plate orientation}
    for i, name in enumerate(pn):
        if i == 0:
            continue
        loc = locus(name)
        if not loc:
            continue
        cand = [(nb - na, sq) for na, nb, sq in by_ctg.get((loc[0], loc[3]), ()) if na <= loc[1] and nb >= loc[2]]
        if not cand:
            continue
        s = min(cand)[1]           # the smallest interval that contains this row's old one
        g5 = s[:add5] if add5 else ""                      # genomic 5' new bases
        g3 = s[len(s) - add3:] if add3 else ""             # genomic 3' new bases
        if name.startswith("_R_"):                          # row is reverse-complemented in the plate
            g5, g3 = g3.translate(COMP)[::-1], g5.translate(COMP)[::-1]
        if g5:
            left[i] = g5
        if g3:
            right[i] = g3
    L, wl = align_block(left, threads)
    R, wr = align_block(right, threads)
    with open(out, "w") as fh:
        for i, (n, s) in enumerate(zip(pn, ps)):
            fh.write(">%s\n%s%s%s\n" % (n, L.get(i, "-" * wl), s, R.get(i, "-" * wr)))
    print("appended 5' block %d cols (%d rows), 3' block %d cols (%d rows)" % (wl, len(L), wr, len(R)))


if __name__ == "__main__":
    main(sys.argv)
