#!/usr/bin/env python3
"""build_composite.py UNITS.tsv GENOME.fa SUBFAMILY "CHAIN" NAME OUT_DIR [--n 60] [--threads 16]

Consensus of a composite element found by composite_scan.py: the copies of SUBFAMILY whose unit
chain is exactly CHAIN (e.g. "r1[he] ~ r3[he]"), best N by summed bitscore, each cut from the start
of its first unit to the end of its last unit (no flank), aligned with MAFFT L-INS-i; majority
consensus over columns where >= half the copies have a base. How r10_groupB was built by hand
(Tal rsi/REFINEMENT.md §3), made repeatable for r1+r3 and r5head+r6 (§11).
Writes OUT_DIR/NAME.fa (consensus), NAME.copies.fa, NAME.aln.fa (consensus as row 1).
"""
import collections
import os
import subprocess
import sys


def arg(name, default, cast=str):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def main():
    units_tsv, genome, sf, chain, name, out = sys.argv[1:7]
    n = arg("--n", 60, int)
    T = arg("--threads", 16, int)
    os.makedirs(out, exist_ok=True)
    cand = []
    for l in open(units_tsv):
        f = l.rstrip("\n").split("\t")
        if f[0] == "locus" or f[1] != sf or f[2] != chain:
            continue
        c, ws, we, st = f[3], int(f[4]), int(f[5]), f[6]
        us = [u.split(":") for u in f[7].split(";")]
        lo = min(int(u[2].split("-")[0]) for u in us)
        hi = max(int(u[2].split("-")[1]) for u in us)
        bits = sum(float(u[4]) for u in us)
        g0, g1 = (ws + lo, ws + hi) if st == "+" else (we - hi, we - lo)
        cand.append((bits, c, g0, g1, st, f[0]))
    cand.sort(reverse=True)
    seen, pick = collections.defaultdict(list), []
    for b, c, g0, g1, st, lid in cand:           # one copy per genomic element
        if any(min(g1, e1) - max(g0, e0) > 0 for e0, e1 in seen[c]):
            continue
        seen[c].append((g0, g1)); pick.append((c, g0, g1, st, lid))
        if len(pick) == n:
            break
    print("%s: %d copies with chain '%s', using %d" % (sf, len(cand), chain, len(pick)))
    bed = os.path.join(out, name + ".bed")
    with open(bed, "w") as fh:
        for c, g0, g1, st, lid in pick:
            fh.write("%s\t%d\t%d\t%s\t0\t%s\n" % (c, g0, g1, lid, st))
    cfa = os.path.join(out, name + ".copies.fa")
    subprocess.run("bedtools getfasta -fi %s -bed %s -s | sed '/^>/s/@U@/_/g' > %s" % (genome, bed, cfa), shell=True, check=True)
    aln = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--nuc", "--quiet",
                          "--thread", str(T), cfa], capture_output=True, text=True, check=True).stdout
    names, seqs = [], []
    for l in aln.splitlines():
        if l.startswith(">"):
            names.append(l[1:]); seqs.append("")
        else:
            seqs[-1] += l.strip().upper()
    cons = ""
    for j in range(len(seqs[0])):
        col = [s[j] for s in seqs]
        k = collections.Counter(x for x in col if x != "-")
        if k and sum(k.values()) >= len(seqs) / 2.0:
            cons += k.most_common(1)[0][0]
    lens = sorted(g1 - g0 for c, g0, g1, st, lid in pick)
    print("  copy length p10/median/p90: %d/%d/%d; consensus %d bp" % (
        lens[len(lens) // 10], lens[len(lens) // 2], lens[9 * len(lens) // 10], len(cons)))
    open(os.path.join(out, name + ".fa"), "w").write(">%s\n%s\n" % (name, cons))
    with open(os.path.join(out, name + ".aln.fa"), "w") as fh:
        fh.write(">%s\n%s\n" % (name, cons))
        for a, s in zip(names, seqs):
            fh.write(">%s\n%s\n" % (a, s))
    print(cons)


if __name__ == "__main__":
    main()
