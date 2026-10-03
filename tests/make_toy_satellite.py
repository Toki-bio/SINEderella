#!/usr/bin/env python3
"""make_toy_satellite.py OUTDIR: a synthetic genome with planted cases for the satellite tools (docs/SATELLITES.md, test plan section 5).

OUTDIR/genome.fa, cons.fa (a 250 bp toy SINE), truth.tsv. Planted (contig, start, end, kind):
  A1  a satellite of 30 monomers = the 3' 120 bp of the SINE (positions 131-250), 8 % diverged, with one full SINE copy at its 5' end
  A2  a satellite of 6 monomers (positions 41-120), 5 % diverged
  B1  an array of 20 units of 1 800 bp, each holding one full SINE copy (units 99 % identical): kind B, not seen by TRF with period <= 300
  S   200 dispersed full-length SINE copies, 10 % diverged
  D   a dimer (two full copies head to tail, 30 bp spacer) x 5 and a trimer of full copies x 3 (must stay in the SINE analysis)
  X   a segmental duplication: a 5 kb region with 3 SINE copies, present twice
"""
import os
import random
import sys

R = random.Random(21)


def rnd(n):
    return "".join(R.choice("ACGT") for _ in range(n))


def mut(s, rate):
    out = []
    for c in s:
        if R.random() < rate:
            r = R.random()
            if r < 0.7:
                out.append(R.choice([x for x in "ACGT" if x != c]))
            elif r < 0.85:
                continue
            else:
                out.append(c + R.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def main():
    d = sys.argv[1]
    os.makedirs(d, exist_ok=True)
    sine = rnd(250)
    contigs = {"ctg%02d" % i: [rnd(200000)] for i in range(40)}
    truth = []

    def plant(c, piece, kind):
        seq = "".join(contigs[c])
        pos = R.randint(2000, len(seq) - 2000)
        contigs[c] = [seq[:pos] + piece + seq[pos:]]
        truth.append((c, pos, pos + len(piece), kind))
        return pos

    # A1
    mono = sine[130:250]
    arr = mut(sine, 0.05) + "".join(mut(mono, 0.08) for _ in range(30))
    plant("ctg00", arr, "A1")
    plant("ctg01", "".join(mut(sine[40:120], 0.05) for _ in range(6)), "A2")
    # B1
    unit = rnd(800) + sine + rnd(750)
    plant("ctg02", "".join(mut(unit, 0.01) for _ in range(20)), "B1")
    # dispersed
    for i in range(200):
        plant("ctg%02d" % R.randint(3, 39), mut(sine, 0.10), "S")
    for i in range(5):
        plant("ctg%02d" % R.randint(3, 39), mut(sine, 0.08) + rnd(30) + mut(sine, 0.08), "D2")
    for i in range(3):
        plant("ctg%02d" % R.randint(3, 39), mut(sine, 0.08) + rnd(25) + mut(sine, 0.08) + rnd(25) + mut(sine, 0.08), "D3")
    region = rnd(1500) + mut(sine, 0.06) + rnd(1000) + mut(sine, 0.06) + rnd(900) + mut(sine, 0.06) + rnd(1000)
    plant("ctg05", region, "X")
    plant("ctg30", mut(region, 0.01), "X")
    with open(os.path.join(d, "genome.fa"), "w") as fh:
        for c, parts in contigs.items():
            s = "".join(parts)
            fh.write(">%s\n" % c + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    open(os.path.join(d, "cons.fa"), "w").write(">toySINE\n" + sine + "\n")
    with open(os.path.join(d, "truth.tsv"), "w") as fh:
        fh.write("contig\tstart\tend\tkind\n")
        for t in truth:
            fh.write("%s\t%d\t%d\t%s\n" % t)
    print("toy written to", d, "| planted", len(truth))


if __name__ == "__main__":
    main()
