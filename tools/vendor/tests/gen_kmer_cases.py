#!/usr/bin/env python3
"""Generate FASTA test cases for kmer_order (test code only; needs Python >= 3.6).

usage: gen_kmer_cases.py OUTDIR [SCALE]
Writes OUTDIR/<case>.fa. SCALE (default 1) multiplies the sizes of the big cases.
"""
import os
import random
import sys

COMP = str.maketrans("ACGTUacgtuRYKMBVDHrykmbvdh", "TGCAAtgcaaYRMKVBHDyrmkvbhd")


def rc(s):
    return s.translate(COMP)[::-1]


def mut(rng, s, r, indel=0.0):
    out = []
    for c in s:
        x = rng.random()
        if x < r:
            out.append(rng.choice("ACGT"))
        elif x < r + indel / 2:
            continue
        elif x < r + indel:
            out.append(c + rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write(path, recs, width=0, crlf=False):
    nl = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as fh:
        for name, seq in recs:
            fh.write(">" + name + nl)
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + nl)
            else:
                fh.write(seq + nl)


def families(rng, nfam, per, length, rate, indel=0.0, trunc=0.0):
    base = rand_seq(rng, length)
    fams = [base]
    for _ in range(nfam - 1):
        p = rng.randrange(len(fams))
        fams.append(mut(rng, fams[p], 0.04, 0.01))
    recs = []
    for f, cons in enumerate(fams):
        for j in range(per):
            s = mut(rng, cons, rate, indel)
            if trunc and rng.random() < trunc:
                s = s[rng.randrange(0, len(s) // 2):]
            recs.append(("f%d_c%d" % (f, j), s))
    rng.shuffle(recs)
    return recs


def main():
    out = sys.argv[1]
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    os.makedirs(out, exist_ok=True)
    rng = random.Random(20261008)

    write(out + "/toy400.fa", families(rng, 2, 200, 300, 0.05))
    write(out + "/fam8_indel_trunc.fa", families(rng, 8, 60, 280, 0.08, 0.03, 0.3), width=60)

    # exact duplicates and near-duplicates: ties everywhere
    recs = families(rng, 3, 30, 200, 0.0)
    recs += [("dup%d" % i, recs[i % 5][1]) for i in range(40)]
    write(out + "/duplicates.fa", recs)
    write(out + "/all_identical.fa", [("s%d" % i, "ACGTACGTAGCTAGCTAGGATCGATTT" * 3) for i in range(25)])

    # mixed case, N, IUPAC, U, gaps, spaces
    recs = []
    for name, s in families(rng, 3, 25, 150, 0.06):
        s = list(s)
        for _ in range(rng.randrange(0, 6)):
            i = rng.randrange(len(s))
            s[i] = rng.choice("NnRYKMSWryk-.u")
        s = "".join(s)
        if rng.random() < 0.3:
            s = s.lower()
        if rng.random() < 0.1:
            s = s.replace("T", "U")
        recs.append((name, s))
    write(out + "/odd_letters.fa", recs, width=70)

    # reverse-complemented members (orientation)
    recs = []
    for name, s in families(rng, 4, 30, 220, 0.05, 0.01):
        recs.append((name, rc(s) if rng.random() < 0.5 else s))
    write(out + "/revcomp.fa", recs, width=80)
    recs = []
    for name, s in families(rng, 3, 25, 150, 0.06):
        s = s.lower() if rng.random() < 0.3 else s
        recs.append((name, rc(s) if rng.random() < 0.5 else s))
    write(out + "/revcomp_mixed_case.fa", recs)

    # tiny and degenerate inputs
    write(out + "/n1.fa", [("only", rand_seq(rng, 100))])
    write(out + "/n2.fa", [("a", rand_seq(rng, 100)), ("b", rand_seq(rng, 100))])
    write(out + "/n3.fa", [("a", rand_seq(rng, 100)), ("b", rand_seq(rng, 90)), ("c", rand_seq(rng, 110))])
    recs = [("s%d" % i, rand_seq(rng, rng.choice([0, 1, 2, 5, 8, 30, 120]))) for i in range(60)]
    write(out + "/short_and_empty.fa", recs)
    write(out + "/crlf.fa", families(rng, 2, 40, 200, 0.05), width=60, crlf=True)
    with open(out + "/crlf.fa", "a", newline="") as fh:
        pass
    write(out + "/names_spaces.fa", [("  id %d  with spaces " % i, s) for i, (_, s) in enumerate(families(rng, 2, 20, 120, 0.05))])

    # random unrelated sequences (distances near 1)
    write(out + "/random_unrelated.fa", [("r%d" % i, rand_seq(rng, rng.randrange(150, 400))) for i in range(300)])

    # larger: sizes scaled
    n = int(2000 * scale)
    write(out + "/big2k.fa", families(rng, 10, n // 10, 300, 0.07, 0.02, 0.2))
    n = int(8000 * scale)
    write(out + "/big8k.fa", families(rng, 20, n // 20, 300, 0.07, 0.02, 0.2))


if __name__ == "__main__":
    main()
