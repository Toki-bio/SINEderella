#!/usr/bin/env python3
"""Simulate a SINE family with a known subfamily tree.

usage: sim_family.py OUTPREFIX [options]
Writes OUTPREFIX.fa (copies), OUTPREFIX.labels (copy id <TAB> subfamily), OUTPREFIX.tree.tsv
(subfamily, parent, diagnostic changes, copies, mean private divergence).

Model (deliberately simple, every knob is a parameter so a benchmark can sweep it):
  * root consensus of --length random bases, GC-rich if --gc is set
  * subfamily 0 is the root; subfamily i descends from a parent chosen among 0..i-1
    (--chain P: probability that the parent is i-1, the nested chain of Yb8 -> Yb8a1 -> Yb10 -> Yb11)
    and carries all its ancestors' diagnostic changes plus --diag new ones
    (substitutions at positions at least --spacing apart from every other diagnostic position;
    a fraction --indel-frac of subfamilies use a 3-base indel as one of them)
  * every copy: the subfamily consensus plus private substitutions at a per-copy divergence drawn from a
    gamma distribution whose mean falls from --div-old (subfamily 0) to --div-young (youngest),
    plus rare private indels (--indel-rate per base)
  * --trunc: fraction of copies cut at the 5' end by up to --trunc-max of their length
  * --revcomp: fraction of copies written as reverse complement
  * --conv: fraction of copies that carry the diagnostic changes of one OTHER subfamily on a short tract
    (a gene-conversion style mosaic); their label stays that of the host subfamily
Sizes: --sizes N (every subfamily N copies), or a comma list, or lognorm:MEAN,SIGMA.
"""
import argparse
import math
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--length", type=int, default=250)
    ap.add_argument("--nsub", type=int, default=8)
    ap.add_argument("--sizes", default="200")
    ap.add_argument("--diag", type=int, default=2, help="new diagnostic changes per subfamily")
    ap.add_argument("--indel-frac", type=float, default=0.0)
    ap.add_argument("--spacing", type=int, default=10)
    ap.add_argument("--chain", type=float, default=0.5)
    ap.add_argument("--div-old", type=float, default=0.12)
    ap.add_argument("--div-young", type=float, default=0.03)
    ap.add_argument("--indel-rate", type=float, default=0.002)
    ap.add_argument("--trunc", type=float, default=0.0)
    ap.add_argument("--trunc-max", type=float, default=0.4)
    ap.add_argument("--revcomp", type=float, default=0.0)
    ap.add_argument("--conv", type=float, default=0.0)
    ap.add_argument("--gc", type=float, default=0.5)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    bases = ["G", "C"] * int(round(a.gc * 10)) + ["A", "T"] * int(round((1 - a.gc) * 10))
    root = "".join(rng.choice(bases) for _ in range(a.length))

    # sizes
    if a.sizes.startswith("lognorm:"):
        mu, sg = map(float, a.sizes[8:].split(","))
        sizes = [max(5, int(rng.lognormvariate(math.log(mu), sg))) for _ in range(a.nsub)]
    elif "," in a.sizes:
        sizes = [int(x) for x in a.sizes.split(",")]
        sizes = (sizes + [sizes[-1]] * a.nsub)[:a.nsub]
    else:
        sizes = [int(a.sizes)] * a.nsub

    # subfamily consensuses; a consensus is a list so indels can be applied in place
    cons = [list(root)]
    parent = [-1]
    diag = [[]]                       # human-readable diagnostic changes of each subfamily
    used_pos = []                     # positions (in root coordinates) already diagnostic
    # indels shift coordinates; keep it simple by working on substitutions in root coordinates and
    # applying a single indel as the LAST change, recorded by its root position
    indel_at = [None]
    for i in range(1, a.nsub):
        p = i - 1 if rng.random() < a.chain else rng.randrange(i)
        c = list(cons[p])
        d = []
        use_indel = rng.random() < a.indel_frac and indel_at[p] is None
        n_sub = a.diag - (1 if use_indel else 0)
        tries = 0
        while len(d) < n_sub and tries < 10000:
            tries += 1
            pos = rng.randrange(5, a.length - 5)
            if any(abs(pos - q) < a.spacing for q in used_pos):
                continue
            if c[pos] in "-":
                continue
            old = c[pos]
            new = rng.choice([b for b in "ACGT" if b != old])
            c[pos] = new
            used_pos.append(pos)
            d.append("%d%s>%s" % (pos + 1, old, new))
        if use_indel:
            for _ in range(10000):
                pos = rng.randrange(10, a.length - 12)
                if all(abs(pos - q) >= a.spacing for q in used_pos) and all(x != "" for x in c[pos:pos + 3]):
                    seg = "".join(c[pos:pos + 3])
                    for k in range(3):
                        c[pos + k] = ""            # 3-base deletion, kept as empty strings so coordinates hold
                    used_pos.append(pos)
                    d.append("%ddel%s" % (pos + 1, seg))
                    indel_at.append(pos)
                    break
            else:
                indel_at.append(None)
        else:
            indel_at.append(indel_at[p])
        cons.append(c)
        parent.append(p)
        diag.append(d)
    cons_str = ["".join(x) for x in cons]

    div = [a.div_old - (a.div_old - a.div_young) * (i / max(1, a.nsub - 1)) for i in range(a.nsub)]

    def private(seq, mean_div):
        k = 4.0                                           # gamma shape: spread of per-copy divergence
        rate = rng.gammavariate(k, mean_div / k)
        out = []
        for ch in seq:
            x = rng.random()
            if x < rate:
                out.append(rng.choice([b for b in "ACGT" if b != ch]))
            elif x < rate + a.indel_rate / 2:
                continue
            elif x < rate + a.indel_rate:
                out.append(ch + rng.choice("ACGT"))
            else:
                out.append(ch)
        return "".join(out)

    n_total = 0
    with open(a.out + ".fa", "w") as fa, open(a.out + ".labels", "w") as lb:
        for i in range(a.nsub):
            for j in range(sizes[i]):
                s = cons_str[i]
                if a.conv and a.nsub > 1 and rng.random() < a.conv:
                    other = rng.choice([k for k in range(a.nsub) if k != i])
                    # copy a 20-60 base tract of the other subfamily's consensus onto this copy where lengths agree
                    L0 = min(len(s), len(cons_str[other]))
                    t0 = rng.randrange(0, max(1, L0 - 60))
                    t1 = t0 + rng.randrange(20, 60)
                    s = s[:t0] + cons_str[other][t0:t1] + s[t1:]
                s = private(s, div[i])
                if a.trunc and rng.random() < a.trunc:
                    s = s[int(len(s) * rng.uniform(0.02, a.trunc_max)):]
                if a.revcomp and rng.random() < a.revcomp:
                    s = rc(s)
                cid = "sf%d|%d" % (i, j)
                fa.write(">%s\n%s\n" % (cid, s)); lb.write("%s\tsf%d\n" % (cid, i)); n_total += 1
    with open(a.out + ".tree.tsv", "w") as tr:
        tr.write("subfamily\tparent\tdiagnostic_changes\tcopies\tmean_private_divergence\n")
        for i in range(a.nsub):
            tr.write("sf%d\t%s\t%s\t%d\t%.3f\n" % (i, "-" if parent[i] < 0 else "sf%d" % parent[i], ",".join(diag[i]) or "-", sizes[i], div[i]))
    print("%d copies, %d subfamilies" % (n_total, a.nsub), file=sys.stderr)


if __name__ == "__main__":
    main()
