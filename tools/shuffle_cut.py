#!/usr/bin/env python3
"""Per-consensus score cuts from a shuffled copy of the library, made with a fixed seed (docs/REPRODUCIBILITY.md).

ssearch36 -z 11 gave every consensus its own threshold: it scored shuffled library sequences and fitted its E-values to
them, so a consensus that scores high against random sequence of the same composition (GC-rich, repeated motifs) needed
a higher score to be reported. Its shuffles were seeded by the clock, so the result changed from run to run. Here the same
is done outside ssearch36, the same way in every run:

  shuffle_cut.py shuffle IN.fa OUT.fa [--max 2000] [--seed 1]
      the first --max sequences of IN, each with its residues shuffled (uniform shuffle of the sequence, as ssearch36
      did; own generator, so the output does not depend on the python version). Lowercase is kept as it is.
  shuffle_cut.py cut HITS.m8 N_SHUFFLED N_LIBRARY [--expect 2]
      HITS.m8: ssearch36 -z 3 -m 8 of the consensuses against the shuffled file with -b =N_SHUFFLED -d 0 (the best
      alignment of every shuffled sequence). Per consensus: the best bit score of each shuffled sequence, an extreme-value
      (Gumbel) fit through their 80 % and 95 % points, and the bit score that --expect shuffled sequences are expected to reach in a library of
      N_LIBRARY sequences: cut = mu + beta * ln(N_LIBRARY / expect). That is what "E <= 2" meant under -z 11.
      Output: consensus <TAB> cut <TAB> mu <TAB> beta <TAB> shuffled sequences.
      A consensus with fewer than 30 scored shuffled sequences, or fewer than 30 % of them, gets no line (the caller keeps its fixed cut).
"""
import math
import sys


def lcg(seed):
    """64-bit linear congruential generator (Knuth MMIX constants); yields the high 32 bits"""
    x = (seed * 2862933555777941757 + 3037000493) & 0xFFFFFFFFFFFFFFFF
    while True:
        x = (x * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        yield x >> 32


def shuffle(src, dst, nmax, seed):
    rnd = lcg(seed)
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        name, seq = None, []

        def flush():
            if name is None:
                return
            s = list("".join(seq))
            for i in range(len(s) - 1, 0, -1):          # Fisher-Yates
                j = next(rnd) % (i + 1)
                s[i], s[j] = s[j], s[i]
            fo.write(">shuf%d\n%s\n" % (n, "".join(s)))

        for line in fi:
            if line.startswith(">"):
                flush()
                if n >= nmax:
                    name = None
                    break
                n += 1
                name, seq = line, []
            elif name is not None:
                seq.append(line.strip())
        else:
            flush()
    print(n)


def cut(m8, nshuf, nlib, expect):
    best = {}
    with open(m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or line.startswith("#"):
                continue
            d = best.setdefault(f[0], {})
            b = float(f[11])
            if b > d.get(f[1], -1.0):
                d[f[1]] = b
    for q in sorted(best):
        v = sorted(best[q].values())
        # ssearch36 does not report the lowest-scoring shuffled sequences, and with several threads which of them it leaves
        # out can change (59-74 % were reported in the test). So only the top is read: the 80 % and 95 % points of all n shuffled sequences (the
        # unreported ones lie below the reported ones), and the Gumbel parameters from these two quantiles.
        n = max(nshuf, len(v))
        if len(v) < 30 or len(v) < 0.3 * n:
            continue
        q80 = v[len(v) - 1 - int(0.2 * n)]
        q95 = v[len(v) - 1 - int(0.05 * n)]
        beta = max((q95 - q80) / 1.47026, 1e-6)      # -ln(-ln 0.95) + ln(-ln 0.8)
        mu = q80 - 1.49994 * beta                    # -ln(-ln 0.8)
        c = mu + beta * math.log(max(nlib, 1) / expect) if nlib > expect else mu
        print("%s\t%.1f\t%.2f\t%.3f\t%d" % (q, c, mu, beta, n))


def main():
    a = sys.argv[1:]
    opt = lambda k, d: type(d)(a[a.index(k) + 1]) if k in a else d
    if len(a) >= 3 and a[0] == "shuffle":
        shuffle(a[1], a[2], opt("--max", 2000), opt("--seed", 1))
    elif len(a) >= 4 and a[0] == "cut":
        cut(a[1], int(a[2]), int(a[3]), opt("--expect", 2.0))
    else:
        raise SystemExit(__doc__)


main()
