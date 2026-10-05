#!/usr/bin/env python3
"""satellite_screen.py: indication of SINE-derived satellites / SINE-containing tandem arrays from hit coordinates only.

Design: docs/SATELLITES.md section 4.1. Input is the hits of one consensus (BED: contig, start, end, name/identity, score, strand; the
hits BEFORE the 80 % length rule, so that partial monomers are there). No alignment, no TRF: seconds even for hundreds of thousands of hits.

Two tests per consensus:
  A  monomer runs   >= MIN_MONO consecutive hits on one contig and strand, each gap <= MONO_GAP (100 bp, the rule of Vassetzky et al. 2023),
                    run spanning > MIN_SPAN (500 bp). Partial hits count: they are the monomers of a SINE-derived satellite.
  B  regular spacing  tools/array_order.regular_runs on the hits (>= 5 copies, gaps <= 6 kb, within a factor 5 of the run's median gap),
                    for arrays whose unit is longer than the SINE (rsi MEG-RS: 2 155 bp unit holding a 107 bp MEG piece).
A hit can be in both kinds of run; the A runs are reported first and the B runs only for hits not already in an A run.
Chance runs: in a hit-rich family (gja Squam3A: 216 000 hits, one per 12 kb) regular spacing occurs by chance, so kind B is calibrated
against a null made of the same hits placed at random positions within each contig's extent (NPERM permutations); the flag uses the
EXCESS over that null, not the raw share. Kind A needs no null (gaps <= 100 bp between >= 4 hits do not occur by chance).

Writes PREFIX.runs.tsv (one row per run: kind, contig, start, end, hits, strand, median hit length, median gap) and PREFIX.summary.tsv
(hits, hits in A runs, hits in B runs, share, loci, largest locus, flag SAT_A / SAT_B / -; flag at >= FLAG_PCT % of the hits).
Nothing is removed or changed: this is the indication step; exclusion is a separate, optional step (docs/SATELLITES.md 4.2).

Usage: satellite_screen.py HITS.bed [HITS2.bed ...] --out PREFIX [--min-mono 4] [--mono-gap 100] [--min-span 500] [--flag-pct 20]
"""
import argparse
import collections
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import array_order as ao  # noqa: E402

MIN_MONO = 4
MONO_GAP = 100
MIN_SPAN = 500
FLAG_PCT = 20.0
NPERM = 5


def read_bed(path):
    hits = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3 or line.startswith("#"):
            continue
        try:
            s, e = int(f[1]), int(f[2])
        except ValueError:
            continue
        hits.append((f[0], s, e, f[5] if len(f) > 5 and f[5] in "+-" else "."))
    return hits


def monomer_runs(hits, min_mono=MIN_MONO, mono_gap=MONO_GAP, min_span=MIN_SPAN):
    """hits: list of (contig, start, end, strand). Returns ([run: list of hit indices], {hit index: run id})."""
    order = sorted(range(len(hits)), key=lambda i: (hits[i][0], hits[i][3], hits[i][1]))
    runs, member = [], {}
    k = 0
    while k < len(order):
        j = k
        while j + 1 < len(order):
            a, b = hits[order[j]], hits[order[j + 1]]
            if a[0] != b[0] or a[3] != b[3] or b[1] - a[2] > mono_gap:
                break
            j += 1
        idx = order[k:j + 1]
        if len(idx) >= min_mono and hits[idx[-1]][2] - hits[idx[0]][1] > min_span:
            for i in idx:
                member[i] = len(runs)
            runs.append(idx)
        k = j + 1
    return runs, member


def null_b(hits, nperm=NPERM, seed=1):
    """share of hits in regular runs when the same number of hits is spread at random over the contigs in proportion to the contig extent
    (genome-average density; the contig extent is the span of its hits, so the null is slightly dense = conservative). A real array
    concentrates hits and must not inflate its own null, so the per-contig counts are NOT kept."""
    import random
    rnd = random.Random(seed)
    ext = {}
    for c, s, e, st in hits:
        lo, hi = ext.get(c, (10 ** 12, 0))
        ext[c] = (min(lo, s), max(hi, e))
    names = sorted(ext)
    span = [max(ext[c][1] - ext[c][0], 1) for c in names]
    total = float(sum(span))
    n = len(hits)
    tot, lens = [], collections.Counter()
    for _ in range(nperm):
        loci = []
        for c, sp in zip(names, span):
            k = int(round(n * sp / total))
            lo = ext[c][0]
            for _ in range(k):
                loci.append((c, lo + rnd.randint(0, sp)))
        runs = ao.regular_runs_wide(loci)
        tot.append(100.0 * len(runs) / max(len(loci), 1))
        for l in collections.Counter(runs.values()).values():
            lens[l] += 1
    return sum(tot) / len(tot), {l: v / float(nperm) for l, v in lens.items()}


def long_run_min(obs_lens, null_lens, max_false=0.05):
    """the smallest run length L such that chance explains fewer than max(1, max_false x observed) of the observed runs with >= L
    units; runs at least that long are satellite loci on their own (the same calibration idea as the TSD minimum of flankscan stage 8).
    None when no length qualifies."""
    for L in sorted(set(obs_lens)):
        o = sum(v for l, v in obs_lens.items() if l >= L)
        e = sum(v for l, v in null_lens.items() if l >= L)
        if o and e < max(1.0, max_false * o):
            return L
    return None


def screen(hits, min_mono=MIN_MONO, mono_gap=MONO_GAP, min_span=MIN_SPAN, nperm=NPERM):
    """returns rows of runs and the summary dict for one consensus"""
    a_runs, a_member = monomer_runs(hits, min_mono, mono_gap, min_span)
    rest = [i for i in range(len(hits)) if i not in a_member]
    reg = ao.regular_runs_wide([(hits[i][0], hits[i][1]) for i in rest])
    b_groups = collections.defaultdict(list)
    for pos, rid in reg.items():
        b_groups[rid].append(rest[pos])
    rows = []

    def row(kind, idx):
        idx = sorted(idx, key=lambda i: hits[i][1])
        lens = [hits[i][2] - hits[i][1] for i in idx]
        gaps = [hits[idx[t + 1]][1] - hits[idx[t]][1] for t in range(len(idx) - 1)]
        return (kind, hits[idx[0]][0], hits[idx[0]][1], hits[idx[-1]][2], len(idx), hits[idx[0]][3],
                int(statistics.median(lens)), int(statistics.median(gaps)) if gaps else 0)

    for idx in a_runs:
        rows.append(row("A", idx))
    for idx in b_groups.values():
        rows.append(row("B", idx))
    na = sum(len(r) for r in a_runs)
    nb = sum(len(v) for v in b_groups.values())
    n = len(hits)
    pa, pb = 100.0 * na / max(n, 1), 100.0 * nb / max(n, 1)
    big = max((r[4] for r in rows), default=0)
    null, null_lens = (null_b(hits, nperm) if (nperm and nb) else (0.0, {}))
    obs_lens = collections.Counter(len(v) for v in b_groups.values())
    lmin = long_run_min(obs_lens, null_lens) if nb else None
    nlong = sum(1 for v in b_groups.values() if lmin is not None and len(v) >= lmin)
    return rows, {"hits": n, "hits_in_A": na, "hits_in_B": nb, "pct_A": pa, "pct_B": pb, "null_B": null, "excess_B": pb - null,
                  "loci_A": len(a_runs), "loci_B": len(b_groups), "largest": big, "long_min": lmin, "long_runs": nlong,
                  "null_lens": null_lens, "obs_lens": dict(obs_lens)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("beds", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mono", type=int, default=MIN_MONO)
    ap.add_argument("--mono-gap", type=int, default=MONO_GAP)
    ap.add_argument("--min-span", type=int, default=MIN_SPAN)
    ap.add_argument("--flag-pct", type=float, default=FLAG_PCT)
    a = ap.parse_args()
    with open(a.out + ".runs.tsv", "w") as ro, open(a.out + ".summary.tsv", "w") as so:
        ro.write("consensus\tkind\tcontig\tstart\tend\thits\tstrand\tmedian_hit_len\tmedian_gap\n")
        so.write("consensus\thits\thits_in_A\tpct_A\tloci_A\thits_in_B\tpct_B\tnull_B\texcess_B\tloci_B\tlargest_locus\tflag\n")
        for p in a.beds:
            name = os.path.basename(p).rsplit(".", 1)[0]
            hits = read_bed(p)
            rows, s = screen(hits, a.min_mono, a.mono_gap, a.min_span)
            for r in rows:
                ro.write(name + "\t" + "\t".join(str(x) for x in r) + "\n")
            flag = "SAT_A" if s["pct_A"] >= a.flag_pct else ("SAT_B" if s["excess_B"] >= a.flag_pct else "-")
            so.write("%s\t%d\t%d\t%.1f\t%d\t%d\t%.1f\t%.1f\t%.1f\t%d\t%d\t%s\n" % (
                name, s["hits"], s["hits_in_A"], s["pct_A"], s["loci_A"], s["hits_in_B"], s["pct_B"], s["null_B"], s["excess_B"],
                s["loci_B"], s["largest"], flag))
            print("%s: %d hits; A runs %d (%d hits, %.1f %%); B runs %d (%d hits, %.1f %%, chance %.1f %%, excess %.1f %%) -> %s" % (
                name, s["hits"], s["loci_A"], s["hits_in_A"], s["pct_A"], s["loci_B"], s["hits_in_B"], s["pct_B"], s["null_B"],
                s["excess_B"], flag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
