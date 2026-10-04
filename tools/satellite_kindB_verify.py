#!/usr/bin/env python3
"""satellite_kindB_verify.py: are the regularly spaced runs (kind B of satellite_screen / satellite_stage) real tandem arrays?

A kind-B run is geometric: >= 5 full-length hits of one consensus at regular spacing. That is what a tandem array of a unit longer than the
SINE looks like, but also what a local cluster of ordinary copies can look like. The test is the sequence: in an array the stretch between
one hit and the next is a copy of the unit, so units are near-identical (rsi MEG-RS: 99.9 %); between ordinary copies it is unrelated
flank. Arrays can be dimeric: the rsi MEG-RS array alternates 2 155 and 1 460 bp units that share only ~800 bp, so a unit is compared with
the next one AND the one after (lags 1 and 2) and the better identity counts. For every run the units are cut from hit start to hit start
(at most --max-units), aligned pairwise (ssearch36, both strands), and the run gets the median identity and the share of units whose best
neighbour identity is >= --min-id. Verdict ARRAY when the median identity is >= --min-id (85: the oldest rsi MEG-RS arrays sit at 89 %,
ordinary copies at 25-55 %).

Input: loci.bed of satellite_stage (kind B rows) and the per-consensus hit beds (gen-*.bed, the pre-exclusion .before_satellites if present).
Output: PREFIX.kindB.tsv (consensus, locus, unit_bp, hits, units_tested, median_id, frac_pairs_ge_min, verdict) and a summary per consensus.

Usage: satellite_kindB_verify.py LOCI.bed --searches DIR --genome G.fa --out PREFIX [--max-units 12] [--min-id 90] [--threads 8]
"""
import argparse
import collections
import glob
import os
import re
import statistics
import subprocess
import sys
import tempfile

MAX_UNITS = 12
MIN_ID = 85.0


def read_fasta_stream(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def bed_for(searches, q):
    for pat in ("gen-*part_%s.bed.before_satellites" % q, "gen-*part_%s.bed" % q, "gen-%s.bed.before_satellites" % q, "gen-%s.bed" % q):
        f = glob.glob(os.path.join(searches, pat))
        if f:
            return f[0]
    return None


def verify_runs(runs, starts, genome, workdir, max_units=MAX_UNITS, min_id=MIN_ID, threads=8):
    """runs: list of (consensus, contig, start, end, unit, hits); starts: {consensus: {contig: sorted [(hit start, hit end)]}}.
    Returns {run index: (units_tested, median_identity or None, frac_ge_min or None, verdict)}."""
    sel = runs
    # regions to cut: unit i = [start_i, start_{i+1})
    regions = []     # (run index, unit index, contig, s, e)
    for i, (q, c, s, e, unit, nh) in enumerate(sel):
        hs = [x for x in starts.get(q, {}).get(c, []) if s <= x[0] <= e]
        hs.sort()
        for k in range(min(len(hs) - 1, max_units)):
            regions.append((i, k, c, hs[k][0], hs[k + 1][0]))
    by_c = collections.defaultdict(list)
    for r in regions:
        by_c[r[2]].append(r)
    tmp = tempfile.mkdtemp(prefix="kindB_", dir=workdir)
    seqs = {}
    for name, seq in read_fasta_stream(genome):
        for i, k, c, s, e in by_c.get(name, ()):
            seqs[(i, k)] = seq[s:e]
    results = {}
    for i, run in enumerate(sel):
        units = [(k, seqs[(i, k)]) for k in range(max_units) if (i, k) in seqs and len(seqs[(i, k)]) >= 50]
        if len(units) < 2:
            results[i] = (len(units), None, None, "-")
            continue
        qf = os.path.join(tmp, "q%d.fa" % i)
        lf = os.path.join(tmp, "l%d.fa" % i)
        with open(qf, "w") as fq, open(lf, "w") as fl:
            for k, sq in units[:-1]:
                fq.write(">u%d\n%s\n" % (k, sq))
            for k, sq in units:
                fl.write(">u%d\n%s\n" % (k, sq))
        res = subprocess.run(["ssearch36", "-m", "8", "-E", "10", "-z", "11", "-T", str(threads), qf, lf], capture_output=True, text=True).stdout
        best = {}
        for line in res.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            k1, k2 = int(f[0][1:]), int(f[1][1:])
            if k2 - k1 not in (1, 2):                                # neighbour or the one after (dimeric arrays)
                continue
            pid, aln = float(f[2]), int(f[3])
            L = min(len(dict(units)[k1]), len(dict(units)[k2]))
            cov = aln / float(L)
            score = pid * min(cov, 1.0)                       # identity over the shorter unit; a short local match counts little
            if k1 not in best or score > best[k1]:
                best[k1] = score
        ids = [best.get(k, 0.0) for k, _ in units[:-1]]
        med = statistics.median(ids)
        frac = sum(1 for x in ids if x >= min_id) / float(len(ids))
        results[i] = (len(units), med, frac, "ARRAY" if med >= min_id else "COPIES")
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    return results


def load_starts(searches, consensuses):
    starts = {}
    for q in consensuses:
        b = bed_for(searches, q)
        d = collections.defaultdict(list)
        if b:
            for line in open(b):
                f = line.split("\t")
                if len(f) >= 3:
                    d[f[0]].append((int(f[1]), int(f[2])))
        for c in d:
            d[c].sort()
        starts[q] = d
    return starts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("loci")
    ap.add_argument("--searches", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-units", type=int, default=MAX_UNITS)
    ap.add_argument("--min-id", type=float, default=MIN_ID)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-runs-per-consensus", type=int, default=400)
    a = ap.parse_args()
    runs = []
    for line in open(a.loci):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "B":
            runs.append((f[1], f[2], int(f[3]), int(f[4]), int(f[5]), int(float(f[6]))))
    per = collections.Counter(r[0] for r in runs)
    keep = collections.Counter()
    sel = []
    for r in sorted(runs, key=lambda r: -r[5]):          # the longest runs of each consensus first
        if keep[r[0]] < a.max_runs_per_consensus:
            keep[r[0]] += 1
            sel.append(r)
    starts = load_starts(a.searches, per)
    results = verify_runs(sel, starts, a.genome, os.path.dirname(os.path.abspath(a.out)) or ".", a.max_units, a.min_id, a.threads)
    summ = collections.defaultdict(lambda: [0, 0, 0])
    with open(a.out + ".kindB.tsv", "w") as o:
        o.write("consensus\tcontig\tstart\tend\tunit_bp\thits\tunits_tested\tmedian_unit_identity\tfrac_pairs_ge_%d\tverdict\n" % int(a.min_id))
        for i, (q, c, s, e, unit, nh) in enumerate(sel):
            n, med, frac, v = results[i]
            o.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" % (q, c, s, e, unit, nh, n, "%.1f" % med if med is not None else "-", "%.2f" % frac if frac is not None else "-", v))
            summ[q][0] += 1
            summ[q][1] += v == "ARRAY"
            summ[q][2] += v == "COPIES"
    with open(a.out + ".kindB.summary.tsv", "w") as o:
        o.write("consensus\truns_tested\tARRAY\tCOPIES\n")
        for q, (n, ar, cp) in sorted(summ.items()):
            o.write("%s\t%d\t%d\t%d\n" % (q, n, ar, cp))
            print("kindB_verify: %s: %d runs tested, %d ARRAY, %d COPIES" % (q, n, ar, cp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
