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

Usage: satellite_kindB_verify.py LOCI.bed --searches DIR --genome G.fa --out PREFIX [--max-units 12] [--min-id 85] [--threads 8] [--serial]

Speed (2026-10-06, docs/SATELLITES.md 5i): only the lag-1/lag-2 pairs are aligned, each once, by --threads ssearch36 processes at once,
without statistics (-z -1); --serial runs the earlier all-against-all implementation (13 h 24 min vs 9.5 min on rsi, same verdicts).
"""
import argparse
import bisect
import collections
import concurrent.futures
import glob
import os
import re
import statistics
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))     # same directory in the repo and in a run's tools/
import satellite_trf_verify as sv  # noqa: E402  (faidx_regions)

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


def verify_runs_serial(runs, starts, genome, workdir, max_units=MAX_UNITS, min_id=MIN_ID, threads=8):
    """The reference implementation (until 2026-10-06): one all-against-all ssearch36 call per run, runs one after the other.
    Kept for validation (--serial, SATELLITE_KINDB_SERIAL=1); verify_runs gives the same verdicts and identities.
    runs: list of (consensus, contig, start, end, unit, hits); starts: {consensus: {contig: sorted [(hit start, hit end)]}}.
    Returns {run index: (units_tested, median_identity or None, frac_ge_min or None, verdict)}."""
    sel = runs
    # regions to cut: unit i = [start_i, start_{i+1})
    regions = []     # (run index, unit index, contig, s, e)
    for i, (q, c, s, e, unit, nh) in enumerate(sel):
        hs = [x for x in starts.get(q, {}).get(c, []) if s <= x[0] <= e]
        hs.sort()
        for k in range(min(len(hs) - 1, max_units)):
            regions.append((i, k, c, hs[k][0], hs[k + 1][0]))
    tmp = tempfile.mkdtemp(prefix="kindB_", dir=workdir)
    seqs = {}
    # random access through samtools faidx when available (one genome pass per consensus was the cost of the stage, audit D10);
    # the streaming path below gives the same sequences
    got = sv.faidx_regions(genome, [(c, s, e) for i, k, c, s, e in regions])
    if got is not None:
        for (i, k, c, s, e), sq in zip(regions, got):
            seqs[(i, k)] = sq
    else:
        by_c = collections.defaultdict(list)
        for r in regions:
            by_c[r[2]].append(r)
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


# --- fast path (2026-10-06) ------------------------------------------------------------------------------------------------------
# The verdict needs only the identity of each unit with the next one and the one after (lags 1 and 2). The serial path aligned every
# unit with every unit of its run (132 alignments for 12 units, 21 used) in one ssearch36 call per run, one run after the other. With
# the 30 kb array tier (13093c9) units grew to 30 kb and the stage took 13 h 24 min on rsi (2 739 runs; 23 min before the tier).
# Here each needed pair is aligned once (a pair shared by a narrow run and the wide run that contains it, or met again for another
# consensus at the same coordinates, is not aligned again), each query unit against its lag-1/lag-2 partners in its own ssearch36
# call, the calls in parallel, with -z -1 (no statistics: the verdict never reads an E-value). Without statistics ssearch36 reports the
# best alignment of each strand and stops; with them it went on searching alternative, non-overlapping local alignments of the pair
# (24 lines for one 26 x 20 kb pair of a real rsi array: 56 s; -z -1: 2 lines, 5.5 s, the same best alignment 99.92 % over 13 445 bp).
# The best alignment is what the identity of a real array comes from, so ARRAY identities are unchanged; a pair of unrelated units can
# score a little lower (an alternative local match no longer listed), which only lowers COPIES numbers far below the threshold.
SSEARCH_OPTS = ["-m", "8", "-E", "10", "-z", "-1"]
_PAIR_CACHE = {}          # (genome, query region, partner region) -> best score of the pair, None when ssearch36 reported no line


def _hits_in(lst, s, e):
    """the (start, end) hits of a sorted list whose start lies in [s, e]; the serial path's linear scan, by bisection"""
    return lst[bisect.bisect_left(lst, (s,)):bisect.bisect_right(lst, (e, float("inf")))]


def _align_query(qseq, partners, tmpdir, ssearch="ssearch36"):
    """partners: [(region, sequence)]. One ssearch36 call of the query against its partners.
    Returns ({partner region: best score or None}, return code)."""
    fd, qf = tempfile.mkstemp(suffix=".q.fa", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        fh.write(">q\n%s\n" % qseq)
    fd, lf = tempfile.mkstemp(suffix=".l.fa", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        for j, (pr, ps) in enumerate(partners):
            fh.write(">p%d\n%s\n" % (j, ps))
    try:
        r = subprocess.run([ssearch] + SSEARCH_OPTS + ["-T", "1", qf, lf], capture_output=True, text=True)
    finally:
        os.unlink(qf)
        os.unlink(lf)
    out = {pr: None for pr, _ in partners}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12 or not f[1].startswith("p"):
            continue
        pr, ps = partners[int(f[1][1:])]
        pid, aln = float(f[2]), int(f[3])
        L = min(len(qseq), len(ps))
        score = pid * min(aln / float(L), 1.0)                # identity over the shorter unit (as in the serial path)
        if out[pr] is None or score > out[pr]:
            out[pr] = score
    return out, r.returncode


def verify_runs(runs, starts, genome, workdir, max_units=MAX_UNITS, min_id=MIN_ID, threads=8):
    """runs: list of (consensus, contig, start, end, unit, hits); starts: {consensus: {contig: sorted [(hit start, hit end)]}}.
    Returns {run index: (units_tested, median_identity or None, frac_ge_min or None, verdict)}, as verify_runs_serial does (same
    verdicts; COPIES identities can be lower, see SSEARCH_OPTS). threads = ssearch36 processes at once, each single-threaded."""
    if os.environ.get("SATELLITE_KINDB_SERIAL") == "1":
        return verify_runs_serial(runs, starts, genome, workdir, max_units, min_id, threads)
    cuts = []
    for q, c, s, e, unit, nh in runs:
        hs = _hits_in(starts.get(q, {}).get(c, []), s, e)
        cuts.append([(c, hs[k][0], hs[k + 1][0]) for k in range(min(len(hs) - 1, max_units))])
    uniq = list(dict.fromkeys(r for cut in cuts for r in cut))
    seq = {}
    got = sv.faidx_regions(genome, uniq)
    if got is not None:
        seq = dict(zip(uniq, got))
    else:
        by_c = collections.defaultdict(list)
        for r in uniq:
            by_c[r[0]].append(r)
        for name, sq in read_fasta_stream(genome):
            for r in by_c.get(name, ()):
                seq[r] = sq[r[1]:r[2]]
    gk = os.path.abspath(genome)
    kept = []                                                # per run: [(unit index, region)] of the units >= 50 bp
    need = collections.defaultdict(set)                      # query region -> partner regions still to align
    for cut in cuts:
        units = [(k, r) for k, r in enumerate(cut) if len(seq.get(r, "")) >= 50]
        kept.append(units)
        byk = dict(units)
        for k1, r1 in units[:-1]:
            for lag in (1, 2):                               # neighbour or the one after (dimeric arrays)
                r2 = byk.get(k1 + lag)
                if r2 is not None and (gk, r1, r2) not in _PAIR_CACHE:
                    need[r1].add(r2)
    failed = 0
    if need:
        tmp = tempfile.mkdtemp(prefix="kindB_", dir=workdir)
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
                # the most expensive calls first, so that the short ones fill the cores at the end
                order = sorted(need, key=lambda r: -(r[2] - r[1]) * sum(p[2] - p[1] for p in need[r]))
                futs = {ex.submit(_align_query, seq[r1], [(r2, seq[r2]) for r2 in sorted(need[r1])], tmp): r1 for r1 in order}
                for fu in concurrent.futures.as_completed(futs):
                    r1 = futs[fu]
                    res, rc = fu.result()
                    failed += rc != 0
                    for r2, sc in res.items():
                        _PAIR_CACHE[(gk, r1, r2)] = sc
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)
    if failed:
        print("kindB_verify: WARNING: %d ssearch36 calls exited with an error (their pairs count as unaligned)" % failed, file=sys.stderr)
    results = {}
    for i, units in enumerate(kept):
        if len(units) < 2:
            results[i] = (len(units), None, None, "-")
            continue
        byk = dict(units)
        ids = []
        for k1, r1 in units[:-1]:
            sc = [_PAIR_CACHE.get((gk, r1, byk[k1 + lag])) for lag in (1, 2) if (k1 + lag) in byk]
            sc = [x for x in sc if x is not None]
            ids.append(max(sc) if sc else 0.0)
        med = statistics.median(ids)
        frac = sum(1 for x in ids if x >= min_id) / float(len(ids))
        results[i] = (len(units), med, frac, "ARRAY" if med >= min_id else "COPIES")
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
    ap.add_argument("--serial", action="store_true", help="the reference implementation (one all-against-all ssearch36 call per run, one run at a time)")
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
    results = (verify_runs_serial if a.serial else verify_runs)(sel, starts, a.genome, os.path.dirname(os.path.abspath(a.out)) or ".", a.max_units, a.min_id, a.threads)
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
