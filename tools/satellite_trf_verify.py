#!/usr/bin/env python3
"""satellite_trf_verify.py: SINE-derived tandem loci (kind A satellites) from hits + TRF on the hit clusters only.

Design: docs/SATELLITES.md section 5c. Hit spacing alone misses most short arrays (the cobra Squam3C hits cover two monomers each and
most loci have 4-9 monomers), so the hits are only a GATE:
  1. gate     clusters of >= MIN_HITS hits (any length, both strands) whose neighbours are <= GAP bp apart become windows (cluster +- PAD);
              a consensus without such clusters costs nothing more.
  2. TRF      Tandem Repeats Finder on the windows only (not the genome), split over --threads processes.
  3. verify   a TRF record is a locus when period >= MIN_PERIOD, copies >= MIN_COPIES and its unit (doubled, to cover rotation) aligns to the
              consensus with ssearch36 (>= MIN_ALN bp, E <= MAX_E). Overlapping records of one locus (period 67 and its dimer 134) are
              reduced to the best by TRF score.
  4. classify the aligned part of the consensus (cons_start-cons_end) is the part of the SINE the monomer holds.

Writes PREFIX.loci.tsv (contig start end period copies pmatch score aln_len pct_id cons_start cons_end), PREFIX.units.fa (monomer consensus per
locus) and PREFIX.summary.tsv. Needs trf, ssearch36 on PATH (or --trf / --ssearch). Nothing is removed or changed.

Usage: satellite_trf_verify.py HITS.bed --genome G.fa --cons CONS.fa --out PREFIX [--threads 16] [--trf trf] [--ssearch ssearch36]
"""
import argparse
import collections
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

GAP = 300
PAD = 1000         # cobra Squam3C: 100 bp padding recovered 60 % of the TRF loci, 1 kb padding 90 % (partial windows cut arrays)
MIN_HITS = 1       # a window around every hit; 2 hits recovered 61 % at 1 kb padding
MAX_WINDOW_MB = 600
MIN_PERIOD = 55
MAX_PERIOD = 300
MIN_COPIES = 4.0
MIN_ALN = 45
MAX_E = 1e-2
TRF_ARGS = ["2", "5", "7", "80", "10", "40", str(MAX_PERIOD), "-h", "-ngs"]


def read_bed(path):
    hits = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and not line.startswith("#"):
            try:
                hits.append((f[0], int(f[1]), int(f[2])))
            except ValueError:
                pass
    return hits


def gate(hits, gap=GAP, pad=PAD, min_hits=MIN_HITS):
    """clusters of >= min_hits hits with neighbours <= gap apart -> windows (contig, start, end, n_hits), start 0-based"""
    hs = sorted(hits)
    out = []
    k = 0
    while k < len(hs):
        j, end = k, hs[k][2]
        while j + 1 < len(hs) and hs[j + 1][0] == hs[k][0] and hs[j + 1][1] - end <= gap:
            j += 1
            end = max(end, hs[j][2])
        if j - k + 1 >= min_hits:
            out.append((hs[k][0], max(0, hs[k][1] - pad), end + pad, j - k + 1))
        k = j + 1
    return out


def merge_windows(wins):
    """merge windows that touch or overlap (padding makes neighbours overlap)"""
    out = []
    for c, s, e, n in sorted(wins):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1] = (c, out[-1][1], max(e, out[-1][2]), out[-1][3] + n)
        else:
            out.append((c, s, e, n))
    return out


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


def faidx_regions(genome, regions, samtools=None):
    """[(contig, start0, end)] -> [sequence] through `samtools faidx -r` (random access through the .fai index; the index is made
    when missing, next to the genome, as step 1 does). Returns None when samtools is not available, so that callers fall back to
    streaming the FASTA. Sequences come back in the order of the regions; an end past the contig is clamped by samtools, like the
    streaming path's seq[s:min(e, len)]. Reading the whole genome once per consensus was the cost of the stage (audit D10, 2026-10-05)."""
    sam = samtools or shutil.which("samtools")
    if not sam or not regions:
        return None
    if not os.path.exists(genome + ".fai"):
        try:
            subprocess.run([sam, "faidx", genome], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, OSError):
            return None
    tmp = tempfile.NamedTemporaryFile("w", suffix=".regions", delete=False)
    try:
        for c, s, e in regions:
            tmp.write("%s:%d-%d\n" % (c, s + 1, e))
        tmp.close()
        try:
            r = subprocess.run([sam, "faidx", "-r", tmp.name, genome], capture_output=True, text=True)
        except OSError:          # samtools path given but not runnable: streaming fallback
            return None
    finally:
        os.unlink(tmp.name)
    if r.returncode != 0:
        return None
    out, cur = [], None
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            cur = []
            out.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    seqs = ["".join(x) for x in out]
    return seqs if len(seqs) == len(regions) else None


def extract_windows(genome, wins, out_fa, samtools=None):
    n = 0
    seqs = faidx_regions(genome, [(c, s, e) for c, s, e, k in wins], samtools)
    with open(out_fa, "w") as fh:
        if seqs is not None:
            for (c, s, e, k), sub in zip(wins, seqs):
                if len(sub) >= 100:
                    fh.write(">%s:%d-%d\n%s\n" % (c, s, s + len(sub), sub))
                    n += 1
            return n
        by = collections.defaultdict(list)
        for w in wins:
            by[w[0]].append(w)
        for name, seq in read_fasta_stream(genome):
            for c, s, e, k in by.get(name, ()):
                sub = seq[s:min(e, len(seq))]
                if len(sub) >= 100:
                    fh.write(">%s:%d-%d\n%s\n" % (c, s, s + len(sub), sub))
                    n += 1
    return n


def split_fasta(path, parts, outdir):
    recs = []
    cur = None
    for line in open(path):
        if line.startswith(">"):
            cur = [line, []]
            recs.append(cur)
        else:
            cur[1].append(line)
    files = [os.path.join(outdir, "w%03d.fa" % i) for i in range(parts)]
    fhs = [open(f, "w") for f in files]
    sizes = [0] * parts
    for hdr, body in sorted(recs, key=lambda r: -sum(len(x) for x in r[1])):
        i = sizes.index(min(sizes))
        fhs[i].write(hdr + "".join(body))
        sizes[i] += sum(len(x) for x in body)
    for fh in fhs:
        fh.close()
    return [f for f, s in zip(files, sizes) if s]


def run_trf(trf, files, workdir):
    procs = []
    for f in files:
        out = open(f + ".trf", "w")
        procs.append((subprocess.Popen([trf, os.path.basename(f)] + TRF_ARGS, cwd=workdir, stdout=out, stderr=subprocess.DEVNULL), out))
    for p, out in procs:
        p.wait()
        out.close()
    return [f + ".trf" for f in files]


def parse_trf(path):
    """records: (contig, start, end, period, copies, pmatch, score, unit); coordinates mapped back from the window name contig:start-end"""
    recs = []
    wc = ws = None
    for line in open(path):
        if line.startswith("@"):
            m = re.match(r"@(.+):(\d+)-(\d+)$", line.strip())
            wc, ws = (m.group(1), int(m.group(2))) if m else (None, 0)
            continue
        f = line.split()
        if wc is None or len(f) < 14:
            continue
        recs.append((wc, ws + int(f[0]), ws + int(f[1]), int(f[2]), float(f[3]), int(f[5]), int(f[7]), f[13]))
    return recs


def best_per_locus(recs):
    """overlapping records of one locus (period p and its multiples): keep the best TRF score, drop records overlapping an accepted one by > 50 %"""
    out = []
    by = collections.defaultdict(list)
    for r in sorted(recs, key=lambda r: -r[6]):
        keep = True
        for a in by[r[0]]:
            ov = min(a[2], r[2]) - max(a[1], r[1])
            if ov > 0.5 * min(a[2] - a[1], r[2] - r[1]):
                keep = False
                break
        if keep:
            by[r[0]].append(r)
            out.append(r)
    return sorted(out)


def verify(recs, cons, ssearch, workdir, threads):
    q = os.path.join(workdir, "units2.fa")
    with open(q, "w") as fh:
        for i, r in enumerate(recs):
            fh.write(">u%d\n%s\n" % (i, r[7] + r[7]))
    # ssearch36 cannot open a file whose path is longer than ~120 characters (FASTA36 file-name limit): every kind-A call failed
    # in the Sicista runs (2026-10-06). The consensus is copied next to the units and both are given by bare name inside workdir.
    shutil.copyfile(cons, os.path.join(workdir, "cons.fa"))
    cmd = [ssearch, "-m", "8", "-E", str(MAX_E), "-z", "11", "-T", str(threads), "units2.fa", "cons.fa"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        # a failed call used to read as "no SINE-derived locus" (Sicista 2026-10-06: 0 loci in every consensus); run it once more,
        # and if it fails again say so instead of returning an empty answer silently
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
        if r.returncode != 0:
            print("satellite_trf_verify: WARNING: ssearch36 exited with %d twice; the TRF units of this consensus are untested: %s"
                  % (r.returncode, " ".join(r.stderr.split())[:300]), file=sys.stderr)
    res = r.stdout
    best = {}
    for line in res.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        i, pid, aln, s0, s1, bits = int(f[0][1:]), float(f[2]), int(f[3]), int(f[8]), int(f[9]), float(f[11])
        if aln >= MIN_ALN and (i not in best or bits > best[i][3]):
            best[i] = (pid, aln, (min(s0, s1), max(s0, s1)), bits)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hits")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--cons", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-hits", type=int, default=MIN_HITS, help="hits in a cluster to open a window (1 = every hit, the default; 2 is faster and recovers less)")
    ap.add_argument("--max-window-mb", type=int, default=MAX_WINDOW_MB)
    ap.add_argument("--gap", type=int, default=GAP)
    ap.add_argument("--pad", type=int, default=PAD)
    ap.add_argument("--trf", default=shutil.which("trf") or "trf")
    ap.add_argument("--ssearch", default=shutil.which("ssearch36") or "ssearch36")
    a = ap.parse_args()
    name = os.path.basename(a.hits).rsplit(".", 1)[0]
    hits = read_bed(a.hits)
    wins = merge_windows(gate(hits, a.gap, a.pad, a.min_hits))
    mb = sum(w[2] - w[1] for w in wins) / 1e6
    if mb > a.max_window_mb and a.min_hits < 2:      # a hit-rich family: fall back to clusters of >= 2 hits and say so
        print("%s: windows %.0f Mb exceed --max-window-mb %d with --min-hits 1; using --min-hits 2 (less sensitive)" % (name, mb, a.max_window_mb), flush=True)
        a.min_hits = 2
        wins = merge_windows(gate(hits, a.gap, a.pad, a.min_hits))
        mb = sum(w[2] - w[1] for w in wins) / 1e6
    print("%s: %d hits -> %d windows (%.1f Mb)" % (name, len(hits), len(wins), mb), flush=True)
    rows, units = [], []
    if wins:
        tmp = tempfile.mkdtemp(prefix="satverify_", dir=os.path.dirname(os.path.abspath(a.out)) or ".")
        try:
            wfa = os.path.join(tmp, "windows.fa")
            n = extract_windows(a.genome, wins, wfa)
            files = split_fasta(wfa, min(a.threads, max(n, 1)), tmp)
            recs = []
            for t in run_trf(a.trf, files, tmp):
                recs += [r for r in parse_trf(t) if r[3] >= MIN_PERIOD and r[4] >= MIN_COPIES]
            recs = best_per_locus(recs)
            print("%s: %d TRF loci (period >= %d, copies >= %d)" % (name, len(recs), MIN_PERIOD, MIN_COPIES), flush=True)
            if recs:
                best = verify(recs, a.cons, a.ssearch, tmp, a.threads)
                for i, r in enumerate(recs):
                    if i in best:
                        rows.append((r, best[i]))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    with open(a.out + ".loci.tsv", "w") as lo, open(a.out + ".units.fa", "w") as uo:
        lo.write("contig\tstart\tend\tperiod\tcopies\tpmatch\tscore\taln_len\tpct_id\tcons_start\tcons_end\tspan_ok\n")
        for r, b in rows:
            span_ok = "yes" if (b[2][1] - b[2][0] + 1) <= 1.15 * r[3] else "no"      # the aligned part of the SINE cannot exceed the monomer
            lo.write("%s\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%.1f\t%d\t%d\t%s\n" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], b[1], b[0], b[2][0], b[2][1], span_ok))
            uo.write(">%s:%d-%d|p%d|c%.1f\n%s\n" % (r[0], r[1], r[2], r[3], r[4], r[7]))
    cop = [r[4] for r, b in rows]
    per = collections.Counter(r[3] for r, b in rows)
    cons_part = collections.Counter((b[2][0] // 10 * 10, b[2][1] // 10 * 10) for r, b in rows)
    with open(a.out + ".summary.tsv", "w") as so:
        so.write("consensus\thits\twindows\twindows_mb\tloci\tmonomers\tmedian_copies\tmax_copies\tmode_period\tcons_part_mode\n")
        so.write("%s\t%d\t%d\t%.1f\t%d\t%d\t%s\t%s\t%s\t%s\n" % (
            name, len(hits), len(wins), mb, len(rows), int(sum(cop)), ("%.1f" % statistics.median(cop)) if cop else "-",
            ("%.0f" % max(cop)) if cop else "-", per.most_common(1)[0][0] if per else "-",
            ("%d-%d" % cons_part.most_common(1)[0][0]) if cons_part else "-"))
    print("%s: %d SINE-derived tandem loci, %d monomers; period mode %s; consensus part mode %s" % (
        name, len(rows), int(sum(cop)), per.most_common(3), cons_part.most_common(2)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
