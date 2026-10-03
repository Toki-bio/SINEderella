#!/usr/bin/env python3
"""satellite_stage.py: the satellite screen of SINEderella, run inside step 1 after the searches and before extraction / SubFam.

Design: docs/SATELLITES.md (sections 4 and 5). Per consensus:
  kind A  SINE-derived satellites (monomer = a part of the SINE, or the SINE plus a bit): tools/satellite_trf_verify on the raw hits of
          sear (searches/gen-<q>.rawhits.tsv: merged hits at the homology cut, any length; falls back to the length-filtered gen-<q>.bed),
          i.e. a window around every hit, TRF on the windows, unit aligned to the consensus. Every verified locus is a satellite locus.
  kind B  arrays whose unit is longer than the SINE (rsi MEG-RS): tools/satellite_screen regular-spacing runs with a chance null on the
          full-length hits (gen-<q>.bed); the family is SAT_B when the excess over chance is >= 20 % of its hits.
Writes OUT/indication.tsv (per consensus), OUT/loci.bed (kind, consensus, locus, monomers, SINE part), OUT/units.fa, OUT/<q>.kindA.loci.tsv.
Exclusion (default on; --no-exclude writes the tables only): hits of gen-<q>.bed that overlap a kind-A locus, or a kind-B run of a SAT_B
consensus (--exclude-b all: every kind-B run), are removed from gen-<q>.bed (the original is kept as gen-<q>.bed.before_satellites and the
removed hits in OUT/excluded_hits.bed), so that extraction, SubFam (the peel input) and assignment never see them. Nothing is deleted.

Usage: satellite_stage.py --searches DIR --genome G.fa --cons CONSENSUSES.fa --out OUT [--threads 16] [--no-exclude] [--exclude-b flagged|all|none]
"""
import argparse
import collections
import glob
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import satellite_screen as ss  # noqa: E402

FLAG_PCT = 20.0


def read_fa(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n is not None:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


def query_name(bed):
    """searches/gen-consensuses.clean.part_<NAME>.bed -> NAME"""
    b = os.path.basename(bed)
    b = re.sub(r"\.(bed|rawhits\.tsv)$", "", b)
    return b.split(".part_", 1)[1] if ".part_" in b else re.sub(r"^gen-", "", b)


def overlaps(lst, s, e):
    """lst sorted by start: any interval overlapping [s, e)"""
    import bisect
    i = bisect.bisect_left(lst, (s, -1))
    j = i - 1
    while j >= 0 and lst[j][1] > s:
        if lst[j][0] < e:
            return True
        j -= 1
    while i < len(lst) and lst[i][0] < e:
        if lst[i][1] > s:
            return True
        i += 1
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--searches", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--cons", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--no-exclude", action="store_true")
    ap.add_argument("--exclude-b", choices=["flagged", "all", "none"], default="flagged")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cons = read_fa(a.cons)
    beds = sorted(glob.glob(os.path.join(a.searches, "gen-*.bed")))
    beds = [b for b in beds if not b.endswith(".before_satellites")]
    if not beds:
        print("satellite_stage: no gen-*.bed in %s, nothing to do" % a.searches)
        return 0
    verify = os.path.join(HERE, "satellite_trf_verify.py")
    rows, loci, units, excluded = [], [], [], []
    for bed in beds:
        q = query_name(bed)
        if q not in cons:
            # the bed may carry the file-part name while the bank carries the id; try both
            cand = [k for k in cons if k == q or k.replace("|", "_") == q]
            if not cand:
                print("satellite_stage: %s: consensus %s not in %s, skipped" % (os.path.basename(bed), q, a.cons))
                continue
            q = cand[0]
        raw = re.sub(r"\.bed$", ".rawhits.tsv", bed)
        hits_src = raw if os.path.exists(raw) and os.path.getsize(raw) else bed
        cfa = os.path.join(a.out, q + ".cons.fa")
        open(cfa, "w").write(">%s\n%s\n" % (q, cons[q]))
        pref = os.path.join(a.out, q + ".kindA")
        subprocess.run([sys.executable, verify, hits_src, "--genome", a.genome, "--cons", cfa, "--out", pref, "--threads", str(a.threads)],
                       check=False, stdout=open(pref + ".log", "w"), stderr=subprocess.STDOUT)
        A = []
        if os.path.exists(pref + ".loci.tsv"):
            for l in list(open(pref + ".loci.tsv"))[1:]:
                f = l.rstrip("\n").split("\t")
                A.append((f[0], int(f[1]), int(f[2]), int(f[3]), float(f[4]), int(f[9]), int(f[10])))
        if os.path.exists(pref + ".units.fa"):
            units.append(open(pref + ".units.fa").read())
        # kind B on the full-length hits
        full = ss.read_bed(bed)
        brows, bsum = ss.screen(full)
        flagB = bsum["excess_B"] >= FLAG_PCT
        B = [(r[1], r[2], r[3], r[4], r[7]) for r in brows if r[0] == "B"]      # contig start end hits median_gap
        mono = int(sum(x[4] for x in A))
        rows.append((q, len(full), len(A), mono, max((x[4] for x in A), default=0), len(B), bsum["pct_B"], bsum["excess_B"],
                     "SAT_A" if A else "-", "SAT_B" if flagB else "-"))
        for c, s, e, per, cop, cs, ce in A:
            loci.append(("A", q, c, s, e, per, cop, "%d-%d" % (cs, ce)))
        for c, s, e, n, gap in B:
            loci.append(("B", q, c, s, e, gap, n, "-"))
        print("satellite_stage: %s: %d full hits; kind A %d loci (%d monomers); kind B %d runs, excess %.1f %% %s%s" % (
            q, len(full), len(A), mono, len(B), bsum["excess_B"], "SAT_A " if A else "", "SAT_B" if flagB else ""), flush=True)
        # exclusion
        if a.no_exclude:
            continue
        ex = collections.defaultdict(list)
        for c, s, e, per, cop, cs, ce in A:
            ex[c].append((s, e))
        if a.exclude_b == "all" or (a.exclude_b == "flagged" and flagB):
            for c, s, e, n, gap in B:
                ex[c].append((s, e))
        if not ex:
            continue
        for c in ex:
            ex[c].sort()
        keep, drop = [], []
        for line in open(bed):
            f = line.rstrip("\n").split("\t")
            try:
                c, s, e = f[0], int(f[1]), int(f[2])
            except (ValueError, IndexError):
                keep.append(line)
                continue
            (drop if overlaps(ex.get(c, []), s, e) else keep).append(line)
        if drop:
            if not os.path.exists(bed + ".before_satellites"):
                shutil.copy2(bed, bed + ".before_satellites")
            open(bed, "w").writelines(keep)
            excluded += [q + "\t" + l for l in drop]
            print("satellite_stage: %s: %d of %d full hits removed from %s (kept in .before_satellites)" % (q, len(drop), len(keep) + len(drop), os.path.basename(bed)), flush=True)
    with open(os.path.join(a.out, "indication.tsv"), "w") as o:
        o.write("consensus\tfull_hits\tkindA_loci\tkindA_monomers\tkindA_largest\tkindB_runs\tkindB_pct\tkindB_excess_pct\tflag_A\tflag_B\n")
        for r in rows:
            o.write("%s\t%d\t%d\t%d\t%.0f\t%d\t%.1f\t%.1f\t%s\t%s\n" % r)
    with open(os.path.join(a.out, "loci.bed"), "w") as o:
        o.write("#kind\tconsensus\tcontig\tstart\tend\tperiod_or_unit\tmonomers_or_hits\tsine_part\n")
        for l in sorted(loci, key=lambda x: (x[2], x[3])):
            o.write("\t".join(str(x) for x in l) + "\n")
    open(os.path.join(a.out, "units.fa"), "w").write("".join(units))
    with open(os.path.join(a.out, "excluded_hits.bed"), "w") as o:
        o.writelines(excluded)
    na = sum(r[2] for r in rows)
    print("satellite_stage: %d consensuses, %d kind-A loci, %d consensuses SAT_B, %d hits excluded -> %s" % (
        len(rows), na, sum(1 for r in rows if r[9] == "SAT_B"), len(excluded), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
