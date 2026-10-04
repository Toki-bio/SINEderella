#!/usr/bin/env python3
"""Run tools/length_variants.py on every candidate pair of a finished SINEderella run.

Candidate pairs: the shorter consensus of the bank is the 5' part of a longer one (consensus_bank_lib.find_length_variant_pairs,
>= 85 % identity along the shorter, the longer 10-120 bp longer). For each pair whose two families both have >= MIN_ASSIGNED assigned
copies the copies of both are pooled, aligned to the longer consensus, and length_variants.analyse decides
TWO_VERSIONS / UNLINKED_ENDS / SINGLE_MODE / UNRESOLVED. Nothing is merged or renamed: results/length_variants/summary.tsv is a decision input.

Usage: length_variants_run.py RUN_DIR [--threads 16] [--max-pairs 12] [--bank FILE]
RUN_DIR needs genome.clean.fa, results/assignment_full.tsv, results/assignment_stats.tsv and the bank (consensuses.clean.fa).
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import length_variants as lv  # noqa: E402
from consensus_bank_lib import find_length_variant_pairs, read_fa  # noqa: E402

MIN_ASSIGNED = 100


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--max-pairs", type=int, default=12)
    ap.add_argument("--min-id", type=float, default=0.85)
    ap.add_argument("--bank", default=None)
    a = ap.parse_args()
    run = os.path.abspath(a.run_dir)
    bank = a.bank or os.path.join(run, "consensuses.clean.fa")
    genome = os.path.join(run, "genome.clean.fa")
    hits = os.path.join(run, "results", "assignment_full.tsv")
    stats = os.path.join(run, "results", "assignment_stats.tsv")
    for f in (bank, genome, hits, stats):
        if not os.path.exists(f):
            print("length_variants_run: missing %s, skipped" % f)
            return 0
    cons = read_fa(bank)
    n_assigned = {r["Subfamily"]: int(r["Assigned"]) for r in csv.DictReader(open(stats), delimiter="\t")}
    cands = find_length_variant_pairs(cons, min_id=a.min_id)
    pairs = [p for p in cands if n_assigned.get(p["short"], 0) >= MIN_ASSIGNED and n_assigned.get(p["long"], 0) >= MIN_ASSIGNED]
    untested = [p for p in cands if p not in pairs]
    capped = []
    out = os.path.join(run, "results", "length_variants")
    os.makedirs(out, exist_ok=True)
    if len(pairs) > a.max_pairs:
        print("length_variants_run: %d candidate pairs, testing the first %d (--max-pairs); the rest are listed as NOT_TESTED" % (len(pairs), a.max_pairs))
        pairs, capped = pairs[:a.max_pairs], pairs[a.max_pairs:]
    rows = []
    for p in pairs:
        tag = "%s__%s" % (p["short"], p["long"])
        print("length_variants_run: %s (%.1f %% identity, +%d bp)" % (tag, p["identity"], p["extra"]), flush=True)
        copies, lcons, skipped, nrows = lv.collect(genome, hits, {p["short"], p["long"]}, bank, p["long"], os.path.join(out, tag + ".work"), a.threads)
        rep = lv.analyse(copies, lcons)
        txt = lv.format_report(rep, len(lcons)) + "\n(%d hits used, %d with mixed strand skipped)\n" % (nrows, skipped)
        open(os.path.join(out, tag + ".report.txt"), "w").write(txt)
        with open(os.path.join(out, tag + ".ends.tsv"), "w") as o:
            o.write("end\tcopies_5prime_complete\n")
            for i, v in enumerate(rep["hist"]):
                if v:
                    o.write("%d\t%d\n" % (i, v))
        with open(os.path.join(out, tag + ".copies.tsv"), "w") as o:
            o.write("fam\tqs\tqe\n")
            for c in copies:
                o.write("%s\t%d\t%d\n" % (c.fam, c.qs, c.qe))
        t = rep["tests"]
        lk = t.get("linkage") or {}
        rows.append({
            "short": p["short"], "long": p["long"], "consensus_identity": p["identity"], "extra_bp": p["extra"],
            "verdict": rep["verdict"], "copies_5prime_complete": rep["n_5prime_complete"],
            "mode_ends": ",".join(str(m["peak"]) for m in rep["modes"]),
            "valley_ratio": "%.3f" % t["valley_ratio_max"] if "valley_ratio_max" in t else "",
            "linkage_index": "" if lk.get("index") is None else "%.2f" % lk["index"],
            "diagnostic_columns": len(lk.get("diag_columns", {})),
            "tsd_excess_points": ",".join("%.0f" % x["excess"] for x in t.get("tsd", []) if x["excess"] is not None),
            "notes": "; ".join(rep["why"]),
        })
        import shutil
        shutil.rmtree(os.path.join(out, tag + ".work"), ignore_errors=True)
    empty = {"copies_5prime_complete": "", "mode_ends": "", "valley_ratio": "", "linkage_index": "", "diagnostic_columns": "", "tsd_excess_points": ""}
    for p in capped:     # beyond --max-pairs: listed, so that no candidate pair disappears from the table (2026-10-05)
        rows.append(dict(empty, short=p["short"], long=p["long"], consensus_identity=p["identity"], extra_bp=p["extra"], verdict="NOT_TESTED",
                         notes="beyond --max-pairs %d: rerun with LENGTH_VARIANTS_MAX_PAIRS raised, or by hand (tools/length_variants.py)" % a.max_pairs))
    for p in untested:   # a candidate pair with too few copies: flagged, to be tested on a species that has more
        rows.append(dict(empty, short=p["short"], long=p["long"], consensus_identity=p["identity"], extra_bp=p["extra"], verdict="NOT_TESTED",
                         notes="too few assigned copies (%s: %d, %s: %d; need %d each): test on another species with more copies" % (
                             p["short"], n_assigned.get(p["short"], 0), p["long"], n_assigned.get(p["long"], 0), MIN_ASSIGNED)))
    cols = ["short", "long", "consensus_identity", "extra_bp", "verdict", "copies_5prime_complete", "mode_ends", "valley_ratio",
            "linkage_index", "diagnostic_columns", "tsd_excess_points", "notes"]
    with open(os.path.join(out, "summary.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("length_variants_run: %d pair(s) tested, results/length_variants/summary.tsv" % len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
