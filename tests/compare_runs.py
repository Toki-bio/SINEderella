#!/usr/bin/env python3
"""compare_runs.py RUN_A RUN_B [--loci OUT.tsv]   exit 0 only when two runs of the same input concluded the same.

Compared, whole lines and bytes (docs/REPRODUCIBILITY.md):
  loci          genome.clean_step1/extracted.fasta            the loci and their sequences
  hits          genome.clean_step1/searches/all_hits.labeled.bed
  assignment    results/assignment_full.tsv                   subfamily, summed bit score, votes, status, threshold, runner-up
  verdict       <step2 output>/loci_verdict.tsv               SINE / grey / not_SINE of every copy
  plates        results/alignments/*.fa                       every published alignment, byte for byte
  cuts          the per-consensus score cuts (step 1: "cut for" lines of step1.stdout.log; step 2: consensus_cuts.tsv)
  consensuses   consensuses.publish.fa
  audit         results/consensus_audit/summary.tsv           the rebuild of every consensus and the identity between its copies
                                                              (both are drawn with fixed seeds, so they must repeat exactly)
A missing file in both runs is skipped (a run without publish has no plates); missing in one run is a difference.
--loci writes one row per locus: subfamily, status and verdict in A and in B (to see which loci changed).
Used by tests/test_reproducible_run.sh.
"""
import glob
import os
import sys


def lines(p):
    with open(p, encoding="utf-8", errors="replace") as fh:
        return [l.rstrip("\n") for l in fh]


def step2_dir(run):
    d = sorted(glob.glob(os.path.join(run, "step2", "step2_output*")), key=os.path.getmtime)
    return d[-1] if d else os.path.join(run, "step2", "step2_output")


def table(p, cols):
    """first column -> tuple of the named columns (header names); {} when the file is missing"""
    if not os.path.isfile(p):
        return {}
    rows = [l.split("\t") for l in lines(p)]
    if not rows:
        return {}
    idx = [rows[0].index(c) if c in rows[0] else None for c in cols]
    return {r[0]: tuple(r[i] if i is not None and i < len(r) else "" for i in idx) for r in rows[1:] if r and r[0]}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    loci_out = sys.argv[sys.argv.index("--loci") + 1] if "--loci" in sys.argv else None
    if loci_out in args:
        args.remove(loci_out)
    if len(args) != 2:
        raise SystemExit(__doc__)
    A, B = (os.path.abspath(x) for x in args)
    sa, sb = step2_dir(A), step2_dir(B)
    files = [("loci", "genome.clean_step1/extracted.fasta", "genome.clean_step1/extracted.fasta"),
             ("hits", "genome.clean_step1/searches/all_hits.labeled.bed", "genome.clean_step1/searches/all_hits.labeled.bed"),
             ("assignment", "results/assignment_full.tsv", "results/assignment_full.tsv"),
             ("verdict", os.path.relpath(os.path.join(sa, "loci_verdict.tsv"), A), os.path.relpath(os.path.join(sb, "loci_verdict.tsv"), B)),
             ("consensuses", "consensuses.publish.fa", "consensuses.publish.fa"),
             ("cuts step 2", os.path.relpath(os.path.join(sa, "consensus_cuts.tsv"), A), os.path.relpath(os.path.join(sb, "consensus_cuts.tsv"), B))]
    # the tables of step 3, of the main run and of the family view (--families): similarity scores, the table of all loci,
    # the summary per subfamily (2026-11-23: the family view's sim_scores.tsv differed between two full-genome runs)
    for sub, da, db in (("", os.path.relpath(sa, A), os.path.relpath(sb, B)),
                        ("family view ", "family_view/step2/step2_output", "family_view/step2/step2_output")):
        for f in ("sim_scores.tsv", "all_sines.bedlike.ALL.tsv", "summary.by_subfam.tsv", "unassigned.tsv", "assignment_full.tsv", "ties.tsv"):
            files.append((sub + f, da + "/" + f, db + "/" + f))
    for f in ("family_view/family_view.tsv", "family_view/family_view.summary.tsv",
              "results/consensus_audit/summary.tsv"):
        files.append((f, f, f))
    bad = 0
    for name, ra, rb in files:
        pa, pb = os.path.join(A, ra), os.path.join(B, rb)
        ea, eb = os.path.isfile(pa), os.path.isfile(pb)
        if not ea and not eb:
            continue
        if ea != eb:
            print("%-12s DIFFERENT: present in A %s, in B %s" % (name, ea, eb))
            bad += 1
            continue
        la, lb = lines(pa), lines(pb)
        if la == lb:
            print("%-12s identical (%d lines)" % (name, len(la)))
        else:
            xa, xb = set(la), set(lb)
            print("%-12s DIFFERENT: lines A %d, B %d, only in A %d, only in B %d%s" % (
                name, len(la), len(lb), len(xa - xb), len(xb - xa), "" if xa != xb else " (same lines, other order)"))
            bad += 1
    # the cut of every query of step 1 (sear prints it; tools/shuffle_cut.py)
    cuts = [[l for l in lines(os.path.join(r, "step1.stdout.log")) if l.startswith("cut for ")] if os.path.isfile(os.path.join(r, "step1.stdout.log")) else [] for r in (A, B)]
    print("%-12s %s (%d and %d queries)" % ("cuts step 1", "identical" if cuts[0] == cuts[1] else "DIFFERENT", len(cuts[0]), len(cuts[1])))
    bad += cuts[0] != cuts[1]
    pa = {os.path.basename(f) for f in glob.glob(os.path.join(A, "results", "alignments", "*.fa"))}
    pb = {os.path.basename(f) for f in glob.glob(os.path.join(B, "results", "alignments", "*.fa"))}
    if pa or pb:
        same = rows = 0
        for f in sorted(pa & pb):
            fa, fb = os.path.join(A, "results", "alignments", f), os.path.join(B, "results", "alignments", f)
            same += open(fa, "rb").read() == open(fb, "rb").read()
            rows += {l for l in lines(fa) if l.startswith(">")} == {l for l in lines(fb) if l.startswith(">")}
        n = len(pa | pb)
        ok = same == n
        print("%-12s %s: %d plates in A, %d in B, in both %d; identical files %d, same set of rows %d" % (
            "plates", "identical" if ok else "DIFFERENT", len(pa), len(pb), len(pa & pb), same, rows))
        bad += not ok
    else:
        print("%-12s in neither run (skipped)" % "plates")

    # per locus: what changed
    cols = ("Subfamily", "Status")
    ta, tb = table(os.path.join(A, "results/assignment_full.tsv"), cols), table(os.path.join(B, "results/assignment_full.tsv"), cols)
    va, vb = table(os.path.join(sa, "loci_verdict.tsv"), ("verdict",)), table(os.path.join(sb, "loci_verdict.tsv"), ("verdict",))
    loci = sorted(set(ta) | set(tb))
    both = [k for k in loci if k in ta and k in tb]
    n_sf = sum(ta[k][0] != tb[k][0] for k in both)
    n_st = sum(ta[k][1] != tb[k][1] for k in both)
    n_v = sum(va.get(k) != vb.get(k) for k in both)
    print("loci: A %d, B %d, in both %d, only A %d, only B %d; of those in both: subfamily differs %d, status differs %d, verdict differs %d" % (
        len(ta), len(tb), len(both), len(set(ta) - set(tb)), len(set(tb) - set(ta)), n_sf, n_st, n_v))
    tot = lambda v: " ".join("%s %d" % (x, sum(1 for y in v.values() if y[0] == x)) for x in sorted({y[0] for y in v.values()}))
    print("verdicts A: %s | B: %s" % (tot(va) or "-", tot(vb) or "-"))
    if loci_out:
        with open(loci_out, "w") as fh:
            fh.write("locus\tsubfamily_A\tstatus_A\tverdict_A\tsubfamily_B\tstatus_B\tverdict_B\tchanged\n")
            for k in loci:
                a = ta.get(k, ("absent", "absent")) + va.get(k, ("",))
                b = tb.get(k, ("absent", "absent")) + vb.get(k, ("",))
                fh.write("\t".join((k,) + a + b + ("yes" if a != b else "no",)) + "\n")
    print("RESULT: %s" % ("the two runs are identical" if not bad else "the two runs DIFFER (%d groups)" % bad))
    sys.exit(1 if bad else 0)


main()
