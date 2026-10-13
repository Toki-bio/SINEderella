#!/usr/bin/env python3
"""Separation test for two subfamilies (his "50 best copies" criterion, two tiers).

    separation_test.py --assignment assignment_full.tsv --copies extracted.fasta --a SUBA --b SUBB [--n 50] [--draws 20] [--out DIR]

Question: when the copies of two subfamilies are aligned TOGETHER, do two groups still hold up?
    tier 1  the N best copies (highest bitscore, status 'assigned') of each subfamily
    tier 2  DRAWS random samples of N firm copies of each subfamily (resampled); reported as the fraction of draws that pass
Per set (2N copies, aligned together with MAFFT):
    knn_purity   mean over copies of the share of their k=5 nearest neighbours (p-distance on columns where both have a base) with the same label
    runs         label changes along the UPGMA leaf order (2 = perfectly two blocks)
    null         the same two numbers with the labels shuffled (NULL shuffles); p = share of shuffles at least as extreme
    control      one subfamily split at random into two halves (A1 vs A2): purity must sit at the null, otherwise the test sees structure that is
                 not a subfamily difference
A set "holds" if purity p <= 0.01 and runs p <= 0.01. Output: TSV line per tier (and per draw with --verbose), plus a one-line verdict.
No threshold here is calibrated; his calls decide where the line goes.
"""
from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile

import numpy as np

K = 5
NULL = 500


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def pdist_matrix(aligned):
    """p-distance over columns where both sequences have a base (gaps and N ignored); needs >= 30 shared columns else 1.0"""
    arr = np.array([list(s.upper()) for s in aligned])
    isb = np.isin(arr, list("ACGT"))
    n = len(aligned)
    d = np.ones((n, n))
    for i in range(n):
        both = isb[i] & isb
        diff = (arr[i] != arr) & both
        shared = both.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            row = np.where(shared >= 30, diff.sum(axis=1) / np.maximum(shared, 1), 1.0)
        d[i] = row
    np.fill_diagonal(d, 0.0)
    return (d + d.T) / 2


def knn_purity(d, labels, k=K):
    labels = np.asarray(labels)
    n = len(labels)
    k = min(k, n - 1)
    idx = np.argsort(d + np.eye(n) * 1e9, axis=1)[:, :k]
    return float((labels[idx] == labels[:, None]).mean())


def leaf_order(d):
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    z = linkage(squareform(d, checks=False), method="average")
    return list(leaves_list(z))


def runs_along(order, labels):
    lab = [labels[i] for i in order]
    return 1 + sum(1 for a, b in zip(lab, lab[1:]) if a != b)


def score_set(d, labels, rng, null=NULL):
    labels = list(labels)
    order = leaf_order(d)
    pur, runs = knn_purity(d, labels), runs_along(order, labels)
    pn, rn = [], []
    lab = labels[:]
    for _ in range(null):
        rng.shuffle(lab)
        pn.append(knn_purity(d, lab))
        rn.append(runs_along(order, lab))
    p_pur = (1 + sum(1 for x in pn if x >= pur)) / (null + 1)
    p_runs = (1 + sum(1 for x in rn if x <= runs)) / (null + 1)
    return {"purity": pur, "runs": runs, "null_purity": float(np.mean(pn)), "null_runs": float(np.mean(rn)),
            "p_purity": p_pur, "p_runs": p_runs, "holds": p_pur <= 0.01 and p_runs <= 0.01}


def mafft_align(named, mafft="mafft", threads=1):
    if shutil.which(mafft) is None:
        raise SystemExit(f"{mafft} not found")
    with tempfile.TemporaryDirectory() as t:
        p = os.path.join(t, "in.fa")
        with open(p, "w") as fh:
            for k, v in named:
                fh.write(f">{k}\n{v}\n")
        r = subprocess.run([mafft, "--auto", "--quiet", "--thread", str(threads), p], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode:
            raise SystemExit("mafft failed: " + r.stderr[-300:])
    out, name = {}, None
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def test_set(copies_a, copies_b, seqs, rng, aligned=False, mafft="mafft", null=NULL):
    """copies_*: lists of sequence ids; returns scores for A vs B and for the control (A split in two)."""
    ids = list(copies_a) + list(copies_b)
    labels = ["A"] * len(copies_a) + ["B"] * len(copies_b)
    if aligned:
        al = [seqs[i] for i in ids]
    else:
        m = mafft_align([(i, seqs[i]) for i in ids], mafft)
        al = [m[i] for i in ids]
    d = pdist_matrix(al)
    res = score_set(d, labels, rng, null)
    # control: A split at random into halves, B dropped
    na = len(copies_a)
    perm = list(range(na))
    rng.shuffle(perm)
    half = set(perm[: na // 2])
    sub = list(range(na))
    dc = d[np.ix_(sub, sub)]
    clab = ["A1" if i in half else "A2" for i in sub]
    ctrl = score_set(dc, clab, rng, null)
    res["control_purity"], res["control_p_purity"], res["control_runs"] = ctrl["purity"], ctrl["p_purity"], ctrl["runs"]
    return res


def read_assignment(path, sub):
    rows = []
    for i, line in enumerate(open(path)):
        f = line.rstrip("\r\n").split("\t")
        if i == 0 and f[0] == "Sequence":
            continue
        if len(f) >= 5 and f[1] == sub and f[4] == "assigned":
            rows.append((f[0], float(f[2])))
    return rows


def run(assignment, copies, a, b, n=50, draws=20, seed=1, mafft="mafft", aligned=False, null=NULL, verbose=False):
    rng = random.Random(seed)
    seqs = read_fasta(copies)
    ra, rb = read_assignment(assignment, a), read_assignment(assignment, b)
    ra = [x for x in ra if x[0] in seqs]
    rb = [x for x in rb if x[0] in seqs]
    if len(ra) < n or len(rb) < n:
        n = min(len(ra), len(rb))
    if n < 10:
        raise SystemExit(f"too few firm copies: {a}={len(ra)} {b}={len(rb)}")
    best_a = [x[0] for x in sorted(ra, key=lambda x: -x[1])[:n]]
    best_b = [x[0] for x in sorted(rb, key=lambda x: -x[1])[:n]]
    out = [("best", 0, test_set(best_a, best_b, seqs, rng, aligned, mafft, null))]
    for j in range(draws):
        sa = [x[0] for x in rng.sample(ra, n)]
        sb = [x[0] for x in rng.sample(rb, n)]
        out.append(("random", j + 1, test_set(sa, sb, seqs, rng, aligned, mafft, null)))
    return n, len(ra), len(rb), out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--draws", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--null", type=int, default=NULL)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--out", default=None, help="write separation_<A>_<B>.tsv here")
    a = ap.parse_args(argv)
    n, na, nb, out = run(a.assignment, a.copies, a.a, a.b, a.n, a.draws, a.seed, a.mafft, False, a.null, a.verbose)
    cols = ["tier", "draw", "purity", "null_purity", "p_purity", "runs", "null_runs", "p_runs", "holds", "control_purity", "control_p_purity", "control_runs"]
    lines = ["\t".join(cols)]
    for tier, j, r in out:
        lines.append("\t".join([tier, str(j)] + [f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols[2:]]))
    rnd = [r for t, _, r in out if t == "random"]
    held = sum(1 for r in rnd if r["holds"])
    ctrl_bad = sum(1 for r in rnd if r["control_p_purity"] <= 0.01)
    verdict = (f"{a.a} vs {a.b}: N={n} (firm {na} / {nb}); best-{n}: {'HOLDS' if out[0][2]['holds'] else 'does not hold'} "
               f"(purity {out[0][2]['purity']:.2f}, runs {out[0][2]['runs']}); random-{n}: holds in {held}/{len(rnd)} draws; "
               f"control (one subfamily split in two) looks like a real split in {ctrl_bad}/{len(rnd)} draws")
    if a.verbose or a.out is None:
        print("\n".join(lines))
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, f"separation_{a.a}_{a.b}.tsv"), "w") as fh:
            fh.write("\n".join(lines) + "\n# " + verdict + "\n")
    print(verdict, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
