#!/usr/bin/env python3
"""Test 2 - measurements for the three properties nothing currently captures.

From PLAN.md section 4:

  1. a subset of copies differing   - everything current is per-column; he was
                                      describing groups of ROWS
  2. where the right edge actually is - the naive rule fails because copy
                                      identity collapses in the poly-A tail
  3. internal instability in a region - a local measure, not a global one

Each is written to be checked against his 72 judgements, not against anything
generated.
"""
import io
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, ".")
import measure_c as M
import profiles as P


# ---------------------------------------------------------------- property 1
def row_groups(path, min_grp=5):
    """Does a SUBSET OF COPIES behave differently from the rest?

    His words: "11 lower sequences need manual attention"; "top proper about
    half longer similarity sequences ... then bottom more discordant shorter
    ones". Both describe groups of ROWS. Everything the tool measures today is
    per-column, which is why it cannot see this.

    Method: describe every copy by three row-level numbers - its identity to
    the consensus, its element length, and how much of the element it covers -
    then ask whether the copies fall into two groups rather than one spread.
    Reported as the separation between the two halves relative to the spread
    inside them, so it does not simply track divergence.
    """
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return {}
    k = ci[0]
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    if len(nz) < 60:
        return {}
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, k, axis=0)[:, lo:hi + 1]
    ref = cons[lo:hi + 1]
    keep = ref != M.GAP
    C, ref = C[:, keep], ref[keep]
    n = C.shape[0]
    if n < 2 * min_grp:
        return {}
    pres = C != M.GAP
    agree = (C == ref[None, :]) & pres
    with np.errstate(invalid="ignore"):
        ident = agree.sum(axis=1) / np.maximum(pres.sum(axis=1), 1)
    elen = pres.sum(axis=1).astype(float)
    cover = elen / C.shape[1]

    F = np.vstack([ident, elen / max(1.0, elen.max()), cover]).T
    F = (F - F.mean(axis=0)) / (F.std(axis=0) + 1e-9)

    # split on the leading direction of variation
    u, s, vt = np.linalg.svd(F - F.mean(axis=0), full_matrices=False)
    proj = (F - F.mean(axis=0)) @ vt[0]
    order = np.argsort(proj)
    best = None
    for cut in range(min_grp, n - min_grp):
        a, b = proj[order[:cut]], proj[order[cut:]]
        sep = abs(b.mean() - a.mean())
        spread = (a.std() + b.std()) / 2 + 1e-9
        val = sep / spread
        if best is None or val > best[0]:
            best = (val, cut)
    val, cut = best
    grp = np.zeros(n, bool)
    grp[order[:cut]] = True
    return {"rowsplit": round(float(val), 3),
            "rowsplit_small": int(min(cut, n - cut)),
            "rowsplit_frac": round(float(min(cut, n - cut) / n), 3),
            "rowsplit_dident": round(float(abs(ident[grp].mean() - ident[~grp].mean())), 3),
            "rowsplit_dlen": round(float(abs(elen[grp].mean() - elen[~grp].mean())), 1)}


# ---------------------------------------------------------------- property 2
def right_edge(path, win=9):
    """Where does the element actually end on the 3' side?

    Identity alone fails: it collapses in the poly-A tail because A-run length
    differs between copies (AluY: 0.93 in the body, 0.27 in the last 90bp,
    while A+T rises 0.26 -> 0.63). An identity rule amputates the tail from
    every SINE.

    So: scan outward from the element centre and stop at the first position
    where identity has fallen to background AND composition is no longer
    A/T-rich. Report that position relative to the consensus end - negative
    means the element ends before the consensus does.
    """
    pr = P.profile(path)
    if pr is None:
        return {}
    xs = np.asarray(pr["x"], float)
    pid = np.asarray(pr["pair_id"], float)
    at = np.asarray(pr["at"], float)
    cov = np.asarray(pr["cover"], float)

    def sm(y):
        ok = np.isfinite(y)
        if ok.sum() < 3:
            return y
        f = np.interp(np.arange(len(y)), np.flatnonzero(ok), y[ok])
        return np.convolve(f, np.ones(win) / win, mode="same")

    pid, at = sm(pid), sm(at)
    inside = xs >= 0
    if inside.sum() < 60:
        return {}
    L = int(xs[inside].max()) + 1
    body = inside & (xs > L * 0.2) & (xs < L * 0.6)
    if body.sum() < 10:
        return {}
    plateau = float(np.nanpercentile(pid[body], 75))
    at_body = float(np.nanmedian(at[body]))
    far = ((xs <= -35) | (xs >= L + 34)) & (cov >= 0.5)
    bg = float(np.nanmedian(pid[far])) if far.sum() >= 5 else 0.25
    if plateau - bg < 0.12:
        return {"redge_usable": 0}

    half = (plateau + bg) / 2.0
    at_hi = at_body + 0.12          # "still A/T-rich" means clearly above body

    idx = np.argsort(xs)
    xo, po, ao, co = xs[idx], pid[idx], at[idx], cov[idx]
    start = int(np.searchsorted(xo, L * 0.6))
    end = None
    for i in range(start, len(xo)):
        if co[i] < 0.25:
            continue
        if po[i] < half and not (ao[i] > at_hi):
            end = xo[i]
            break
    if end is None:
        end = xo[-1]
    naive = None
    for i in range(start, len(xo)):
        if co[i] < 0.25:
            continue
        if po[i] < half:
            naive = xo[i]
            break
    return {"redge_usable": 1,
            "redge_pos": float(end),
            "redge_shift": round(float(end - (L - 1)), 1),
            "redge_naive_shift": round(float((naive if naive is not None else end) - (L - 1)), 1),
            "redge_tail_len": round(float(end - (naive if naive is not None else end)), 1)}


# ---------------------------------------------------------------- property 3
def instability(path, win=21):
    """Is there a LOCAL region of instability inside an otherwise good element?

    His words on the worst mosaic case: "pre-tail region containing an island of
    instability". A global spread cannot see this - the element is fine
    everywhere else.

    Method: identity along the element, smoothed; find the deepest sustained dip
    relative to the element's own plateau, and report its depth, its width and
    where it sits as a fraction of element length.
    """
    pr = P.profile(path)
    if pr is None:
        return {}
    xs = np.asarray(pr["x"], float)
    pid = np.asarray(pr["pair_id"], float)
    inside = xs >= 0
    if inside.sum() < 80:
        return {}
    x = xs[inside]
    y = pid[inside]
    ok = np.isfinite(y)
    if ok.sum() < 40:
        return {}
    y = np.interp(np.arange(len(y)), np.flatnonzero(ok), y[ok])
    ys = np.convolve(y, np.ones(win) / win, mode="same")
    L = len(ys)
    core = ys[int(L * 0.1):int(L * 0.9)]
    if len(core) < 30:
        return {}
    plateau = float(np.percentile(core, 75))
    dip = plateau - core
    j = int(np.argmax(dip))
    depth = float(dip[j])
    thr = depth * 0.5
    a = j
    while a > 0 and dip[a] > thr:
        a -= 1
    b = j
    while b < len(dip) - 1 and dip[b] > thr:
        b += 1
    return {"dip_depth": round(depth, 3),
            "dip_width": int(b - a),
            "dip_pos": round(float((int(L * 0.1) + j) / max(1, L)), 3)}


def main():
    rows = [l.split("\t") for l in
            io.open("calls.tsv", encoding="utf-8").read().rstrip("\n").split("\n")[1:]]
    out = {}
    for st, corp, call, txt in rows:
        p = None
        for d in (corp, "aln_v2", "bench", "tim_bench", "aln_ext", "aln_c"):
            q = os.path.join(d, st + ".aln.fa")
            if os.path.exists(q):
                p = q
                break
        if p is None:
            continue
        d = {}
        for fn in (row_groups, right_edge, instability):
            try:
                d.update(fn(p) or {})
            except Exception:
                pass
        out[st] = (call, d)

    keys = sorted({k for _, d in out.values() for k in d})
    io.open("test2_vars.tsv", "w", encoding="utf-8").write(
        "set\tcall\t" + "\t".join(keys) + "\n" +
        "\n".join(st + "\t" + c + "\t" +
                  "\t".join("" if k not in d else str(d[k]) for k in keys)
                  for st, (c, d) in out.items()) + "\n")

    g = defaultdict(list)
    for st, (c, d) in out.items():
        g[c].append(d)
    show = ["rowsplit", "rowsplit_frac", "rowsplit_dlen",
            "redge_shift", "redge_naive_shift", "redge_tail_len",
            "dip_depth", "dip_width"]
    print("alignments: %d   new measurements: %d" % (len(out), len(keys)))
    print()
    print("%-24s %3s " % ("what you said", "n") + " ".join("%>10s" % s for s in show).replace(">", ""))
    for c in sorted(g, key=lambda c: -len(g[c])):
        vals = []
        for k in show:
            v = [d[k] for d in g[c] if k in d and d[k] == d[k]]
            vals.append("%10.2f" % np.mean(v) if v else "%10s" % "-")
        print("%-24s %3d " % (c, len(g[c])) + " ".join(vals))


if __name__ == "__main__":
    main()
