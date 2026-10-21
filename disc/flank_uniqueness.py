#!/usr/bin/env python3
"""Per-side flank uniqueness: do copies share genomic context?

Each flank (left and right separately) is scanned for groups of similar
sequences among the copies.  Shared flanks mean the loci are not independent
insertions — satellite, duplication, or nested repeat.

Severity (by copy set tier):
  rand100 / random  — high flag if a large fraction share a flank
  top100            — medium flag (may be one real subgroup)
  other             — medium by default

Method: extract ungapped flank per copy (element-adjacent first), pairwise
identity on the overlapping prefix, single-linkage clustering at share_thr.
Report cluster structure and the fraction of copies with a genuinely unique
flank (no partner above threshold).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fix_alignments import consensus_index, is_seed, judged_span, read_fa

SHARE_THR = 0.55          # match verdict.FLANK_SHARE
MIN_FLANK_BP = 25         # copies with shorter flanks skipped in matrix
MIN_CLUSTER = 2           # shared group size floor


def tier_from_path(path):
    base = os.path.basename(path).replace(".aln.fa", "")
    if "rand100" in base or "random" in base:
        return "rand100"
    if "top100" in base or "best50" in base:
        return "top100"
    return "other"


def element_bounds(seq):
    upper = [i for i, c in enumerate(seq) if c.isupper()]
    if upper:
        return upper[0], upper[-1]
    nz = [i for i, c in enumerate(seq) if c != "-"]
    return (nz[0], nz[-1]) if nz else (0, len(seq) - 1)


def extract_flanks(seqs, ci, names=None, path=None):
    cons = seqs[ci]
    # judged_span: the element as rebuilt, without the lowercase trim proposals (fix_alignments)
    lo, hi = judged_span(names, seqs, ci, path) if names else element_bounds(cons)
    lefts, rights, idx = [], [], []
    for i, s in enumerate(seqs):
        if i == ci or (names and is_seed(names[i], names[ci])):
            continue
        l = "".join(c for c in s[:lo] if c != "-")
        r = "".join(c for c in s[hi + 1:] if c != "-")
        lefts.append(l[::-1])   # index 0 = adjacent to element
        rights.append(r)
        idx.append(i)
    return lefts, rights, lo, hi


def pair_identity(a, b, maxlen=120):
    """Identity on equal-length prefix, ungapped."""
    m = min(len(a), len(b), maxlen)
    if m < MIN_FLANK_BP:
        return 0.0
    return sum(1 for i in range(m) if a[i] == b[i]) / float(m)


def identity_matrix(flanks, maxlen=120):
    n = len(flanks)
    M = np.eye(n)
    for i in range(n):
        if len(flanks[i]) < MIN_FLANK_BP:
            continue
        for j in range(i + 1, n):
            if len(flanks[j]) < MIN_FLANK_BP:
                continue
            v = pair_identity(flanks[i], flanks[j], maxlen=maxlen)
            M[i, j] = M[j, i] = v
    return M


def cluster_labels(M, thr=SHARE_THR):
    """Single-linkage components at thr."""
    n = M.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] >= thr:
                union(i, j)
    roots = [find(i) for i in range(n)]
    uniq = {}
    labels = []
    for r in roots:
        if r not in uniq:
            uniq[r] = len(uniq)
        labels.append(uniq[r])
    return np.array(labels, int)


def side_report(flanks, usable_mask=None):
    n = len(flanks)
    if n < 4:
        return {"n": n, "measured": False, "reason": "too_few_copies"}

    lens = np.array([len(f) for f in flanks])
    ok = lens >= MIN_FLANK_BP
    if usable_mask is not None:
        ok &= usable_mask
    n_ok = int(ok.sum())
    if n_ok < 4:
        return {"n": n, "n_measured": n_ok, "measured": False,
                "reason": "too_few_flanks"}

    sub = [flanks[i] for i in range(n) if ok[i]]
    M = identity_matrix(sub)
    labels = cluster_labels(M)
    sizes = {}
    for lb in labels:
        sizes[lb] = sizes.get(lb, 0) + 1
    cluster_list = sorted(sizes.values(), reverse=True)
    n_shared = sum(s for s in cluster_list if s >= MIN_CLUSTER)
    n_unique = sum(1 for s in cluster_list if s == 1)
    largest = cluster_list[0] if cluster_list else 1
    n_multi = sum(1 for s in cluster_list if s >= MIN_CLUSTER)

    # weighted mean identity inside multi-member clusters
    within = []
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            if labels[i] == labels[j] and labels[i] >= 0:
                if M[i, j] > 0:
                    within.append(M[i, j])
    mean_within = float(np.mean(within)) if within else 0.0

    return {
        "n": n,
        "n_measured": n_ok,
        "measured": True,
        "n_clusters": len(cluster_list),
        "n_multi_clusters": n_multi,
        "largest_cluster": int(largest),
        "largest_cluster_frac": round(largest / float(n_ok), 3),
        "unique_frac": round(n_unique / float(n_ok), 3),
        "shared_copy_frac": round(n_shared / float(n_ok), 3),
        "mean_within_shared": round(mean_within, 3),
        "cluster_sizes": cluster_list[:8],
        "median_flank_bp": int(np.median(lens[ok])),
    }


def flag_text(side_name, side, tier, severity):
    base = (
        "%s flank: %.0f%% of copies sit in shared groups "
        "(largest cluster %.0f%%, mean within %.2f)."
        % (side_name.capitalize(),
           100 * side.get("shared_copy_frac", 0),
           100 * side.get("largest_cluster_frac", 0),
           side.get("mean_within_shared", 0))
    )
    if tier == "rand100" and severity == "high":
        return base + " Random copies should not share flanks — strong evidence of non-independent loci."
    if tier == "top100":
        return base + " Top copies may include one genomic-context subgroup; check other loci for unique flanks."
    return base + " Shared flanks warrant inspection."


def flag_side(side, tier):
    if not side.get("measured"):
        return None
    frac = side.get("shared_copy_frac", 0.0)
    largest = side.get("largest_cluster_frac", 0.0)
    # a "group" of one copy is not shared: with few measured copies a singleton alone passed the
    # 15 % largest-group bar (cse MEG-TR rand100, 6 copies, all unique -> flagged "high", 2026-09-28)
    if side.get("largest_cluster", 0) < 2:
        return None
    if largest < 0.10 and frac < 0.15:
        return None
    if tier == "rand100":
        if largest >= 0.15 or frac >= 0.25:
            return "high"
        if largest >= 0.08 or frac >= 0.12:
            return "medium"
    else:
        if largest >= 0.25 or frac >= 0.40:
            return "high"
        if largest >= 0.10 or frac >= 0.15:
            return "medium"
    return None


def scan(path, share_thr=SHARE_THR, maxlen=120):
    names, seqs = read_fa(path)
    if len(seqs) < 5:
        return {"set": os.path.basename(path).replace(".aln.fa", ""),
                "error": "too_few_sequences"}
    ci = consensus_index(names)
    lefts, rights, lo, hi = extract_flanks(seqs, ci, names, path)
    tier = tier_from_path(path)
    L = side_report(lefts)
    R = side_report(rights)
    out = {
        "set": os.path.basename(path).replace(".aln.fa", ""),
        "tier": tier,
        "share_thr": share_thr,
        "element_cols": [lo, hi],
        "left": L,
        "right": R,
        "flags": [],
    }
    for side_name, side in ("left", L), ("right", R):
        fl = flag_side(side, tier)
        if fl:
            out["flags"].append({
                "side": side_name,
                "severity": fl,
                "largest_cluster_frac": side.get("largest_cluster_frac"),
                "shared_copy_frac": side.get("shared_copy_frac"),
                "text": flag_text(side_name, side, tier, fl),
            })
    worst = None
    for f in out["flags"]:
        if f["severity"] == "high":
            worst = "high"
        elif f["severity"] == "medium" and worst != "high":
            worst = "medium"
    out["worst_flag"] = worst
    return out


def main():
    import json
    paths = sys.argv[1:] if len(sys.argv) > 1 else []
    for p in paths:
        r = scan(p)
        print(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()
