#!/usr/bin/env python3
"""Trim abandoned alignment columns from the display MSA (ends only).

Two modes (agreed 2026-09-08):
  occupancy  — variant A for oma-style publish: drop prefix/suffix columns where
               too few copies have any base (MAFFT padding for long-flank outliers).
  hybrid     — variant C for the general corpus: occupancy outer bound plus
               trim_flanks inner cap on gap fraction per side.

Rules (universal, no per-set tuning):
  * slice all rows to [col_first, col_last]; never remove columns inside element
  * min_occ = max(10, ceil(0.10 * n_rows))
  * refuse trim if median retained ungapped flank < min_flank_bp (default 25)
  * element [lo, hi] must lie fully inside the window

Ground truth: oma_SINE10_top100 → cols 430–969 after trim (occ >= 10).
"""
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fix_alignments import consensus_index, is_seed, read_fa
import trim_flanks as TF

MIN_FLANK_BP = 25
MIN_OCC_ABS = 10
MIN_OCC_FRAC = 0.10


def min_occ_threshold(n_rows):
    return max(MIN_OCC_ABS, int(math.ceil(MIN_OCC_FRAC * n_rows)))


def element_bounds(seqs, ci):
    cons = seqs[ci]
    upper = [i for i, c in enumerate(cons) if c.isupper()]
    if not upper:
        nz = [i for i, c in enumerate(cons) if c != "-"]
        if not nz:
            return None
        return nz[0], nz[-1]
    return upper[0], upper[-1]


def column_occupancy(seqs):
    L = max(len(s) for s in seqs)
    return [sum(1 for s in seqs if j < len(s) and s[j] != "-") for j in range(L)]


def occupancy_edges(occ, min_occ):
    hit = [j for j, o in enumerate(occ) if o >= min_occ]
    if not hit:
        return None
    return hit[0], hit[-1]


def median_flank_bp(seqs, ci, lo, hi, left, right):
    others = [s for i, s in enumerate(seqs) if i != ci]
    lb = [sum(1 for c in s[left:lo] if c != "-") for s in others]
    rb = [sum(1 for c in s[hi + 1:right + 1] if c != "-") for s in others]
    return float(np.median(lb)), float(np.median(rb))


def hybrid_edges(seqs, ci, lo, hi, occ_left, occ_right):
    others = [s for i, s in enumerate(seqs) if i != ci]
    lb = [sum(1 for c in s[:lo] if c != "-") for s in others]
    rb = [sum(1 for c in s[hi + 1:] if c != "-") for s in others]
    keepL = TF.width(lb, lo)
    keepR = TF.width(rb, len(seqs[ci]) - hi - 1)
    left = max(occ_left, lo - keepL)
    right = min(occ_right, hi + keepR)
    return left, right


def compute_window(seqs, ci, mode="occupancy", min_flank_bp=MIN_FLANK_BP):
    """Return (left, right, diag) or (None, None, reason) if trim refused."""
    bounds = element_bounds(seqs, ci)
    if bounds is None:
        return None, None, {"reason": "no_element"}
    lo, hi = bounds
    n = len(seqs)
    min_occ = min_occ_threshold(n)
    occ = column_occupancy(seqs)
    occ_edge = occupancy_edges(occ, min_occ)
    if occ_edge is None:
        return None, None, {"reason": "no_occupied_columns", "min_occ": min_occ}
    occ_left, occ_right = occ_edge

    if mode == "hybrid":
        left, right = hybrid_edges(seqs, ci, lo, hi, occ_left, occ_right)
    else:
        left, right = occ_left, occ_right

    if left > lo or right < hi:
        return None, None, {
            "reason": "element_outside_window",
            "lo": lo, "hi": hi, "left": left, "right": right,
        }

    med_l, med_r = median_flank_bp(seqs, ci, lo, hi, left, right)
    if med_l < min_flank_bp or med_r < min_flank_bp:
        return None, None, {
            "reason": "flank_too_short",
            "median_left_bp": med_l,
            "median_right_bp": med_r,
            "min_flank_bp": min_flank_bp,
        }

    if left == 0 and right == len(seqs[ci]) - 1:
        return None, None, {"reason": "no_trim_needed"}

    return left, right, {
        "min_occ": min_occ,
        "lo": lo, "hi": hi,
        "median_left_bp": med_l,
        "median_right_bp": med_r,
        "mode": mode,
        "width_before": len(seqs[ci]),
        "width_after": right - left + 1,
    }


def slice_rows(names, seqs, left, right):
    out = []
    for n, s in zip(names, seqs):
        s = s.ljust(right + 1, "-")
        out.append((n, s[left:right + 1]))
    return out


def trim_inplace(path, mode="occupancy", min_flank_bp=MIN_FLANK_BP):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    # the window is measured on the consensus and the copies; the original row (row 2) has no flank
    keep = [i for i in range(len(seqs)) if i == ci or not is_seed(names[i], names[ci])]
    ci_k = keep.index(ci)
    left, right, diag = compute_window([seqs[i] for i in keep], ci_k, mode=mode, min_flank_bp=min_flank_bp)
    if left is None:
        return diag

    old_elem = "".join(c for c in seqs[ci] if c.isupper())
    sliced = slice_rows(names, seqs, left, right)
    new_elem = "".join(c for c in sliced[ci][1] if c.isupper())

    with io.open(path, "w", encoding="utf-8") as fh:
        for n, s in sliced:
            fh.write(">%s\n" % n)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + "\n")

    diag["elem_len_before"] = len(old_elem)
    diag["elem_len_after"] = len(new_elem)
    diag["elem_unchanged"] = old_elem == new_elem
    diag["left"] = left
    diag["right"] = right
    return diag


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paths", nargs="+", help="alignment files (.aln.fa)")
    ap.add_argument("--mode", choices=("occupancy", "hybrid"), default="occupancy")
    ap.add_argument("--min-flank-bp", type=int, default=MIN_FLANK_BP)
    args = ap.parse_args()
    n = 0
    for p in args.paths:
        d = trim_inplace(p, mode=args.mode, min_flank_bp=args.min_flank_bp)
        if not d:
            continue
        if d.get("reason"):
            print("%-40s SKIP %s" % (os.path.basename(p), d["reason"]))
            continue
        n += 1
        print("%-40s %d->%d cols  flank med L/R %.0f/%.0f"
              % (os.path.basename(p), d["width_before"], d["width_after"],
                 d["median_left_bp"], d["median_right_bp"]))
    print("trim_display_flanks: %d trimmed" % n)


if __name__ == "__main__":
    main()
