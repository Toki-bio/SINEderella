#!/usr/bin/env python3
"""Derive the element window from the copies, not from the consensus.

The defect this replaces, in measure_c.py:91 and justify_all.py:24:

    nz = np.where(cons != GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])

The window is the consensus's first and last non-gap column, taken on faith.
Nothing checks it against the copies, so an over-extended consensus puts real
flank inside the element (where justify never degaps it) and an under-extended
one puts real element outside (where justify degaps it away).

The rule here is step7's, applied per column of an existing alignment instead of
by re-extracting windows from the genome:

  occupancy(col)    fraction of copies with a base there
  conservation(col) fraction of those bases that are the majority base

step7 calls a window background when its ELEVATED fraction - pairs above 45 %
identity - falls to what random genomic regions give. The per-column analogue of
"elevated" is conservation > 0.45; random DNA sits near 0.25-0.30, and measured
on soft-masked Timema the true genomic background is 0.272.

A column is part of the element when occupancy >= OCC and conservation >= CONS.
Starting from the consensus's own supported core, walk outward on each side and
stop after MISS consecutive failing columns, so a single ragged column does not
end the element. That both trims an over-extended consensus and extends an
under-extended one, with the same rule and no per-set tuning.

CONS_EDGE (0.50) is stricter: used when assigning letters to the consensus row
after step8a (`rebuild_consensus_row.py`).  Columns at 0.45-0.49 are elevated
but ambiguous vs random DNA (~0.25-0.30) and must not become consensus bases.
Display justify keeps the 0.45 walk; consensus rebuild does not re-derive span
from copy walks (span stays from border-loop geometry in row 0).
"""
import os

import numpy as np

# BOUNDARY_MODE=seed (opt-in, 2026-09-27; default stays walk until tested): the element is the seed consensus row's
# own span, first to last nucleotide - the project's consensus-anchored rule. No outward
# column walk; growth past the seed is done per copy on ungapped flanks by
# correct_published_aln.py (extend_unique_flank). The column walk took A-rich columns that
# pass 45 % agreement easily (rsi r7 rand100 5' TATAAATAAAAA) and loses shared sequence
# that MAFFT scatters over columns. BOUNDARY_MODE=walk (default) is the old column walk.
MODE = os.environ.get("BOUNDARY_MODE", "walk")

OCC = 0.50      # half the copies must reach the column
CONS = 0.45     # step7's elevated cutoff — element window walk / display
CONS_EDGE = 0.50  # consensus row: letter only when clearly above noise
MISS = 8        # consecutive failures that end the element
GAPS = set("-.")


def column_supported(rows, j, cons_min=CONS):
    """True when copies treat column j as part of the element at cons_min."""
    o, c = column_stats(rows, j)
    return o >= OCC and c >= cons_min


def trim_consensus_edges(letters, rows, lo, hi, cons_min=CONS_EDGE):
    """Drop outermost consensus letters that fail the stricter edge cutoff.

    Skips internal gaps so a trailing gap cannot block trimming spurious
    letters on the far side of it (oma_SINE18 3 prime ATT case).
    """
    def leftmost():
        for j in range(lo, hi + 1):
            if letters[j] != "-":
                return j
        return None

    def rightmost():
        for j in range(hi, lo - 1, -1):
            if letters[j] != "-":
                return j
        return None

    while True:
        j = leftmost()
        if j is None or column_supported(rows, j, cons_min):
            break
        letters[j] = "-"
    while True:
        j = rightmost()
        if j is None or column_supported(rows, j, cons_min):
            break
        letters[j] = "-"
    j = leftmost()
    k = rightmost()
    if j is None or k is None:
        return lo, hi
    return j, k


def column_stats(rows, j):
    b = [r[j] for r in rows if j < len(r) and r[j] not in GAPS]
    if not b:
        return 0.0, 0.0
    occ = len(b) / float(len(rows))
    u = [x.upper() for x in b]
    top = max(u.count(c) for c in set(u))
    return occ, top / float(len(u))


def seed_window(cons):
    """The seed consensus row's first and last non-gap column."""
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz:
        return 0, len(cons) - 1, {"mode": "seed", "note": "empty consensus row"}
    return nz[0], nz[-1], {"mode": "seed", "old_span": nz[-1] - nz[0] + 1}


def element_window(cons, rows, cons_min=CONS):
    """Return (lo, hi, diagnostics) for the element, measured on the copies.

    cons_min controls how strict the per-column test is.  Display justify uses
    the default CONS=0.45 (step7 elevated cutoff).  Consensus rebuild uses
    CONS_EDGE=0.50 so borderline columns at 0.45-0.49 never become letters.
    """
    L = len(cons)
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not rows:
        return (0, L - 1, {})
    old_lo, old_hi = nz[0], nz[-1]

    ok = {}

    def good(j):
        if j not in ok:
            o, c = column_stats(rows, j)
            ok[j] = (o >= OCC and c >= cons_min)
        return ok[j]

    # start from the consensus columns that are themselves supported; if none
    # are, fall back to the consensus's own span so we never return nothing
    core = [j for j in nz if good(j)]
    if not core:
        return (old_lo, old_hi, {"note": "no supported consensus column"})
    lo, hi = core[0], core[-1]

    miss = 0
    j = lo - 1
    while j >= 0 and miss < MISS:
        if good(j):
            lo = j
            miss = 0
        else:
            miss += 1
        j -= 1

    miss = 0
    j = hi + 1
    while j < L and miss < MISS:
        if good(j):
            hi = j
            miss = 0
        else:
            miss += 1
        j += 1

    d = {
        "old_lo": old_lo, "old_hi": old_hi, "new_lo": lo, "new_hi": hi,
        "trimmed_left": max(0, lo - old_lo), "trimmed_right": max(0, old_hi - hi),
        "extended_left": max(0, old_lo - lo), "extended_right": max(0, hi - old_hi),
        "old_span": old_hi - old_lo + 1, "new_span": hi - lo + 1,
        "cons_bp": len(nz),
        "hit_left_edge": lo == 0, "hit_right_edge": hi == L - 1,
    }
    return (lo, hi, d)
