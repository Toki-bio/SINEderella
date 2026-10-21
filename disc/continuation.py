#!/usr/bin/env python3
"""How far past the original consensus the copies stay similar to one another.

His review of rsi r1_9seqs (2026-09-28): past the right end the copies continue with the same
sequence (an A-run of variable length, then gtcctggaagtacacactgttccccaataaagtcctgttcccc...) for
~170 bp. That region must be shown ALIGNED until similarity is lost; the chain packed it as flank,
and because the A-runs differ in length the shared sequence landed at different offsets.

The column walk cannot see it: MAFFT splits such a stretch into blocks that different copies occupy,
so column occupancy jumps between 0.1 and 0.9 and the walk's half-the-copies rule fails. This
measures it per copy instead. For every copy, its own bases are read outward from the original's
edge; a base scores 1 when it equals the majority base of its column (columns reached by at least
MIN_PRESENT copies). A WIN-base window slides along the copy; the median over copies at each offset
is the similarity profile. Measured on the bat corpus (raw step8a alignments):

    rsi r1_9seqs 3'   0.95 to offset ~160, 0.68 at 180, 0.45-0.50 from 200 on
    cse MEG-T2   3'   1.00, but only 12 % of copies still have sequence at offset 60
    rda Rhin-1   3'   0.55 at the edge, 0.40 from offset 20: the element ends there
    rsi MEG-RL   5'   0.38-0.50 from the edge: no continuation (the TGGGGGAAATA case)

Unrelated flank sits at 0.40-0.47 (MAFFT pulls similar bases into shared columns even there); shared
sequence at 0.9-1.0. SIM_LOST = 0.60 is between the two.

Status per side:
    none        the profile is below SIM_LOST at the edge: nothing shared
    ends        it falls below SIM_LOST while >= MIN_COVER of the copies still have sequence
    unresolved  the copies run out (< MIN_COVER left) while still similar: extract more flank

Usage: continuation.py ALIGNMENT.aln.fa [--need]
  --need  print "<extra_5p_bp> <extra_3p_bp>" to extract (0 = resolved); step8a uses this.
"""
import statistics
import sys
from collections import Counter

WIN = 20
SIM_LOST = 0.60
MIN_COVER = 0.50
MIN_PRESENT_FRAC = 0.05
MIN_PRESENT_ABS = 3
STEP_BP = 150          # extra flank asked for per unresolved round
FAR_QUANTILE = 0.75    # the aligned stretch reaches the column 75 % of copies reach at its end
GAPS = "-."


def read_fa(path):
    names, seqs = [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
        elif seqs:
            seqs[-1].append(line.strip())
    return names, ["".join(s) for s in seqs]


def column_majority(rows, width):
    need = max(MIN_PRESENT_ABS, int(MIN_PRESENT_FRAC * len(rows)))
    maj = {}
    for j in range(width):
        b = [r[j].upper() for r in rows if j < len(r) and r[j] not in GAPS]
        if len(b) >= need:
            maj[j] = Counter(b).most_common(1)[0][0]
    return maj


STRAY_GAP, STRAY_MAX, STRAY_OCC, STRAY_ID = 5, 12, 0.50, 0.40   # as SINEderella tools/add_seed_row.py


def main_block(nz, rows, ref):
    """The reference's letter columns without stray END blocks: <= STRAY_MAX letters, > STRAY_GAP
    empty columns from the next letter, and NOT carried by the copies - fewer than STRAY_OCC of them
    reach it, or those that do match its letters below STRAY_ID (median; flank background ~0.25).
    cth Rhin-1: the original's last 6 letters sat ~50 columns past its body where 10 % of copies
    reach, so the 3' edge was taken there, cover 0.11, and ~80 bp the copies share were missed.
    Identity matters on packed plates, where flank letters fill a stray block's columns."""
    n = float(max(1, len(rows)))

    def stray(block):
        if len(block) > STRAY_MAX:
            return False
        occ = sum(sum(1 for r in rows if j < len(r) and r[j] not in GAPS) for j in block) / (n * len(block))
        if occ < STRAY_OCC:
            return True
        ids = []
        for r in rows:
            pr = [(r[j].upper(), ref[j].upper()) for j in block if j < len(r) and r[j] not in GAPS]
            if pr and len(pr) >= len(block) / 2.0:
                ids.append(sum(a == b for a, b in pr) / float(len(pr)))
        return not ids or statistics.median(ids) < STRAY_ID
    while True:
        runs, cur = [], [nz[0]]
        for j in nz[1:]:
            if j - cur[-1] > STRAY_GAP:
                runs.append(cur)
                cur = [j]
            else:
                cur.append(j)
        runs.append(cur)
        if len(runs) > 1 and stray(runs[0]):
            nz = [j for r in runs[1:] for j in r]
        elif len(runs) > 1 and stray(runs[-1]):
            nz = [j for r in runs[:-1] for j in r]
        else:
            return nz


def shared_extent(ref, rows, side, maj=None):
    """ref: the row whose first/last base is the edge (the original, row 0 of a raw step8a plate).
    -> {"status", "bp", "col", "cover_at_end"}; col = the last alignment column (outward) that
    belongs to the shared stretch in any copy, None when nothing is shared."""
    width = len(ref)
    nz = [j for j, c in enumerate(ref) if c not in GAPS]
    if not nz or not rows:
        return {"status": "none", "bp": 0, "col": None, "cover_at_end": 0.0}
    nz = main_block(nz, rows, ref)
    edge = nz[-1] if side == "3" else nz[0]
    maj = maj if maj is not None else column_majority(rows, width)
    outward = range(edge + 1, width) if side == "3" else range(edge - 1, -1, -1)
    per_copy = []
    for r in rows:
        cols = [j for j in outward if j < len(r) and r[j] not in GAPS]
        per_copy.append((cols, [1 if (j in maj and r[j].upper() == maj[j]) else 0 for j in cols]))
    n = float(len(rows))
    k, end_k, status = 0, None, None
    while True:
        vals = [sum(m[k:k + WIN]) / float(WIN) for cols, m in per_copy if len(m) >= k + WIN]
        cover = len(vals) / n
        if cover < MIN_COVER:
            status = "unresolved" if k > 0 else ("unresolved" if vals and statistics.median(vals) >= SIM_LOST else "none")
            end_k = k
            break
        if statistics.median(vals) < SIM_LOST:
            status = "ends" if k > 0 else "none"
            end_k = k
            break
        k += 1
    if status == "none":
        return {"status": "none", "bp": 0, "col": None, "cover_at_end": round(cover, 2)}
    # the shared stretch is offsets [0, end_k): its outermost column in any copy
    # (end_k == 0 with status unresolved: most copies stop at the edge itself - nothing to keep aligned)
    # The column where the shared stretch ends differs per copy (MAFFT places blocks differently).
    # Not the outermost copy: one far-placed copy dragged rsi r1's window to column 1102 of 1147.
    # The column FAR_QUANTILE of those copies have reached.
    far = sorted(cols[end_k - 1] for cols, m in per_copy if end_k > 0 and len(cols) >= end_k)
    if far:
        q = int(round(FAR_QUANTILE * (len(far) - 1)))
        col = far[q] if side == "3" else far[len(far) - 1 - q]
    else:
        col = None
    return {"status": status, "bp": end_k, "col": col, "cover_at_end": round(cover, 2)}


def extents(names, seqs, ci=0, skip=()):
    rows = [s for i, s in enumerate(seqs) if i != ci and i not in skip]
    maj = column_majority(rows, len(seqs[ci]))
    return {side: shared_extent(seqs[ci], rows, side, maj) for side in ("5", "3")}


def main(argv):
    path = argv[1]
    names, seqs = read_fa(path)
    if "--need" in argv:
        # Tandem-array units (" [array]", SINEderella tools/array_order.py) share their flank far past
        # the element; if they alone drive "unresolved", re-extraction loops to the +600 bp cap for
        # nothing (cse MEG-RS rand100: 4 rounds, > 30 min). Decide on the independent copies only.
        skip = {i for i, n in enumerate(names) if i != 0 and "[array]" in n}
        ex = extents(names, seqs, 0, skip) if len(names) - 1 - len(skip) >= 3 else {s: {"status": "none"} for s in ("5", "3")}
        print(" ".join(str(STEP_BP if ex[s]["status"] == "unresolved" else 0) for s in ("5", "3")))
        return
    ex = extents(names, seqs)
    for s in ("5", "3"):
        print("%s'  %s" % (s, ex[s]))


if __name__ == "__main__":
    main(sys.argv)
