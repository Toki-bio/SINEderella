#!/usr/bin/env python3
"""Copy-supported element window + display justify for publish alignments.

Imports boundary.element_window and fix_alignments.justify — does not
reimplement either. Replaces extract_alignments.sh postprocess_flanks on
MAFFT output (must run on the raw alignment, before any consensus-only
degapping).
"""
import glob
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boundary as B
import continuation as CT
from fix_alignments import consensus_index, is_seed, justify, read_fa

CONT_FIELDS = ["plate", "side", "status", "bp", "cover_at_end"]


def record_continuation(path, ex):
    """continuation.tsv next to the plate: one row per plate and side, replaced on re-run."""
    tsv = os.path.join(os.path.dirname(os.path.abspath(path)), "continuation.tsv")
    plate = os.path.basename(path)
    # plates of one directory may be processed in parallel: serialise the read-modify-write
    # (a parallel test run on the bat corpus interleaved rows without this)
    lock = io.open(tsv + ".lock", "w")
    try:
        import fcntl
        fcntl.flock(lock, fcntl.LOCK_EX)
    except ImportError:          # Windows: no fcntl; the chain there runs one plate at a time
        pass
    try:
        rows = []
        if os.path.isfile(tsv):
            with io.open(tsv, encoding="utf-8") as fh:
                rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()][1:]
        rows = [r for r in rows if len(r) == len(CONT_FIELDS) and r[0] != plate
                and not r[0].endswith("_subfam.aln.fa")]
        for side in ("5", "3"):
            e = ex[side]
            rows.append([plate, side, e["status"], str(e["bp"]), str(e["cover_at_end"])])
        tmp = "%s.%d.tmp" % (tsv, os.getpid())
        with io.open(tmp, "w", encoding="utf-8") as fh:
            fh.write("\t".join(CONT_FIELDS) + "\n")
            for r in sorted(rows):
                fh.write("\t".join(r) + "\n")
        os.replace(tmp, tsv)
    finally:
        lock.close()


def continuation_sides(path):
    """Sides of this plate with sequence shared past the original ('5', '3') - from continuation.tsv."""
    tsv = os.path.join(os.path.dirname(os.path.abspath(path)), "continuation.tsv")
    if not os.path.isfile(tsv):
        return set()
    plate = os.path.basename(path)
    with io.open(tsv, encoding="utf-8") as fh:
        return {r[1] for r in (l.rstrip("\n").split("\t") for l in fh)
                if r[0] == plate and r[2] in ("ends", "unresolved") and r[3] != "0"}


def flanked_element_window(cons, rows):
    """Copy-supported window for 50L/70R publish alignments.

    boundary.element_window walks outward on copy support. Border loop + rebuilt
    consensus must fix the anchor row before publish; do not cap 3' extension
    back to the seed consensus span (that blocked SINE10 GAC recovery).
    """
    if B.MODE == "walk":
        lo_b, hi_b, diag = B.element_window(cons, rows)
    else:
        lo_b, hi_b, diag = B.seed_window(cons)
    lo = lo_b
    hi = hi_b
    diag["display_lo"] = lo
    diag["display_hi"] = hi
    diag["final_span"] = hi - lo + 1
    return lo, hi, diag


def apply(path):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    cons = seqs[ci]
    others = [s for i, s in enumerate(seqs) if i != ci and not is_seed(names[i], names[ci])]
    lo, hi, diag = flanked_element_window(cons, others)
    # Sequence the copies still share past the original's ends stays ALIGNED: the window is widened
    # to it, so justify() packs only what lies beyond (continuation.py; rsi r1_9seqs 3' +189 bp).
    # Not for SubFam plates (his call 2026-09-28): their rows are chunk consensuses, element only - no
    # copy runs past the element, so "shared past the end" has no meaning and read "unresolved 0 bp".
    if os.environ.get("CONTINUATION", "1") == "1" and not path.endswith("_subfam.aln.fa"):
        # measured from the element WINDOW, not from row 0's own first/last letter: a consensus the
        # border loop widened into the flanks (cth Rhin-1: row 0 from column 0) put the edge where few
        # copies reach, and the ~80 bp the copies share past the element went unseen (cover 0.26)
        ref = "".join(c if lo <= j <= hi else "-" for j, c in enumerate(cons))
        ex = CT.extents(names, [ref] + others, 0)
        if ex["5"]["col"] is not None and ex["5"]["col"] < lo:
            diag["continuation_left"] = lo - ex["5"]["col"]
            lo = ex["5"]["col"]
        if ex["3"]["col"] is not None and ex["3"]["col"] > hi:
            diag["continuation_right"] = ex["3"]["col"] - hi
            hi = ex["3"]["col"]
        record_continuation(path, ex)

    out = []
    for i, s in enumerate(seqs):
        s = s.ljust(len(cons), "-")
        # the original consensus (row 2) is carried unchanged: never packed or recased
        out.append(s if (i != ci and is_seed(names[i], names[ci])) else justify(s, lo, hi))

    order = [ci] + [i for i in range(len(names)) if i != ci]
    with io.open(path, "w", encoding="utf-8") as fh:
        for i in order:
            h, s = names[i], out[i]
            fh.write(u">%s\n" % h)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + u"\n")
    return diag


def main():
    paths = sys.argv[1:] if len(sys.argv) > 1 else sorted(
        glob.glob(os.path.join("alignments", "*.aln.fa")))
    n = 0
    ext_l = ext_r = 0
    for p in paths:
        d = apply(p)
        if not d:
            continue
        n += 1
        ext_l += d.get("extended_left", 0)
        ext_r += d.get("extended_right", 0)
        el = d.get("extended_left", 0)
        er = d.get("extended_right", 0)
        cap = d.get("capped_right", 0)
        if el or er or d.get("trimmed_left") or d.get("trimmed_right") or cap:
            print("%-40s extend L=%d R=%d  trim L=%d R=%d  capR=%d  span %d->%d"
                  % (os.path.basename(p), el, er,
                     d.get("trimmed_left", 0), d.get("trimmed_right", 0),
                     cap, d.get("old_span", 0), d.get("final_span", 0)))
    print("boundary_justify: %d files" % n)


if __name__ == "__main__":
    main()
