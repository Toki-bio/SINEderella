#!/usr/bin/env python3
"""Apply his display convention to every alignment: flanks degapped, not aligned.

Two defects he pointed at on NEGTRUNC5__saq__s5_5seqs:

  1. the consensus is over-extended on the left and not supported by the copies
  2. the left flank is aligned, with gaps through it, instead of degapped

They are one defect. His postprocess_flanks() in extract_alignments.sh defines
the element as the consensus row's first and last non-gap column, and degaps
everything outside it. That is correct only when the consensus is right. On
s5_5seqs the consensus's leftmost 12 bases sit in columns where only 34 of 100
copies have any base at all, so the body boundary lands 40+ columns too far
left, and the true left flank ends up INSIDE the body - where his degapping
never reaches. Measured across the corpus: the consensus spans 1.75x its own
length even in clean POS sets (43 % of the body is interior gap columns), so
this is systematic, not one bad file.

So: find the element edges by copy SUPPORT, then apply his rule.

  support(col) = fraction of copies with a non-gap base in that column

Walk in from each end over the consensus's own non-gap columns and stop at the
first one reaching MIN_SUPPORT. Everything outside becomes flank and is
degapped, lowercased and butted against the element exactly as his awk does:
internal gaps removed, bases pushed to the element edge, outer side padded so
the column count never changes.

The consensus bases that get trimmed are NOT deleted - they move into the flank
and are lowercased with everything else, so an over-extended consensus stays
visible as lowercase leader rather than being silently hidden.
"""
import glob
import io
import os
import sys

MIN_SUPPORT = 0.50      # a consensus column is real if half the copies reach it
GAPS = "-."


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


SEED_TAG = "_seed_as_searched"   # row-2 name before 2026-09-28; now the original keeps its plain name
EXT_SUFFIX = "_extended"         # row 1, the consensus rebuilt from the copies, since 2026-09-28


def base_name(name):
    """<subfamily> from a consensus row name (drops _R_, the description and the _extended suffix)."""
    n = name.split()[0] if name.split() else name
    n = n[3:] if n.startswith("_R_") else n
    return n[:-len(EXT_SUFFIX)] if n.endswith(EXT_SUFFIX) else n


def is_seed(name, cons_name=None):
    """The original consensus row, not a copy: every per-copy measurement must skip it.

    Row 2 of a published plate is the consensus exactly as searched (SINEderella add_seed_row.py):
    named <subfamily>, with row 1 named <subfamily>_extended; before 2026-09-28 it was
    <subfamily>_seed_as_searched. Pass the consensus row's name to recognise the new form."""
    if SEED_TAG in name:
        return True
    return bool(cons_name) and base_name(name) == base_name(cons_name) and name != cons_name


ADD_SUPPORT_MIN = 0.50   # an addition counts toward the judged element only if the copies carry it


def _addition_support(path):
    """{'5': support, '3': support} for this plate from proposals.tsv next to it (ungapped identity
    of the copies' own flanks to the proposed bases; unrelated DNA ~0.25), or {}."""
    if not path:
        return {}
    tsv = os.path.join(os.path.dirname(os.path.abspath(path)), "proposals.tsv")
    if not os.path.isfile(tsv):
        return {}
    plate = os.path.basename(path)
    with io.open(tsv, encoding="utf-8") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if not lines:
        return {}
    head = lines[0]
    for r in lines[1:]:
        if r and r[0] == plate:
            d = dict(zip(head, r))
            out = {}
            for s in ("5", "3"):
                try:
                    out[s] = float(d.get("add%s_ungapped" % s) or "nan")
                except ValueError:
                    pass
            return out
    return {}


def _continuation_sides(path):
    """{'5','3'} sides with shared sequence past the original (continuation.tsv next to the plate)."""
    if not path:
        return set()
    tsv = os.path.join(os.path.dirname(os.path.abspath(path)), "continuation.tsv")
    if not os.path.isfile(tsv):
        return set()
    plate = os.path.basename(path)
    with io.open(tsv, encoding="utf-8") as fh:
        return {r[1] for r in (l.rstrip("\n").split("\t") for l in fh)
                if len(r) >= 4 and r[0] == plate and r[2] in ("ends", "unresolved") and r[3] != "0"}


def judged_span(names, seqs, ci, path=None):
    """(lo, hi) of the element as the verdict judges it: the consensus row's letters without its
    lowercase TRIM proposals (lowercase inside the original's span, row 2 since 2026-09-28), and
    with its lowercase ADDITIONS only where the copies carry them (proposals.tsv ungapped support
    >= ADD_SUPPORT_MIN; without a proposals.tsv all additions count, as before the marking).
    Judging only the uppercase span pushed A-tails and TSDs into the flank (rsi r3 core 82 -> 1);
    judging every addition let background junk decide the call (msc MEG-RS: +205/+84 bp of a
    tandem copy group's shared flank at support 0.28/0.31 -> NO_ELEMENT).
    Plates without an original row: uppercase span, else all letters (the old rule)."""
    row = seqs[ci]
    orig = [i for i in range(len(names)) if i != ci and is_seed(names[i], names[ci])]
    if orig:
        oc = [j for j, c in enumerate(seqs[orig[0]]) if c not in GAPS]
        if oc:
            # the original's MAIN block (stray end blocks the copies do not reach are trim proposals,
            # not the span - SINEderella add_seed_row.mark does the same)
            from continuation import main_block
            oc = main_block(oc, [s for i, s in enumerate(seqs) if i != ci and i not in orig], seqs[orig[0]])
        sup = _addition_support(path)
        # A side where continuation.py found sequence the copies share past the original keeps its
        # additions whatever their position-by-position support: that measure collapses at the
        # first indel (rsi r10 3' 0.39, r2 3' 0.37 - shared tails), and dropping them turned the
        # shared tail into "shared flank" and the core to SMALL_CORE (Not SINE 45).
        cont = _continuation_sides(path)
        drop5 = sup.get("5", 1.0) < ADD_SUPPORT_MIN and "5" not in cont
        drop3 = sup.get("3", 1.0) < ADD_SUPPORT_MIN and "3" not in cont
        o_row = seqs[orig[0]]
        keep = [j for j, c in enumerate(row) if c not in GAPS
                and not (c.islower() and oc and oc[0] <= j <= oc[-1])
                and not (c.islower() and o_row[j] not in GAPS)      # a restored stray original letter
                and not (c.islower() and oc and j < oc[0] and drop5)
                and not (c.islower() and oc and j > oc[-1] and drop3)]
    else:
        keep = [j for j, c in enumerate(row) if c.isupper()] or [j for j, c in enumerate(row) if c not in GAPS]
    return (keep[0], keep[-1]) if keep else (0, len(row) - 1)


def consensus_index(names):
    for i, h in enumerate(names):
        if "CONSENSUS" in h.upper():
            return i
    return 0


def element_bounds(cons, others):
    """His rule (first/last non-gap consensus base), corrected by copy support."""
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not others:
        return (nz[0], nz[-1]) if nz else (0, len(cons) - 1)
    n = float(len(others))
    sup = {}
    for j in nz:
        sup[j] = sum(1 for s in others if j < len(s) and s[j] not in GAPS) / n

    lo = next((j for j in nz if sup[j] >= MIN_SUPPORT), nz[0])
    hi = next((j for j in reversed(nz) if sup[j] >= MIN_SUPPORT), nz[-1])
    if hi <= lo:                       # nothing supported - keep his plain rule
        return nz[0], nz[-1]
    return lo, hi


def justify(seq, lo, hi):
    """His postprocess_flanks, verbatim in behaviour.

    left flank  : gaps dropped, lowercased, right-justified against the element
    right flank : gaps dropped, lowercased, left-justified against the element
    body        : forced UPPERCASE - see below
    width       : unchanged

    His awk lowercases the flanks and leaves the body alone, which is enough in
    his pipeline because his extracted copies arrive uppercase. This corpus does
    not: NEGLINEORF__teu__r00 is 98788 lowercase against 19674 uppercase because
    the genome's soft-masking was carried straight through, while
    NEGTRUNC5__saq__s5_5seqs is uppercase throughout. If the body is left as-is,
    lowercase stops meaning "flank" and the viewer's case cue is noise. So the
    body is uppercased explicitly and only the flank is lowered.
    """
    left_len = lo
    lf = "".join(c.lower() for c in seq[:lo] if c not in GAPS)
    lf = "-" * (left_len - len(lf)) + lf

    body = seq[lo:hi + 1].upper()

    right_len = len(seq) - hi - 1
    rf = "".join(c.lower() for c in seq[hi + 1:] if c not in GAPS)
    rf = rf + "-" * (right_len - len(rf))

    return lf + body + rf


def fix(path, out_path):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    cons = seqs[ci]
    others = [s for i, s in enumerate(seqs) if i != ci]

    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz:
        return None
    lo, hi = element_bounds(cons, others)
    trimmed_l = sum(1 for j in nz if j < lo)
    trimmed_r = sum(1 for j in nz if j > hi)

    out = []
    for i, s in enumerate(seqs):
        s = s.ljust(len(cons), "-")
        # the consensus is justified too: its trimmed ends become lowercase
        # flank, so an over-extended consensus stays visible instead of vanishing
        out.append(justify(s, lo, hi))

    with io.open(out_path, "w", encoding="utf-8") as fh:
        for h, s in zip(names, out):
            fh.write(u">%s\n" % h)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + u"\n")
    return trimmed_l, trimmed_r, len(nz)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "alignments"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    n = 0
    tl = tr = 0
    worst = []
    for p in files:
        r = fix(p, p)
        if not r:
            continue
        n += 1
        a, b, ln = r
        tl += a
        tr += b
        if a + b:
            worst.append((a + b, a, b, ln, os.path.basename(p)))
    worst.sort(reverse=True)
    print("rewrote %d alignments" % n)
    print("consensus bases moved into the flank: %d left, %d right" % (tl, tr))
    print("alignments with an over-extended consensus: %d" % len(worst))
    print("\nworst 15:")
    print("  %-46s %5s %5s %6s" % ("set", "left", "right", "conslen"))
    for t, a, b, ln, nm in worst[:15]:
        print("  %-46s %5d %5d %6d" % (nm[:46], a, b, ln))


if __name__ == "__main__":
    main()
