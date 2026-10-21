#!/usr/bin/env python3
"""Correct a published alignment's orientation and element edge.

Two failures this repairs:

1. MAFFT --adjustdirection reverse-complemented every copy, and the consensus
   was then rebuilt from those copies. The seed consensus is the authority:
   if the alignment matches its reverse complement better than itself, the
   whole alignment is turned back.

2. The element edge stops where flanks become unique, and does not take a
   polyA run with it. A run already sitting at the 5' end of the consensus
   is moved back into the flank. The 3' A-tail that is already inside the
   seed is not stripped.

Flanks are the ungapped sequence outside the current consensus, position 0
touching the edge. Run this on an alignment whose flanks are already packed
against the element (the published file), so a distance in bases is a column.
"""
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extend_unique_flank as U
import boundary_justify as BJ
from fix_alignments import base_name, is_seed, justify, read_fa as read_fa_full

COMP = str.maketrans("ACGTNacgtn", "TGCANTGCAN")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def ungapped(seq):
    return "".join(c for c in seq if c not in "-.").upper()


def best_identity(query, target):
    """Identity of the shorter sequence placed on the longer one."""
    if not query or not target:
        return 0.0
    q, t = query, target
    if len(q) > len(t):
        q, t = t, q
    if len(t) < len(q):
        return 0.0
    best = 0.0
    n = len(q)
    for i in range(0, len(t) - n + 1):
        same = 0
        for k in range(n):
            if q[k] == t[i + k]:
                same += 1
        ident = same / float(n)
        if ident > best:
            best = ident
            if best == 1.0:
                return best
    return best


K = 8


def kmer_share(ref, body, k=K):
    """Fraction of ref's k-mers that occur in body. An indel costs k k-mers, not the
    rest of an ungapped placement, so old and indel-rich consensuses still score."""
    ks = {ref[i:i + k] for i in range(len(ref) - k + 1)}
    if not ks:
        return 0.0
    kb = {body[i:i + k] for i in range(len(body) - k + 1)}
    return len(ks & kb) / float(len(ks))


def should_flip(cons, ref):
    """Turn the file when the consensus shares clearly more k-mers with the
    seed's reverse complement than with the seed. The ungapped identity test
    (rev >= 0.65) left the diverged hla Rhin-1 SubFam file backwards
    (fwd 0.338, rev 0.486)."""
    body = ungapped(cons)
    fwd = kmer_share(ref.upper(), body)
    rev = kmer_share(revcomp(ref.upper()), body)
    return rev >= 0.10 and rev > 2 * fwd, round(fwd, 3), round(rev, 3)


def uppercase_bounds(seq):
    upper = [i for i, c in enumerate(seq) if c.isupper()]
    if upper:
        return upper[0], upper[-1]
    nz = [i for i, c in enumerate(seq) if c not in "-."]
    return (nz[0], nz[-1]) if nz else (0, max(0, len(seq) - 1))


def trim_leading_polya(cons, lo, hi):
    """Move lo past a 5' polyA. A base or two in front of that run goes too.

    The 3' A-tail is not touched. A real 5' end that merely contains an A
    is left alone: the polyA has to sit at the start, after at most a stub
    shorter than the run itself.
    """
    cols = [k for k in range(lo, hi + 1) if cons[k] not in "-."]
    i = 0
    stub = 0
    while i < len(cols) and cons[cols[i]].upper() != "A":
        stub += 1
        i += 1
        if stub >= U.POLYA:
            return lo
    run = 0
    while i < len(cols) and cons[cols[i]].upper() == "A":
        run += 1
        i += 1
    if run >= U.POLYA and i < len(cols):
        return cols[i]
    return lo


def fill_consensus(cons, lo, hi, lefts, rights, left_n, right_n):
    row = list(cons)
    maj_l = U.majority(lefts, left_n)
    for k in range(left_n):
        col = lo - 1 - k
        if 0 <= col < len(row) and maj_l[k] in "ACGT":
            row[col] = maj_l[k]
    maj_r = U.majority(rights, right_n)
    for k in range(right_n):
        col = hi + 1 + k
        if 0 <= col < len(row) and maj_r[k] in "ACGT":
            row[col] = maj_r[k]
    return "".join(row)


def copies_upper_span(names, seqs, ci):
    """First and last uppercase column over the copies: boundary_justify's window, which already
    includes any sequence the copies share past the original (it is left aligned there)."""
    ups = [[j for j, c in enumerate(s) if c.isupper()] for i, s in enumerate(seqs)
           if i != ci and not is_seed(names[i], names[ci])]
    ups = [u for u in ups if u]
    if not ups:
        return None
    return min(u[0] for u in ups), max(u[-1] for u in ups)


def correct(names, seqs, ref, cont_sides=()):
    """cont_sides: '5'/'3' sides where the copies share sequence past the original
    (boundary_justify.continuation_sides). There the aligned continuation is kept: no ungapped
    extension (it assumes packed flanks, where a distance in bases is a column) and no repacking."""
    ci = U.consensus_index(names)
    flipped = False
    fwd = rev = None
    if ref:
        flipped, fwd, rev = should_flip(seqs[ci], ref)
        if flipped:
            seqs = [revcomp(s) for s in seqs]
            names = [n[3:] if n.startswith("_R_") else n for n in names]
    cons = seqs[ci]
    lo, hi = uppercase_bounds(cons)
    lo2 = trim_leading_polya(cons, lo, hi)
    if lo2 != lo:
        row = list(cons)
        for j in range(lo, lo2):
            if row[j] not in "-.":
                row[j] = row[j].lower()
        cons = "".join(row)
        seqs[ci] = cons
        lo = lo2
    # flanks against the (possibly trimmed) edge
    seqs[ci] = cons
    skip = {i for i, n in enumerate(names) if i != ci and is_seed(n, names[ci])}
    lefts, rights, lo, hi = U.flanks_of(seqs, ci, skip)
    # Both sides stop at a polyA run. The A-tail already inside the seed
    # stays; a run out in the flank is not added to the consensus.
    left_n = U.extension_edge(lefts, polya="stop")
    right_n = U.extension_edge(rights, polya="stop")
    left_n = min(left_n, lo)
    right_n = min(right_n, len(cons) - hi - 1)
    # a flipped plate swaps the sides the continuation was recorded on
    sides = {{"5": "3", "3": "5"}[x] for x in cont_sides} if flipped else set(cont_sides)
    if "5" in sides:
        left_n = 0
    if "3" in sides:
        right_n = 0
    cons = fill_consensus(cons, lo, hi, lefts, rights, left_n, right_n)
    seqs[ci] = cons
    new_lo = lo - left_n
    new_hi = hi + right_n
    w = copies_upper_span(names, seqs, ci) if sides else None
    if w and "5" in sides:
        new_lo = min(new_lo, w[0])
    if w and "3" in sides:
        new_hi = max(new_hi, w[1])
    width = len(cons)
    out = []
    for i, s in enumerate(seqs):
        s = s.ljust(width, "-")[:width]
        # the original consensus (row 2) is carried unchanged: never packed or recased
        out.append(s if (i != ci and is_seed(names[i], names[ci])) else justify(s, new_lo, new_hi))
    info = {
        "flipped": flipped,
        "fwd": fwd,
        "rev": rev,
        "left": left_n,
        "right": right_n,
        "consensus": ungapped(out[ci]),
    }
    return names, out, info


def write_fa(path, names, seqs):
    ci = U.consensus_index(names)
    order = [ci] + [i for i in range(len(names)) if i != ci]
    with io.open(path, "w", encoding="utf-8") as fh:
        for i in order:
            fh.write(">%s\n" % names[i])
            s = seqs[i]
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + "\n")


def subfam_of(path):
    base = os.path.basename(path)
    # rsi_r1_9seqs_top100.aln.fa -> r1_9seqs
    parts = base.split("_")
    # rsi, r1, 9seqs, top100...
    if len(parts) >= 3 and parts[1].startswith("r") and parts[2].endswith("seqs"):
        return parts[1] + "_" + parts[2]
    return None


def load_bank(path):
    names, seqs = U.read_fa(path)
    return {n.split()[0]: ungapped(s) for n, s in zip(names, seqs)}


def ref_for(bank, path, names):
    """The forward seed for this alignment.

    The consensus header is the name SINEderella put in the file. The rsi
    filename pattern is only a fallback for a header that no longer matches.
    """
    ci = U.consensus_index(names)
    key = base_name(names[ci])
    if key in bank:
        return bank[key]
    sf = subfam_of(path)
    if sf and sf in bank:
        return bank[sf]
    return None


def main(argv):
    write = "--write" in argv
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        raise SystemExit("usage: correct_published_aln.py BANK.fa ALIGNMENT.fa... [--write]")
    bank = load_bank(args[0])
    for path in args[1:]:
        # full headers: U.read_fa keeps only the first word, and write_fa below then dropped the plate
        # marks " [soft]" / " [array]" from every row (found on the cse republish, 2026-09-28)
        names, seqs = read_fa_full(path)
        ref = ref_for(bank, path, names)
        names2, seqs2, info = correct(names, seqs, ref, BJ.continuation_sides(path))
        # info['consensus'] includes lowercase flank; show the element only
        ci = U.consensus_index(names2)
        body = "".join(c for c in seqs2[ci] if c.isupper())
        print("%s  flip=%s fwd=%s rev=%s  L+%d R+%d" % (
            os.path.basename(path), info["flipped"],
            None if info["fwd"] is None else round(info["fwd"], 3),
            None if info["rev"] is None else round(info["rev"], 3),
            info["left"], info["right"]))
        print("  5': %s" % body[:50])
        print("  3': %s" % body[-40:])
        if write:
            write_fa(path, names2, seqs2)


if __name__ == "__main__":
    main(sys.argv)
