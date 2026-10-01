#!/usr/bin/env python3
"""CpG-adjusted divergence for one SINEderella plate alignment.

Usage:
    cpg_div.py PLATE.aln.fa [--out TSV]

Plate format (SINEderella plates):
    - FASTA. The first record is the consensus, every other record is a copy.
    - All rows are aligned to the same length; '-' is a gap.
    - Lower case = flank, upper case = element.

For every copy, over the element columns (consensus upper-case and not a
gap), ignoring columns where the copy has a gap:

    raw_p        mismatches / aligned non-gap positions
    cpg_adj_p    the same ratio after RepeatMasker's CpG-adjusted rule: at
                 consensus CpG sites (the C or the G of a CG dinucleotide) a
                 transition C->T or G->A counts 1/10 of a substitution, and
                 every other change counts 1
    k2p_adj      Kimura 2-parameter distance with the transition proportion P
                 and the transversion proportion Q taken after the same CpG
                 rule: K = -0.5*ln(1-2P-Q) - 0.25*ln(1-2Q)
    n_aligned    element columns where the copy is not a gap
    n_cpg_sites  CpG sites among those aligned columns

CpG sites are consensus positions that are the C or the G of a CG
dinucleotide; consensus gaps are removed first (a gap in the consensus row is
an insertion in a copy, not a consensus base), and case is ignored.

NA is written whenever a value cannot be computed: no aligned columns, or a
saturated K2P logarithm (non-positive argument).

Output TSV (stdout, or the file given with --out), one row per copy:
    copy  raw_p  cpg_adj_p  k2p_adj  n_aligned  n_cpg_sites
"""

import argparse
import math
import sys

# RepeatMasker CpG-adjusted rule: a C->T or G->A transition at a consensus
# CpG site counts this fraction of a substitution.
CPG_TRANSITION_WEIGHT = 0.1

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

HEADER = ("copy", "raw_p", "cpg_adj_p", "k2p_adj", "n_aligned", "n_cpg_sites")


def parse_fasta(path):
    """Parse FASTA; yield (name, sequence). name = first word of the header."""
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                words = line[1:].split()
                name = words[0] if words else ""
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def element_columns(cons):
    """Element columns of the consensus: upper-case and not a gap."""
    return [i for i, ch in enumerate(cons) if ch != "-" and ch.isupper()]


def cpg_columns(cons):
    """Columns that are the C or the G of a CG dinucleotide in the consensus.

    Consensus gaps are removed first, so a C and a G separated only by a
    consensus gap still form a CG. Case is ignored (a lower-case flank base
    can still complete a CG).
    """
    sites = set()
    prev_col = None
    prev_base = ""
    for col, ch in enumerate(cons):
        if ch == "-":
            continue
        base = ch.upper()
        if prev_base == "C" and base == "G":
            sites.add(prev_col)
            sites.add(col)
        prev_base = base
        prev_col = col
    return sites


def k2p_distance(p, q):
    """Kimura 2-parameter distance; None when saturated (non-positive log)."""
    d1 = 1.0 - 2.0 * p - q
    d2 = 1.0 - 2.0 * q
    if d1 <= 0.0 or d2 <= 0.0:
        return None
    return -0.5 * math.log(d1) - 0.25 * math.log(d2)


def copy_stats(cons, copy, elem, cpg_sites):
    """Divergence of one copy against the consensus (see module docstring).

    cons: consensus row of the alignment
    copy: copy row of the alignment (same length as cons)
    elem: element columns (consensus upper-case, not a gap)
    cpg_sites: set of consensus CpG columns (cpg_columns())
    """
    n_aligned = 0
    mismatches = 0
    n_cpg_sites = 0
    n_cpg_transition = 0    # C->T / G->A at CpG sites (weight 1/10)
    n_other_transition = 0
    n_transversion = 0
    for col in elem:
        if copy[col] == "-":
            continue  # gaps ignored
        n_aligned += 1
        cons_base = cons[col].upper()
        copy_base = copy[col].upper()
        at_cpg = col in cpg_sites
        if at_cpg:
            n_cpg_sites += 1
        if cons_base == copy_base:
            continue
        mismatches += 1
        if (cons_base, copy_base) in TRANSITIONS:
            # At a CpG site the consensus is C or G, so a transition there
            # is exactly C->T or G->A.
            if at_cpg:
                n_cpg_transition += 1
            else:
                n_other_transition += 1
        else:
            n_transversion += 1
    if n_aligned == 0:
        return {"raw_p": None, "cpg_adj_p": None, "k2p_adj": None,
                "n_aligned": 0, "n_cpg_sites": 0}
    p = (CPG_TRANSITION_WEIGHT * n_cpg_transition +
         n_other_transition) / n_aligned
    q = n_transversion / n_aligned
    weighted = (CPG_TRANSITION_WEIGHT * n_cpg_transition +
                n_other_transition + n_transversion)
    return {
        "raw_p": mismatches / n_aligned,
        "cpg_adj_p": weighted / n_aligned,
        "k2p_adj": k2p_distance(p, q),
        "n_aligned": n_aligned,
        "n_cpg_sites": n_cpg_sites,
    }


def compute_plate(records):
    """Compute stats for every copy of a plate.

    records: [(name, sequence), ...]; records[0] is the consensus row.
    Returns (consensus_name, [(copy_name, stats), ...],
             n_element_columns, n_cpg_element_columns).
    """
    if not records:
        raise ValueError("no FASTA records")
    cons_name, cons = records[0]
    elem = element_columns(cons)
    cpg_sites = cpg_columns(cons)
    n_cpg_elem = sum(1 for col in elem if col in cpg_sites)
    copies = []
    for name, seq in records[1:]:
        if len(seq) != len(cons):
            raise ValueError(f"row {name!r} has length {len(seq)} but the "
                             f"consensus has {len(cons)}")
        copies.append((name, copy_stats(cons, seq, elem, cpg_sites)))
    return cons_name, copies, len(elem), n_cpg_elem


def fmt(value):
    """Format one TSV value; NA for None; floats rounded to 12 decimals."""
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    value = round(float(value), 12)
    if value == 0.0:
        value = 0.0  # never print "-0.0"
    return str(value)


def plate_tsv(records):
    """Full TSV report for one plate (as a string)."""
    _, copies, _, _ = compute_plate(records)
    lines = ["\t".join(HEADER)]
    for name, stats in copies:
        lines.append("\t".join(
            [name, fmt(stats["raw_p"]), fmt(stats["cpg_adj_p"]),
             fmt(stats["k2p_adj"]), fmt(stats["n_aligned"]),
             fmt(stats["n_cpg_sites"])]))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="CpG-adjusted divergence for one plate alignment")
    parser.add_argument("plate",
                        help="plate alignment FASTA (row 1 = consensus)")
    parser.add_argument("--out", metavar="TSV",
                        help="write the TSV here instead of stdout")
    args = parser.parse_args(argv)

    try:
        records = list(parse_fasta(args.plate))
    except OSError as exc:
        print(f"cpg_div: cannot read {args.plate}: {exc}", file=sys.stderr)
        return 1
    if not records:
        print(f"cpg_div: no FASTA records in {args.plate}", file=sys.stderr)
        return 1
    try:
        report = plate_tsv(records)
    except ValueError as exc:
        print(f"cpg_div: {args.plate}: {exc}", file=sys.stderr)
        return 1

    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(report)
    else:
        sys.stdout.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
