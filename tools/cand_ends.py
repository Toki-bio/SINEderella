#!/usr/bin/env python3
"""End evidence for the open composite candidates.

Reads one or more candidate alignments in FASTA format:

  row 1     = candidate consensus, dashes where there is only flank,
  row 2..n  = copies, lower case = flank, upper case = element.

For the 5' and the 3' end of every candidate the tool reports:

  - the consensus end column (1-based),
  - the first column outside the end whose top-base share falls under
    60 %; the top-base share of a column is the share of all copies
    that show the commonest base at that column,
  - the share of copies that carry bases beyond the consensus end for
    10, 30 and 60 columns (the word "carrying"); a copy carries N
    columns beyond the end when its own sequence still reaches N
    columns past the end,
  - the commonest 12 bases outside the end, per copy the first 12
    bases past the end, with the share of copies that show them,
  - whether the outside is dominated by an A run (a run of >= 8 A in
    >= 40 % of the copies) or by an oligo-G run (a run of >= 8 G in
    >= 40 % of the copies),
  - one verdict: closed, open <share> (the end continues; <share> is
    the share of copies that show the commonest base over the whole
    continuation), poly-A tail or oligo-G tail.

Usage: python tools/cand_ends.py ALN.fa [ALN.fa ...] [--out TSV]

Arguments may contain wildcards, for example
  python tools/cand_ends.py composites/v2/cand/*P9*.aln.fa --out ends.tsv
"""

import argparse
import glob
import os
import sys
from collections import Counter

GAPS = {"-", "."}
TOP_BASE_MIN = 0.60         # a column outside the end counts as part of
                            # the continuation until its top-base share
                            # falls under this
CARRY_STEPS = (10, 30, 60)  # columns beyond the end for the carrying shares
WORD_LEN = 12               # length of the commonest outside word
RUN_MIN = 8                 # A / oligo-G run length
RUN_SHARE = 0.40            # share of copies for A / oligo-G dominance


def read_fasta(path):
    """Read aligned FASTA, return a list of (header, sequence) in file order."""
    records = []
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def candidate_name(path):
    """File name without directories and without a FASTA suffix."""
    name = os.path.basename(path)
    for suffix in (".aln.fa", ".aln.fasta", ".fasta", ".fa"):
        if name.lower().endswith(suffix) and len(name) > len(suffix):
            return name[: -len(suffix)]
    return name


def walk_outside(end_col, side, width):
    """Columns outside the end, walking away from the element."""
    if side == 5:
        return range(end_col - 1, -1, -1)
    return range(end_col + 1, width)


def outside_width(end_col, side, width):
    """Number of columns between the end and the alignment border."""
    return end_col if side == 5 else width - 1 - end_col


def top_base(copies, col):
    """Commonest base at a column and its share over all copies."""
    counts = Counter(seq[col].upper() for seq in copies if seq[col] not in GAPS)
    if not counts:
        return None, 0.0
    base, count = counts.most_common(1)[0]
    return base, count / len(copies)


def copy_reaches(copies, end_col, side):
    """How many columns beyond the end each copy still carries bases."""
    reaches = []
    for seq in copies:
        cols = [col for col in walk_outside(end_col, side, len(seq))
                if seq[col] not in GAPS]
        if not cols:
            reaches.append(0)
        elif side == 5:
            reaches.append(end_col - cols[-1])
        else:
            reaches.append(cols[-1] - end_col)
    return reaches


def first_word(seq, end_col, side, length):
    """The copy's first `length` bases outside the end, upper cased."""
    letters = []
    for col in walk_outside(end_col, side, len(seq)):
        if seq[col] not in GAPS:
            letters.append(seq[col].upper())
            if len(letters) == length:
                return "".join(letters)
    return None


def run_share(copies, end_col, side, base):
    """Share of copies with a run of >= RUN_MIN `base` outside the end."""
    hits = 0
    for seq in copies:
        run = 0
        for col in walk_outside(end_col, side, len(seq)):
            if seq[col] not in GAPS and seq[col].upper() == base:
                run += 1
                if run >= RUN_MIN:
                    hits += 1
                    break
            else:
                run = 0
    return hits / len(copies)


def end_evidence(copies, end_col, side, width):
    """Collect the end evidence for one end (side is 5 or 3)."""
    evidence = {
        "side": side,
        "consensus_end_col": end_col + 1,  # 1-based column
        "outside_cols": outside_width(end_col, side, width),
    }

    # First column outside the end whose top-base share is under TOP_BASE_MIN.
    first_col = None
    walked = 0
    for col in walk_outside(end_col, side, width):
        walked += 1
        if top_base(copies, col)[1] < TOP_BASE_MIN:
            first_col = col
            break
    if first_col is None:
        evidence["first_under_col"] = None
        evidence["first_under_beyond"] = None
        evidence["continues_cols"] = evidence["outside_cols"]
    else:
        evidence["first_under_col"] = first_col + 1  # 1-based column
        evidence["first_under_beyond"] = walked      # columns beyond the end
        evidence["continues_cols"] = walked - 1     # columns >= TOP_BASE_MIN

    # Shares of copies that carry bases beyond the end.
    reaches = copy_reaches(copies, end_col, side)
    evidence["carrying"] = {
        step: sum(reach >= step for reach in reaches) / len(copies)
        for step in CARRY_STEPS
    }

    # Commonest 12 bases outside the end.
    words = [w for w in (first_word(seq, end_col, side, WORD_LEN) for seq in copies) if w]
    if words:
        word, count = Counter(words).most_common(1)[0]
        evidence["word"] = word
        evidence["word_share"] = count / len(copies)
    else:
        evidence["word"] = None
        evidence["word_share"] = 0.0

    # A-run / oligo-G dominance.
    evidence["a_run_share"] = run_share(copies, end_col, side, "A")
    evidence["g_run_share"] = run_share(copies, end_col, side, "G")

    # Share of copies that continue the conserved sequence beyond the end.
    continues = evidence["continues_cols"]
    if continues > 0:
        tops = []
        for col in walk_outside(end_col, side, width):
            tops.append((col, top_base(copies, col)[0]))
            if len(tops) == continues:
                break
        matches = sum(
            all(seq[col] not in GAPS and seq[col].upper() == base
                for col, base in tops)
            for seq in copies
        )
        evidence["open_share"] = matches / len(copies)
    else:
        evidence["open_share"] = 0.0

    # One verdict per end.
    if evidence["a_run_share"] >= RUN_SHARE:
        evidence["verdict"] = "poly-A tail"
    elif evidence["g_run_share"] >= RUN_SHARE:
        evidence["verdict"] = "oligo-G tail"
    elif continues > 0:
        evidence["verdict"] = f"open {evidence['open_share']:.2f}"
    else:
        evidence["verdict"] = "closed"
    return evidence


def analyse(path):
    """Analyse one candidate alignment and return its report."""
    records = read_fasta(path)
    if len(records) < 2:
        raise SystemExit(f"{path}: expected a consensus row and at least one copy")
    width = max(len(seq) for _, seq in records)
    rows = [seq.ljust(width, "-") for _, seq in records]
    consensus, copies = rows[0], rows[1:]
    cols = [col for col, ch in enumerate(consensus) if ch not in GAPS]
    if not cols:
        raise SystemExit(f"{path}: the consensus row has no bases")
    return {
        "name": candidate_name(path),
        "path": path,
        "consensus_header": records[0][0],
        "copies": len(copies),
        "width": width,
        "ends": [
            end_evidence(copies, cols[0], 5, width),
            end_evidence(copies, cols[-1], 3, width),
        ],
    }


def fmt_share(value):
    """Format a share with two decimals, e.g. 0.75."""
    return f"{value:.2f}"


def print_report(report):
    """Print the human-readable report for one candidate."""
    print(f"candidate: {report['name']}")
    print(f"file: {report['path']}")
    print(f"consensus: {report['consensus_header']}")
    print(f"copies: {report['copies']}, alignment columns: {report['width']}")
    for evidence in report["ends"]:
        print(f"{evidence['side']}' end:")
        print(f"  consensus end column: {evidence['consensus_end_col']}")
        if evidence["first_under_col"] is None:
            if evidence["outside_cols"]:
                reason = f"all {evidence['outside_cols']} outside columns >= 60 %"
            else:
                reason = "no columns outside the end"
            print(f"  first column outside the end with top-base share under 60 %: none ({reason})")
        else:
            print(f"  first column outside the end with top-base share under 60 %: "
                  f"{evidence['first_under_col']} "
                  f"({evidence['first_under_beyond']} columns beyond the end)")
        print(f"  conserved continuation beyond the end: {evidence['continues_cols']} columns")
        for step in CARRY_STEPS:
            print(f"  carrying {step} columns beyond: {fmt_share(evidence['carrying'][step])}")
        if evidence["word"] is None:
            print(f"  commonest {WORD_LEN} bases outside the end: none")
        else:
            print(f"  commonest {WORD_LEN} bases outside the end: {evidence['word']} "
                  f"(share {fmt_share(evidence['word_share'])})")
        a_run = evidence["a_run_share"]
        g_run = evidence["g_run_share"]
        print(f"  A run >= {RUN_MIN} in {fmt_share(a_run)} of the copies "
              f"(>= {fmt_share(RUN_SHARE)}): {'yes' if a_run >= RUN_SHARE else 'no'}")
        print(f"  oligo-G run >= {RUN_MIN} in {fmt_share(g_run)} of the copies "
              f"(>= {fmt_share(RUN_SHARE)}): {'yes' if g_run >= RUN_SHARE else 'no'}")
        print(f"  verdict: {evidence['verdict']}")
    print()


TSV_COLUMNS = (
    "candidate",
    "end",
    "consensus_end_col",
    "first_under60_col",
    "first_under60_beyond",
    "continues_cols",
) + tuple(f"carrying_{step}" for step in CARRY_STEPS) + (
    "word12",
    "word12_share",
    "a_run_share",
    "oligo_g_run_share",
    "verdict",
)


def tsv_row(report, evidence):
    """One TSV row for one end."""
    first_col = evidence["first_under_col"]
    beyond = evidence["first_under_beyond"]
    row = [
        report["name"],
        str(evidence["side"]),
        str(evidence["consensus_end_col"]),
        "" if first_col is None else str(first_col),
        "" if beyond is None else str(beyond),
        str(evidence["continues_cols"]),
    ]
    row += [fmt_share(evidence["carrying"][step]) for step in CARRY_STEPS]
    row += [
        evidence["word"] if evidence["word"] is not None else "none",
        fmt_share(evidence["word_share"]),
        fmt_share(evidence["a_run_share"]),
        fmt_share(evidence["g_run_share"]),
        evidence["verdict"],
    ]
    return row


def expand_args(patterns):
    """Expand wildcard arguments (cmd.exe does not do it for us)."""
    paths = []
    for pattern in patterns:
        if "*" in pattern or "?" in pattern:
            matches = sorted(glob.glob(pattern))
            if not matches:
                raise SystemExit(f"no file matches {pattern}")
            paths.extend(matches)
        else:
            paths.append(pattern)
    return paths


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="End evidence for the open composite candidates.")
    parser.add_argument("aln", nargs="+", metavar="ALN.fa",
                        help="candidate alignment, row 1 = consensus, "
                             "other rows = copies")
    parser.add_argument("--out", metavar="TSV",
                        help="also write the evidence as a TSV table")
    args = parser.parse_args(argv)

    rows = []
    for path in expand_args(args.aln):
        report = analyse(path)
        print_report(report)
        for evidence in report["ends"]:
            rows.append(tsv_row(report, evidence))

    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(TSV_COLUMNS) + "\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
