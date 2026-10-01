#!/usr/bin/env python3
"""CpG-adjusted divergence across all species: summarize every plate.

Finds every plate at

    <root>/*/alignments/*top100.aln.fa        (default root: C:/work/Tal)

(one species per folder) and writes summary.tsv next to this script, one row
per plate:

    species              species folder name
    family               plate file name without the "top100.aln.fa" suffix
    n_copies             copies on the plate
    consensus_cpg_share  consensus CpG sites / consensus element columns
    median_raw_p         median raw_p over the plate's copies
    median_cpg_adj_p     median cpg_adj_p over the plate's copies
    spearman_raw_vs_adj  Spearman rank correlation of raw_p vs cpg_adj_p
                         over the plate's copies
    bin_change_share     share of copies that fall in a different bin when
                         the copies are split into 5 equal-size divergence
                         bins by rank

Copies with no aligned element columns (NA in cpg_div.py) still count in
n_copies but are left out of the medians, the correlation and the bin share.

Usage:
    run_all.py [ROOT] [--out TSV]
"""

import argparse
import glob
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cpg_div

DEFAULT_ROOT = "C:/work/Tal"
N_BINS = 5
PLATE_SUFFIX = "top100.aln.fa"

HEADER = ("species", "family", "n_copies", "consensus_cpg_share",
          "median_raw_p", "median_cpg_adj_p", "spearman_raw_vs_adj",
          "bin_change_share")


def average_ranks(values):
    """1-based ranks; tied values get the average of their ranks."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # average of the 1-based ranks i+1 .. j+1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    """Spearman rank correlation; None if undefined (n < 2 or no variance)."""
    n = len(xs)
    if n < 2:
        return None
    ranks_x = average_ranks(xs)
    ranks_y = average_ranks(ys)
    mean_x = sum(ranks_x) / n
    mean_y = sum(ranks_y) / n
    cov = sum((a - mean_x) * (b - mean_y) for a, b in zip(ranks_x, ranks_y))
    var_x = sum((a - mean_x) ** 2 for a in ranks_x)
    var_y = sum((b - mean_y) ** 2 for b in ranks_y)
    if var_x == 0.0 or var_y == 0.0:
        return None
    return cov / math.sqrt(var_x * var_y)


def rank_bin(rank, n):
    """Bin 0..N_BINS-1 of a 1-based (average) rank; N_BINS equal-size bins."""
    return min(N_BINS - 1, int((rank - 1.0) * N_BINS / n))


def bin_change_share(xs, ys):
    """Share of values that land in a different equal-size rank bin."""
    n = len(xs)
    if n == 0:
        return None
    changed = sum(1 for a, b in zip(average_ranks(xs), average_ranks(ys))
                  if rank_bin(a, n) != rank_bin(b, n))
    return changed / n


def species_of(path):
    """.../SPECIES/alignments/plate.top100.aln.fa -> SPECIES"""
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def family_of(path):
    """Family name: plate file name without the top100.aln.fa suffix."""
    name = os.path.basename(path)
    if name.endswith(PLATE_SUFFIX):
        name = name[: -len(PLATE_SUFFIX)]
    return name.rstrip("._-") or os.path.basename(path)


def summarize_plate(path):
    """One summary row (dict) for one plate."""
    records = list(cpg_div.parse_fasta(path))
    _, copies, n_elem, n_cpg_elem = cpg_div.compute_plate(records)
    raws = []
    adjs = []
    for _, stats in copies:
        if stats["raw_p"] is None or stats["cpg_adj_p"] is None:
            continue  # no aligned columns: keep it out of the summaries
        raws.append(stats["raw_p"])
        adjs.append(stats["cpg_adj_p"])
    return {
        "species": species_of(path),
        "family": family_of(path),
        "n_copies": len(copies),
        "consensus_cpg_share": (n_cpg_elem / n_elem) if n_elem else None,
        "median_raw_p": statistics.median(raws) if raws else None,
        "median_cpg_adj_p": statistics.median(adjs) if adjs else None,
        "spearman_raw_vs_adj": spearman(raws, adjs),
        "bin_change_share": bin_change_share(raws, adjs),
    }


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="CpG-adjusted divergence summary over every plate")
    parser.add_argument("root", nargs="?", default=DEFAULT_ROOT,
                        help="species root (default: %(default)s)")
    parser.add_argument("--out", default=os.path.join(here, "summary.tsv"),
                        help="output TSV (default: %(default)s)")
    args = parser.parse_args(argv)

    pattern = os.path.join(args.root, "*", "alignments", "*" + PLATE_SUFFIX)
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f"run_all: no plates found: {pattern}", file=sys.stderr)
        return 1

    rows = []
    skipped = 0
    for path in paths:
        try:
            row = summarize_plate(path)
        except Exception as exc:  # one bad plate must not stop the survey
            skipped += 1
            print(f"run_all: WARNING: skipping {path}: {exc}", file=sys.stderr)
            continue
        rows.append(row)
        print(f"run_all: {row['species']} / {row['family']}: "
              f"{row['n_copies']} copies", file=sys.stderr)

    lines = ["\t".join(HEADER)]
    for row in rows:
        lines.append("\t".join([
            row["species"],
            row["family"],
            str(row["n_copies"]),
            cpg_div.fmt(row["consensus_cpg_share"]),
            cpg_div.fmt(row["median_raw_p"]),
            cpg_div.fmt(row["median_cpg_adj_p"]),
            cpg_div.fmt(row["spearman_raw_vs_adj"]),
            cpg_div.fmt(row["bin_change_share"]),
        ]))
    text = "\n".join(lines) + "\n"

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(text)
    sys.stdout.write(text)
    print(f"run_all: wrote {args.out} ({len(rows)} plates, {skipped} skipped)",
          file=sys.stderr)
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
