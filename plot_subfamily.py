#!/usr/bin/env python3
"""
plot_subfamily.py — Generate per-subfamily diagnostic plots for SINEderella.

Plot 1: Histogram of % divergence from consensus (100 − %identity).
Plot 2: Per-position nucleotide frequency stacked bar chart.

Usage:
  python3 plot_subfamily.py \
    --subfamily <name> \
    --pctid <pctid.tsv> \
    --msa <msa.fasta> \
    --outdir <dir> \
    --consensus-name <name>
"""

import argparse
import sys
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_pctid(path):
    """Read tab-separated seqID  pctid file."""
    vals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                vals.append(float(parts[1]))
    return np.array(vals)


def read_msa(path):
    """Read FASTA alignment, return list of (name, seq) tuples."""
    records = []
    name, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                # strip _R_ prefix added by mafft --adjustdirection
                if name.startswith("_R_"):
                    name = name[3:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records


# ---------------------------------------------------------------------------
# Plot 1: Divergence histogram
# ---------------------------------------------------------------------------

def plot_divergence(pctid_values, subfamily, outdir):
    """Histogram of % divergence = 100 − %identity."""
    divergence = 100.0 - pctid_values

    fig, ax = plt.subplots(figsize=(8, 5))

    # Determine bin range
    lo = max(0, np.floor(divergence.min()) - 1)
    hi = min(100, np.ceil(divergence.max()) + 1)
    nbins = min(50, max(10, int((hi - lo) / 0.5)))

    ax.hist(divergence, bins=nbins, range=(lo, hi),
            color="#4C72B0", edgecolor="white", linewidth=0.5, alpha=0.85)

    mean_div = np.mean(divergence)
    median_div = np.median(divergence)
    ax.axvline(mean_div, color="#C44E52", linestyle="--", linewidth=1.5,
               label=f"Mean: {mean_div:.1f}%")
    ax.axvline(median_div, color="#DD8452", linestyle=":", linewidth=1.5,
               label=f"Median: {median_div:.1f}%")

    ax.set_xlabel("Divergence from consensus (%)", fontsize=12)
    ax.set_ylabel("Number of copies", fontsize=12)
    ax.set_title(f"{subfamily}  —  Divergence Distribution  (n={len(divergence)})",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{subfamily}_divergence.{ext}"),
                    dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2: Per-position nucleotide frequency stacked bar chart
# ---------------------------------------------------------------------------

NUC_COLORS = {
    "A": "#2ca02c",   # green
    "T": "#d62728",   # red
    "C": "#1f77b4",   # blue
    "G": "#ff7f0e",   # orange
    "gap": "#cccccc", # grey
}

NUCS = ["A", "T", "C", "G", "gap"]


def compute_nuc_freqs(msa_records, consensus_name):
    """Compute per-position nucleotide fractions from MSA.

    Returns (positions, freqs_dict) where freqs_dict[nuc] is an array
    of fractions at each consensus-aligned position.
    The consensus row is used to define positions (non-gap columns).
    """
    # Find consensus sequence
    cons_seq = None
    copy_seqs = []
    for name, seq in msa_records:
        if name == consensus_name:
            cons_seq = seq
        else:
            copy_seqs.append(seq)

    if cons_seq is None:
        # Fallback: first sequence is consensus
        cons_seq = msa_records[0][1]
        copy_seqs = [s for _, s in msa_records[1:]]

    if not copy_seqs:
        return None, None

    aln_len = len(cons_seq)

    # Find consensus non-gap positions
    cons_positions = [i for i in range(aln_len) if cons_seq[i] not in "-. "]
    n_pos = len(cons_positions)

    if n_pos == 0:
        return None, None

    # Count nucleotides at each consensus position across copies
    counts = {nuc: np.zeros(n_pos) for nuc in NUCS}
    n_copies = len(copy_seqs)

    for seq in copy_seqs:
        for pi, ai in enumerate(cons_positions):
            if ai < len(seq):
                base = seq[ai]
            else:
                base = "-"
            if base in "ATCG":
                counts[base][pi] += 1
            else:
                counts["gap"][pi] += 1

    # Convert to fractions
    freqs = {}
    for nuc in NUCS:
        freqs[nuc] = counts[nuc] / max(n_copies, 1)

    # Consensus bases for annotation
    cons_bases = [cons_seq[i] for i in cons_positions]

    return cons_bases, freqs, n_copies


def plot_nucfreq(msa_records, consensus_name, subfamily, outdir):
    """Per-position stacked bar chart of nucleotide frequencies."""
    result = compute_nuc_freqs(msa_records, consensus_name)
    if result is None or result[0] is None:
        return

    cons_bases, freqs, n_copies = result
    n_pos = len(cons_bases)

    # For very long consensuses, use a wider figure
    fig_width = max(10, min(60, n_pos * 0.12))
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x = np.arange(n_pos)
    bar_width = 0.9

    bottom = np.zeros(n_pos)
    for nuc in NUCS:
        label = nuc if nuc != "gap" else "Gap/N"
        ax.bar(x, freqs[nuc], bar_width, bottom=bottom,
               color=NUC_COLORS[nuc], label=label, linewidth=0)
        bottom += freqs[nuc]

    # X-axis: show consensus nucleotide at each position
    tick_step = 1 if n_pos <= 80 else (5 if n_pos <= 400 else 10)
    tick_positions = list(range(0, n_pos, tick_step))
    tick_labels = [cons_bases[i] if i < n_pos else "" for i in tick_positions]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=max(4, min(8, 500 // n_pos)),
                       fontfamily="monospace")

    # Add position numbers on top
    if n_pos <= 400:
        pos_step = 10 if n_pos <= 200 else 20
        for i in range(0, n_pos, pos_step):
            ax.annotate(str(i + 1), xy=(i, 1.01), fontsize=6,
                        ha="center", va="bottom", color="gray")

    ax.set_xlim(-0.5, n_pos - 0.5)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Consensus position (bp)", fontsize=11)
    ax.set_ylabel("Nucleotide frequency", fontsize=11)
    ax.set_title(
        f"{subfamily}  —  Nucleotide Composition  "
        f"(n={n_copies}, L={n_pos} bp)",
        fontsize=13, fontweight="bold")

    # Legend outside plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper right", fontsize=9,
              ncol=5, framealpha=0.8)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{subfamily}_nucfreq.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate SINEderella per-subfamily diagnostic plots")
    parser.add_argument("--subfamily", required=True, help="Subfamily name")
    parser.add_argument("--pctid", required=True,
                        help="Tab-separated file: seqID <tab> %%identity")
    parser.add_argument("--msa", required=True,
                        help="MAFFT MSA FASTA (consensus + copies)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--consensus-name", required=True,
                        help="Name of consensus sequence in MSA")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Plot 1: Divergence histogram
    pctid_vals = read_pctid(args.pctid)
    if len(pctid_vals) > 0:
        plot_divergence(pctid_vals, args.subfamily, args.outdir)
        print(f"  [OK] {args.subfamily}_divergence.png/pdf "
              f"({len(pctid_vals)} copies)")
    else:
        print(f"  [SKIP] {args.subfamily}: no %identity data", file=sys.stderr)

    # Plot 2: Nucleotide frequency
    msa_records = read_msa(args.msa)
    if len(msa_records) >= 2:
        plot_nucfreq(msa_records, args.consensus_name,
                     args.subfamily, args.outdir)
        print(f"  [OK] {args.subfamily}_nucfreq.png/pdf "
              f"({len(msa_records)-1} copies)")
    else:
        print(f"  [SKIP] {args.subfamily}: MSA has <2 sequences",
              file=sys.stderr)


if __name__ == "__main__":
    main()
