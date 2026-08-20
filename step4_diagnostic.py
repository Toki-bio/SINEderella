#!/usr/bin/env python3
"""
step4_diagnostic.py â€” Subfamily diagnostic analysis for SINEderella.

Integrates into SINEderella's post-processing pipeline. Reads step2 output,
computes position-level diagnostic weights, PCA validation, per-copy diagnostic
state, volcano plots, and integrates with SINEplot.

Usage:
    python3 step4_diagnostic.py RUN_ROOT [--threads N] [--top-k 20]

Inputs (auto-detected from RUN_ROOT):
    step2/step2_output*/all_votes.tsv
    step2/step2_output*/assignment_full.tsv
    step2/step2_output*/subfamilies/*.fasta
    consensuses.clean.fa

Outputs (written to step2/step2_output*/diagnostic/):
    pca_subfamilies.png/pdf         â€” PCA of bitscore matrix
    position_weights.tsv             â€” per-position MI, RF, KL weights
    diagnostic_positions.tsv         â€” selected top-K diagnostic positions
    synergy_pairs.tsv                â€” synergistic position pairs
    copy_diagnostic_state.tsv        â€” per-copy nucleotide at diagnostic positions
    diagnostic_scores.tsv            â€” per-copy diagnostic scores vs each subfamily
    diagnostic_flags.tsv             â€” discordance flags (fp_risk, fn_risk, chimera)
    diagnostic_heatmap.png/pdf       â€” copies Ã— positions heatmap
    diagnostic_vs_bitscore.png/pdf   â€” scatter with discordance quadrants
    volcano_{A}_vs_{B}.png/pdf       â€” per-pair volcano plots
    sineplot_input.txt               â€” all-vs-all ssearch36 output for SINEplot
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
from collections import defaultdict
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from scipy.spatial.distance import jensenshannon
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings("ignore")

# â”€â”€ Colour palettes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
NUC_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728", "-": "#cccccc"}
SUBFAM_COLORS = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
    "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494"
]

# â”€â”€ I/O helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def log(msg):
    print(f"[step4_diagnostic] {msg}", file=sys.stderr)

def run_cmd(cmd, shell=True):
    """Run command, return stdout."""
    try:
        return subprocess.check_output(cmd, shell=True, text=True,
                                       stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        log(f"WARNING: command failed: {e}")
        return ""

def read_fasta(path):
    """Read FASTA, return dict {name: seq}."""
    records = {}
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name:
        records[name] = "".join(seq).upper()
    return records

def read_fasta_list(path):
    """Read FASTA, return list of (name, seq)."""
    records = []
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(seq).upper()))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name:
        records.append((name, "".join(seq).upper()))
    return records

def load_tsv(path):
    """Load TSV, skip comments."""
    return pd.read_csv(path, sep="\t", comment="#")

# â”€â”€ Step 1: Extract bitscore matrix â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def build_bitscore_matrix(assignment_path, consensus_path, subfamily_fastas,
                           run_root, threads):
    """
    Build (copies Ã— subfamilies) bitscore matrix.

    Strategy: run ssearch36 with ALL consensus sequences as queries
    against ALL assigned copies as subjects (single M-query Ã— N-subject run).
    Parse bitscores to build the complete matrix.

    Falls back to using assignment_full.tsv if ssearch36 fails.
    """
    log("Building bitscore matrix ...")

    # Collect all assigned copy sequences
    all_copies_fa = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False).name
    n_copies = 0
    with open(all_copies_fa, "w") as out:
        for sf, fasta_path in sorted(subfamily_fastas.items()):
            if os.path.exists(fasta_path):
                recs = read_fasta_list(fasta_path)
                for name, seq in recs:
                    out.write(f">{name}\n{seq}\n")
                    n_copies += 1

    if n_copies == 0:
        log("  WARNING: No assigned copies found, using assignment data only")
        os.unlink(all_copies_fa)
        return _fallback_matrix(assignment_path)

    # Run ssearch36: consensus queries vs copy subjects
    ssearch_out = tempfile.NamedTemporaryFile(
        mode="w", suffix=".m8", delete=False).name
    cmd = (f"ssearch36 -Q -n -z 11 -E 100 -T {threads} -m 8 "
           f"{consensus_path} {all_copies_fa} > {ssearch_out}")
    log(f"  Running ssearch36 ({n_copies} copies Ã— query consensuses) ...")
    ret = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if ret.returncode != 0 or not os.path.exists(ssearch_out) or os.path.getsize(ssearch_out) == 0:
        log(f"  WARNING: ssearch36 failed, falling back to assignment data")
        os.unlink(all_copies_fa)
        if os.path.exists(ssearch_out):
            os.unlink(ssearch_out)
        return _fallback_matrix(assignment_path)

    # Parse: query=$1, subject=$2, bitscore=$12
    scores = defaultdict(dict)
    subfam_set = set()
    copy_set = set()
    with open(ssearch_out) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                query = parts[0]
                subject = parts[1]
                bs = float(parts[11])
                if subject not in scores or query not in scores[subject] or bs > scores[subject][query]:
                    scores[subject][query] = bs
                subfam_set.add(query)
                copy_set.add(subject)
            except (ValueError, IndexError):
                continue

    os.unlink(all_copies_fa)
    os.unlink(ssearch_out)

    if not scores:
        return _fallback_matrix(assignment_path)

    # Build matrix DataFrame
    subfamilies = sorted(subfam_set)
    rows = []
    for copy_id in sorted(copy_set):
        row = {"seqID": copy_id}
        for sf in subfamilies:
            row[sf] = scores[copy_id].get(sf, 0.0)
        rows.append(row)

    matrix = pd.DataFrame(rows).set_index("seqID")
    matrix = matrix.fillna(0)

    # Labels from assignment
    labels = _get_labels(assignment_path, matrix.index)

    log(f"  Matrix: {matrix.shape[0]} copies Ã— {matrix.shape[1]} subfamilies")
    assigned_n = sum(1 for v in labels.values() if v != "unassigned")
    log(f"  Assigned: {assigned_n}, Unassigned: {len(labels) - assigned_n}")
    return matrix, labels, subfamilies


def _fallback_matrix(assignment_path):
    """Fallback: build matrix from assignment_full.tsv only."""
    log("  Building from assignment data (single column per assigned subfamily)")
    assign = load_tsv(assignment_path)
    assign.columns = ["Sequence", "Subfamily", "Bitscore", "Votes", "Status", "Threshold"]
    assigned = assign[assign["Status"] == "assigned"]

    subfamilies = sorted(assigned["Subfamily"].unique())
    rows = []
    for _, row in assigned.iterrows():
        r = {"seqID": row["Sequence"]}
        for sf in subfamilies:
            r[sf] = row["Bitscore"] if row["Subfamily"] == sf else 0.0
        rows.append(r)

    matrix = pd.DataFrame(rows).set_index("seqID")
    labels = _get_labels(assignment_path, matrix.index)
    log(f"  Fallback matrix: {matrix.shape[0]} copies Ã— {matrix.shape[1]} subfamilies")
    return matrix, labels, subfamilies


def _get_labels(assignment_path, copy_ids):
    """Get subfamily labels for each copy from assignment."""
    assign = load_tsv(assignment_path)
    assign.columns = ["Sequence", "Subfamily", "Bitscore", "Votes", "Status", "Threshold"]
    assigned = assign[assign["Status"] == "assigned"].set_index("Sequence")
    labels = {}
    for cid in copy_ids:
        labels[cid] = assigned.loc[cid, "Subfamily"] if cid in assigned.index else "unassigned"
    return labels


# â”€â”€ Step 2: PCA â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def run_pca(matrix, labels, subfamilies, outdir):
    """PCA on bitscore matrix, produce scatter plot."""
    log("Running PCA on bitscore matrix ...")
    if len(subfamilies) < 2:
        log("  SKIP: need â‰¥2 subfamilies for PCA")
        return

    X = matrix.values
    pca = PCA(n_components=min(3, len(subfamilies)))
    projected = pca.fit_transform(X)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # â”€â”€ Plot 1: copies coloured by assigned subfamily â”€â”€
    ax = axes[0]
    unique_labels = sorted(set(labels.values()))
    for i, sf in enumerate(unique_labels):
        mask = [labels[idx] == sf for idx in matrix.index]
        color = SUBFAM_COLORS[i % len(SUBFAM_COLORS)] if sf != "unassigned" else "#cccccc"
        alpha = 0.4 if sf == "unassigned" else 0.6
        zorder = 1 if sf == "unassigned" else 2
        ax.scatter(projected[mask, 0], projected[mask, 1],
                   c=color, alpha=alpha, s=8, label=f"{sf} ({sum(mask)})",
                   zorder=zorder, edgecolors="none")

    # Loadings arrows
    for i, sf in enumerate(subfamilies):
        ax.arrow(0, 0, pca.components_[0, i] * 5, pca.components_[1, i] * 5,
                 head_width=0.15, head_length=0.15, fc="black", ec="black", alpha=0.5, lw=0.8)
        ax.text(pca.components_[0, i] * 5.3, pca.components_[1, i] * 5.3,
                sf, fontsize=8, alpha=0.7)

    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100
    ax.set_xlabel(f"PC1 ({var1:.1f}%)")
    ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    ax.set_title(f"PCA of SINE subfamily bitscores\n(n={matrix.shape[0]} copies, {len(subfamilies)} subfamilies)")
    ax.legend(fontsize=7, loc="upper right", ncol=2, markerscale=2)

    # â”€â”€ Plot 2: cumulative explained variance â”€â”€
    ax = axes[1]
    cumvar = np.cumsum(pca.explained_variance_ratio_) * 100
    ax.bar(range(1, len(cumvar) + 1), np.diff([0] + list(cumvar)),
           color="#4C72B0", alpha=0.7)
    ax.plot(range(1, len(cumvar) + 1), cumvar, "ro-", markersize=5)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)")
    ax.set_title("PCA explained variance")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, 105)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"pca_subfamilies.{ext}"), dpi=200)
    plt.close(fig)
    log(f"  PCA plots saved (PC1={var1:.1f}%, PC2={var2:.1f}%)")

    # Save PC coordinates
    pc_df = pd.DataFrame({
        "seqID": matrix.index,
        "PC1": projected[:, 0],
        "PC2": projected[:, 1],
        "assigned_subfamily": [labels[idx] for idx in matrix.index]
    })
    if projected.shape[1] >= 3:
        pc_df["PC3"] = projected[:, 2]
    pc_df.to_csv(os.path.join(outdir, "pca_coordinates.tsv"), sep="\t", index=False)

    return pca


# â”€â”€ Step 3: Diagnostic position weights â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def compute_position_weights(consensus_seqs, outdir, top_k=30):
    """
    Compute per-position diagnostic weights using three methods:
    1. Mutual Information between position and subfamily label
    2. Random Forest feature importance (captures interactions)
    3. KL divergence (pairwise discrimination)
    """
    log("Computing diagnostic position weights ...")
    if len(consensus_seqs) < 2:
        log("  SKIP: need â‰¥2 subfamilies with assigned copies")
        return None, None

    # â”€â”€ Build alignment of all copies â”€â”€
    # We align each subfamily's copies to the consensus bank
    subfamilies = sorted(consensus_seqs.keys())
    subfam_list = sorted(subfamily_assignments.keys())

    # First, find the reference length (longest consensus)
    ref_len = max(len(seq) for seq in consensus_seqs.values())
    ref_name = max(consensus_seqs, key=lambda k: len(consensus_seqs[k]))
    log(f"  Reference: {ref_name} ({ref_len} bp)")

    # Build per-position frequency tables per subfamily
    # For simplicity, we use the consensus sequences aligned to each other
    # via MAFFT, then treat the alignment columns as positions.
    # Alternative: align copies to reference using mafft --add.

    # â”€â”€ Simplified approach: use consensus alignment â”€â”€
    # Write consensuses to temp file, align, use alignment columns
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for sf in subfamilies:
            tmp.write(f">{sf}\n{consensus_seqs[sf]}\n")
        cons_file = tmp.name

    aln_file = cons_file + ".aln"
    run_cmd(f"mafft --quiet --auto {cons_file} > {aln_file}")
    if not os.path.exists(aln_file) or os.path.getsize(aln_file) == 0:
        log("  WARNING: MAFFT alignment failed, using raw consensus positions")
        os.unlink(cons_file)
        return None, None

    cons_aln = read_fasta(aln_file)
    os.unlink(cons_file)
    os.unlink(aln_file)

    # Extract alignment matrix: positions Ã— subfamilies
    aln_matrix = {}
    for sf in subfamilies:
        if sf in cons_aln:
            aln_matrix[sf] = cons_aln[sf]

    if not aln_matrix:
        return None, None

    # Find non-gap positions in reference (informative columns)
    ref_seq = aln_matrix[ref_name]
    valid_positions = [i for i, c in enumerate(ref_seq) if c not in "-."]
    n_pos = len(valid_positions)
    log(f"  Alignment: {n_pos} non-gap positions (of {len(ref_seq)} total columns)")

    # For each valid position, collect nucleotide counts per subfamily
    # Here we use the consensus itself as the representative.
    # For a full implementation, you'd align all copies and count per subfamily.
    # We use a hybrid: the consensus represents the "ideal" profile.

    # â”€â”€ MI computation â”€â”€
    mi_scores = np.zeros(n_pos)
    for pi, aln_col in enumerate(valid_positions):
        # Build contingency table: subfamilies Ã— nucleotides
        counts = defaultdict(lambda: defaultdict(int))
        for sf in subfamilies:
            if sf in aln_matrix and aln_col < len(aln_matrix[sf]):
                nuc = aln_matrix[sf][aln_col]
                if nuc in "ACGT":
                    counts[nuc][sf] += 1
                else:
                    counts["-"][sf] += 1

        # Flatten to contingency table
        nucs = sorted(set(k for d in counts.values() for k in d))
        if len(nucs) < 2:
            mi_scores[pi] = 0.0
            continue

        # Compute MI: I(Position; Subfamily) = H(Subfamily) + H(Position) - H(Position, Subfamily)
        # Using empirical frequencies.
        total = len(subfamilies)
        # Marginal entropy H(Subfamily): uniform by design
        h_subfam = np.log2(len(subfamilies))

        # H(Position)
        col_counts = {n: sum(counts[n].values()) for n in nucs}
        h_pos = -sum((c / total) * np.log2(c / total) for c in col_counts.values() if c > 0)

        # Joint entropy H(Position, Subfamily)
        h_joint = 0
        for n in nucs:
            for sf in subfamilies:
                c = counts[n].get(sf, 0)
                if c > 0:
                    h_joint -= (c / total) * np.log2(c / total)

        # MI = H(X) + H(Y) - H(X,Y)
        mi = h_subfam + h_pos - h_joint
        # Normalize by min(H_subfam, H_pos) to [0,1]
        norm = min(h_subfam, h_pos) if min(h_subfam, h_pos) > 0 else 1
        mi_scores[pi] = mi / norm if norm > 0 else 0

    # â”€â”€ KL divergence (max over all subfamily pairs) â”€â”€
    kl_scores = np.zeros(n_pos)
    for pi, aln_col in enumerate(valid_positions):
        max_kl = 0
        for sf_a, sf_b in combinations(subfamilies, 2):
            if sf_a not in aln_matrix or sf_b not in aln_matrix:
                continue
            nuc_a = aln_matrix[sf_a][aln_col] if aln_col < len(aln_matrix[sf_a]) else "-"
            nuc_b = aln_matrix[sf_b][aln_col] if aln_col < len(aln_matrix[sf_b]) else "-"
            if nuc_a != nuc_b:
                # Binary KL: 1 bit if different, scaled by conservation
                max_kl = max(max_kl, 1.0)
        kl_scores[pi] = max_kl

    # â”€â”€ RF importance â”€â”€
    rf_scores = np.zeros(n_pos)
    try:
        # One-hot encode alignment positions
        X_rf = []
        y_rf = []
        for sf in subfamilies:
            if sf in aln_matrix:
                row = []
                for aln_col in valid_positions:
                    nuc = aln_matrix[sf][aln_col] if aln_col < len(aln_matrix[sf]) else "-"
                    for base in "ACGT-":
                        row.append(1 if nuc == base else 0)
                X_rf.append(row)
                y_rf.append(sf)

        X_rf = np.array(X_rf)
        if len(set(y_rf)) >= 2 and X_rf.shape[0] >= 3:
            rf = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
            rf.fit(X_rf, y_rf)
            # Aggregate per-position from per-base importances
            imp = rf.feature_importances_
            for pi in range(n_pos):
                rf_scores[pi] = np.sum(imp[pi*5:(pi+1)*5])
    except Exception:
        log("  WARNING: RF importance computation failed, using zeros")

    # â”€â”€ Combined weight â”€â”€
    mi_norm = mi_scores / max(mi_scores.max(), 1e-10)
    kl_norm = kl_scores / max(kl_scores.max(), 1e-10)
    rf_norm = rf_scores / max(rf_scores.max(), 1e-10)

    combined = 0.4 * mi_norm + 0.3 * kl_norm + 0.3 * rf_norm

    # â”€â”€ Save position_weights.tsv â”€â”€
    rows = []
    for pi, aln_col in enumerate(valid_positions):
        ref_nuc = ref_seq[aln_col]
        rows.append({
            "position": pi + 1,
            "ref_nuc": ref_nuc,
            "MI": round(mi_scores[pi], 4),
            "MI_norm": round(mi_norm[pi], 4),
            "KL_max": round(kl_scores[pi], 4),
            "KL_norm": round(kl_norm[pi], 4),
            "RF_imp": round(rf_scores[pi], 6),
            "RF_norm": round(rf_norm[pi], 4),
            "combined_weight": round(combined[pi], 4)
        })
    wt_df = pd.DataFrame(rows)
    wt_df = wt_df.sort_values("combined_weight", ascending=False)
    wt_df.to_csv(os.path.join(outdir, "position_weights.tsv"), sep="\t", index=False)

    # Select top-K diagnostic positions
    top_positions = wt_df.head(top_k)["position"].tolist()
    diag_df = wt_df[wt_df["position"].isin(top_positions)].copy()
    diag_df.to_csv(os.path.join(outdir, "diagnostic_positions.tsv"), sep="\t", index=False)
    log(f"  Selected {len(top_positions)} diagnostic positions (top-{top_k})")

    # â”€â”€ Synergy detection â”€â”€
    log("  Computing pairwise synergies ...")
    synergy_rows = []
    n_top = min(30, len(wt_df))
    top_indices = [valid_positions[wt_df.iloc[i]["position"] - 1] for i in range(n_top)]
    # Use MI on consensus alone for synergy (full copy-based would be expensive)
    for i in range(min(n_top, 15)):
        for j in range(i + 1, min(n_top, 15)):
            # Quick heuristic: positions with different nucleotides across subfamilies
            # that co-vary are synergistic
            pi, pj = wt_df.iloc[i]["position"] - 1, wt_df.iloc[j]["position"] - 1
            ai, aj = valid_positions[pi], valid_positions[pj]

            # Collect (nuc_i, nuc_j, subfamily) triples
            patterns = defaultdict(int)
            for sf in subfamilies:
                if sf in aln_matrix:
                    ni = aln_matrix[sf][ai] if ai < len(aln_matrix[sf]) else "-"
                    nj = aln_matrix[sf][aj] if aj < len(aln_matrix[sf]) else "-"
                    patterns[(ni, nj)] += 1

            if len(patterns) == len(subfamilies):
                # Each subfamily has a unique (ni, nj) pair â€” synergistic!
                synergy = 1.0 - 1.0 / len(patterns)
            elif len(patterns) == 1:
                synergy = 0.0
            else:
                # Partial synergy: some pairs shared, some unique
                max_p = max(patterns.values())
                synergy = (len(patterns) - 1) / (len(subfamilies) - 1) * (1 - max_p / len(subfamilies))

            if synergy > 0.1:
                synergy_rows.append({
                    "pos1": pi + 1,
                    "pos2": pj + 1,
                    "pos1_weight": round(combined[pi], 4),
                    "pos2_weight": round(combined[pj], 4),
                    "synergy": round(synergy, 4)
                })

    if synergy_rows:
        syn_df = pd.DataFrame(synergy_rows).sort_values("synergy", ascending=False)
        syn_df.to_csv(os.path.join(outdir, "synergy_pairs.tsv"), sep="\t", index=False)
        log(f"  Found {len(syn_df)} synergistic pairs")
    else:
        log("  No synergistic pairs found")

    return wt_df, top_positions


# â”€â”€ Step 4: Per-copy diagnostic state â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def compute_copy_diagnostic_state(subfamily_fastas, consensus_seqs,
                                   top_positions, position_weights_df,
                                   outdir, threads):
    """
    For each copy, extract nucleotide at each diagnostic position.
    Align each subfamily's copies to the reference consensus.
    """
    log("Computing per-copy diagnostic state ...")
    if not top_positions:
        log("  SKIP: no diagnostic positions identified")
        return None, None

    # Build reference-indexed position map
    # We need to map "position numbers" (1-based, from consensus alignment)
    # back to actual alignment columns.
    ref_name = max(consensus_seqs, key=lambda k: len(consensus_seqs[k]))
    ref_seq = consensus_seqs[ref_name]

    # Write reference consensus
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp.write(f">{ref_name}\n{ref_seq}\n")
        ref_file = tmp.name

    all_states = {}  # seqID â†’ {pos: nucleotide}
    subfam_profiles = {}  # subfamily â†’ {pos: expected_nucleotide}

    for sf, fasta_path in sorted(subfamily_fastas.items()):
        if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
            continue

        records = read_fasta_list(fasta_path)
        if not records:
            continue
        copies = [r for r in records if r[0] != sf]
        if not copies:
            continue

        # Write copies to temp file
        copies_file = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False).name
        with open(copies_file, "w") as fh:
            for name, seq in copies:
                fh.write(f">{name}\n{seq}\n")

        # MAFFT --add to align copies to reference alignment
        combined = ref_file + ".combined.fa"
        with open(combined, "w") as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")
            for name, seq in copies:
                fh.write(f">{name}\n{seq}\n")

        aln_out = copies_file + ".aln"
        cmd = (f"mafft --thread {threads} --quiet --auto "
               f"--preservecase {combined} > {aln_out}")
        run_cmd(cmd)

        if os.path.exists(aln_out) and os.path.getsize(aln_out) > 0:
            aln_records = read_fasta(aln_out)
            # Find reference in alignment
            if ref_name in aln_records:
                aln_ref = aln_records[ref_name]
                # Map consensus 1-based positions â†’ alignment columns
                pos_to_col = {}
                rp = 0
                for ci, c in enumerate(aln_ref):
                    if c not in "-.":
                        rp += 1
                        pos_to_col[rp] = ci

                # Extract per-copy nucleotides at diagnostic positions
                for name, aln_seq in aln_records.items():
                    if name == ref_name:
                        continue
                    state = {}
                    for pos in top_positions:
                        if pos in pos_to_col:
                            col = pos_to_col[pos]
                            nuc = aln_seq[col] if col < len(aln_seq) else "-"
                            state[pos] = nuc if nuc in "ACGT" else "-"
                        else:
                            state[pos] = "-"
                    all_states[name] = state

                # Build subfamily profile from alignment
                profile = {}
                for pos in top_positions:
                    if pos in pos_to_col:
                        col = pos_to_col[pos]
                        counts = defaultdict(int)
                        for name, aln_seq in aln_records.items():
                            if name == ref_name:
                                continue
                            nuc = aln_seq[col] if col < len(aln_seq) else "-"
                            counts[nuc if nuc in "ACGT" else "-"] += 1
                        profile[pos] = max(counts, key=counts.get)
                    else:
                        profile[pos] = "-"
                subfam_profiles[sf] = profile

        # Cleanup
        for f in [copies_file, combined, aln_out]:
            if os.path.exists(f):
                os.unlink(f)

    os.unlink(ref_file)

    if not all_states:
        return None, None

    # Save copy_diagnostic_state.tsv
    all_pos = sorted(top_positions)
    rows = []
    for seqid in sorted(all_states):
        row = {"seqID": seqid}
        for pos in all_pos:
            row[f"pos_{pos}"] = all_states[seqid].get(pos, "-")
        rows.append(row)

    state_df = pd.DataFrame(rows)
    state_df.to_csv(os.path.join(outdir, "copy_diagnostic_state.tsv"), sep="\t", index=False)
    log(f"  Extracted diagnostic states for {len(all_states)} copies")

    return all_states, subfam_profiles


# â”€â”€ Step 5: Diagnostic scoring and discordance â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def compute_diagnostic_scores(all_states, subfam_profiles, position_weights_df,
                               assignment_path, matrix, labels, outdir):
    """Score each copy against each subfamily using diagnostic positions."""
    log("Computing diagnostic scores and discordance flags ...")
    if not all_states or not subfam_profiles:
        return

    wt_df = position_weights_df.set_index("position")
    all_pos = sorted(next(iter(all_states.values())).keys())

    # Compute diagnostic scores
    scores = {}  # (seqID, subfamily) â†’ score
    for seqid, state in all_states.items():
        for sf, profile in subfam_profiles.items():
            w_sum = 0
            w_match = 0
            for pos in all_pos:
                w = wt_df.loc[pos, "combined_weight"] if pos in wt_df.index else 0
                w_sum += w
                if state.get(pos, "-") == profile.get(pos, "-"):
                    w_match += w
            scores[(seqid, sf)] = w_match / w_sum if w_sum > 0 else 0

    # Save diagnostic_scores.tsv
    rows = []
    for (seqid, sf), score in sorted(scores.items()):
        rows.append({"seqID": seqid, "subfamily": sf, "diagnostic_score": round(score, 4)})
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "diagnostic_scores.tsv"), sep="\t", index=False)

    # â”€â”€ Discordance detection â”€â”€
    # Read assignment
    assign = load_tsv(assignment_path)
    assign.columns = ["Sequence", "Subfamily", "Bitscore", "Votes", "Status", "Threshold"]

    sim_path = os.path.join(os.path.dirname(assignment_path), "sim_scores.tsv")
    sim_scores = {}
    if os.path.exists(sim_path):
        sim = load_tsv(sim_path)
        sim_scores = dict(zip(sim.iloc[:, 0], sim.iloc[:, 3]))  # seqID â†’ sim_ratio

    flag_rows = []
    for seqid in sorted(all_states):
        assigned_sf = None
        if seqid in assign["Sequence"].values:
            row = assign[assign["Sequence"] == seqid].iloc[0]
            assigned_sf = row["Subfamily"]
            status = row["Status"]
        else:
            assigned_sf = "unassigned"
            status = "unassigned"

        # Best diagnostic subfamily
        best_sf = max(subfam_profiles.keys(),
                      key=lambda sf: scores.get((seqid, sf), 0))
        best_diag = scores.get((seqid, assigned_sf), 0) if assigned_sf != "unassigned" else 0

        # Find runner-up
        sf_scores = {sf: scores.get((seqid, sf), 0) for sf in subfam_profiles}
        sorted_sf = sorted(sf_scores.items(), key=lambda x: x[1], reverse=True)
        runner_up = sorted_sf[1][0] if len(sorted_sf) > 1 else None
        runner_score = sorted_sf[1][1] if len(sorted_sf) > 1 else 0

        # Overall similarity
        sim = sim_scores.get(seqid, None)

        flag = "concordant"
        detail = ""

        if status == "assigned" and assigned_sf:
            if best_diag < 0.5 and runner_score > best_diag:
                flag = "fp_risk"
                detail = f"overall_to_{assigned_sf}_but_diag_matches_{best_sf}"
            elif best_diag < 0.3:
                flag = "fp_risk"
                detail = "low_diagnostic_confidence"
        elif status != "assigned":
            if best_diag > 0.7:
                flag = "fn_risk"
                detail = f"diagnostic_matches_{best_sf}"

        # Chimera check: mixed pattern across positions
        if status == "assigned" and assigned_sf:
            # Check if 5' half matches one subfamily and 3' half matches another
            mid = len(all_pos) // 2
            first_half = {sf: 0 for sf in subfam_profiles}
            second_half = {sf: 0 for sf in subfam_profiles}
            for i, pos in enumerate(all_pos):
                nuc = all_states[seqid].get(pos, "-")
                for sf, prof in subfam_profiles.items():
                    if nuc == prof.get(pos, "-"):
                        if i < mid:
                            first_half[sf] += 1
                        else:
                            second_half[sf] += 1
            best_first = max(first_half, key=first_half.get)
            best_second = max(second_half, key=second_half.get)
            if (best_first != best_second and
                first_half[best_first] >= 0.6 * mid and
                second_half[best_second] >= 0.6 * (len(all_pos) - mid)):
                flag = "chimera_suspected"
                detail = f"5prime_{best_first}_3prime_{best_second}"

        flag_rows.append({
            "seqID": seqid,
            "assigned_subfamily": assigned_sf,
            "best_diagnostic_subfamily": best_sf,
            "diagnostic_score": round(best_diag, 4),
            "runner_up_subfamily": runner_up,
            "runner_up_score": round(runner_score, 4),
            "overall_sim_ratio": round(sim, 4) if sim is not None else ".",
            "flag": flag,
            "detail": detail
        })

    flag_df = pd.DataFrame(flag_rows)
    flag_df.to_csv(os.path.join(outdir, "diagnostic_flags.tsv"), sep="\t", index=False)

    n_fp = sum(1 for r in flag_rows if r["flag"] == "fp_risk")
    n_fn = sum(1 for r in flag_rows if r["flag"] == "fn_risk")
    n_ch = sum(1 for r in flag_rows if r["flag"] == "chimera_suspected")
    log(f"  Flagged: {n_fp} fp_risk, {n_fn} fn_risk, {n_ch} chimera_suspected")

    return flag_rows


# â”€â”€ Step 6: Diagnostic heatmap â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def plot_diagnostic_heatmap(all_states, subfam_profiles, position_weights_df,
                             labels, assignment_path, outdir, max_copies=200):
    """Plot copies Ã— diagnostic positions heatmap."""
    log("Plotting diagnostic heatmap ...")
    if not all_states or not subfam_profiles:
        return

    wt_df = position_weights_df.sort_values("combined_weight", ascending=False)
    top_positions = wt_df["position"].head(30).tolist()

    # Group copies by subfamily
    assign = load_tsv(assignment_path)
    assign.columns = ["Sequence", "Subfamily", "Bitscore", "Votes", "Status", "Threshold"]
    assigned_map = dict(zip(assign["Sequence"], assign["Subfamily"]))

    by_sf = defaultdict(list)
    for seqid in sorted(all_states):
        sf = assigned_map.get(seqid, "unassigned")
        by_sf[sf].append(seqid)

    # Sample if too many
    sampled = {}
    for sf, seqs in sorted(by_sf.items()):
        if len(seqs) > max_copies // max(len(by_sf), 1):
            np.random.seed(42)
            sampled[sf] = list(np.random.choice(seqs, max_copies // len(by_sf), replace=False))
        else:
            sampled[sf] = seqs

    # Build ordered list
    order = []
    for sf in sorted(sampled.keys()):
        order.extend([(seqid, sf) for seqid in sampled[sf]])

    # Build color matrix
    nuc_to_color = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
    cmap_colors = ["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728", "#eeeeee"]
    color_matrix = np.zeros((len(order), len(top_positions)))
    for ri, (seqid, sf) in enumerate(order):
        for ci, pos in enumerate(top_positions):
            nuc = all_states[seqid].get(pos, "-")
            color_matrix[ri, ci] = nuc_to_color.get(nuc, 4)

    fig, ax = plt.subplots(figsize=(max(12, len(top_positions) * 0.3),
                                    max(6, len(order) * 0.08)))
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(cmap_colors)
    ax.imshow(color_matrix, aspect="auto", cmap=cmap, interpolation="nearest")

    # Subfamily separator lines
    y_pos = 0
    for sf in sorted(sampled.keys()):
        n = len(sampled[sf])
        ax.axhline(y=y_pos + n - 0.5, color="black", linewidth=1.5)
        ax.text(-1, y_pos + n / 2 - 0.5, sf, fontsize=8, ha="right", va="center",
                fontweight="bold")
        y_pos += n

    # X-axis labels: position numbers
    ax.set_xticks(range(len(top_positions)))
    ax.set_xticklabels([str(p) for p in top_positions], fontsize=6,
                       rotation=90, fontfamily="monospace")
    ax.set_xlabel("Diagnostic position", fontsize=11)
    ax.set_ylabel(f"Copies (n={len(order)})", fontsize=11)
    ax.set_title("Diagnostic Position Heatmap\n(colour = nucleotide at position)",
                 fontsize=13, fontweight="bold")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=n) for n, c in
                       zip(["A", "C", "G", "T", "Gap/N"], cmap_colors)]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8, ncol=5)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"diagnostic_heatmap.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    log(f"  Heatmap saved ({len(order)} copies Ã— {len(top_positions)} positions)")


# â”€â”€ Step 7: Diagnostic vs bitscore scatter â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def plot_diagnostic_vs_bitscore(flag_rows, outdir):
    """Scatter plot of diagnostic_score vs sim_ratio with discordance quadrants."""
    log("Plotting diagnostic vs bitscore scatter ...")
    if not flag_rows:
        return

    df = pd.DataFrame(flag_rows)
    # Only points with valid sim_ratio
    plot_data = df[df["overall_sim_ratio"] != "."].copy()
    if len(plot_data) == 0:
        return
    plot_data["overall_sim_ratio"] = plot_data["overall_sim_ratio"].astype(float)

    fig, ax = plt.subplots(figsize=(8, 7))

    colors = {"concordant": "#2ecc71", "fp_risk": "#e74c3c",
              "fn_risk": "#f39c12", "chimera_suspected": "#9b59b6"}
    for flag in ["concordant", "fp_risk", "fn_risk", "chimera_suspected"]:
        subset = plot_data[plot_data["flag"] == flag]
        if len(subset) == 0:
            continue
        ax.scatter(subset["overall_sim_ratio"], subset["diagnostic_score"],
                   c=colors[flag], alpha=0.5, s=10, label=f"{flag} ({len(subset)})",
                   edgecolors="none")

    # Quadrant lines
    ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0.5, color="gray", linestyle="--", alpha=0.5)
    # Label quadrants
    ax.text(0.25, 0.85, "fp_risk", fontsize=9, color="#e74c3c", ha="center",
            fontweight="bold", alpha=0.6)
    ax.text(0.75, 0.15, "fn_risk", fontsize=9, color="#f39c12", ha="center",
            fontweight="bold", alpha=0.6)
    ax.text(0.25, 0.15, "Weak (both low)", fontsize=8, color="gray", ha="center", alpha=0.5)

    ax.set_xlabel("Overall similarity (sim_ratio)", fontsize=12)
    ax.set_ylabel("Diagnostic score", fontsize=12)
    ax.set_title("Diagnostic vs Overall Similarity\n(discordance detection)",
                 fontsize=13, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8, markerscale=2)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"diagnostic_vs_bitscore.{ext}"), dpi=200)
    plt.close(fig)
    log(f"  Scatter saved ({len(plot_data)} copies)")


# â”€â”€ Step 8: Volcano plots per subfamily pair â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def plot_volcano_plots(consensus_seqs, position_weights_df,
                        outdir, alpha=0.01):
    """Volcano plot for each subfamily pair."""
    log("Generating volcano plots ...")
    subfamilies = sorted(consensus_seqs.keys())
    if len(subfamilies) < 2:
        return

    # Align all consensuses
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for sf in subfamilies:
            tmp.write(f">{sf}\n{consensus_seqs[sf]}\n")
        cons_file = tmp.name

    aln_file = cons_file + ".aln"
    run_cmd(f"mafft --quiet --auto {cons_file} > {aln_file}")
    if not os.path.exists(aln_file) or os.path.getsize(aln_file) == 0:
        os.unlink(cons_file)
        return

    aln = read_fasta(aln_file)
    os.unlink(cons_file)
    os.unlink(aln_file)

    # Build per-position statistical test for each pair
    for sf_a, sf_b in combinations(subfamilies, 2):
        if sf_a not in aln or sf_b not in aln:
            continue

        seq_a = aln[sf_a]
        seq_b = aln[sf_b]
        positions = []
        pvals = []
        effects = []

        for i in range(len(seq_a)):
            nuc_a = seq_a[i] if seq_a[i] in "ACGT" else "-"
            nuc_b = seq_b[i] if seq_b[i] in "ACGT" else "-"
            if nuc_a == nuc_b:
                positions.append(i + 1)
                pvals.append(1.0)
                effects.append(0.0)
                continue

            # Fisher's exact test (simplified: binary same/different)
            # In full implementation: use aligned copy counts
            table = [[1, 0], [0, 1]]  # placeholder for single-sequence comparison
            try:
                _, pval = fisher_exact(table)
            except Exception:
                pval = 1.0
            positions.append(i + 1)
            pvals.append(pval)
            effects.append(1.0)  # different nucleotide = effect

        # Multiple testing correction (Bonferroni)
        n_tests = len(pvals)
        adj_pvals = [min(p * n_tests, 1.0) for p in pvals]

        # Volcano plot
        fig, ax = plt.subplots(figsize=(8, 6))

        neg_log_p = [-np.log10(max(p, 1e-300)) for p in adj_pvals]
        sig_threshold = -np.log10(alpha)

        colors = []
        for pos, eff, p in zip(positions, effects, adj_pvals):
            if p < alpha and eff > 0.5:
                colors.append("#e74c3c")  # significant + large effect
            elif p < alpha:
                colors.append("#f39c12")  # significant, small effect
            else:
                colors.append("#bdc3c7")  # not significant

        ax.scatter(effects, neg_log_p, c=colors, s=40, alpha=0.7, edgecolors="none")

        # Label top hits
        sig_positions = [(pos, eff, p) for pos, eff, p in
                         zip(positions, effects, adj_pvals) if p < alpha and eff > 0.5]
        sig_positions.sort(key=lambda x: x[2])
        for pos, eff, p in sig_positions[:10]:
            ax.annotate(str(pos), (eff, -np.log10(max(p, 1e-300))),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")

        ax.axhline(y=sig_threshold, color="gray", linestyle="--", alpha=0.5,
                   label=f"p_adj = {alpha}")
        ax.set_xlabel("Effect size (nucleotide difference)", fontsize=12)
        ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
        ax.set_title(f"Volcano Plot: {sf_a} vs {sf_b}\n"
                     f"(diagnostic positions for subfamily discrimination)",
                     fontsize=13, fontweight="bold")
        ax.legend(fontsize=8)

        plt.tight_layout()
        safe_name = f"{sf_a}_vs_{sf_b}".replace("/", "_").replace("\\", "_")
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(outdir, f"volcano_{safe_name}.{ext}"), dpi=200)
        plt.close(fig)

        n_sig = sum(1 for p in adj_pvals if p < alpha)
        log(f"  Volcano {sf_a} vs {sf_b}: {n_sig} significant positions")


# â”€â”€ Step 9: Generate SINEplot input â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def generate_sineplot_input(consensus_seqs, subfamily_fastas, outdir, threads):
    """Generate all-vs-all ssearch36 output for SINEplot."""
    log("Generating SINEplot input ...")

    # Combine consensuses + all copies into one FASTA
    combined = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False).name
    with open(combined, "w") as fh:
        for sf in sorted(consensus_seqs):
            fh.write(f">{sf}\n{consensus_seqs[sf]}\n")
        for sf, fasta_path in sorted(subfamily_fastas.items()):
            if os.path.exists(fasta_path):
                recs = read_fasta_list(fasta_path)
                for name, seq in recs:
                    fh.write(f">{name}\n{seq}\n")

    scores_out = os.path.join(outdir, "sineplot_input.txt")
    cmd = (f"ssearch36 -m 8 -T {threads} {combined} {combined} > {scores_out}")
    log(f"  Running ssearch36 all-vs-all (this may take a while) ...")
    ret = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if ret.returncode != 0:
        log(f"  WARNING: ssearch36 failed: {ret.stderr[:200]}")
        log(f"  You can generate SINEplot input manually with:")
        log(f"    cat consensuses.fa subfamilies/*.fasta > all.fa")
        log(f"    ssearch36 -m 8 -T 4 all.fa all.fa > {scores_out}")
    else:
        log(f"  SINEplot input: {scores_out}")

    os.unlink(combined)
    return scores_out if os.path.exists(scores_out) else None


# â”€â”€ Main â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

def main():
    parser = argparse.ArgumentParser(
        description="SINEderella step4: diagnostic analysis pipeline")
    parser.add_argument("run_root", help="Path to SINEderella run directory")
    parser.add_argument("--threads", type=int, default=None,
                        help="CPU threads (default: auto-detect)")
    parser.add_argument("--top-k", type=int, default=20,
                        help="Number of diagnostic positions to select (default: 20)")
    parser.add_argument("--skip-sineplot", action="store_true",
                        help="Skip SINEplot input generation (slow)")
    args = parser.parse_args()

    run_root = os.path.abspath(args.run_root)
    threads = args.threads or int(subprocess.check_output(
        "nproc 2>/dev/null || echo 4", shell=True, text=True).strip())

    # â”€â”€ Locate step2 output â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    step2_dirs = sorted([
        d for d in os.listdir(os.path.join(run_root, "step2"))
        if d.startswith("step2_output")
    ], reverse=True)
    if not step2_dirs:
        log(f"ERROR: No step2_output found in {run_root}/step2/")
        sys.exit(1)
    step2_out = os.path.join(run_root, "step2", step2_dirs[0])
    log(f"Using step2 output: {step2_out}")

    # Create diagnostic output directory
    outdir = os.path.join(step2_out, "diagnostic")
    os.makedirs(outdir, exist_ok=True)

    # Required inputs
    assign_path = os.path.join(step2_out, "assignment_full.tsv")
    sim_path = os.path.join(step2_out, "sim_scores.tsv")
    subfam_dir = os.path.join(step2_out, "subfamilies")
    cons_path = os.path.join(run_root, "consensuses.clean.fa")

    for p, name in [(assign_path, "assignment_full.tsv"),
                     (cons_path, "consensuses.clean.fa")]:
        if not os.path.exists(p):
            log(f"ERROR: Missing required file: {p}")
            sys.exit(1)

    # Read consensus sequences
    consensus_seqs = read_fasta(cons_path)
    log(f"Loaded {len(consensus_seqs)} consensus sequences")

    # Read subfamily FASTA files
    subfamily_fastas = {}
    if os.path.isdir(subfam_dir):
        for fname in sorted(os.listdir(subfam_dir)):
            if fname.endswith(".fasta") or fname.endswith(".fa"):
                sf_name = fname.rsplit(".", 1)[0]
                subfamily_fastas[sf_name] = os.path.join(subfam_dir, fname)
    log(f"Found {len(subfamily_fastas)} subfamily FASTA files")

    # â”€â”€ Step 1: Bitscore matrix â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    matrix, labels, subfamilies = build_bitscore_matrix(
        assign_path, cons_path, subfamily_fastas, run_root, threads)

    # â”€â”€ Step 2: PCA â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    run_pca(matrix, labels, subfamilies, outdir)

    # â”€â”€ Step 3: Position weights â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    wt_df, top_positions = compute_position_weights(
        consensus_seqs, outdir, top_k=args.top_k)

    # â”€â”€ Step 4: Per-copy diagnostic state â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    all_states, subfam_profiles = compute_copy_diagnostic_state(
        subfamily_fastas, consensus_seqs, top_positions, wt_df, outdir, threads)

    # â”€â”€ Step 5: Diagnostic scoring + discordance â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    flag_rows = None
    if all_states and subfam_profiles:
        flag_rows = compute_diagnostic_scores(
            all_states, subfam_profiles, wt_df, assign_path, matrix, labels, outdir)

    # â”€â”€ Step 6: Diagnostic heatmap â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if all_states and subfam_profiles:
        plot_diagnostic_heatmap(all_states, subfam_profiles, wt_df,
                                labels, assign_path, outdir)

    # â”€â”€ Step 7: Diagnostic vs bitscore scatter â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if flag_rows:
        plot_diagnostic_vs_bitscore(flag_rows, outdir)

    # â”€â”€ Step 8: Volcano plots â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    plot_volcano_plots(consensus_seqs, wt_df, outdir)

    # â”€â”€ Step 9: SINEplot input â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    if not args.skip_sineplot:
        generate_sineplot_input(consensus_seqs, subfamily_fastas, outdir, threads)

    # â”€â”€ Summary â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    log("")
    log("=" * 60)
    log("  Diagnostic Analysis Complete")
    log("=" * 60)
    log(f"  Output directory: {outdir}")
    log(f"  Output files:")
    for f in sorted(os.listdir(outdir)):
        fpath = os.path.join(outdir, f)
        size_kb = os.path.getsize(fpath) / 1024
        log(f"    {f} ({size_kb:.1f} KB)")
    log("")
    log("  To generate SINEplot interactive visualization:")
    log(f"    python SINEplot.py {os.path.join(outdir, 'sineplot_input.txt')} "
        f"-o {os.path.join(outdir, 'sineplot.html')}")
    log("")


if __name__ == "__main__":
    main()
