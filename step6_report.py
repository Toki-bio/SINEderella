#!/usr/bin/env python3
"""
step6_report.py - Build a single self-contained interactive HTML report
                  for a completed SINEderella run.

Design:
* Numeric tables (with column legends) for things that compress well into
  small tables: step1 hits, assignment stats, per-subfamily composition,
  threshold vs self-bits, quality flags.
* Interactive Plotly only where it adds insight: pipeline funnel,
  per-position consensus-base conservation, similarity histogram,
  similarity violins.
* Per-subfamily PNG gallery from step4 is embedded inline (base64).
* Optional SINEplot (https://github.com/Toki-bio/SINEplot) PCA: if
  SINEplot.py is reachable on PATH (or under $SCRIPT_DIR/SINEplot/),
  step6 builds a downsampled all-vs-all ssearch36 score file and embeds
  the resulting standalone HTML via <iframe srcdoc>.

Stdlib only (urllib used to fetch plotly.js once and cache).
"""
from __future__ import annotations

import argparse
import base64
import csv
import glob
import html
import json
import math
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from urllib.parse import quote
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PLOTLY_VERSION = "2.35.2"
PLOTLY_URL = f"https://cdn.plot.ly/plotly-{PLOTLY_VERSION}.min.js"
CACHE_DIR = Path(os.environ.get("XDG_CACHE_HOME",
                                Path.home() / ".cache")) / "sinederella"

LOG = logging.getLogger("step6_report")


# ===========================================================================
# Discovery
# ===========================================================================

def find_step2_out(run_root: Path) -> Path:
    candidates = sorted(
        glob.glob(str(run_root / "step2" / "step2_output*")),
        key=os.path.getmtime, reverse=True,
    )
    for c in candidates:
        if (Path(c) / "assignment_full.tsv").is_file():
            return Path(c)
    raise SystemExit(f"No step2_output* with assignment_full.tsv under {run_root}")


# ===========================================================================
# File parsers (stdlib only)
# ===========================================================================

def read_kv_manifest(path: Path) -> Dict[str, str]:
    out = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            out[k.strip()] = v.strip()
    return out


def read_tsv(path: Path, has_header: bool = True,
             max_rows: Optional[int] = None
             ) -> Tuple[List[str], List[List[str]]]:
    if not path.is_file():
        return [], []
    rows: List[List[str]] = []
    header: List[str] = []
    with path.open(newline="", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for i, row in enumerate(reader):
            if i == 0 and has_header:
                header = row
                continue
            rows.append(row)
            if max_rows is not None and len(rows) >= max_rows:
                break
    return header, rows


def parse_step1_hits(stderr_log: Path, stdout_log: Path
                     ) -> Tuple[Dict[str, int], Optional[int]]:
    hits: Dict[str, int] = {}
    total: Optional[int] = None
    rx_hit = re.compile(r"->\s*([^\s:]+):\s*(\d+)\s*hits", re.I)
    rx_per = re.compile(r"^\s*([^\s:]+):\s*(\d+)\s+hits\s*$", re.I)
    rx_total = re.compile(r"Total\s+merged\s+hits:\s*(\d+)", re.I)
    for path in (stderr_log, stdout_log):
        if not path.is_file():
            continue
        for line in path.read_text(errors="replace").splitlines():
            m = rx_hit.search(line) or rx_per.match(line)
            if m and m.group(1) and m.group(2):
                hits[m.group(1)] = int(m.group(2))
            m2 = rx_total.search(line)
            if m2:
                total = int(m2.group(1))
    return hits, total


def parse_step2_summary(summary_txt: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not summary_txt.is_file():
        return out
    txt = summary_txt.read_text(errors="replace")
    for key, pat in [
        ("total",      r"Total sequences:\s*(\d+)"),
        ("unanimous",  r"Unanimous\s+\d+/\d+:\s*(\d+)"),
        ("assigned",   r"Assigned\s*\(passed threshold\):\s*(\d+)"),
        ("unassigned", r"Unassigned:\s*(\d+)"),
        ("date",       r"Date:\s*(.+)"),
    ]:
        m = re.search(pat, txt)
        if m:
            out[key] = m.group(1).strip()
    return out


def stratified_sample_sim(sim_path: Path, assign_path: Path,
                          per_group: int = 3000
                          ) -> Dict[str, List[float]]:
    """Return {subfam: [sim_ratio*100, ...]} sampled per subfamily."""
    if not sim_path.is_file() or not assign_path.is_file():
        return {}
    seq_to_sf: Dict[str, str] = {}
    with assign_path.open(errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[4] == "assigned":
                seq_to_sf[parts[0]] = parts[1]

    rng = random.Random(42)
    by_sf: Dict[str, List[float]] = {}
    counts: Dict[str, int] = {}
    with sim_path.open(errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            sid = parts[0]
            sf = seq_to_sf.get(sid)
            if not sf:
                continue
            try:
                v = float(parts[3]) * 100.0
            except ValueError:
                continue
            counts[sf] = counts.get(sf, 0) + 1
            buf = by_sf.setdefault(sf, [])
            if len(buf) < per_group:
                buf.append(v)
            else:
                j = rng.randrange(counts[sf])
                if j < per_group:
                    buf[j] = v
    return by_sf


def count_flags_per_subfam(bedlike_path: Path
                           ) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    if not bedlike_path.is_file():
        return out
    with bedlike_path.open(errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            sf = parts[4] if len(parts) > 4 else ""
            flag = parts[8] if len(parts) > 8 else ""
            d = out.setdefault(sf, {"CONFLICT": 0, "LEAK": 0,
                                    "OK": 0, "total": 0})
            d["total"] += 1
            if "CONFLICT" in flag:
                d["CONFLICT"] += 1
            elif "LEAK" in flag:
                d["LEAK"] += 1
            else:
                d["OK"] += 1
    return out


def read_nucfreq_tsv(path: Path) -> Optional[Tuple[List[int], List[str],
                                                   Dict[str, List[float]]]]:
    """Read step4 companion TSV: pos, cons_base, A, T, C, G, gap.

    Returns (positions, cons_bases, freqs_by_nuc) or None.
    """
    if not path.is_file():
        return None
    pos: List[int] = []
    cb: List[str] = []
    cols = {"A": [], "T": [], "C": [], "G": [], "gap": []}
    with path.open(errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("pos\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                pos.append(int(parts[0]))
                cb.append(parts[1])
                cols["A"].append(float(parts[2]))
                cols["T"].append(float(parts[3]))
                cols["C"].append(float(parts[4]))
                cols["G"].append(float(parts[5]))
                cols["gap"].append(float(parts[6]))
            except ValueError:
                continue
    if not pos:
        return None
    return pos, cb, cols


def conservation_curve(nucfreq: Tuple[List[int], List[str],
                                      Dict[str, List[float]]]
                       ) -> Tuple[List[int], List[float]]:
    """Per-position frequency of the consensus base."""
    pos, cb, cols = nucfreq
    y = []
    for i, base in enumerate(cb):
        b = base.upper()
        if b in cols and b != "gap":
            y.append(cols[b][i])
        else:
            y.append(0.0)
    return pos, y


# ===========================================================================
# Plotly figures (returned as dicts; no plotly Python lib used)
# ===========================================================================

COLORS = {
    "primary":  "#4C72B0",
    "ok":       "#55A868",
    "warn":     "#DD8452",
    "danger":   "#C44E52",
    "muted":    "#8C8C8C",
    "soft":     "#8172B2",
}

# Qualitative palette for per-subfamily traces (colour-blind friendly).
SF_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
    "#aec7e8", "#ffbb78",
]


def fig_funnel(stats: Dict[str, str]) -> dict:
    labels = ["Hits found (step1)", "Unanimous 10/10",
              "Assigned", "Unassigned"]
    values = [
        int(stats.get("total", 0) or 0),
        int(stats.get("unanimous", 0) or 0),
        int(stats.get("assigned", 0) or 0),
        int(stats.get("unassigned", 0) or 0),
    ]
    return {
        "data": [{
            "type": "funnel",
            "y": labels,
            "x": values,
            "textinfo": "value+percent initial",
            "marker": {"color": [COLORS["primary"], COLORS["soft"],
                                 COLORS["ok"], COLORS["warn"]]},
        }],
        "layout": {
            "title": "Pipeline funnel: candidates &rarr; assignments",
            "margin": {"l": 180, "r": 20, "t": 50, "b": 30},
            "height": 360,
        },
    }


def funnel_html(stats: Dict[str, str]) -> str:
    """Simple two-branch HTML summary of pipeline counts (no chart)."""
    def fmt(n: str) -> str:
        try:
            return f"{int(n):,}"
        except Exception:
            return str(n) if n else "?"

    def pct(num: str, den: str) -> str:
        try:
            n_, d_ = int(num), int(den)
            if d_ <= 0:
                return ""
            return (f" <span class='muted small'>({100.0 * n_ / d_:.1f}%)</span>")
        except Exception:
            return ""

    total    = stats.get("total", "")
    unan     = stats.get("unanimous", "")
    assigned = stats.get("assigned", "")
    unassgn  = stats.get("unassigned", "")
    return (
        "<table class='tbl' style='max-width:540px'>"
        "<colgroup><col style='width:62%'><col></colgroup>"
        f"<tr><td>Total sequences (step1 merged hits entering step2)</td>"
        f"<td><b>{fmt(total)}</b></td></tr>"
        f"<tr><td style='padding-left:1.4em'>&#9492; Unanimous"
        f" (10&#47;10 sub-sample votes)</td>"
        f"<td>{fmt(unan)}{pct(unan, total)}</td></tr>"
        f"<tr><td><b>&#10003;&thinsp;Assigned</b>"
        f" (passed bitscore threshold)</td>"
        f"<td><b>{fmt(assigned)}</b>{pct(assigned, total)}</td></tr>"
        f"<tr><td>&#10007;&thinsp;Unassigned</td>"
        f"<td>{fmt(unassgn)}{pct(unassgn, total)}</td></tr>"
        "</table>"
        "<p class='small muted' style='margin-top:6px'>"
        "Assigned + Unassigned &asymp; Total. "
        "Unanimous is a strict subset; assigned includes unanimous + "
        "soft-assigned copies that passed the bitscore threshold.</p>"
    )


def kde_curve(vals: List[float], n_pts: int = 300,
              ) -> Tuple[List[float], List[float]]:
    """Gaussian KDE (no scipy). Scott's bandwidth: h = std * n^(-1/5)."""
    n = len(vals)
    if n < 2:
        return [], []
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / (n - 1)
    std = math.sqrt(max(var, 1e-12))
    h = std * n ** (-0.2)
    lo = min(vals) - 2.5 * h
    hi = max(vals) + 2.5 * h
    step = (hi - lo) / (n_pts - 1)
    x_arr = [lo + i * step for i in range(n_pts)]
    c = 1.0 / (n * h * math.sqrt(2 * math.pi))
    y_arr = [
        c * sum(math.exp(-0.5 * ((xi - v) / h) ** 2) for v in vals)
        for xi in x_arr
    ]
    return x_arr, y_arr


def fig_similarity_kde(by_sf: Dict[str, List[float]]) -> dict:
    """Per-subfamily KDE density lines of similarity-to-consensus (%)."""
    sf_sorted = sorted(by_sf.keys())
    traces = []
    for i, sf in enumerate(sf_sorted):
        vals = by_sf[sf]
        if not vals:
            continue
        x, y = kde_curve(vals)
        if not x:
            continue
        traces.append({
            "type": "scatter", "mode": "lines",
            "x": [round(v, 3) for v in x],
            "y": [round(v, 6) for v in y],
            "name": sf,
            "line": {"color": SF_PALETTE[i % len(SF_PALETTE)], "width": 2},
            "hovertemplate": (
                "%{fullData.name}<br>sim %{x:.1f}%"
                "<br>density %{y:.5f}<extra></extra>"),
        })
    return {
        "data": traces,
        "layout": {
            "title": "Similarity to consensus \u2014 per-copy distribution (KDE)",
            "xaxis": {"title": "Similarity to subfamily consensus (%)"},
            "yaxis": {"title": "Density"},
            "legend": {"title": {"text": "Subfamily (click to toggle)"}},
            "height": 460,
            "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
        },
    }


def fig_divergence_kde(by_sf: Dict[str, List[float]]) -> dict:
    """Per-subfamily KDE density lines of divergence = max(0, 100 - similarity)."""
    sf_sorted = sorted(by_sf.keys())
    traces = []
    for i, sf in enumerate(sf_sorted):
        # clamp to [0, 100]: bitscore-based similarity can exceed 100% due to
        # local alignment scoring (copy bitscore > consensus self-bitscore)
        vals = [max(0.0, 100.0 - v) for v in by_sf[sf]]
        if not vals:
            continue
        x, y = kde_curve(vals)
        if not x:
            continue
        traces.append({
            "type": "scatter", "mode": "lines",
            "x": [round(max(0.0, v), 3) for v in x],
            "y": [round(v, 6) for v in y],
            "name": sf,
            "line": {"color": SF_PALETTE[i % len(SF_PALETTE)], "width": 2},
            "hovertemplate": (
                "%{fullData.name}<br>divergence %{x:.1f}%"
                "<br>density %{y:.5f}<extra></extra>"),
        })
    return {
        "data": traces,
        "layout": {
            "title": "Bitscore-based divergence from consensus \u2014 per-copy KDE",
            "xaxis": {"title": "Divergence (100 \u2212 bitscore/self-bits \u00d7 100%)",
                      "rangemode": "nonnegative"},
            "yaxis": {"title": "Density"},
            "legend": {"title": {"text": "Subfamily (click to toggle)"}},
            "height": 460,
            "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
        },
    }


def fig_sim_violins(by_sf: Dict[str, List[float]]) -> dict:
    sf_sorted = sorted(by_sf.keys())
    traces = []
    for sf in sf_sorted:
        vals = by_sf[sf]
        if not vals:
            continue
        traces.append({
            "type": "violin",
            "y": vals,
            "name": sf,
            "box": {"visible": True},
            "meanline": {"visible": True},
            "points": False,
        })
    return {
        "data": traces,
        "layout": {
            "title": "Similarity to consensus -- distribution per subfamily",
            "yaxis": {"title": "Similarity (%)"},
            "xaxis": {"title": "Subfamily", "tickangle": -30},
            "height": 460,
            "showlegend": False,
            "margin": {"t": 50, "r": 20, "b": 100, "l": 60},
        },
    }


def fig_conservation(curves: Dict[str, Tuple[List[int], List[float]]]
                     ) -> dict:
    """One trace per subfamily: x=consensus position (bp), y=cons-base freq."""
    sf_sorted = sorted(curves.keys())
    traces = []
    for i, sf in enumerate(sf_sorted):
        x, y = curves[sf]
        traces.append({
            "type": "scatter", "mode": "lines",
            "x": x, "y": y, "name": sf,
            "line": {"width": 1.4,
                     "color": SF_PALETTE[i % len(SF_PALETTE)]},
            "hovertemplate": (
                "%{fullData.name}<br>pos %{x} bp<br>"
                "cons-base freq %{y:.3f}<extra></extra>"),
        })
    return {
        "data": traces,
        "layout": {
            "title": "Per-position conservation "
                     "(fraction of copies matching the consensus base)",
            "xaxis": {"title": "Position on consensus (bp)"},
            "yaxis": {"title": "Fraction of copies = consensus base",
                      "range": [0, 1.02]},
            "legend": {"title": {"text": "Subfamily (click to toggle)"}},
            "height": 460,
            "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
        },
    }


# ===========================================================================
# SINEplot integration
# ===========================================================================

def find_sineplot(script_dir: Optional[Path]) -> Optional[Path]:
    env = os.environ.get("SINEPLOT_PY")
    if env and Path(env).is_file():
        return Path(env)
    if script_dir:
        for cand in [script_dir / "SINEplot.py",
                     script_dir / "SINEplot" / "SINEplot.py"]:
            if cand.is_file():
                return cand
    on_path = shutil.which("SINEplot.py") or shutil.which("SINEplot")
    if on_path:
        return Path(on_path)
    return None


def _read_fasta(path: Path) -> List[Tuple[str, str]]:
    out, name, seq = [], None, []
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def _write_fasta(records: List[Tuple[str, str]], path: Path) -> None:
    with path.open("w") as fh:
        for n, s in records:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def build_sineplot_iframe(s2_out: Path,
                          script_dir: Optional[Path],
                          max_copies_per_sf: int,
                          threads: int) -> Optional[str]:
    """Run SINEplot, return its HTML content for <iframe srcdoc>.

    Returns None (and logs the reason) if SINEplot or its inputs are missing.
    """
    sp = find_sineplot(script_dir)
    if sp is None:
        LOG.info("SINEplot not found; skipping PCA panel")
        return None
    if not shutil.which("ssearch36"):
        LOG.info("ssearch36 not in PATH; skipping SINEplot PCA panel")
        return None
    cons = s2_out.parent.parent / "consensuses.clean.fa"
    subfam_dir = s2_out / "subfamilies"
    if not cons.is_file() or not subfam_dir.is_dir():
        LOG.info("No consensus / subfamily fastas; skipping SINEplot")
        return None

    LOG.info("Running SINEplot (%s, max %d copies/sf)",
             sp, max_copies_per_sf)
    rng = random.Random(42)
    work = Path(tempfile.mkdtemp(prefix="sineplot_"))
    try:
        all_fa = work / "all.fa"
        records = _read_fasta(cons)
        # Try per-subfamily fastas first; fall back to assigned.fasta when
        # those are empty (some pipeline versions only populate assigned.fasta).
        sf_buckets: Dict[str, List[Tuple[str, str]]] = {}
        for fa in sorted(subfam_dir.glob("*.fasta")):
            sf_records = _read_fasta(fa)
            if sf_records:
                sf_buckets[fa.stem] = sf_records
        if not sf_buckets:
            assigned = s2_out / "assigned.fasta"
            if assigned.is_file():
                LOG.info("  subfamily fastas empty; reading %s", assigned.name)
                for n, s in _read_fasta(assigned):
                    parts = n.split("|")
                    if len(parts) < 2:
                        continue
                    sf = parts[1]
                    sf_buckets.setdefault(sf, []).append((parts[0], s))
        for sf, sf_records in sorted(sf_buckets.items()):
            if len(sf_records) > max_copies_per_sf:
                sf_records = rng.sample(sf_records, max_copies_per_sf)
            sf_records = [(r[0].split("|")[0], r[1]) for r in sf_records]
            records.extend(sf_records)
        _write_fasta(records, all_fa)

        scores = work / "scores.txt"
        LOG.info("  ssearch36 all-vs-all on %d sequences", len(records))
        with scores.open("w") as fh:
            r = subprocess.run(
                ["ssearch36", "-m", "8", "-T", str(threads),
                 str(all_fa), str(all_fa)],
                stdout=fh, stderr=subprocess.PIPE, check=False)
        if r.returncode != 0:
            LOG.warning("ssearch36 failed (rc=%d): %s",
                        r.returncode, r.stderr.decode(errors="replace")[:300])
            return None

        out_html = work / "sineplot.html"
        # SINEplot internally samples 2*max_points for positioning, then
        # downsamples to max_points for display. Keep the cap conservative
        # so the positioning loop stays fast.
        sineplot_max_display = min(max_copies_per_sf, 300)
        LOG.info("  invoking SINEplot.py (--max-points %d, timeout 1800s)",
                 sineplot_max_display)
        try:
            r2 = subprocess.run(
                [sys.executable, str(sp), str(scores),
                 "-o", str(out_html),
                 "-t", "SINEplot PCA -- bitscore-based subfamily layout",
                 "--max-points", str(sineplot_max_display)],
                capture_output=True, check=False, timeout=1800)
        except subprocess.TimeoutExpired:
            LOG.warning("SINEplot.py timed out (>1800s); skipping PCA panel")
            return None
        if r2.returncode != 0 or not out_html.is_file():
            LOG.warning("SINEplot.py failed (rc=%d): %s",
                        r2.returncode,
                        (r2.stderr or b"").decode(errors="replace")[:300])
            return None
        return out_html.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        LOG.warning("SINEplot embedding failed: %s", exc)
        return None
    finally:
        shutil.rmtree(work, ignore_errors=True)


def build_pca_fig(run_root: Path, s2_out: Path,
                  n_per_sf: int = 200,
                  threads: int = 8) -> Optional[dict]:
    """Mutation-space PCA of assigned SINE copies.

    Each copy is represented as a binary vector: 1 at alignment positions
    where it differs from its assigned-subfamily consensus, 0 elsewhere.
    mafft --add adds copies into the consensus alignment (keeplength),
    then SVD gives PC1/PC2.  Returns a Plotly figure dict + metadata dict,
    or None if prerequisites are missing.
    """
    try:
        import numpy as _np
    except ImportError:
        LOG.warning("numpy not available; skipping mutation-space PCA")
        return None
    if not shutil.which("mafft"):
        LOG.warning("mafft not in PATH; skipping PCA")
        return None

    # Locate consensus FASTA
    cons: Optional[Path] = None
    for cand in [run_root / "consensuses.clean.fa",
                 run_root / "results" / "consensuses.fa",
                 s2_out.parent.parent / "consensuses.clean.fa"]:
        if cand.is_file():
            cons = cand
            break
    if cons is None:
        LOG.warning("No consensuses FASTA found; skipping PCA")
        return None

    assigned = s2_out / "assigned.fasta"
    if not assigned.is_file():
        LOG.warning("No assigned.fasta; skipping PCA")
        return None

    cons_records = _read_fasta(cons)
    cons_name_set = {n for n, _ in cons_records}
    if len(cons_name_set) < 2:
        LOG.warning("Fewer than 2 consensuses; PCA skipped")
        return None

    # Stratified sample: up to n_per_sf copies per subfamily.
    # Subfamily parsed from assigned.fasta header: >id|sf|score
    LOG.info("PCA: reading assigned.fasta for stratified sample ...")
    all_copies = _read_fasta(assigned)
    seq_to_sf: Dict[str, str] = {}
    by_sf: Dict[str, list] = {}
    for rec in all_copies:
        parts = rec[0].split("|")
        sf = parts[1] if len(parts) >= 2 else "unknown"
        seq_to_sf[rec[0]] = sf
        by_sf.setdefault(sf, []).append(rec)

    rng = random.Random(42)
    sampled: list = []
    for sf in sorted(by_sf):
        pool = by_sf[sf]
        take = min(len(pool), n_per_sf)
        sampled.extend(rng.sample(pool, take))

    n_copies = len(sampled)
    LOG.info("PCA: %d copies (up to %d per subfamily, %d subfamilies)",
             n_copies, n_per_sf, len(by_sf))

    # Use simplified names so mafft never truncates long headers
    simple_to_sf: Dict[str, str] = {}
    simple_records: list = []
    for i, (orig, seq) in enumerate(sampled):
        sname = f"cp{i}"
        simple_to_sf[sname] = seq_to_sf.get(orig, "unknown")
        simple_records.append((sname, seq))

    work = Path(tempfile.mkdtemp(prefix="pca_"))
    try:
        copies_fa = work / "copies.fa"
        cons_tmp  = work / "cons.fa"
        cons_aln  = work / "cons_aln.fa"
        all_aln   = work / "all_aln.fa"

        _write_fasta(simple_records, copies_fa)
        _write_fasta(cons_records,   cons_tmp)

        # Step 1: align consensus sequences
        LOG.info("PCA: mafft aligning %d consensuses ...", len(cons_records))
        r1 = subprocess.run(
            ["mafft", "--auto", "--quiet",
             "--thread", str(threads), "--nuc", str(cons_tmp)],
            capture_output=True, check=False, timeout=120)
        if r1.returncode != 0:
            LOG.warning("mafft (consensus) failed: %s",
                        r1.stderr.decode(errors="replace")[:300])
            return None
        cons_aln.write_bytes(r1.stdout)

        # Step 2: add copies (fixed length = consensus alignment)
        LOG.info("PCA: mafft --add %d copies ...", n_copies)
        r2 = subprocess.run(
            ["mafft", "--add", str(copies_fa),
             "--keeplength", "--quiet",
             "--thread", str(threads), "--nuc", str(cons_aln)],
            capture_output=True, check=False, timeout=900)
        if r2.returncode != 0:
            LOG.warning("mafft --add failed: %s",
                        r2.stderr.decode(errors="replace")[:300])
            return None
        all_aln.write_bytes(r2.stdout)

        # Parse alignment
        cons_aln_dict: Dict[str, str] = {}
        copy_aln: list = []
        for name, seq in _read_fasta(all_aln):
            if name in cons_name_set:
                cons_aln_dict[name] = seq.upper()
            else:
                copy_aln.append((name, seq.upper()))

        if not cons_aln_dict or not copy_aln:
            LOG.warning("PCA: alignment parsing yielded no data")
            return None

        aln_len = len(next(iter(cons_aln_dict.values())))
        copy_names = [n for n, _ in copy_aln]
        assigned_sfs = [simple_to_sf.get(n, "unknown") for n in copy_names]
        n_copies_actual = len(copy_aln)

        # Binary mutation matrix: 1 where copy ≠ assigned-sf consensus
        LOG.info("PCA: building mutation matrix (%d × %d) ...",
                 n_copies_actual, aln_len)
        X = _np.zeros((n_copies_actual, aln_len), dtype=_np.float32)
        for i, (cname, cseq) in enumerate(copy_aln):
            sf = assigned_sfs[i]
            ref = cons_aln_dict.get(sf) or next(iter(cons_aln_dict.values()))
            for j in range(aln_len):
                cc = cseq[j]   if j < len(cseq) else "-"
                rc = ref[j]    if j < len(ref)  else "-"
                if cc not in "-Nn" and rc not in "-Nn" and cc != rc:
                    X[i, j] = 1.0

        # Keep variable columns (mutation rate 3 %–97 %)
        col_rate = X.mean(axis=0)
        keep = (col_rate >= 0.03) & (col_rate <= 0.97)
        X = X[:, keep]
        n_var = int(keep.sum())
        LOG.info("PCA: %d variable columns retained", n_var)
        if n_var < 2:
            LOG.warning("PCA: fewer than 2 variable columns; skipping")
            return None

        # SVD
        X_c = (X - X.mean(axis=0)).astype(_np.float64)
        U, S, Vt = _np.linalg.svd(X_c, full_matrices=False)
        pc1 = (X_c @ Vt[0]).tolist()
        pc2 = (X_c @ Vt[1]).tolist()
        var_total = float((S ** 2).sum()) or 1.0
        pct1 = round(float(S[0] ** 2) / var_total * 100, 1)
        pct2 = round(float(S[1] ** 2) / var_total * 100, 1)

        # eta² — how well PC1 separates subfamilies
        sf_uniq = sorted(set(assigned_sfs))
        pc1_arr = _np.array(pc1)
        grand_mean = float(pc1_arr.mean())
        ss_total = float(((pc1_arr - grand_mean) ** 2).sum())
        ss_between = sum(
            float((pc1_arr[_np.array(assigned_sfs) == sf] - grand_mean).sum()) ** 2
            / max(1, int((_np.array(assigned_sfs) == sf).sum()))
            for sf in sf_uniq
        )
        eta2 = round(ss_between / ss_total, 3) if ss_total > 0 else 0.0

        subfams_ordered = sf_uniq
        traces = []
        sf_arr = _np.array(assigned_sfs)
        for idx, sf in enumerate(subfams_ordered):
            mask = sf_arr == sf
            xi = [round(float(v), 4) for v in pc1_arr[mask]]
            yi = [round(float(v), 4) for v in _np.array(pc2)[mask]]
            ti = [copy_names[i] for i in range(n_copies_actual) if assigned_sfs[i] == sf]
            if not xi:
                continue
            traces.append({
                "type": "scatter", "mode": "markers",
                "x": xi, "y": yi, "name": sf, "text": ti,
                "marker": {
                    "size": 5, "opacity": 0.65,
                    "color": SF_PALETTE[idx % len(SF_PALETTE)],
                },
                "hovertemplate": (
                    "<b>%{fullData.name}</b><br>%{text}<br>"
                    "PC1&thinsp;%{x:.3f} &nbsp;PC2&thinsp;%{y:.3f}"
                    "<extra></extra>"),
            })
        return {
            "data": traces,
            "layout": {
                "title": (
                    f"Mutation-space PCA \u2014 per-copy alignment differences "
                    f"from assigned-subfamily consensus"
                    f"<br><sub>n\u2009=\u2009{n_copies_actual:,} copies &middot; "
                    f"{n_var} variable cols &middot; "
                    f"PC1\u2009{pct1}\u2009% &middot; "
                    f"PC2\u2009{pct2}\u2009% variance &middot; "
                    f"\u03b7\u00b2(SF\u2192PC1)\u2009=\u2009{eta2}</sub>"
                ),
                "xaxis": {"title": f"PC1 ({pct1}\u2009% variance)",
                          "zeroline": True},
                "yaxis": {"title": f"PC2 ({pct2}\u2009% variance)",
                          "zeroline": True},
                "legend": {"title": {"text": "Subfamily"}},
                "height": 600,
                "margin": {"t": 80, "r": 20, "b": 60, "l": 70},
                "_meta": {"n_copies": n_copies_actual, "n_var": n_var,
                          "pct1": pct1, "pct2": pct2, "eta2": eta2},
            },
        }
    except Exception as exc:
        LOG.warning("PCA failed: %s", exc, exc_info=True)
        return None
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ===========================================================================
# HTML helpers
# ===========================================================================

def get_plotly_js(inline: bool) -> Tuple[str, str]:
    if not inline:
        return (f'<script src="{PLOTLY_URL}" charset="utf-8"></script>',
                "CDN")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / f"plotly-{PLOTLY_VERSION}.min.js"
    if not cache.is_file():
        LOG.info("Downloading Plotly.js to cache: %s", cache)
        with urllib.request.urlopen(PLOTLY_URL, timeout=60) as r:
            cache.write_bytes(r.read())
    js = cache.read_text(encoding="utf-8")
    return f"<script>{js}</script>", "inline"


def img_to_data_uri(p: Path) -> Optional[str]:
    if not p.is_file():
        return None
    try:
        b = p.read_bytes()
    except OSError:
        return None
    return "data:image/png;base64," + base64.b64encode(b).decode("ascii")


def render_table(header: List[str], rows: List[List[str]],
                 max_rows: int = 100,
                 col_titles: Optional[Dict[str, str]] = None,
                 escape: bool = True,
                 highlight_cols: Optional[List[str]] = None) -> str:
    if not header and not rows:
        return "<p><em>(no data)</em></p>"
    titles = col_titles or {}
    hl_set = set(highlight_cols or [])
    hl_idx = {i for i, h in enumerate(header) if h in hl_set}
    head_cells = []
    for h in header:
        title = titles.get(h, "")
        attr = f' title="{html.escape(title)}"' if title else ""
        hl = ' class="hl"' if h in hl_set else ""
        head_cells.append(f"<th{attr}{hl}>{html.escape(h)}</th>")
    head = "".join(head_cells)
    body_rows = rows[:max_rows]

    def cell(c, idx):
        s = str(c)
        v = html.escape(s) if escape else s
        cls = ' class="hl"' if idx in hl_idx else ""
        return f"<td{cls}>{v}</td>"

    body = "".join(
        "<tr>" + "".join(cell(c, i) for i, c in enumerate(r)) + "</tr>"
        for r in body_rows
    )
    extra = ""
    if len(rows) > max_rows:
        extra = (f"<p class='small muted'>... showing first {max_rows} "
                 f"of {len(rows)} rows.</p>")
    return (f"<table class='tbl'><thead><tr>{head}</tr></thead>"
            f"<tbody>{body}</tbody></table>{extra}")


def render_legend(items: List[Tuple[str, str]]) -> str:
    if not items:
        return ""
    lis = "".join(
        f"<li><code>{html.escape(k)}</code> &mdash; {v}</li>"
        for k, v in items
    )
    return (f"<details class='legend'><summary>Column legend</summary>"
            f"<ul>{lis}</ul></details>")


CSS = """
:root { --fg:#222; --bg:#fafafa; --card:#fff; --accent:#4C72B0;
        --muted:#666; --border:#e2e2e2; }
* { box-sizing: border-box; }
body { font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial,
       sans-serif; background: var(--bg); color: var(--fg);
       margin: 0; padding: 0; }
header { background: linear-gradient(135deg, #2c3e50, #4C72B0);
         color: #fff; padding: 24px 32px; }
header h1 { margin: 0 0 6px 0; font-size: 1.6rem; }
header .sub { opacity: 0.85; font-size: 0.95rem; }
main { max-width: 1280px; margin: 0 auto; padding: 24px; }
section.card { background: var(--card); border: 1px solid var(--border);
               border-radius: 8px; padding: 18px 22px; margin: 18px 0;
               box-shadow: 0 1px 3px rgba(0,0,0,0.04); }
section.card h2 { margin-top: 0; font-size: 1.15rem; color: #2c3e50;
                  border-bottom: 1px solid var(--border);
                  padding-bottom: 8px; }
section.card h3 { font-size: 1rem; margin-top: 20px; color: #444; }
section.card p.intro { color: var(--muted); font-size: 0.9rem;
                       margin-top: 4px; margin-bottom: 14px; }
.kv { display: grid; grid-template-columns: max-content 1fr;
      gap: 4px 16px; font-size: 0.92rem; }
.kv .k { color: var(--muted); }
.kv .v { font-family: ui-monospace, "Cascadia Mono", Menlo, monospace;
         word-break: break-all; }
.tbl { border-collapse: collapse; width: 100%; font-size: 0.88rem;
       margin: 8px 0; }
.tbl th, .tbl td { border: 1px solid var(--border); padding: 5px 8px;
                   text-align: left; }
.tbl th { background: #f0f3f7; cursor: help; }
.tbl tbody tr:nth-child(even) { background: #fafbfd; }
.tbl td:not(:first-child) { font-variant-numeric: tabular-nums; }
.metrics { display: flex; flex-wrap: wrap; gap: 14px; margin: 8px 0 16px 0; }
.metric { background: #f0f3f7; border-left: 4px solid var(--accent);
          padding: 10px 14px; border-radius: 4px; min-width: 160px; }
.metric .num { font-size: 1.4rem; font-weight: 600; color: #2c3e50; }
.metric .lbl { font-size: 0.78rem; color: var(--muted);
               text-transform: uppercase; letter-spacing: 0.04em; }
.subfam-grid { display: grid; gap: 18px;
               grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); }
.subfam-block { border: 1px solid var(--border); border-radius: 6px;
                padding: 10px 12px; background: #fff; }
.subfam-block h3 { margin: 0 0 8px 0; font-size: 0.95rem; color: #4C72B0; }
.subfam-block img { max-width: 100%; height: auto; display: block;
                    margin: 6px 0; border: 1px solid #eee;
                    cursor: zoom-in; transition: opacity .15s; }
.subfam-block img:hover { opacity: .85; }
/* Lightbox */
#lightbox { display: none; position: fixed; inset: 0;
            background: rgba(0,0,0,.88); z-index: 9999;
            align-items: center; justify-content: center;
            cursor: zoom-out; }
#lightbox.active { display: flex; }
#lightbox img { max-width: 94vw; max-height: 94vh; border-radius: 4px;
                box-shadow: 0 8px 40px rgba(0,0,0,.6); }
/* Alignment links */
.aln-link { display: inline-block; background: var(--accent); color: #fff;
            padding: 3px 10px; border-radius: 4px; font-size: .83rem;
            text-decoration: none; margin: 1px 0; }
.aln-link:hover { background: #2c3e50; }
.aln-link.orange { background: #e07b39; }
.aln-link.green  { background: #28a745; }
.aln-link.green:hover { background: #1e7e34; }
.small { font-size: 0.8rem; }
.muted { color: var(--muted); }
/* Highlighted table columns */
.tbl th.hl { background: #d6e8ff; }
.tbl td.hl { font-weight: 600; color: #1a3a6e; background: #f4f8ff; }
/* Collapsible card (details element styled like section.card) */
details.card { background: var(--card); border: 1px solid var(--border);
               border-radius: 8px; margin: 18px 0;
               box-shadow: 0 1px 3px rgba(0,0,0,0.04); }
details.card > summary { padding: 14px 22px; cursor: pointer;
               list-style: none; display: flex; align-items: center;
               user-select: none; }
details.card > summary::-webkit-details-marker { display: none; }
details.card > summary h2 { margin: 0; font-size: 1.15rem; color: #2c3e50;
               padding-bottom: 0; border-bottom: none; flex: 1; }
details.card > summary::before { content: '\25B6'; margin-right: 10px;
               color: var(--accent); font-size: .75rem;
               transition: transform .18s; }
details.card[open] > summary::before { transform: rotate(90deg); }
details.card[open] > summary { border-bottom: 1px solid var(--border); }
.card-body { padding: 4px 22px 18px; }
nav.toc { background: #fff; border: 1px solid var(--border);
          border-radius: 8px; padding: 12px 18px; margin: 18px 0;
          font-size: 0.92rem; }
nav.toc a { color: var(--accent); text-decoration: none;
            margin-right: 14px; }
nav.toc a:hover { text-decoration: underline; }
footer { padding: 18px 32px; color: var(--muted); font-size: 0.8rem;
         text-align: center; }
.plot { width: 100%; }
details.legend { margin: 6px 0 12px 0; font-size: 0.85rem; }
details.legend summary { color: var(--accent); cursor: pointer; }
details.legend ul { margin: 6px 0 0 0; padding-left: 20px;
                    color: var(--muted); }
details.legend code { background: #f0f3f7; padding: 1px 4px;
                      border-radius: 3px; color: #2c3e50; }
.skip-note { background: #fff8e1; border-left: 4px solid #DD8452;
             padding: 10px 14px; border-radius: 4px;
             color: #6a5028; font-size: 0.9rem; }
iframe.embed { width: 100%; height: 800px; border: 1px solid var(--border);
               border-radius: 6px; background: #fff; }
"""


# ===========================================================================
# Section builders
# ===========================================================================

LEG_ASSIGN_STATS = [
    ("Subfamily",      "Subfamily name from the consensus FASTA."),
    ("Assigned",       "Copies that passed all assignment criteria "
                       "(unanimous 10/10 vote AND bitscore &ge; threshold)."),
    ("TopN_Bitscore",  "Bitscore of the N-th best per-subfamily hit "
                       "(N = min(10, count)); reference for the threshold."),
    ("Threshold",      "Cutoff applied: 0.45 &times; TopN_Bitscore "
                       "&times; 100. Copies below this are unassigned."),
]

LEG_BYSUBFAM = [
    ("subfam",         "Subfamily name."),
    ("firm_assigned",  "Copies whose primary assignment uses 10/10 votes "
                       "and passed the threshold (high-confidence)."),
    ("soft_assigned",  "Copies promoted via soft-rules "
                       "(ssearch36 tie-breaking on otherwise-unassigned)."),
    ("total_assigned", "firm_assigned + soft_assigned."),
    ("leak_n",         "Copies of OTHER subfamilies that, in step3 sanity "
                       "checks, leaked into this subfamily's intervals."),
    ("conf_alt_n",     "Genomic intervals where this subfamily is the "
                       "winning name but at least one other subfamily "
                       "also had hits (CONFLICT-flagged)."),
    ("firm_pct",       "firm_assigned as % of all firm assignments."),
    ("total_pct",      "total_assigned as % of all assignments."),
    ("leak_pct",       "leak_n as % of total_assigned."),
    ("conf_alt_pct",   "conf_alt_n as % of total_assigned."),
    ("sim_mean",       "Mean per-copy bitscore / consensus self-bits."),
    ("sim_median",     "Median per-copy bitscore / consensus self-bits."),
]

LEG_THRESHOLDS = [
    ("Subfamily",        "Subfamily name."),
    ("Threshold (x100)", "Bitscore cutoff used by step2 "
                         "(stored as integer = bits &times; 100)."),
    ("RealSelfBits",     "ssearch36 bitscore of the consensus aligned "
                         "against itself (raw bitscore)."),
    ("Threshold/Self",   "Threshold (in bits, divided back by 100) "
                         "/ RealSelfBits. Roughly the minimum fractional "
                         "similarity to consensus a copy must reach."),
]

LEG_FLAGS = [
    ("Subfamily", "Subfamily name (winner of the merged interval)."),
    ("Total",     "Number of merged genomic intervals labelled with "
                  "this subfamily."),
    ("OK",        "Intervals with no conflict and no leak flag."),
    ("LEAK",      "Intervals where this subfamily 'leaked' a hit into "
                  "an interval owned by a different subfamily."),
    ("CONFLICT",  "Intervals where multiple subfamilies had hits in "
                  "the same merged region."),
    ("%CONFLICT", "100 &times; CONFLICT / Total."),
]

LEG_STEP1_HITS = [
    ("Query",    "Subfamily consensus FASTA used as ssearch36 query."),
    ("RawHits",  "Number of ssearch36 hits returned for this query "
                 "BEFORE bedtools-merge collapses overlapping intervals."),
]

# Human-readable column names for summary.by_subfam.tsv
_BYSF_RENAME: Dict[str, str] = {
    "subfam":         "Subfamily",
    "firm_assigned":  "Firm",
    "soft_assigned":  "Soft",
    "total_assigned": "Total copies",
    "leak_n":         "Leaks",
    "conf_alt_n":     "Conflicts",
    "firm_pct":       "Firm %",
    "total_pct":      "Total %",
    "leak_pct":       "Leak %",
    "conf_alt_pct":   "Conflict %",
    "sim_mean":       "Sim mean",
    "sim_median":     "Sim median",
}
_BYSF_HIGHLIGHT = {"Total copies", "Sim mean", "Sim median", "Conflict %"}


def render_bysf_table(header: List[str], rows: List[List[str]],
                      max_rows: int = 200) -> str:
    """Render summary.by_subfam table with renamed columns, legend, highlights."""
    disp_header = [_BYSF_RENAME.get(h, h) for h in header]
    ren_leg = [(_BYSF_RENAME.get(k, k), v) for k, v in LEG_BYSUBFAM]
    leg_html = render_legend(ren_leg)
    tbl_html = render_table(
        disp_header, rows, max_rows,
        col_titles={_BYSF_RENAME.get(k, k): _strip_html(v)
                    for k, v in LEG_BYSUBFAM},
        highlight_cols=list(_BYSF_HIGHLIGHT),
    )
    return leg_html + tbl_html


def _strip_html(s: str) -> str:
    return re.sub("<[^>]+>", "",
                  s.replace("&mdash;", "-")
                   .replace("&times;", "x")
                   .replace("&ge;", ">="))


def step1_hits_table(hits: Dict[str, int]) -> str:
    rows = sorted([(k, v) for k, v in hits.items()],
                  key=lambda x: -x[1])
    total = sum(v for _, v in rows)
    body = [[k, f"{v:,}"] for k, v in rows]
    body.append(["<b>TOTAL (sum of raw hits)</b>", f"<b>{total:,}</b>"])
    return render_table(["Query", "RawHits"], body,
                        max_rows=10000, escape=False,
                        col_titles={k: _strip_html(v)
                                    for k, v in LEG_STEP1_HITS})


def thresholds_table(stats_rows: List[List[str]],
                     self_bits_real_rows: List[List[str]]) -> str:
    real = {r[0]: float(r[1]) for r in self_bits_real_rows
            if len(r) >= 2}
    rows = []
    for r in stats_rows:
        if len(r) < 4:
            continue
        sf = r[0]
        thr = int(r[3]) if r[3].isdigit() else None
        rb = real.get(sf)
        ratio = ""
        if thr is not None and rb:
            ratio = f"{(thr / 100.0) / rb:.3f}"
        rows.append([sf,
                     f"{thr:,}" if thr is not None else "",
                     f"{rb:.1f}" if rb is not None else "",
                     ratio])
    return render_table(
        ["Subfamily", "Threshold (x100)", "RealSelfBits", "Threshold/Self"],
        rows, max_rows=10000,
        col_titles={k: _strip_html(v) for k, v in LEG_THRESHOLDS})


def flags_table(flags: Dict[str, Dict[str, int]]) -> str:
    rows = []
    for sf, d in sorted(flags.items(), key=lambda kv: -kv[1]["total"]):
        tot = d["total"] or 1
        rows.append([sf, f"{d['total']:,}", f"{d['OK']:,}",
                     f"{d['LEAK']:,}", f"{d['CONFLICT']:,}",
                     f"{100*d['CONFLICT']/tot:.1f}%"])
    return render_table(
        ["Subfamily", "Total", "OK", "LEAK", "CONFLICT", "%CONFLICT"],
        rows, max_rows=10000,
        col_titles={k: _strip_html(v) for k, v in LEG_FLAGS})


def build_alignment_section(
    species_code: str,
    subfams: List[str],
    msa_url: str = "https://toki-bio.github.io/MSA-viewer/",
    raw_base: str = "https://raw.githubusercontent.com/Toki-bio/Tal/main/",
) -> str:
    raw_aln    = f"{raw_base}{species_code}/alignments/"
    raw_subfam = f"{raw_base}{species_code}/subfam/"

    def msa_href(url: str, title: str) -> str:
        return (f"{msa_url}?url={quote(url, safe='')}"
                f"&title={quote(title, safe='')}")

    consi_href = msa_href(raw_aln + f"{species_code}_consensuses.fa",
                          f"{species_code} all consensi")
    subfam_input_href = msa_href(raw_aln + f"{species_code}_subfam_input.aln.fa",
                                 f"{species_code} SubFam input (all families)")

    rows_html = ""
    for sf in sorted(subfams):
        t100  = msa_href(raw_aln    + f"{species_code}_{sf}_top100.aln.fa",
                         f"{species_code} {sf} top100")
        r100  = msa_href(raw_aln    + f"{species_code}_{sf}_rand100.aln.fa",
                         f"{species_code} {sf} rand100")
        subfam_href = msa_href(raw_subfam + f"{sf}.al",
                               f"{species_code} {sf} SubFam")
        rows_html += (
            f"<tr><td><code>{html.escape(sf)}</code></td>"
            f"<td><a class='aln-link' href='{t100}' target='_blank'>"
            f"top 100 by score</a></td>"
            f"<td><a class='aln-link orange' href='{r100}' target='_blank'>"
            f"100 random</a></td>"
            f"<td><a class='aln-link green' href='{subfam_href}' target='_blank'>"
            f"SubFam</a></td></tr>"
        )
    return (
        "<section class='card' id='alignments'>"
        "<h2>Subfamily Alignments &mdash; open in MSA Viewer</h2>"
        "<p class='intro'>Copies re-extracted with "
        "<strong>50&thinsp;bp upstream + 70&thinsp;bp downstream</strong> "
        "genomic flanks (strand-aware).&nbsp;&nbsp;"
        f"<a class='aln-link' href='{consi_href}' target='_blank'>"
        "All consensi</a>"
        "&nbsp;&nbsp;"
        f"<a class='aln-link green' href='{subfam_input_href}' target='_blank'>"
        "SubFam input alignment (all families combined)</a></p>"
        "<table class='tbl'>"
        "<thead><tr><th>Subfamily</th>"
        "<th>Top 100 by bitscore</th>"
        "<th>100 random copies</th>"
        "<th>SubFam (chunk consensuses)</th></tr></thead>"
        f"<tbody>{rows_html}</tbody>"
        "</table></section>"
    )


# ===========================================================================
# Build
# ===========================================================================

def build_html(run_root: Path,
               out_path: Path,
               inline_plotly: bool,
               max_table_rows: int,
               embed_images: bool,
               sineplot: bool,
               sineplot_max: int,
               threads: int,
               tal_species_code: Optional[str] = None) -> None:
    LOG.info("Building report for %s", run_root)
    s2 = find_step2_out(run_root)
    LOG.info("step2 output: %s", s2)

    manifest = read_kv_manifest(run_root / "manifest.txt")
    summary_stats = parse_step2_summary(s2 / "summary.txt")
    step1_hits, step1_total = parse_step1_hits(
        run_root / "step1.stderr.log", run_root / "step1.stdout.log")
    if step1_total and "total" not in summary_stats:
        summary_stats["total"] = str(step1_total)

    stats_hdr, stats_rows = read_tsv(s2 / "assignment_stats.tsv")
    bysf_hdr, bysf_rows   = read_tsv(s2 / "summary.by_subfam.tsv")
    _, sbr_rows           = read_tsv(s2 / "self_bits_real.tsv",
                                     has_header=False)

    sim_by_sf = stratified_sample_sim(
        s2 / "sim_scores.tsv", s2 / "assignment_full.tsv", per_group=3000)

    flags = count_flags_per_subfam(s2 / "all_sines.bedlike.ALL.tsv")

    # Conservation curves from step4 companion TSVs.
    data_dir = s2 / "plots" / "data"
    curves: Dict[str, Tuple[List[int], List[float]]] = {}
    if data_dir.is_dir():
        for tsv in sorted(data_dir.glob("*_nucfreq.tsv")):
            sf = tsv.stem.replace("_nucfreq", "")
            nf = read_nucfreq_tsv(tsv)
            if nf is not None:
                curves[sf] = conservation_curve(nf)

    figs = {
        "div_kde":     fig_divergence_kde(sim_by_sf),
        "sim_violins": fig_sim_violins(sim_by_sf),
    }
    if curves:
        figs["conservation"] = fig_conservation(curves)

    # Per-subfamily PNG gallery
    plots_dir = s2 / "plots"
    image_blocks: List[str] = []
    if embed_images and plots_dir.is_dir():
        subfams = sorted({p.name.rsplit("_divergence.", 1)[0]
                          for p in plots_dir.glob("*_divergence.png")})
        for sf in subfams:
            div_uri = img_to_data_uri(plots_dir / f"{sf}_divergence.png")
            nuc_uri = img_to_data_uri(plots_dir / f"{sf}_nucfreq.png")
            block = [f"<div class='subfam-block'><h3>{html.escape(sf)}</h3>"]
            if div_uri:
                block.append(f"<img alt='{sf} divergence' src='{div_uri}'>")
            if nuc_uri:
                block.append(f"<img alt='{sf} nucfreq' src='{nuc_uri}'>")
            block.append("</div>")
            image_blocks.append("".join(block))

    # SINEplot iframe
    sineplot_html = None
    if sineplot:
        sd = manifest.get("SCRIPT_DIR", "")
        script_dir = Path(sd) if sd else None
        sineplot_html = build_sineplot_iframe(
            s2, script_dir, sineplot_max, threads)

    # Built-in PCA (mutation-space: mafft alignment columns)
    pca_fig = build_pca_fig(run_root, s2, n_per_sf=200, threads=threads)
    if pca_fig:
        # strip the internal _meta key before JSON serialisation
        pca_fig["layout"].pop("_meta", None)
        figs["pca"] = pca_fig
    pca_section_html = (
        "<p class='intro'>Each point is one assigned SINE copy, coloured by "
        "subfamily. Axes are PC1&thinsp;/&thinsp;PC2 of a binary mutation "
        "matrix: each alignment column where the copy nucleotide differs from "
        "its assigned-subfamily consensus is encoded as&thinsp;1, matches as&thinsp;0. "
        "Only variable columns (mutation rate 3&ndash;97&thinsp;%) are used. "
        "Copies with distinct mutation patterns cluster together &mdash; "
        "subfamily clouds are expected when subfamilies have different "
        "characteristic mutations.</p>"
        "<div class='plot' id='plot_pca'></div>"
        "<p class='small muted'>Method: <code>mafft --add --keeplength</code> "
        "adds sampled copies into the consensus alignment (up to 200 per "
        "subfamily). Binary mutation matrix &rarr; SVD. "
        "&eta;&sup2;(subfamily&rarr;PC1) = fraction of PC1 variance explained "
        "by subfamily label (0&thinsp;=&thinsp;no separation, "
        "1&thinsp;=&thinsp;perfect).</p>"
    ) if pca_fig else (
        "<div class='skip-note'><b>PCA not available.</b> "
        "Requires <code>mafft</code> in PATH and "
        "<code>numpy</code>. Run with <code>--verbose</code> "
        "for details.</div>"
    )

    # Plotly script
    plotly_tag, plotly_mode = get_plotly_js(inline_plotly)
    fig_init = "\n".join(
        f"Plotly.newPlot('plot_{name}', "
        f"{json.dumps(spec['data'])}, {json.dumps(spec['layout'])}, "
        f"{{responsive: true, displaylogo: false}});"
        for name, spec in figs.items()
    )

    # Overview prose
    def _fmt(n):
        try:
            return f"{int(n):,}"
        except Exception:
            return str(n)

    def _ipct(num, den):
        try:
            n_, d_ = int(num), int(den)
            if d_ <= 0:
                return ""
            return f" ({100.0 * n_ / d_:.1f}%)"
        except Exception:
            return ""

    n_subfam = len(stats_rows)
    genome_in_full = manifest.get("GENOME_IN", "")
    genome_basename = (Path(genome_in_full).name
                       if genome_in_full else "")
    species = manifest.get("SPECIES") or manifest.get("GENOME_LABEL") \
        or genome_basename or "this genome"
    cons_basename = (Path(manifest.get("CONS_IN", "")).name
                     if manifest.get("CONS_IN") else "")
    total_v      = summary_stats.get("total", "")
    unan_v       = summary_stats.get("unanimous", "")
    assigned_v   = summary_stats.get("assigned", "")
    unassigned_v = summary_stats.get("unassigned", "")

    # Try to find the merged-hit total (after bedtools merge)
    merged_hits = None
    try:
        if int(total_v) > 0:
            merged_hits = int(total_v)
    except Exception:
        pass

    overview_paragraph = (
        f"<p>This run searched <b>{html.escape(species)}</b> with "
        f"<b>{n_subfam}</b> subfamily "
        f"consensus{'es' if n_subfam != 1 else ''}"
        + (f" from <code>{html.escape(cons_basename)}</code>"
           if cons_basename else "")
        + ". "
        + (f"After merging overlapping hits across queries, "
           f"<b>{_fmt(merged_hits)}</b> candidate SINE copies entered "
           f"step2. " if merged_hits is not None else "")
        + (f"Of those, <b>{_fmt(unan_v)}</b>{_ipct(unan_v, total_v)} "
           "were unanimously called by all 10 sub-samples, "
           if unan_v else "")
        + (f"<b>{_fmt(assigned_v)}</b>{_ipct(assigned_v, total_v)} "
           "passed the bitscore threshold and received a final "
           "subfamily label, " if assigned_v else "")
        + (f"and <b>{_fmt(unassigned_v)}</b>"
           f"{_ipct(unassigned_v, total_v)} remained unassigned."
           if unassigned_v else "")
        + "</p>"
    )

    mani_html = "<div class='kv'>" + "".join(
        f"<div class='k'>{html.escape(k)}</div>"
        f"<div class='v'>{html.escape(v)}</div>"
        for k, v in manifest.items()
    ) + "</div>"

    # Tables (raw TSV rows)
    stats_table = render_table(
        stats_hdr, stats_rows, max_table_rows,
        col_titles={k: _strip_html(v) for k, v in LEG_ASSIGN_STATS})
    bysf_table = render_bysf_table(bysf_hdr, bysf_rows, max_table_rows)
    funnel_table = funnel_html(summary_stats)

    # Conservation panel HTML (per-position nucfreq when available)
    if curves:
        conservation_html = (
            "<div class='plot' id='plot_conservation'></div>"
            "<p class='small muted'>Per-position conservation from "
            "<code>plots/data/*_nucfreq.tsv</code>.</p>"
        )
    else:
        conservation_html = ""

    # Alignment section (optional, requires --tal-species-code)
    alignment_section = ""
    if tal_species_code:
        subfams_for_aln = sorted({r[0] for r in stats_rows if r})
        alignment_section = build_alignment_section(
            tal_species_code, subfams_for_aln)

    # SINEplot panel HTML
    if sineplot_html:
        srcdoc = (sineplot_html.replace("&", "&amp;")
                                .replace('"', "&quot;"))
        sineplot_section = (
            f"<iframe class='embed' srcdoc=\"{srcdoc}\" "
            "sandbox='allow-scripts allow-same-origin allow-popups' "
            "title='SINEplot PCA'></iframe>"
        )
    elif sineplot:
        sineplot_section = (
            "<div class='skip-note'>"
            "<b>SINEplot PCA skipped.</b><br>"
            "Could not locate <code>SINEplot.py</code> "
            "(checked <code>$SINEPLOT_PY</code>, "
            "<code>$SCRIPT_DIR/SINEplot/SINEplot.py</code>, and "
            "<code>$PATH</code>) or <code>ssearch36</code> "
            "is unavailable. To enable: "
            "<pre><code>git clone https://github.com/Toki-bio/SINEplot \\\n"
            "  $SCRIPT_DIR/SINEplot\n"
            "pip install pandas numpy plotly scikit-learn</code></pre>"
            "Then re-run <code>step6_report.sh</code>.</div>"
        )
    else:
        sineplot_section = (
            "<p class='small muted'>SINEplot panel disabled "
            "(<code>--no-sineplot</code>).</p>"
        )

    gallery_html = ""
    if image_blocks:
        gallery_html = (
            "<section class='card' id='gallery'>"
            "<h2>Per-subfamily diagnostic plots (step4)</h2>"
            f"<p class='intro'>{len(image_blocks)} subfamilies. "
            "Each block: divergence histogram (top) + nucleotide frequency "
            "stacked bar (bottom). Embedded as base64 PNGs.</p>"
            "<div class='subfam-grid'>"
            + "".join(image_blocks) +
            "</div></section>"
        )

    title = manifest.get("RUN", str(run_root)).rstrip("/").split("/")[-1]
    genome_name = manifest.get("GENOME_IN", "?").rstrip("/").split("/")[-1]
    cons_name = manifest.get("CONS_IN", "?").rstrip("/").split("/")[-1]
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Cross-species nav bar (only when a species code is provided)
    if tal_species_code:
        _others = {
            "saq": ("../ccr/report.html", "ccr"),
            "ccr": ("../saq/report.html", "saq"),
        }
        _other_href, _other_label = _others.get(
            tal_species_code, ("", ""))
        other_link = (
            f"<a href='{_other_href}' "
            f"style='color:rgba(255,255,255,.7);text-decoration:none;"
            f"margin-right:16px;'>{_other_label}</a>"
            if _other_href else ""
        )
        cross_nav = (
            "<div style='background:#1a2634;color:rgba(255,255,255,.8);"
            "padding:6px 32px;font-size:.85rem;'>"
            "<a href='../index.html' style='color:rgba(255,255,255,.85);"
            "text-decoration:none;margin-right:16px;'>&#8592; All species</a>"
            + other_link
            + "<a href='https://github.com/Toki-bio/Tal' target='_blank' "
            "style='color:rgba(255,255,255,.6);text-decoration:none;'>"
            "GitHub</a></div>"
        )
    else:
        cross_nav = ""

    html_doc = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SINEderella report &mdash; {html.escape(title)}</title>
<style>{CSS}</style>
{plotly_tag}
</head><body>
<div id="lightbox" onclick="this.classList.remove('active')">
  <img id="lightbox-img" src="" alt="">
</div>
{cross_nav}
<header>
  <h1>SINEderella run report &mdash; {html.escape(title)}</h1>
  <div class="sub">Genome: <code>{html.escape(genome_name)}</code> &middot;
       Consensus: <code>{html.escape(cons_name)}</code> &middot;
       Generated {generated} &middot; plotly.js: {plotly_mode}</div>
</header>
<main>
  <nav class="toc">
    <strong>Sections:</strong>
    {'<a href="#alignments">Alignments</a>' if alignment_section else ''}
    <a href="#overview">Overview</a>
    <a href="#composition">Composition</a>
    <a href="#divergence">Divergence&thinsp;/&thinsp;Similarity</a>
    <a href="#pca">PCA</a>
    {'<a href="#gallery">Gallery</a>' if image_blocks else ''}
    <span style="opacity:.45">&bull;</span>
    <a href="#run" style="opacity:.7">Run info</a>
    <a href="#step1" style="opacity:.7">Step1</a>
    <a href="#assignment" style="opacity:.7">Step2</a>
    <a href="#thresholds" style="opacity:.7">Thresholds</a>
    <a href="#flags" style="opacity:.7">Flags</a>
  </nav>

  {alignment_section}

  <section class="card" id="overview">
    <h2>Overview</h2>
    {overview_paragraph}
    <h3 style="margin-top:16px">Pipeline counts</h3>
    {funnel_table}
  </section>

  <section class="card" id="composition">
    <h2>Subfamily composition</h2>
    <p class="intro">Source: <code>summary.by_subfam.tsv</code>.
    Highlighted columns (<span style="font-weight:600;color:#1a3a6e">bold</span>):
    total copies, similarity mean&thinsp;/&thinsp;median, conflict rate.
    Hover a column header for its full description.</p>
    {bysf_table}
  </section>

  <section class="card" id="divergence">
    <h2>Divergence from consensus &mdash; per copy</h2>
    <p class="intro"><b>Metric:</b> bitscore-based divergence =
    100&thinsp;&minus;&thinsp;(copy bitscore&thinsp;/&thinsp;consensus self-bitscore &times; 100&thinsp;%).
    This is a proxy for sequence divergence, not a direct nucleotide count.
    Values are clamped to 0 (local alignment can occasionally score a copy
    <i>above</i> the self-bitscore, which would otherwise appear as negative
    divergence). One KDE curve per subfamily (up to 3,000 copies per
    subfamily sampled for display); click legend entries to toggle.</p>
    {conservation_html}
    <div class="plot" id="plot_div_kde"></div>
    <h3>Distributions per subfamily (violin)</h3>
    <p class="intro">Same data as above shown as symmetric violin plots.
    Width encodes density; inner box = IQR; dashed line = mean.</p>
    <div class="plot" id="plot_sim_violins"></div>
    <p class="small muted">Source: <code>sim_scores.tsv</code> joined with
    <code>assignment_full.tsv</code>.</p>
  </section>

  <section class="card" id="pca">
    <h2>PCA &mdash; mutation landscape</h2>
    {pca_section_html}
  </section>

  {gallery_html}

  <details class="card" id="run">
    <summary><h2>Run info (<code>manifest.txt</code>)</h2></summary>
    <div class="card-body">{mani_html}</div>
  </details>

  <details class="card" id="step1">
    <summary><h2>Step1 &mdash; raw hits per query consensus</h2></summary>
    <div class="card-body">
      <p class="intro">Raw <code>ssearch36</code> hit counts per query consensus,
      before <code>bedtools merge</code>. The pre-merge total is the sum across
      all queries; after merging, overlapping intervals are collapsed.</p>
      {render_legend(LEG_STEP1_HITS)}
      {step1_hits_table(step1_hits)}
    </div>
  </details>

  <details class="card" id="assignment">
    <summary><h2>Step2 &mdash; assignment stats per subfamily</h2></summary>
    <div class="card-body">
      <p class="intro">Output of <code>step2_asSINEment.sh</code>
      (<code>assignment_stats.tsv</code>). <i>Assigned</i> = copies that
      passed the bitscore threshold.</p>
      {render_legend(LEG_ASSIGN_STATS)}
      {stats_table}
    </div>
  </details>

  <details class="card" id="thresholds">
    <summary><h2>Bitscore thresholds vs consensus self-bits</h2></summary>
    <div class="card-body">
      <p class="intro">Threshold is stored as integer bits&times;100; real
      self-bits = ssearch36 score of the consensus vs itself.
      Threshold/Self &asymp; minimum fractional similarity to pass.</p>
      {render_legend(LEG_THRESHOLDS)}
      {thresholds_table(stats_rows, sbr_rows)}
    </div>
  </details>

  <details class="card" id="flags">
    <summary><h2>Quality flags per subfamily</h2></summary>
    <div class="card-body">
      <p class="intro">Merged genomic intervals grouped by sanity flag from
      <code>step3_postprocess.sh</code>. High <code>%CONFLICT</code> suggests
      ambiguous subfamily boundaries at those loci.</p>
      {render_legend(LEG_FLAGS)}
      {flags_table(flags)}
    </div>
  </details>

</main>
<footer>
  Generated by <code>step6_report.py</code> &middot;
  SINEderella pipeline &middot; {generated}
</footer>
<script>
{fig_init}
(function(){{
  var lb = document.getElementById('lightbox');
  var lbimg = document.getElementById('lightbox-img');
  document.querySelectorAll('.subfam-block img').forEach(function(img){{
    img.addEventListener('click', function(e){{
      e.stopPropagation();
      lbimg.src = img.src;
      lbimg.alt = img.alt;
      lb.classList.add('active');
    }});
  }});
  document.addEventListener('keydown', function(e){{
    if (e.key === 'Escape') lb.classList.remove('active');
  }});
}})();
</script>
</body></html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_doc, encoding="utf-8")
    size_mb = out_path.stat().st_size / (1024 * 1024)
    LOG.info("Wrote %s (%.2f MB)", out_path, size_mb)


# ===========================================================================
# CLI
# ===========================================================================

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_root", type=Path,
                    help="Path to a completed SINEderella run directory")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output HTML path "
                         "(default: <RUN_ROOT>/results/report.html)")
    ap.add_argument("--inline-plotly", dest="inline_plotly",
                    action="store_true", default=False,
                    help="Inline plotly.js (~4 MB larger but fully self-contained).")
    ap.add_argument("--no-inline-plotly", dest="inline_plotly",
                    action="store_false",
                    help="Use plotly.js from CDN (default; requires internet).")
    ap.add_argument("--no-embed-images", dest="embed_images",
                    action="store_false", default=True,
                    help="Do not embed step4 PNGs (smaller HTML).")
    ap.add_argument("--max-rows", type=int, default=200,
                    help="Max rows shown in tables (default: 200).")
    ap.add_argument("--sineplot", dest="sineplot",
                    action="store_true", default=True,
                    help="Run SINEplot if available (default on).")
    ap.add_argument("--no-sineplot", dest="sineplot",
                    action="store_false",
                    help="Skip the SINEplot PCA panel.")
    ap.add_argument("--sineplot-max", type=int, default=400,
                    help="Max copies per subfamily fed to SINEplot "
                         "(default: 400).")
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("THREADS",
                                               os.cpu_count() or 1)),
                    help="Threads for ssearch36 inside SINEplot stage.")
    ap.add_argument("--tal-species-code", default=None,
                    help="Species code (e.g. 'saq', 'ccr') for Tal GitHub "
                         "Pages alignment links and cross-species nav bar.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%H:%M:%S")

    run_root = args.run_root.resolve()
    if not run_root.is_dir():
        ap.error(f"RUN_ROOT does not exist: {run_root}")

    out_path = args.out or (run_root / "results" / "report.html")
    build_html(run_root, out_path,
               inline_plotly=args.inline_plotly,
               max_table_rows=args.max_rows,
               embed_images=args.embed_images,
               sineplot=args.sineplot,
               sineplot_max=args.sineplot_max,
               threads=args.threads,
               tal_species_code=args.tal_species_code)
    return 0


if __name__ == "__main__":
    sys.exit(main())
