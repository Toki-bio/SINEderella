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

Stdlib only (urllib used to fetch plotly.js once and cache), except the
per-subfamily alignment-composition diagram (--profile-diagrams), which needs
numpy — imported lazily in build_alignment_section() so the rest of the report
still builds on a bare-stdlib install if numpy is missing.
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

SORT_JS = r'''<script>
// every table.tbl with a header row sorts by a clicked column (numbers by value, others as text; click again to reverse)
document.querySelectorAll('table.tbl').forEach(function(t){
  var hr=t.tHead?t.tHead.rows[0]:null; if(!hr)return;
  Array.prototype.forEach.call(hr.cells,function(th,ci){
    th.style.cursor='pointer'; th.title='click to sort';
    th.addEventListener('click',function(){
      var tb=t.tBodies[0],rows=Array.prototype.slice.call(tb.rows),dir=th.dataset.dir==='asc'?-1:1;
      Array.prototype.forEach.call(hr.cells,function(o){delete o.dataset.dir;});
      th.dataset.dir=dir===1?'asc':'desc';
      var val=function(r){var x=r.cells[ci]?r.cells[ci].textContent.trim().replace(/,/g,'').replace(/%$/,''):'';var n=parseFloat(x);return isNaN(n)?x.toLowerCase():n;};
      rows.sort(function(a,b){var x=val(a),y=val(b);if(typeof x==='number'&&typeof y==='number')return dir*(x-y);if(typeof x==='number')return -1;if(typeof y==='number')return 1;return dir*(x<y?-1:x>y?1:0);});
      rows.forEach(function(r){tb.appendChild(r);});
    });
  });
});
</script>
'''

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


def run_inputs(manifest: Dict[str, str]) -> Tuple[str, str]:
    """(genome file name, consensus bank name) of a run. An --add / --exclude run's manifest names only
    SOURCE_RUN (and ADD_FILE), so the header said "Genome: ? - Consensus: ?" on every add run (all 25
    bat pages, 2026-09-28): follow SOURCE_RUN back to the full run and list every added bank."""
    genome, cons, added, m, seen = "", "", [], manifest, set()
    while m:
        if m.get("ADD_FILE"):
            added.insert(0, m["ADD_FILE"].rstrip("/").split("/")[-1])
        genome = genome or m.get("GENOME_IN", "")
        cons = cons or m.get("CONS_IN", "")
        src = m.get("SOURCE_RUN", "")
        if (genome and cons) or not src or src in seen:
            break
        seen.add(src)
        m = read_kv_manifest(Path(src) / "manifest.txt")
    genome = genome.rstrip("/").split("/")[-1] or "?"
    cons = " + ".join([cons.rstrip("/").split("/")[-1] or "?"] + added)
    return genome, cons


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


def load_pctid_by_sf(plots_dir: Path) -> Dict[str, List[float]]:
    """step4 *_pctid.tsv: col2 = ssearch36 %identity to consensus."""
    by_sf: Dict[str, List[float]] = {}
    if not plots_dir.is_dir():
        return by_sf
    for p in sorted(plots_dir.glob("*_pctid.tsv")):
        sf = p.name.replace("_pctid.tsv", "")
        vals = []
        with p.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        vals.append(float(parts[1]))
                    except ValueError:
                        pass
        if vals:
            by_sf[sf] = vals
    return by_sf


def bin_divergence(vals: List[float], bin_width: float = 1.0) -> Dict[float, int]:
    counts: Dict[float, int] = {}
    for v in vals:
        d = max(0.0, 100.0 - float(v))
        bin_start = math.floor((d + 1e-9) / bin_width) * bin_width
        bin_start = round(bin_start, 6)
        counts[bin_start] = counts.get(bin_start, 0) + 1
    return counts


def fig_pctid_spline_divergence(by_sf: Dict[str, List[float]],
                                bin_width: float = 1.0) -> dict:
    """1% bin midpoints, Y = copy count, smooth spline (step4 gallery metric)."""
    sf_sorted = sorted(by_sf.keys())
    traces = []
    max_bin_end = 0.0
    for i, sf in enumerate(sf_sorted):
        counts = bin_divergence(by_sf[sf], bin_width)
        if not counts:
            continue
        bins = sorted(counts)
        max_bin_end = max(max_bin_end, max(bins) + bin_width)
        traces.append({
            "type": "scatter",
            "mode": "lines",
            "x": [round(b + bin_width / 2.0, 1) for b in bins],
            "y": [counts[b] for b in bins],
            "name": sf,
            "line": {
                "color": SF_PALETTE[i % len(SF_PALETTE)],
                "width": 2,
                "shape": "spline",
            },
            "hovertemplate": (
                "%{fullData.name}<br>divergence ~%{x:.0f}%"
                "<br>copies %{y:,d}<extra></extra>"),
        })
    x_range_max = min(100.0, max(5.0, math.ceil(max_bin_end / 5.0) * 5.0))
    return {
        "data": traces,
        "layout": {
            "title": "ssearch36 %identity divergence (step4)",
            "xaxis": {
                "title": "Divergence (100 \u2212 %identity to consensus)",
                "range": [0, x_range_max],
                "dtick": 5,
            },
            "yaxis": {"title": "Copies", "rangemode": "tozero"},
            "legend": {"title": {"text": "Subfamily (click to toggle)"}},
            "height": 460,
            "margin": {"t": 60, "r": 20, "b": 60, "l": 70},
        },
    }


DIAGRAM_JS = r"""
// Alignment-composition diagram + structural-feature track, ported in full
// this time from SINE_discriminator's site/index.html (drawProfile,
// buildTrackLegend, FEATON/trackOn, the SHOW chip list) and annotate.py
// (via report_annotate.py / ANNOTATIONS below), 2026-09-09. Colors, labels,
// checkboxes and chip formatting are the ORIGINAL's, not a paraphrase --
// this report should read identically to the site's per-family view.
const TRACKS = [
  { k: "pair_id", label: "pairwise identity between copies", color: "var(--t1)", on: true },
  { k: "cons_id", label: "identity to known consensus", color: "var(--t2)", on: true },
  { k: "cover", label: "copies present (coverage)", color: "var(--t3)", on: true },
  { k: "at", label: "A+T fraction", color: "var(--t4)", on: false }
];
const BG_LEVEL = 0.25;
let TOTAL_H = 0, MOTIF_H = 0;
const FEATCOLOR = {
  abox: "var(--f1)", bbox: "var(--f2)", trna_region: "var(--f1)",
  conserved_core: "var(--f4)", simple_repeat: "var(--f5)",
  tail_repeat: "var(--f6)", internal_dup: "var(--f7)", tsd: "var(--f8)",
  terminator: "var(--f3)"
};
const FEATALPHA = {trna_region: 0.28, terminator: 0.85};
const FEATLIST = [
  ["abox", "A box"], ["bbox", "B box"], ["trna_region", "tRNA region"],
  ["conserved_core", "conserved block"], ["simple_repeat", "simple repeat"],
  ["tail_repeat", "3' tail repeat"], ["internal_dup", "internal duplication"],
  ["tsd", "TSD"], ["terminator", "Pol III terminator"],
  ["abox_p", "A box match profile"], ["bbox_p", "B box match profile"],
  ["selfsim_p", "similarity to tRNA head"], ["at_heat", "A+T composition"]
];
const SHOW = ["cliff", "cons_identity_med", "frac_supported", "elem_len_cv", "frac_full",
              "tsd_frac", "tsd_len_med", "res_asymmetry", "rank1_excess", "flank_id",
              "elem_len_med", "cons_bp", "n_copies"];

function pathFor(xs, ys, sx, sy) {
  let d = "", pen = false;
  for (let i = 0; i < xs.length; i++) {
    if (ys[i] === null || ys[i] === undefined) { pen = false; continue; }
    const X = sx(xs[i]), Y = sy(ys[i]);
    d += (pen ? "L" : "M") + X.toFixed(1) + " " + Y.toFixed(1) + " ";
    pen = true;
  }
  return d;
}

function heatRow(xs, vals, sx, y, h, opts) {
  const {lo, hi, hue, diverge, mid, label, unit} = opts;
  let out = "";
  for (let i = 0; i < xs.length - 1; i++) {
    const v = vals[i];
    if (v === null || v === undefined) continue;
    const x0 = sx(xs[i]), x1 = sx(xs[i + 1]);
    let col, a;
    if (diverge) {
      const t = (v - mid) / (v >= mid ? (hi - mid) : (mid - lo));
      a = Math.min(1, Math.abs(t));
      col = v >= mid ? "var(--div-hi)" : "var(--div-lo)";
      if (a < 0.08) { col = "var(--div-mid)"; a = 0.85; }
    } else {
      a = Math.min(1, Math.max(0, (v - lo) / (hi - lo)));
      a = Math.pow(a, opts.gamma || 1);
      col = hue;
      if (a <= 0.01) continue;
    }
    out += '<rect x="' + x0.toFixed(1) + '" y="' + y + '" width="' +
           Math.max(0.8, x1 - x0).toFixed(1) + '" height="' + h +
           '" fill="' + col + '" fill-opacity="' + a.toFixed(2) + '">' +
           '<title>' + (label || "") + " at position " + xs[i] + ": " +
           v.toFixed(2) + (unit || "") + '</title></rect>';
  }
  return out;
}

function FEATON(k) {
  const cb = document.querySelector(".ftr[data-k='" + k + "']");
  return !cb || cb.checked;
}
function trackOn(k) {
  const cb = document.querySelector(".trk[data-k='" + k + "']");
  return !cb || cb.checked;
}

function buildDiagLegend(host) {
  host.innerHTML =
    "<div class='diag-legend-group'><span class='diag-legend-title'>Profile tracks</span>" +
    TRACKS.map(t =>
      "<label><input type='checkbox' class='trk' data-k='" + t.k + "'" + (t.on ? " checked" : "") +
      "><span class='sw' style='background:" + t.color + "'></span>" + t.label + "</label>").join("") +
    "<label style='cursor:default'><span class='sw' style='background:none;" +
    "border-top:2px dashed var(--warn)'></span>unrelated DNA, 0.25</label></div>" +
    "<div class='diag-legend-group'><span class='diag-legend-title'>Structural features</span>" +
    FEATLIST.map(([k, lab]) =>
      "<label><input type='checkbox' class='ftr' data-k='" + k + "' checked>" +
      "<span class='swb' style='background:" + (FEATCOLOR[k] ||
        {abox_p: "var(--f1)", bbox_p: "var(--f2)", selfsim_p: "var(--f3)",
         at_heat: "linear-gradient(90deg,var(--div-lo),var(--div-mid),var(--div-hi))"}[k] ||
        "var(--muted)") + "'></span>" + lab + "</label>").join("") + "</div>";
}

function diagChipsHTML(setName) {
  const p = PROFILES[setName], m = (p && p._measure) || {};
  return SHOW.filter(k => k in m).map(k => {
    let v = m[k];
    if (v === null || v === undefined) return "";
    v = (Math.abs(v) >= 100 || Number.isInteger(v)) ? v : v.toFixed(3);
    return "<span class='diag-chip'>" + k + " <b>" + v + "</b></span>";
  }).join("");
}

function drawProfile(setName, host) {
  const p = PROFILES[setName];
  if (!p) { host.innerHTML = "<p class='small muted'>No profile for this set.</p>"; return; }
  const W = 900, H = 190, PADL = 82, PADR = 12, PADT = 12, PADB = 24;
  const MH = 48;
  const xs = p.x, L = p.elem_len;
  const x0 = xs[0], x1 = xs[xs.length - 1];
  const sx = v => PADL + (v - x0) / (x1 - x0) * (W - PADL - PADR);
  const sy = v => PADT + (1 - v) * (H - PADT - PADB);

  let s = "<div class='diag-chips'>" + diagChipsHTML(setName) + "</div>" +
          '<svg viewBox="0 0 ' + W + ' ' + (H + MH) + '" width="100%" role="img" ' +
          'aria-label="positional profile for ' + setName + '">';

  s += '<rect x="' + sx(0) + '" y="' + PADT + '" width="' + (sx(L - 1) - sx(0)) +
       '" height="' + (H - PADT - PADB) + '" fill="var(--band)"/>';
  [0, 0.25, 0.5, 0.75, 1].forEach(g => {
    s += '<line x1="' + PADL + '" y1="' + sy(g) + '" x2="' + (W - PADR) + '" y2="' + sy(g) +
         '" stroke="var(--grid)" stroke-width="1"/>' +
         '<text x="' + (PADL - 6) + '" y="' + (sy(g) + 3.5) + '" text-anchor="end" ' +
         'font-size="9" fill="var(--muted)" font-family="monospace">' + g.toFixed(2) + '</text>';
  });
  s += '<line x1="' + PADL + '" y1="' + sy(BG_LEVEL) + '" x2="' + (W - PADR) + '" y2="' +
       sy(BG_LEVEL) + '" stroke="var(--warn)" stroke-width="1" stroke-dasharray="4 3"/>';
  [0, L - 1].forEach(e => {
    s += '<line x1="' + sx(e) + '" y1="' + PADT + '" x2="' + sx(e) + '" y2="' + (H - PADB) +
         '" stroke="var(--ink)" stroke-width="1.2"/>';
  });
  s += '<text x="' + sx(0) + '" y="' + (H - PADB + 13) + '" font-size="9.5" ' +
       'fill="var(--ink-2)" font-family="monospace">0</text>' +
       '<text x="' + sx(L - 1) + '" y="' + (H - PADB + 13) + '" text-anchor="end" font-size="9.5" ' +
       'fill="var(--ink-2)" font-family="monospace">' + (L - 1) + ' bp</text>' +
       '<text x="' + PADL + '" y="' + (H - PADB + 13) + '" font-size="9.5" ' +
       'fill="var(--muted)" font-family="monospace">5&#8242; flank</text>' +
       '<text x="' + (W - PADR) + '" y="' + (H - PADB + 13) + '" text-anchor="end" font-size="9.5" ' +
       'fill="var(--muted)" font-family="monospace">3&#8242; flank</text>';

  for (let i = 0; i < xs.length; i += Math.max(1, Math.round(xs.length / 220))) {
    const xv = xs[i];
    const parts2 = TRACKS.filter(t => trackOn(t.k))
      .map(t => p[t.k] && p[t.k][i] !== null && p[t.k][i] !== undefined
           ? t.label + " " + p[t.k][i].toFixed(3) : null).filter(Boolean);
    const where = xv < 0 ? (-xv) + " bp into the 5' flank"
                : xv >= L ? (xv - L + 1) + " bp into the 3' flank"
                : "element position " + xv;
    s += '<rect x="' + (sx(xv) - 2).toFixed(1) + '" y="' + PADT + '" width="5" height="' +
         (H - PADT - PADB) + '" fill="transparent"><title>' + where +
         String.fromCharCode(10) + parts2.join(String.fromCharCode(10)) +
         '</title></rect>';
  }

  TRACKS.forEach(t => {
    if (!trackOn(t.k)) return;
    s += '<path d="' + pathFor(xs, p[t.k], sx, sy) + '" fill="none" stroke="' + t.color +
         '" stroke-width="1.6" stroke-linejoin="round"/>';
  });

  // ---- mosaic strip -------------------------------------------------------
  const my0 = H + 6, mh = MH - 22;
  const mvals = (p.mosaic || []).filter(v => v !== null);
  const msort = mvals.slice().sort((a, b) => a - b);
  const q = f => msort.length ? msort[Math.min(msort.length - 1,
                 Math.max(0, Math.round(f * (msort.length - 1))))] : 0;
  let mlo = q(0.02), mhi = q(0.98);
  if (mhi - mlo < 0.05) { mhi = mlo + 0.05; }
  const pad = (mhi - mlo) * 0.15;
  mlo -= pad; mhi += pad;
  const smy = v => my0 + (1 - (v - mlo) / (mhi - mlo)) * mh;
  s += '<rect x="' + sx(0) + '" y="' + my0 + '" width="' + (sx(L - 1) - sx(0)) +
       '" height="' + mh + '" fill="var(--band)"/>';
  s += '<path d="' + pathFor(p.mosaic_x || [], p.mosaic || [], sx, smy) + '" fill="none" ' +
       'stroke="var(--t5)" stroke-width="1.6"/>';
  s += '<text x="' + (PADL - 6) + '" y="' + (my0 + 9) + '" text-anchor="end" font-size="9" ' +
       'fill="var(--muted)" font-family="monospace">' + mhi.toFixed(2) + '</text>' +
       '<text x="' + (PADL - 6) + '" y="' + (my0 + mh) + '" text-anchor="end" font-size="9" ' +
       'fill="var(--muted)" font-family="monospace">' + mlo.toFixed(2) + '</text>';
  s += '<text x="' + PADL + '" y="' + (my0 + mh + 13) + '" font-size="9.5" fill="var(--t5)" ' +
       'font-family="monospace">mosaicism &mdash; residual after rank-1 fit (high = a copy subset ' +
       'carries its own block here)</text>';

  // ---- motif score strip: every position scored, not just the best hit ----
  const MOTIFS = [["abox_p", "A box", "var(--f1)"],
                  ["bbox_p", "B box", "var(--f2)"],
                  ["selfsim_p", "tRNA rpt", "var(--f3)"]];
  const anyMotif = MOTIFS.some(([k]) => p[k] && FEATON(k));
  let motifH = 0;
  if (anyMotif) {
    const y0 = H + MH + 10, mh2 = 46;
    motifH = mh2 + 22;
    const RH = 11, GAP2 = 3;
    let yy = y0;
    MOTIFS.forEach(([k, lab, col]) => {
      if (!p[k] || !FEATON(k)) return;
      s += '<rect x="' + sx(0) + '" y="' + yy + '" width="' + (sx(L - 1) - sx(0)) +
           '" height="' + RH + '" fill="var(--band)"/>';
      s += heatRow(p.motif_x || [], p[k], sx, yy, RH, {lo: 0, hi: 5, hue: col, gamma: 2.2,
                   label: lab + ' match', unit: ' (z vs shuffled consensus)'});
      s += '<text x="' + (PADL - 6) + '" y="' + (yy + RH - 2) + '" text-anchor="end" ' +
           'font-size="8.5" fill="var(--ink-2)" font-family="monospace">' + lab + '</text>';
      yy += RH + GAP2;
    });
    if (p.at && FEATON("at_heat")) {
      const K = 9, sm = [];
      for (let i = 0; i < p.at.length; i++) {
        let t = 0, c2 = 0;
        for (let j = Math.max(0, i - (K >> 1)); j <= Math.min(p.at.length - 1, i + (K >> 1)); j++) {
          if (p.at[j] !== null && p.at[j] !== undefined) { t += p.at[j]; c2++; }
        }
        sm.push(c2 ? t / c2 : null);
      }
      s += heatRow(xs, sm, sx, yy, RH, {lo: 0.35, hi: 0.85, mid: 0.62, diverge: true,
                     label: 'A+T fraction', unit: ' (genomic background 0.62)'});
      s += '<text x="' + (PADL - 6) + '" y="' + (yy + RH - 2) + '" text-anchor="end" ' +
           'font-size="8.5" fill="var(--ink-2)" font-family="monospace">A+T</text>';
      yy += RH + GAP2;
    }
    motifH = (yy - y0) + 14;
    s += '<text x="' + PADL + '" y="' + (yy + 9) + '" font-size="9" ' +
         'fill="var(--muted)" font-family="monospace">motif rows: darker = higher z vs a ' +
         'shuffled consensus (faint = chance). A+T: blue GC-rich, orange AT-rich, ' +
         'grey = background</text>';
  }
  MOTIF_H = motifH;

  // ---- structural feature track (from ANNOTATIONS, ported annotate.py) ----
  const ann = (typeof ANNOTATIONS !== "undefined") ? ANNOTATIONS[setName] : null;
  if (ann && ann.features && ann.features.length) {
    const ay0 = H + MH + 4 + MOTIF_H, ah = 8;
    const spans = ann.features.filter(f => f.start !== undefined && f.end !== undefined)
                              .filter(f => FEATON(f.type))
                              .slice().sort((a, b) => a.start - b.start);
    const CH = 5.5;
    const lanes = [];
    spans.forEach(f => {
      const x0 = sx(f.start);
      const wpx = Math.max(2.5, sx(f.end) - x0);
      const lab = f.type.length * CH + 12;
      f._x = x0; f._w = wpx;
      f._inside = wpx > lab + 6;
      f._hangL = !f._inside && (x0 + wpx + lab > W - PADR);
      f._x0eff = x0 - (f._hangL ? lab : 0);
      f._x1 = x0 + wpx + (f._inside || f._hangL ? 0 : lab);
      let ln = 0;
      while (lanes[ln] !== undefined && lanes[ln] > f._x0eff - 10) ln++;
      lanes[ln] = f._x1;
      f._lane = ln;
    });
    const nlane = Math.max(1, lanes.length);
    if (!spans.length && !FEATON("tsd")) { TOTAL_H = H + MH; }
    spans.forEach(f => {
      const c = FEATCOLOR[f.type] || "var(--muted)";
      const x = f._x, w = f._w;
      const y = ay0 + f._lane * (ah + 3);
      const op = (f.type === "terminator" && f["class"] === "moderate") ? 0.4
               : (FEATALPHA[f.type] !== undefined ? FEATALPHA[f.type] : 0.75);
      s += '<rect x="' + x + '" y="' + y + '" width="' + w + '" height="' + ah +
           '" rx="2" fill="' + c + '" fill-opacity="' + op + '"><title>' +
           f.type + (f["class"] ? " (" + f["class"] + ")" : "") + " " + f.start + "-" + f.end +
           (f.unit ? " (" + f.unit + ")x" + f.units : "") +
           (f.mismatches !== undefined ? " " + f.seq + ", " + f.mismatches +
            "/" + f.max_mm + " mismatches" : "") +
           (f.note ? " &mdash; " + f.note : "") + "</title></rect>";
      let lx = f._inside ? x + 4 : x + w + 4, anchor = "start";
      if (f._hangL) { lx = x - 4; anchor = "end"; }
      s += '<text x="' + lx + '" y="' + (y + ah - 1.2) + '" text-anchor="' + anchor +
           '" font-size="8.5" fill="' + (f._inside ? "#fff" : "var(--ink-2)") +
           '" font-family="monospace">' + f.type + "</text>";
      if (f.type === "internal_dup" && f.partner_start !== undefined) {
        const px = sx(f.partner_start), pw = Math.max(2.5, sx(f.partner_end) - px);
        s += '<rect x="' + px + '" y="' + y + '" width="' + pw + '" height="' + ah +
             '" rx="2" fill="' + c + '" fill-opacity="0.45"/>' +
             '<line x1="' + (x + w / 2) + '" y1="' + (y + ah / 2) + '" x2="' + (px + pw / 2) +
             '" y2="' + (y + ah / 2) + '" stroke="' + c + '" stroke-width="1" ' +
             'stroke-dasharray="2 2"/>';
      }
    });
    const tsd = FEATON("tsd") ? ann.features.find(f => f.type === "tsd") : null;
    if (tsd) {
      const y = ay0 + nlane * (ah + 3);
      [[sx(-12), "5'"], [sx(L - 1 + 2), "3'"]].forEach(([x]) => {
        s += '<rect x="' + x + '" y="' + y + '" width="' + Math.max(6, sx(10) - sx(0)) +
             '" height="' + ah + '" rx="2" fill="' + FEATCOLOR.tsd +
             '" fill-opacity="0.75"><title>TSD in ' + (tsd.frac * 100).toFixed(0) +
             ' % of copies, median ' + tsd.len_med + ' bp</title></rect>';
      });
      s += '<text x="' + (sx(L - 1 + 2) + Math.max(6, sx(10) - sx(0)) + 5) + '" y="' +
           (y + ah - 1.2) + '" font-size="7.5" fill="var(--ink-2)" font-family="monospace">' +
           'TSD ' + (tsd.frac * 100).toFixed(0) + ' % of copies, ' + tsd.len_med + ' bp</text>';
    }
    const tc = FEATON("terminator") && ann.features.find(f => f.type === "terminator_copies");
    let extraLane = 0;
    if (tc) {
      const y2 = ay0 + (nlane + (tsd ? 1 : 0)) * (ah + 3);
      extraLane = 1;
      s += '<rect x="' + sx(L - 1 + 2) + '" y="' + y2 + '" width="' +
           Math.max(6, sx(14) - sx(0)) + '" height="' + ah + '" rx="2" fill="' +
           FEATCOLOR.terminator + '" fill-opacity="0.85"/>' +
           '<text x="' + (sx(L - 1 + 2) - 5) + '" y="' + (y2 + ah - 1.2) +
           '" text-anchor="end" font-size="8.5" fill="var(--ink-2)" ' +
           'font-family="monospace">Pol III term: ' +
           (tc.strong_frac * 100).toFixed(0) + '% strong, ' +
           (tc.moderate_frac * 100).toFixed(0) + '% moderate, ~' +
           tc.dist_med + ' bp out</text>';
    }
    TOTAL_H = H + MH + 4 + MOTIF_H + (nlane + 1 + extraLane) * (ah + 3) + 6;
  } else {
    TOTAL_H = H + MH + MOTIF_H;
  }
  s = s.replace('viewBox="0 0 ' + W + ' ' + (H + MH) + '"',
                'viewBox="0 0 ' + W + ' ' + TOTAL_H + '"');
  s += "</svg>";
  host.innerHTML = s;
}
"""

REPORT_PLOTS_JS = """
(function(){
  function addBar(plotId){
    var plot=document.getElementById(plotId);
    if(!plot)return;
    var bar=document.createElement('div');
    bar.style.cssText='margin:4px 0 6px;font-size:.72rem;';
    bar.innerHTML='<button type="button" style="margin-right:4px;padding:1px 6px;font-size:.72rem;cursor:pointer;" data-a="0">Hide all</button><button type="button" style="padding:1px 6px;font-size:.72rem;cursor:pointer;" data-a="1">Show all</button>';
    plot.parentNode.insertBefore(bar,plot);
    bar.querySelectorAll('button').forEach(function(b){
      b.addEventListener('click',function(){
        var gd=document.getElementById(plotId);
        if(!gd||!gd.data)return;
        var show=b.getAttribute('data-a')==='1';
        Plotly.restyle(gd,{visible:gd.data.map(function(){return show?true:'legendonly';})});
      });
    });
  }
  addBar('plot_div_kde');
  addBar('plot_pca');
})();
"""


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
        --muted:#666; --border:#e2e2e2;
        /* Diagram track/feature colors -- exact values from
           SINE_discriminator's site/index.html :root, so a diagram here reads
           identically to the site's (light-mode values only; this report has
           no dark-mode CSS of its own to match against). */
        --t1:#1baf7a; --t2:#2a78d6; --t3:#eda100; --t4:#4a3aa7; --t5:#e34948;
        --f1:#2a78d6; --f2:#eb6834; --f3:#1baf7a; --f4:#eda100; --f5:#e87ba4;
        --f6:#008300; --f7:#4a3aa7; --f8:#e34948;
        --div-lo:#2a78d6; --div-mid:#b9bcb6; --div-hi:#eb6834;
        --band:#e9efe9; --grid:#e2e7e0; --warn:#a8501d; --ink-2:#3d4b47; }
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
.diag-icon { font-size: 0.75rem; line-height: 1; background: none;
  border: 1px solid var(--muted); border-radius: 3px; width: 1.5em; height: 1.5em;
  padding: 0; cursor: pointer; color: inherit; vertical-align: middle; }
.diag-icon:hover { background: rgba(120,120,120,0.12); }
.diag-icon.diag-open { background: var(--accent); color: #fff; border-color: var(--accent); }
/* Diagram modal -- NOT inline in the table cell. A narrow table column
   (this one sits among 8 others) crushed the SVG to ~300px wide, which made
   every track and label illegible -- found 2026-09-09 by actually looking at
   the rendered page, not just checking that an <svg> element existed. The SVG
   itself scales via viewBox + width:100%, so the fix is giving it a large
   container, not changing the chart's internals. */
#diag-modal { display: none; position: fixed; inset: 0;
  background: rgba(0,0,0,.75); z-index: 9998;
  align-items: flex-start; justify-content: center; padding: 4vh 3vw;
  overflow-y: auto; }
#diag-modal.active { display: flex; }
#diag-modal-inner { background: var(--card); border-radius: 8px;
  padding: 20px 24px; width: min(1100px, 92vw); box-shadow: 0 8px 40px rgba(0,0,0,.5);
  position: relative; }
#diag-modal-title { font-size: 1rem; font-weight: 600; margin: 0 32px 10px 0; }
#diag-modal-close { position: absolute; top: 14px; right: 16px;
  background: none; border: none; font-size: 1.4rem; line-height: 1;
  cursor: pointer; color: var(--muted); padding: 4px; }
#diag-modal-close:hover { color: var(--fg); }
#diag-modal-msalink { display: inline-block; margin-bottom: 10px; }
.diag-legend { display: flex; flex-wrap: wrap; gap: 10px 18px; margin-bottom: 10px;
  padding-bottom: 8px; border-bottom: 1px solid var(--border); font-size: 0.76rem; }
.diag-legend-group { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 10px; }
.diag-legend-title { font-weight: 600; color: var(--muted); margin-right: 2px; }
.diag-legend label { display: inline-flex; align-items: center; gap: 3px; cursor: pointer; }
.diag-legend .sw { display: inline-block; width: 14px; height: 3px; border-radius: 1px; }
.diag-legend .swb { display: inline-block; width: 10px; height: 10px; border-radius: 2px; }
.diag-chips { display: flex; flex-wrap: wrap; gap: 6px 14px; margin-bottom: 6px;
  font-size: 0.78rem; }
.diag-chip b { font-weight: 600; margin-right: 3px; }
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


VCHIP_CSS = """
.vchip { display: inline-block; font-size: .78rem; font-weight: 600;
         padding: 2px 7px; border-radius: 3px; white-space: nowrap; }
.vchip.v-ok { background: #e2efe9; color: #1f6f5c; }
.vchip.v-edge { background: #f5eed8; color: #7a5c12; }
.vchip.v-warn { background: #f6e8de; color: #a8501d; }
.vchip.v-muted { background: #eceee8; color: #68766f; }
"""


def _vchip(label: str, kind: str, title: str = "") -> str:
    t = f" title='{html.escape(title, quote=True)}'" if title else ""
    return f"<span class='vchip v-{kind}'{t}>{html.escape(label)}</span>"


def _vcodes(v):
    return {f["code"] for f in v.get("flags", [])}


def status_flanks(v):
    """Ported from SINE_discriminator's inject_oma_aln_section.py (2026-09-09),
    unmodified logic, against report_verdict.verdict()'s output."""
    if not v or v.get("error"):
        return _vchip("n/a", "muted", v.get("error", "not scored") if v else "")
    codes = _vcodes(v)
    if "NO_FLANKS_PRESENT" in codes:
        return _vchip("No flanks", "muted", "Alignment carries no flank sequence.")
    if "NOT_ISOLATED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "NOT_ISOLATED")
        return _vchip("Satellite / dup", "warn", f.get("text", ""))
    if "FRAGMENT_OF_LONGER" in codes or "ELEMENT_CONTINUES" in codes:
        f = next(x for x in v["flags"]
                 if x["code"] in ("FRAGMENT_OF_LONGER", "ELEMENT_CONTINUES"))
        return _vchip("Fragment / LINE", "warn", f.get("text", ""))
    if "SHARED_FLANKS" in codes:
        f = next(x for x in v["flags"] if x["code"] == "SHARED_FLANKS")
        return _vchip("Shared flanks", "warn", f.get("text", ""))
    if "FLANK_ISLANDS" in codes:
        f = next(x for x in v["flags"] if x["code"] == "FLANK_ISLANDS")
        return _vchip("Flank islands", "edge", f.get("text", ""))
    if "FLANKS_UNMEASURED" in codes:
        return _vchip("Short-flank view", "muted",
                      "50L/70R publish geometry; no 400 bp decay profile.")
    fb = v.get("flank_bg")
    if fb is None:
        return _vchip("Not measured", "muted")
    if fb < 0.32:
        return _vchip("Independent", "ok",
                      "Flank pairwise identity %.2f vs ~0.25 background." % fb)
    if fb < 0.42:
        return _vchip("Borderline", "edge",
                      "Flank background %.2f — check by eye." % fb)
    return _vchip("Raised flanks", "warn",
                  "Flank background %.2f — copies may share context." % fb)


def status_element(v):
    if not v or v.get("error"):
        return _vchip("n/a", "muted", v.get("error", "not scored") if v else "")
    codes = _vcodes(v)
    g = v.get("groups", {}).get("element", 0.0)
    if "NO_ELEMENT" in codes or g < 0.25:
        f = next((x for x in v["flags"] if x["code"] == "NO_ELEMENT"), None)
        return _vchip("No element", "warn", (f or {}).get("text", ""))
    if "MICROSATELLITE_ELEMENT" in codes:
        f = next(x for x in v["flags"] if x["code"] == "MICROSATELLITE_ELEMENT")
        return _vchip("Microsatellite", "warn", f.get("text", ""))
    if "CONSENSUS_OVEREXTENDED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "CONSENSUS_OVEREXTENDED")
        return _vchip("Over-extended", "edge", f.get("text", ""))
    if "CONSENSUS_UNDEREXTENDED" in codes:
        f = next(x for x in v["flags"] if x["code"] == "CONSENSUS_UNDEREXTENDED")
        return _vchip("Under-extended", "edge", f.get("text", ""))
    if "SMALL_CORE" in codes:
        f = next(x for x in v["flags"] if x["code"] == "SMALL_CORE")
        return _vchip("Small core", "edge", f.get("text", ""))
    if g >= 0.85:
        return _vchip("Strong", "ok",
                      "%d/%d copies support the consensus."
                      % (v.get("n_supported", 0), v.get("n", 0)))
    if g >= 0.5:
        return _vchip("Supported", "ok",
                      "%d/%d copies support the consensus."
                      % (v.get("n_supported", 0), v.get("n", 0)))
    return _vchip("Weak", "edge", "Element group score %.2f." % g)


def status_overall(v):
    if not v or v.get("error"):
        return _vchip("n/a", "muted", v.get("error", "not scored") if v else "")
    if v.get("deferred"):
        return _vchip("Deferred", "edge", "Mixture — split before a firm call.")
    if not v.get("assessable", True):
        return _vchip("Cannot assess", "muted", "Too little evidence for a score.")
    s = float(v.get("score", 0))
    note = "%.0f/100" % s
    if v.get("capped_by"):
        note += " (capped: %s)" % ", ".join(v["capped_by"])
    if s >= 90:
        return _vchip("SINE", "ok", note)
    if s >= 75:
        return _vchip("SINE (caveats)", "ok", note)
    if s >= 55:
        return _vchip("Grey zone", "edge", note)
    return _vchip("Not SINE", "warn", note)


def _side_note(side):
    if not side.get("measured"):
        return side.get("reason", "not measured")
    return (
        "%.0f%% unique; largest shared group %.0f%% (%d copies)"
        % (100 * side.get("unique_frac", 0),
           100 * side.get("largest_cluster_frac", 0),
           side.get("largest_cluster", 0))
    )


def status_flank_context(top, rand):
    """Per-side flank clustering: rand100 high, top100 medium."""
    parts = []
    worst = None
    for label, r in (("rand", rand), ("top", top)):
        if not r or r.get("error"):
            continue
        sev = r.get("worst_flag")
        if sev == "high":
            worst = "high"
        elif sev == "medium" and worst != "high":
            worst = "medium"
        for f in r.get("flags", []):
            parts.append("%s %s: %s" % (label, f["side"], f["text"]))
    if not top and not rand:
        return _vchip("n/a", "muted", "no alignment")
    if not parts:
        note = ""
        if top and top.get("left", {}).get("measured"):
            note = "Top100 L/R " + _side_note(top["left"]) + "; " + _side_note(top["right"])
        return _vchip("Independent", "ok", note or "No shared-flank clusters above threshold.")
    title = " ".join(parts)
    if worst == "high":
        return _vchip("Shared context", "warn", title)
    return _vchip("Subgroup", "edge", title)


def compute_verdicts(species_code: str, subfams: List[str], aln_dir: Path) -> Dict[str, dict]:
    """{sf: (verdict_dict_or_None, flank_scan_top100_or_None, flank_scan_rand100_or_None)}
    Ported from inject_oma_aln_section.py's score_top100()/scan_flanks(), which
    scored only the top100 tier for the verdict itself (element/flanks/overall
    columns) while flank-context looks at both tiers."""
    try:
        import report_verdict as V
        import report_flank_uniqueness as FU
    except ImportError as exc:
        LOG.warning("Skipping verdict columns: %s (numpy required)", exc)
        return {}
    out = {}
    for sf in subfams:
        t100 = aln_dir / f"{species_code}_{sf}_top100.aln.fa"
        r100 = aln_dir / f"{species_code}_{sf}_rand100.aln.fa"
        v = None
        if t100.is_file():
            try:
                v = V.verdict(str(t100))
            except Exception as exc:
                v = {"error": str(exc)}
        fu_t = None
        if t100.is_file():
            try:
                fu_t = FU.scan(str(t100))
            except Exception as exc:
                fu_t = {"error": str(exc)}
        fu_r = None
        if r100.is_file():
            try:
                fu_r = FU.scan(str(r100))
            except Exception as exc:
                fu_r = {"error": str(exc)}
        out[sf] = {"verdict": v, "flank_top": fu_t, "flank_rand": fu_r}
    return out


def compute_profiles(species_code: str, subfams: List[str], aln_dir: Path,
                      max_diagrams: int) -> Dict[str, dict]:
    """Positional alignment-composition diagram data (report_profile.py, ported
    from SINE_discriminator's profiles.py/measure_c.py 2026-09-09) for every
    subfamily x {top100, rand100} pair whose alignment file actually exists on
    disk. Computed once here at report-build time and embedded as JSON — a
    static HTML report has no server to compute on click, so "lazy" means only
    computed for tiers that exist, not deferred past build time; the button
    itself only renders the (already-embedded) SVG on expand, which is cheap.

    max_diagrams caps total tiers computed (each is one MAFFT-scale numpy pass)
    so a report with many subfamilies doesn't stall step6 — same pattern as
    --sineplot-max.
    """
    try:
        import report_profile as RP
    except ImportError as exc:
        LOG.warning("Skipping alignment diagrams: %s (numpy required)", exc)
        return {}
    out: Dict[str, dict] = {}
    n = 0
    for sf in subfams:
        for tier in ("top100", "rand100"):
            if n >= max_diagrams:
                LOG.info("Alignment diagrams: stopped at cap (%d)", max_diagrams)
                return out
            fn = aln_dir / f"{species_code}_{sf}_{tier}.aln.fa"
            if not fn.is_file():
                continue
            try:
                p = RP.profile(str(fn))
                m = RP.measure(str(fn))
            except Exception as exc:
                LOG.warning("Alignment diagram failed for %s/%s: %s", sf, tier, exc)
                continue
            if p is None:
                continue
            p["_measure"] = m
            out[f"{sf}_{tier}"] = p
            n += 1
    return out


def compute_annotations(species_code: str, subfams: List[str], aln_dir: Path,
                         profiles: Dict[str, dict], max_diagrams: int) -> Dict[str, dict]:
    """Structural-feature annotations (report_annotate.py, ported from
    SINE_discriminator's annotate.py 2026-09-09) for the same tiers
    compute_profiles already found a profile for -- one alignment read each,
    reusing the already-computed profile for the conserved-block detection."""
    try:
        import report_annotate as RA
    except ImportError as exc:
        LOG.warning("Skipping structural-feature annotations: %s", exc)
        return {}
    out: Dict[str, dict] = {}
    n = 0
    for sf in subfams:
        for tier in ("top100", "rand100"):
            key = f"{sf}_{tier}"
            if key not in profiles or n >= max_diagrams:
                continue
            fn = aln_dir / f"{species_code}_{sf}_{tier}.aln.fa"
            if not fn.is_file():
                continue
            try:
                a = RA.annotate(str(fn), profiles[key])
            except Exception as exc:
                LOG.warning("Annotation failed for %s/%s: %s", sf, tier, exc)
                continue
            if a is None:
                continue
            out[key] = a
            n += 1
    return out


def build_alignment_section(
    species_code: str,
    subfams: List[str],
    msa_url: str = "https://toki-bio.github.io/MSA-viewer/",
    raw_base: Optional[str] = None,
    aln_dir: Optional[Path] = None,
    profiles: Optional[Dict[str, dict]] = None,
    verdicts: Optional[Dict[str, dict]] = None,
) -> str:
    """Alignment table. With raw_base (http URL): MSA-viewer links. Else: relative paths."""
    use_remote = bool(raw_base and raw_base.startswith("http"))
    profiles = profiles or {}
    verdicts = verdicts or {}

    def aln_href(fn: str, title: str, remote_fn: Optional[str] = None) -> str:
        if use_remote:
            base = raw_base.rstrip("/") + "/"
            if "alignments" not in base:
                base = f"{raw_base.rstrip('/')}/{species_code}/alignments/"
            url = base + (remote_fn or fn)
            return (f"{msa_url}?url={quote(url, safe='')}"
                    f"&title={quote(title, safe='')}")
        return f"alignments/{fn}"

    def aln_link(fn: str, title: str, label: str, css: str = "", remote_fn: Optional[str] = None) -> str:
        href = aln_href(fn, title, remote_fn)
        cls = "aln-link" + (f" {css}" if css else "")
        return (f'<a class="{cls}" href="{html.escape(href, quote=True)}" '
                f'target="_blank">{html.escape(label)}</a>')

    intro_links = ""
    if use_remote:
        intro_links = (
            "&nbsp;&nbsp;"
            + aln_link(f"{species_code}_consensuses.fa",
                       f"{species_code} all consensi", "All consensi")
            + "&nbsp;&nbsp;"
            + aln_link(f"{species_code}_subfam_input.aln.fa",
                       f"{species_code} SubFam input", "SubFam input")
        )
    elif aln_dir:
        intro_links = (
            "<span class='small muted'>&nbsp;(Relative links; pass "
            "<code>--aln-base</code> with a published raw URL for MSA viewer.)</span>"
        )

    rows_html = ""
    for sf in sorted(subfams):
        t100_fn = f"{species_code}_{sf}_top100.aln.fa"
        r100_fn = f"{species_code}_{sf}_rand100.aln.fa"
        sub_fn = f"{species_code}_{sf}_subfam.aln.fa"
        if aln_dir and not (aln_dir / t100_fn).is_file():
            continue
        # Some runs' subfamily names already carry the species prefix (this
        # oma run's do -- its own consensus headers are `>oma_SINE10`, so
        # step8a's own naming doubles it to `oma_oma_SINE10_...` on disk),
        # while a raw_base publish typically strips that duplication (verified
        # against the real published copies: alignments/oma/oma_SINE10_top100
        # .aln.fa on GitHub, single-prefixed, 2026-09-09). Local disk lookups
        # must use the actual on-disk (possibly doubled) name; remote links
        # must use the deduplicated one, or they 404 against a real publish.
        remote_sf = sf[len(species_code) + 1:] if sf.startswith(f"{species_code}_") else sf
        t100_remote = f"{species_code}_{remote_sf}_top100.aln.fa"
        r100_remote = f"{species_code}_{remote_sf}_rand100.aln.fa"
        sub_remote = f"{species_code}_{remote_sf}_subfam.aln.fa"
        has_t100 = f"{sf}_top100" in profiles
        has_r100 = f"{sf}_rand100" in profiles

        def _diag_icon(tier_key: str, tier_label: str, remote_fn: str) -> str:
            # tiny icon inline with the link it diagrams -- opens the SHARED
            # full-size modal (#diag-modal), not a per-row box. A per-row div
            # confined to this table's column width crushed the SVG to ~300px,
            # illegible (found 2026-09-09 by actually looking at the render).
            href = aln_href(remote_fn, f"{species_code} {tier_key}", remote_fn=remote_fn)
            return (f"<button type='button' class='diag-icon' "
                    f"data-tier='{tier_key}' data-label='{html.escape(sf)} {tier_label}' "
                    f"data-msa-href='{html.escape(href, quote=True)}' "
                    f"title='{tier_label} alignment diagram'>&#9656;</button>")

        t100_icon = _diag_icon(f"{sf}_top100", "top100", t100_remote) if has_t100 else ""
        r100_icon = _diag_icon(f"{sf}_rand100", "rand100", r100_remote) if has_r100 else ""
        vd = verdicts.get(sf) or {}
        verdict_cells = (
            f"<td>{status_flanks(vd.get('verdict'))}</td>"
            f"<td>{status_flank_context(vd.get('flank_top'), vd.get('flank_rand'))}</td>"
            f"<td>{status_element(vd.get('verdict'))}</td>"
            f"<td>{status_overall(vd.get('verdict'))}</td>"
        ) if vd else (
            "<td class='small muted'>n/a</td>" * 4
        )
        rows_html += (
            f"<tr><td><code>{html.escape(sf)}</code></td>"
            f"<td>{aln_link(t100_fn, f'{species_code} {sf} top100', 'top 100 by score', remote_fn=t100_remote)}"
            f" {t100_icon}</td>"
            f"<td>{aln_link(r100_fn, f'{species_code} {sf} rand100', '100 random', 'orange', remote_fn=r100_remote)}"
            f" {r100_icon}</td>"
            f"<td>{aln_link(sub_fn, f'{species_code} {sf} subfam', 'SubFam', 'green', remote_fn=sub_remote)}</td>"
            f"{verdict_cells}"
            "</tr>"
        )
    if not rows_html:
        return ""
    return (
        "<section class='card' id='alignments'>"
        "<h2>Subfamily alignments</h2>"
        "<p class='intro'>Copies re-extracted with "
        "<strong>50&thinsp;bp upstream + 70&thinsp;bp downstream</strong> "
        "genomic flanks (strand-aware). The alignment diagram plots per-position "
        "pairwise identity, coverage, consensus identity, A+T fraction and "
        "mosaicism across the 5&prime; flank/element/3&prime; flank, plus A "
        "box/B box/tRNA-head self-similarity motif scores — a track reads flat "
        "at z&#8776;0 (\"not detected\"), never disappears, on families that "
        "are not tRNA-derived. The last four columns score the "
        "<strong>top 100</strong> alignment with <code>report_verdict.py</code> "
        "(ported from SINE_discriminator's <code>verdict.py</code>) — each cell "
        "is a coloured label, hover for the measurement behind it."
        + intro_links
        + "</p>"
        "<table class='tbl'>"
        "<thead><tr><th>Subfamily</th>"
        "<th title='Click &#9656; to open the alignment-composition diagram'>"
        "Top 100 by bitscore</th>"
        "<th title='Click &#9656; to open the alignment-composition diagram'>"
        f"100 random copies (seed {html.escape(os.environ.get('RAND_SEED', '42'))})</th>"
        "<th>SubFam (chunk consensuses)</th>"
        "<th>Flanks</th><th>Flank context</th><th>Element</th><th>Overall</th>"
        "</tr></thead>"
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
               tal_species_code: Optional[str] = None,
               aln_base: Optional[str] = None,
               pages_index: Optional[str] = None,
               profile_diagrams: bool = True,
               profile_diagrams_max: int = 200) -> None:
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

    plots_dir = s2 / "plots"
    figs = {}
    pctid_by_sf = load_pctid_by_sf(plots_dir) if plots_dir.is_dir() else {}
    use_pctid = bool(pctid_by_sf)
    if use_pctid:
        figs["div_kde"] = fig_pctid_spline_divergence(pctid_by_sf)
    else:
        figs["div_kde"] = fig_divergence_kde(sim_by_sf)
        figs["sim_violins"] = fig_sim_violins(sim_by_sf)
    if curves:
        figs["conservation"] = fig_conservation(curves)

    # Per-subfamily PNG gallery
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

    # Alignment section (requires --species-code; links need --aln-base for MSA viewer)
    species_code = tal_species_code
    alignment_section = ""
    aln_profiles: Dict[str, dict] = {}
    aln_annotations: Dict[str, dict] = {}
    aln_dir = run_root / "results" / "alignments"
    if not aln_dir.is_dir():
        aln_dir = run_root / "alignments"
    if species_code:
        subfams_for_aln = sorted({r[0] for r in stats_rows if r})
        if aln_dir.is_dir():
            from_disk = sorted({
                p.name.replace("_top100.aln.fa", "").replace(f"{species_code}_", "", 1)
                for p in aln_dir.glob(f"{species_code}_*_top100.aln.fa")
            })
            if from_disk:
                subfams_for_aln = from_disk
        if profile_diagrams and aln_dir.is_dir():
            aln_profiles = compute_profiles(
                species_code, subfams_for_aln, aln_dir, profile_diagrams_max)
            aln_annotations = compute_annotations(
                species_code, subfams_for_aln, aln_dir, aln_profiles, profile_diagrams_max)
        aln_verdicts: Dict[str, dict] = {}
        if profile_diagrams and aln_dir.is_dir():
            # gate on the same flag as the diagrams -- both are the optional,
            # numpy-needing, per-file-read analysis stage
            aln_verdicts = compute_verdicts(species_code, subfams_for_aln, aln_dir)
        alignment_section = build_alignment_section(
            species_code, subfams_for_aln, raw_base=aln_base, aln_dir=aln_dir,
            profiles=aln_profiles, verdicts=aln_verdicts)
    # Element hierarchy (flankscan stage 7: composites drawn as their parts, to scale); "" without it
    hierarchy_section = ""
    try:
        _here = os.path.dirname(os.path.abspath(__file__))       # publish runs step6 from elsewhere
        if _here not in sys.path:
            sys.path.insert(0, _here)
        import report_hierarchy as RH
        hierarchy_section = RH.section(
            run_root, species_code, aln_base,
            subfams_for_aln if species_code else sorted({r[0] for r in stats_rows if r}))
    except Exception as e:  # the report must not fail over an optional panel
        sys.stderr.write("WARNING: element hierarchy skipped: %s\n" % e)
    # Sequence similarity between consensuses (the sequence-part counterpart of the hierarchy); "" when numpy / bank missing
    blocks_section = ""
    try:
        import report_blocks as RB
        blocks_section = RB.section(run_root)
    except Exception as e:  # optional panel
        sys.stderr.write("WARNING: similarity blocks skipped: %s\n" % e)
    profiles_json = json.dumps(aln_profiles, separators=(",", ":"))
    annotations_json = json.dumps(aln_annotations, separators=(",", ":"))

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

    if use_pctid:
        divergence_body = (
            '<p class="intro"><b>Metric:</b> divergence = 100 &minus; ssearch36 '
            '%identity to the subfamily consensus (same as Gallery histograms). '
            'Copies binned at 1% divergence; line connects bin counts (smooth spline). '
            'One line per subfamily; click legend to hide/show.</p>'
            + conservation_html +
            '<div class="plot" id="plot_div_kde"></div>'
            '<p class="small muted">Source: step4 <code>*_pctid.tsv</code> on '
            'assigned copies.</p>'
        )
    else:
        divergence_body = (
            '<p class="intro"><b>Metric:</b> bitscore-based divergence = '
            '100&thinsp;&minus;&thinsp;(copy bitscore&thinsp;/&thinsp;consensus '
            'self-bitscore &times; 100&thinsp;%). One KDE curve per subfamily '
            '(up to 3,000 copies sampled); click legend to toggle.</p>'
            + conservation_html +
            '<div class="plot" id="plot_div_kde"></div>'
            '<h3>Distributions per subfamily (violin)</h3>'
            '<p class="intro">Same data as above shown as violin plots.</p>'
            '<div class="plot" id="plot_sim_violins"></div>'
            '<p class="small muted">Source: <code>sim_scores.tsv</code> joined with '
            '<code>assignment_full.tsv</code>.</p>'
        )

    title = manifest.get("RUN", str(run_root)).rstrip("/").split("/")[-1]
    genome_name, cons_name = run_inputs(manifest)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    cross_nav = ""
    if pages_index:
        cross_nav = (
            "<div style='background:#1a2634;color:rgba(255,255,255,.8);"
            "padding:6px 32px;font-size:.85rem;'>"
            f"<a href='{html.escape(pages_index, quote=True)}' "
            "style='color:rgba(255,255,255,.85);text-decoration:none;'>"
            "&#8592; All species</a></div>"
        )

    html_doc = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SINEderella report &mdash; {html.escape(title)}</title>
<style>{CSS}{VCHIP_CSS}</style>
{plotly_tag}
</head><body>
<div id="lightbox" onclick="this.classList.remove('active')">
  <img id="lightbox-img" src="" alt="">
</div>
<div id="diag-modal">
  <div id="diag-modal-inner">
    <button type="button" id="diag-modal-close" title="Close" aria-label="Close">&times;</button>
    <p id="diag-modal-title"></p>
    <a id="diag-modal-msalink" class="aln-link" target="_blank" rel="noopener" hidden></a>
    <div id="diag-modal-legend" class="diag-legend"></div>
    <div id="diag-modal-body"></div>
  </div>
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
    {'<a href="#hierarchy">Hierarchy</a>' if hierarchy_section else ''}
    {'<a href="#blocks">Similarity</a>' if blocks_section else ''}
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

  {hierarchy_section}

  {blocks_section}

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
    {divergence_body}
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
{REPORT_PLOTS_JS}
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
<script>
// Deliberately a SEPARATE <script> tag from the block above: fig_init calls
// Plotly.newPlot, and if the Plotly CDN is unreachable (blocked host, offline,
// or just a slow load race) that throws and aborts every remaining statement
// in ITS OWN script block -- which silently killed the diagram toggle below
// when both lived in one block (found 2026-09-09, diagrams did not expand on
// a published Artifact because cdn.plot.ly is not on its script allowlist).
// Splitting the block means a Plotly failure can never take this out too.
var PROFILES = {profiles_json};
var ANNOTATIONS = {annotations_json};
{DIAGRAM_JS}
(function(){{
  // Every diag-icon opens the SAME full-size modal (#diag-modal), sized to
  // min(1100px, 92vw) -- not a per-row box confined to a table column, which
  // is what made the chart illegible.
  var modal = document.getElementById('diag-modal');
  var body = document.getElementById('diag-modal-body');
  var title = document.getElementById('diag-modal-title');
  var msalink = document.getElementById('diag-modal-msalink');
  var legend = document.getElementById('diag-modal-legend');
  var closeBtn = document.getElementById('diag-modal-close');
  var legendBuilt = false;
  var currentTier = null;
  function closeModal(){{ modal.classList.remove('active'); body.innerHTML = ''; currentTier = null; }}
  if (closeBtn) closeBtn.addEventListener('click', closeModal);
  if (modal) modal.addEventListener('click', function(e){{
    if (e.target === modal) closeModal();
  }});
  document.addEventListener('keydown', function(e){{
    if (e.key === 'Escape' && modal.classList.contains('active')) closeModal();
  }});
  document.querySelectorAll('.diag-icon').forEach(function(btn){{
    btn.addEventListener('click', function(){{
      if (!legendBuilt) {{
        buildDiagLegend(legend);
        legend.addEventListener('change', function(e){{
          if (e.target.matches('.trk,.ftr') && currentTier) drawProfile(currentTier, body);
        }});
        legendBuilt = true;
      }}
      title.textContent = btn.dataset.label;
      if (btn.dataset.msaHref) {{
        msalink.href = btn.dataset.msaHref;
        msalink.textContent = 'open in MSA-viewer ↗';
        msalink.hidden = false;
      }} else {{
        msalink.hidden = true;
      }}
      currentTier = btn.dataset.tier;
      modal.classList.add('active');
      drawProfile(currentTier, body);
    }});
  }});
}})();
</script>
{SORT_JS}
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
    ap.add_argument("--profile-diagrams", dest="profile_diagrams",
                    action="store_true", default=True,
                    help="Per-subfamily alignment-composition diagram, "
                         "top100/rand100 toggle (default on; needs numpy).")
    ap.add_argument("--no-profile-diagrams", dest="profile_diagrams",
                    action="store_false",
                    help="Skip the alignment-composition diagrams.")
    ap.add_argument("--profile-diagrams-max", type=int, default=200,
                    help="Max subfamily x tier diagrams computed per report "
                         "(default: 200).")
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("THREADS",
                                               os.cpu_count() or 1)),
                    help="Threads for ssearch36 inside SINEplot stage.")
    ap.add_argument("--tal-species-code", default=None,
                    help="Deprecated alias for --species-code.")
    ap.add_argument("--species-code", default=None,
                    help="Species prefix for alignment filenames (e.g. mysp).")
    ap.add_argument("--aln-base", default=None,
                    help="Published raw URL prefix for alignment MSAs (required "
                         "for MSA-viewer links), e.g. "
                         "https://raw.githubusercontent.com/org/repo/main/mysp/alignments/")
    ap.add_argument("--pages-index", default=None,
                    help="Optional URL for 'All species' nav link (multi-species site).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    species = args.species_code or args.tal_species_code
    pages_index = args.pages_index or os.environ.get("PAGES_INDEX") or None

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
               tal_species_code=species,
               aln_base=args.aln_base,
               pages_index=pages_index,
               profile_diagrams=args.profile_diagrams,
               profile_diagrams_max=args.profile_diagrams_max)
    return 0


if __name__ == "__main__":
    sys.exit(main())
