#!/usr/bin/env python3
"""Side-by-side test: conse vs sine_consensus.sh on the same top100 copies.

Same input for both tools: ungapped element sequences cut from the alignment
using copy-supported boundary walk (boundary.py rule). Also reports the peel
loci consensus (conse on rebuild alignment) and the alignment row-0 span.

Writes HTML + JSON under OUT_DIR for Sergei to pick a canonical rebuild method.
"""
from __future__ import annotations

import html
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np

# boundary walk (inlined — same as sinedisc/boundary.py)
OCC, CONS, MISS = 0.50, 0.45, 8
GAPS = set("-.")


def column_stats(rows, j):
    b = [r[j] for r in rows if j < len(r) and r[j] not in GAPS]
    if not b:
        return 0.0, 0.0
    occ = len(b) / float(len(rows))
    u = [x.upper() for x in b]
    top = max(u.count(c) for c in set(u))
    return occ, top / float(len(u))


def element_window(cons, rows):
    L = len(cons)
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not rows:
        return 0, L - 1, {}
    old_lo, old_hi = nz[0], nz[-1]
    ok = {}

    def good(j):
        if j not in ok:
            o, c = column_stats(rows, j)
            ok[j] = (o >= OCC and c >= CONS)
        return ok[j]

    core = [j for j in nz if good(j)]
    if not core:
        return old_lo, old_hi, {"note": "no supported consensus column"}
    lo, hi = core[0], core[-1]

    miss, j = 0, lo - 1
    while j >= 0 and miss < MISS:
        if good(j):
            lo, miss = j, 0
        else:
            miss += 1
        j -= 1

    miss, j = 0, hi + 1
    while j < L and miss < MISS:
        if good(j):
            hi, miss = j, 0
        else:
            miss += 1
        j += 1

    d = {
        "old_lo": old_lo, "old_hi": old_hi, "new_lo": lo, "new_hi": hi,
        "extended_left": max(0, old_lo - lo), "extended_right": max(0, hi - old_hi),
        "old_span": old_hi - old_lo + 1, "new_span": hi - lo + 1,
    }
    return lo, hi, d


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
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


def write_fa(path, pairs):
    with open(path, "w") as fh:
        for n, s in pairs:
            fh.write(">%s\n%s\n" % (n, s))


def ungap(s):
    return re.sub(r"[^ACGTacgt]", "", s).upper()


def sh(cmd, cwd=None, timeout=3600):
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True,
                       text=True, timeout=timeout)
    return r.returncode, r.stdout, r.stderr


def identity(a, b):
    a, b = ungap(a), ungap(b)
    if not a or not b:
        return 0.0
    # global identity after padding shorter to longer (simple)
    n = max(len(a), len(b))
    m = min(len(a), len(b))
    if m == 0:
        return 0.0
    same = sum(1 for i in range(m) if a[i] == b[i])
    return same / float(n)


def nucfreq_profile(consensus, copy_seqs):
    """Per-position frequency of consensus base over ungapped copy strings."""
    L = len(consensus)
    y = []
    for i, base in enumerate(consensus):
        b = base.upper()
        if b not in "ACGT":
            y.append(0.0)
            continue
        n = tot = 0
        for s in copy_seqs:
            if i < len(s) and s[i] in "ACGT":
                tot += 1
                if s[i].upper() == b:
                    n += 1
        y.append(n / tot if tot else 0.0)
    return y


def conservation_from_alignment(cons_row, copy_rows):
    """Consensus-base fraction per column from gapped alignment."""
    cons = list(cons_row.upper())
    y = []
    for j, base in enumerate(cons):
        if base not in "ACGT":
            y.append(None)
            continue
        n = tot = 0
        for r in copy_rows:
            if j >= len(r):
                continue
            c = r[j].upper()
            if c in "ACGT":
                tot += 1
                if c == base:
                    n += 1
        y.append(n / tot if tot else 0.0)
    return y


def run_conse(aln_fa, pct, work):
    out_cons = aln_fa + ".cons"
    code, _, err = sh("bash /usr/bin/conse %s %d" % (aln_fa, pct), cwd=work)
    if code != 0 or not os.path.isfile(out_cons):
        return None, "conse failed: %s" % err[-500:]
    _, seqs = read_fa(out_cons)
    if not seqs:
        return None, "empty conse output"
    return ungap(seqs[0]), None


def run_sine_consensus(copies_fa, work, script):
    code, _, err = sh("bash %s %s 100 50 0.01" % (script, copies_fa), cwd=work,
                      timeout=7200)
    hits = list(Path(work).glob("*_consensus.fasta")) + \
           list(Path(work).glob("final_consensus.fasta"))
    log = Path(work) / "consensus_log.txt"
    iters = None
    if log.is_file():
        m = re.search(r"Converged at iteration (\d+)", log.read_text(errors="replace"))
        if m:
            iters = int(m.group(1))
    if not hits:
        return None, iters, "sine_consensus failed: %s" % err[-500:]
    _, seqs = read_fa(str(hits[0]))
    if not seqs:
        return None, iters, "empty sine_consensus output"
    return ungap(seqs[0]), iters, None


def svg_conservation(curve, title, width=900, height=80):
    pts = [(i, v) for i, v in enumerate(curve) if v is not None]
    if not pts:
        return "<p>No data</p>"
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs), max(xs)
    pad = 8
    w = width - 2 * pad
    h = height - 2 * pad
    path = []
    for i, (x, y) in enumerate(pts):
        px = pad + (x - x0) / max(1, x1 - x0) * w
        py = pad + (1 - y) * h
        path.append(("%.1f,%.1f" % (px, py)) if i else "M %.1f,%.1f" % (px, py))
    d = " ".join(path[1:]) if path else ""
    if path:
        d = path[0] + " L " + " L ".join(path[1:])
    return (
        '<svg viewBox="0 0 %d %d" width="%d" height="%d" role="img" '
        'aria-label="%s"><title>%s</title>'
        '<rect x="0" y="0" width="%d" height="%d" fill="var(--surface-2)"/>'
        '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="var(--rule)"/>'
        '<path d="%s" fill="none" stroke="var(--accent)" stroke-width="1.5"/>'
        '</svg>' % (width, height, width, height, html.escape(title),
                    html.escape(title), width, height,
                    pad, pad + h, width - pad, pad + h, d)
    )


def main():
    aln_path = sys.argv[1]
    out_dir = sys.argv[2]
    peel_aln = sys.argv[3] if len(sys.argv) > 3 else ""
    sine_script = sys.argv[4] if len(sys.argv) > 4 else \
        "/staging/tmp/SINEderella/sine_consensus.sh"
    tag = Path(aln_path).stem.replace(".aln", "")

    os.makedirs(out_dir, exist_ok=True)
    work = os.path.join(out_dir, "work")
    os.makedirs(work, exist_ok=True)

    names, seqs = read_fa(aln_path)
    if len(seqs) < 10:
        sys.exit("too few sequences in %s" % aln_path)

    ci = [i for i, n in enumerate(names) if "CONSENSUS" in n.upper()]
    if not ci:
        ci = [0]
    k = ci[0]
    cons_row = seqs[k].upper()
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]

    lo, hi, bdiag = element_window(cons_row, rows)
    copies = []
    for i, (n, s) in enumerate(zip(names, seqs)):
        if i == k:
            continue
        copies.append((re.sub(r"\s+", "_", n)[:80], ungap(s[lo:hi + 1])))

    copies_fa = os.path.join(work, "copies_element.fa")
    write_fa(copies_fa, copies)
    copy_seqs = [s for _, s in copies]

    row0_ungap = ungap(cons_row[lo:hi + 1])

    results = {}

    # --- conse at 35% and 50% ---
    aln_out = os.path.join(work, "copies.mafft.aln.fa")
    code, _, err = sh(
        "mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 4 "
        "%s > %s 2>/dev/null" % (copies_fa, aln_out), cwd=work)
    if code != 0:
        sys.exit("mafft failed: %s" % err[-800:])

    for pct in (35, 50):
        seq, errm = run_conse(aln_out, pct, work)
        key = "conse_%d" % pct
        results[key] = {"seq": seq, "len": len(seq) if seq else 0, "error": errm,
                        "tool": "conse", "pct": pct}

    # --- sine_consensus.sh ---
    sc_work = os.path.join(work, "sine_consensus")
    os.makedirs(sc_work, exist_ok=True)
    import shutil
    local_copies = os.path.join(sc_work, "copies_element.fa")
    shutil.copy(copies_fa, local_copies)
    seq, iters, errm = run_sine_consensus("copies_element.fa", sc_work, sine_script)
    results["sine_consensus"] = {
        "seq": seq, "len": len(seq) if seq else 0, "error": errm,
        "tool": "sine_consensus.sh", "iters": iters,
    }

    # --- peel loci consensus (conse on rebuild alignment) ---
    if peel_aln and os.path.isfile(peel_aln):
        peel_work = os.path.join(work, "peel")
        os.makedirs(peel_work, exist_ok=True)
        seq, errm = run_conse(peel_aln, 50, peel_work)
        results["peel_loci_conse50"] = {
            "seq": seq, "len": len(seq) if seq else 0, "error": errm,
            "tool": "conse on peel loci MSA", "pct": 50,
        }

    results["boundary_row0"] = {
        "seq": row0_ungap, "len": len(row0_ungap), "error": None,
        "tool": "alignment row-0 span (ungapped)", "boundary": bdiag,
    }

    # --- pairwise identity matrix ---
    keys = [k for k, v in results.items() if v.get("seq")]
    mat = {}
    for a in keys:
        mat[a] = {}
        for b in keys:
            mat[a][b] = round(identity(results[a]["seq"], results[b]["seq"]), 4)

    # --- conservation: align each consensus to copies via profile on ungapped cut ---
    # For display, use mafft profile against copies for conse_35 winner
    profiles = {}
    for key in keys:
        seq = results[key]["seq"]
        if not seq:
            continue
        # project copies onto consensus length by simple left-align (copies already
        # element-cut; lengths differ — use mafft add if needed)
        prof_fa = os.path.join(work, "prof_%s.fa" % key)
        write_fa(prof_fa, [(key, seq)] + copies[:20])  # subsample 20 for speed
        prof_aln = prof_fa + ".aln"
        sh("mafft --retree 1 --maxiterate 0 --quiet %s > %s 2>/dev/null"
           % (prof_fa, prof_aln), cwd=work)
        if os.path.isfile(prof_aln):
            pn, ps = read_fa(prof_aln)
            if ps:
                profiles[key] = conservation_from_alignment(ps[0], ps[1:])
            else:
                profiles[key] = nucfreq_profile(seq, copy_seqs)
        else:
            profiles[key] = nucfreq_profile(seq, copy_seqs)

    payload = {
        "tag": tag, "aln_path": aln_path, "n_copies": len(copies),
        "boundary": bdiag, "results": {k: {kk: vv for kk, vv in v.items()
                                          if kk != "seq"} for k, v in results.items()},
        "identity_matrix": mat,
    }
    # store sequences separately
    seq_path = os.path.join(out_dir, tag + "_consensus_seqs.json")
    with open(seq_path, "w") as fh:
        json.dump({k: v.get("seq") for k, v in results.items()}, fh, indent=2)
    with open(os.path.join(out_dir, tag + "_consensus_test.json"), "w") as fh:
        json.dump(payload, fh, indent=2)

    # --- HTML ---
    html_path = os.path.join(out_dir, tag + "_consensus_test.html")
    rows_html = []
    for key in keys:
        v = results[key]
        extra = ""
        if v.get("pct"):
            extra = " @ %d%%" % v["pct"]
        if v.get("iters"):
            extra += " (%d iters)" % v["iters"]
        rows_html.append(
            "<tr><td class='f'>%s%s</td><td class='n'>%d</td>"
            "<td class='f'>%s</td></tr>" % (
                html.escape(key), html.escape(extra), v["len"],
                html.escape(v.get("error") or "")))

    mat_html = []
    for a in keys:
        cells = ["<td class='f'>%s</td>" % html.escape(a)]
        for b in keys:
            cells.append("<td class='n'>%.3f</td>" % mat[a][b])
        mat_html.append("<tr>" + "".join(cells) + "</tr>")

    prof_html = []
    for key in keys:
        if key not in profiles:
            continue
        prof_html.append("<h3>%s</h3>%s" % (
            html.escape(key), svg_conservation(profiles[key], key)))

    bd = bdiag
    note = (
        "Element window from copy-supported boundary walk on <code>%s</code> "
        "(%d copies). Consensus row: <code>%s</code>. "
        "Old consensus span %d bp &rarr; copy-supported span %d bp "
        "(extended 5′ %d, 3′ %d)." % (
            html.escape(tag), len(copies), html.escape(names[k]),
            bd.get("old_span", 0), bd.get("new_span", 0),
            bd.get("extended_left", 0), bd.get("extended_right", 0)))

    page = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>%s — consensus tool comparison</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@600&display=swap">
<style>
:root{--ground:#f4f6f2;--surface:#fff;--surface-2:#eceee8;--ink:#131c1a;--muted:#68766f;
--rule:#d6dcd4;--accent:#1f6f5c;}
body{margin:0;background:var(--ground);color:var(--ink);font-family:"IBM Plex Sans",sans-serif;
font-size:15px;line-height:1.55}
.wrap{max-width:1100px;margin:0 auto;padding:32px 24px 80px}
h1{font-family:"IBM Plex Serif",serif;font-size:28px;margin:0 0 12px}
h2{font-size:18px;margin:36px 0 10px}
h3{font-size:14px;margin:20px 0 6px;color:var(--muted)}
.note{background:var(--surface);border-left:3px solid var(--accent);padding:12px 16px;margin:20px 0}
table{border-collapse:collapse;width:100%%;font-size:13px;background:var(--surface);
border:1px solid var(--rule);margin:16px 0}
th,td{padding:8px 10px;border-bottom:1px solid var(--rule)}
th{text-align:left;font-family:"IBM Plex Mono",monospace;font-size:10px;text-transform:uppercase;color:var(--muted)}
td.n{font-family:"IBM Plex Mono",monospace;text-align:right}
td.f{font-family:"IBM Plex Mono",monospace;font-size:12px;text-align:left}
code{font-family:"IBM Plex Mono",monospace;background:var(--surface-2);padding:1px 4px}
</style></head><body><div class="wrap">
<h1>%s — consensus rebuild comparison</h1>
<p class="note">%s</p>
<p>Same input for <code>conse</code> and <code>sine_consensus.sh</code>: ungapped element
sequences from top100 copies. Peel reference: conse 50%% on loci rebuild MSA when provided.</p>
<h2>Results</h2>
<table><thead><tr><th>method</th><th>length</th><th>notes</th></tr></thead>
<tbody>%s</tbody></table>
<h2>Pairwise identity (ungapped global)</h2>
<table><thead><tr><th>method</th>%s</tr></thead><tbody>%s</tbody></table>
<h2>Per-column conservation (consensus base fraction)</h2>
%s
<p style="color:var(--muted);font-size:12px;margin-top:40px">Generated on server by
<code>flank_border_consensus_test.py</code>. Sequences in companion JSON.</p>
</div></body></html>""" % (
        html.escape(tag), html.escape(tag), note, "".join(rows_html),
        "".join("<th>%s</th>" % html.escape(k) for k in keys),
        "".join(mat_html), "".join(prof_html))

    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(page)

    print("Wrote", html_path)
    print("Wrote", seq_path)
    for key in keys:
        print("  %-22s %4d bp" % (key, results[key]["len"]))


if __name__ == "__main__":
    main()
