"""report_hierarchy.py - the "Element hierarchy" section of the SINEderella report.

Reads RUN/flankscan/hierarchy.tsv (flankscan stage 7) and draws every composite element to scale as
its parts: part 1, the gap (linker or unexplained middle), part 2, each at its positions in its own
consensus. Each element shows its copies, full-unit share and verdict; the element's name and each part
link to the plates of that family in this report when a plate base URL is given. Returns "" when the
run has no hierarchy.tsv, so runs without composites are unchanged.
"""
import csv
import html
from pathlib import Path
from urllib.parse import quote

PALETTE = ["#3b6ea5", "#7a5aa6", "#2f8a6d", "#b0643a", "#b08a22", "#5e7d3a", "#8a4f6e",
           "#35798a", "#4f6a8f", "#9c5b5b", "#6b7280", "#8a7a3a"]
PX, X0, ROWH, BARH = 1.5, 230, 52, 20


def _plate(fam, code, aln_base, sfnames):
    """MSA-viewer link to the top100 plate of the report subfamily matching the short name fam."""
    if not aln_base or not code:
        return None
    sf = next((s for s in sfnames if s == fam or s.startswith(fam + "_")), None)
    if sf is None:
        return None
    raw = aln_base.rstrip("/") + "/%s_%s_top100.aln.fa" % (code, sf)
    return "https://toki-bio.github.io/MSA-viewer/?url=%s&title=%s" % (quote(raw, safe=""), quote(sf))


def _segs(r):
    """[(unit, start, end) | (None, gap_bp, 0)] from the parts column (any number of units); older two-part
    tables without it fall back to part1 / gap / part2."""
    out = []
    if r.get("parts"):
        for e in r["parts"].split(";"):
            if e.startswith("gap:"):
                out.append((None, int(e[4:]), 0))
            else:
                u, span = e.split(":")
                a, b = span.split("-")
                out.append((u, int(a), int(b)))
        return out
    out.append((r["part1"], int(r["part1_start"]), int(r["part1_end"])))
    if int(r["gap"]) > 0:
        out.append((None, int(r["gap"]), 0))
    out.append((r["part2"], int(r["part2_start"]), int(r["part2_end"])))
    return out



GREY = "#c9c9c9"                     # stretches of a composite that no unit explains
UNIT_PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948", "#b07aa1", "#ff9da7",
                "#9c755f", "#86bcb6", "#d37295", "#a0cbe8", "#8cd17d", "#f1ce63", "#bab0ac"]


def _composition(run_root):
    """step2_output/composite_units.tsv (tools/composite_units.py): consensus -> (kind, [(start, end, unit, id, alts)])."""
    ss = sorted(Path(run_root).glob("step2/step2_output*/composite_units.tsv"), key=lambda q: q.stat().st_mtime)
    if not ss:
        return {}
    out = {}
    for r in csv.DictReader(open(ss[-1], encoding="utf-8"), delimiter="\t"):
        k = out.setdefault(r["consensus"], [r["kind"], []])
        k[1].append((int(r["start"]), int(r["end"]), r["unit"], r["identity"], r["alternatives"]))
    return out


def _full(fam, names):
    """bank name for a short part name (hierarchy.tsv writes r6 for r6_210seqs)"""
    return fam if fam in names else next((n for n in names if n.startswith(fam + "_")), None)


def _short(u):
    return u.split("_")[0] if "_" in u and not u.startswith("MEG") else u


def _track(comp, cons, a, b, x0, y, ucol, h=8):
    """the unit composition of consensus positions a..b, drawn from x0 at row y (units coloured, no-unit grey)"""
    if cons not in comp:
        return ""
    kind, segs = comp[cons]
    out = []
    for s0, e0, u, pid, alts in segs:
        lo, hi = max(a, s0), min(b, e0)
        if lo > hi:
            continue
        x = x0 + (lo - a) * PX
        w = (hi - lo + 1) * PX
        fill = GREY if u == "-" else ucol.get(u, "#888")
        tip = ("%s %d-%d: no unit explains this stretch" % (cons, lo, hi) if u == "-" else
               "%s %d-%d = %s (%s %%)%s" % (cons, lo, hi, u, pid, "" if alts == "-" else "; next: " + alts))
        out.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="%s"><title>%s</title></rect>'
                   % (x, y, w, h, fill, html.escape(tip)))
    return "".join(out)


def _bank_composites(comp, ucol):
    """one bar per composite of the searched bank, drawn as its units (labelled) and its no-unit stretches (grey)"""
    cs = [c for c, (k, _) in comp.items() if k == "composite"]
    if not cs:
        return ""
    cs.sort(key=lambda c: -max(e for _, e, _, _, _ in comp[c][1]))
    maxbp = max(max(e for _, e, _, _, _ in comp[c][1]) for c in cs)
    W = int(X0 + maxbp * PX + 30)
    RH = 40
    H = RH * len(cs) + 26
    s = ['<svg viewBox="0 0 %d %d" width="%d" role="img" aria-label="bank composites as units" '
         'style="font-family:-apple-system,Segoe UI,Roboto,sans-serif">' % (W, H, W)]
    for t in range(0, maxbp + 1, 50):
        x = X0 + t * PX
        s.append('<line x1="%.1f" y1="16" x2="%.1f" y2="%d" stroke="#eee"/><text x="%.1f" y="11" font-size="10" '
                 'fill="#888" text-anchor="middle">%d</text>' % (x, x, H - 4, x, t))
    y = 24
    for c in cs:
        L = max(e for _, e, _, _, _ in comp[c][1])
        s.append('<text x="0" y="%d" font-size="12.5" font-weight="600" fill="#222">%s</text>' % (y + 14, html.escape(c)))
        s.append('<text x="0" y="%d" font-size="10.5" fill="#666">%d bp</text>' % (y + 27, L))
        for s0, e0, u, pid, alts in comp[c][1]:
            x = X0 + (s0 - 1) * PX
            w = (e0 - s0 + 1) * PX
            fill = GREY if u == "-" else ucol.get(u, "#888")
            tip = ("%d-%d: no unit explains this stretch" % (s0, e0) if u == "-" else
                   "%d-%d = %s (%s %%)%s" % (s0, e0, u, pid, "" if alts == "-" else "; next: " + alts))
            s.append('<rect x="%.1f" y="%d" width="%.1f" height="20" rx="2" fill="%s"><title>%s</title></rect>'
                     % (x, y, w, fill, html.escape(tip)))
            if w > 34:
                lab = "?" if u == "-" else "%s %s" % (_short(u), pid if w > 70 else "")
                s.append('<text x="%.1f" y="%d" font-size="10.5" font-weight="600" fill="%s" text-anchor="middle">%s</text>'
                         % (x + w / 2, y + 14, "#555" if u == "-" else "#fff", html.escape(lab.strip())))
        y += RH
    s.append("</svg>")
    leg = " ".join('<span style="display:inline-block;padding:1px 6px;margin:1px;border-radius:3px;background:%s;color:#fff;'
                   'font-size:11px">%s</span>' % (col, html.escape(u)) for u, col in ucol.items())
    leg += (' <span style="display:inline-block;padding:1px 6px;margin:1px;border-radius:3px;background:%s;color:#333;'
            'font-size:11px">no unit</span>' % GREY)
    return ("<h3 style='margin-top:18px'>Composites of the searched bank, as units</h3>"
            "<p class='intro'>Each composite consensus of the bank drawn as the units it is made of: at every position the "
            "unit with the most similar stretch (ssearch36, >= 75 % over >= 50 bp); <b>grey = no unit explains it</b>. The "
            "number is the identity; hover a block for the next most similar units (e.g. r1_r3_P1 201-401 = r3, next r4 at "
            "97.8 %). Composite = holds a near-complete copy of another consensus (>= 70 % of it at >= 90 %) and >= 120 bp "
            "more. Source: <code>composite_units.tsv</code> (tools/composite_units.py).</p><p>" + leg + "</p>"
            "<div style='overflow-x:auto'>" + "\n".join(s) + "</div>")


def section(run_root, species_code=None, aln_base=None, subfamilies=()):
    comp = _composition(run_root)
    names = list(comp)
    units = sorted({u for _, segs in comp.values() for _, _, u, _, _ in segs if u != "-"})
    ucol = {u: UNIT_PALETTE[i % len(UNIT_PALETTE)] for i, u in enumerate(units)}
    bank_block = _bank_composites(comp, ucol)
    p = Path(run_root) / "flankscan" / "hierarchy.tsv"
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")) if p.is_file() else []
    if not rows:
        if not bank_block:
            return ""
        return ("<section class='card' id='hierarchy'><h2>Element hierarchy</h2>" + bank_block + "</section>")
    global ROWH
    ROWH = 64 if comp else 52
    rows.sort(key=lambda r: (r["verdict"] != "accept", -int(r["copies"] or 0)))
    for r in rows:
        r["_segs"] = _segs(r)
    fams = sorted({f for r in rows for f, _, _ in r["_segs"] if f})
    col = {f: PALETTE[i % len(PALETTE)] for i, f in enumerate(fams)}
    maxbp = max(sum((b - a + 1) if f else a for f, a, b in r["_segs"]) for r in rows)
    W = int(X0 + maxbp * PX + 30)
    H = ROWH * len(rows) + 30
    s = ['<svg viewBox="0 0 %d %d" width="%d" role="img" aria-label="element hierarchy, parts to scale" '
         'style="font-family:-apple-system,Segoe UI,Roboto,sans-serif">' % (W, H, W)]
    for t in range(0, maxbp + 1, 50):
        x = X0 + t * PX
        s.append('<line x1="%.1f" y1="16" x2="%.1f" y2="%d" stroke="#ddd"/><text x="%.1f" y="11" font-size="10" '
                 'fill="#888" text-anchor="middle">%d</text>' % (x, x, H - 4, x, t))
    y = 28
    for r in rows:
        name = r["bank_name"]
        link = _plate(name, species_code, aln_base, subfamilies)
        lab = html.escape(name)
        if link:
            lab = '<a href="%s" target="_blank">%s</a>' % (html.escape(link), lab)
        chip = "#2f7d4f" if r["verdict"] == "accept" else ("#a23b3b" if r["verdict"] == "open" else "#a8631a")
        s.append('<text x="0" y="%d" font-size="12.5" font-weight="600" fill="#222">%s</text>' % (y + 13, lab))
        s.append('<text x="0" y="%d" font-size="10.5" fill="%s">%s &#183; %s copies &#183; %s%% full &#183; %s bp</text>'
                 % (y + 27, chip, html.escape(r["verdict"]), r["copies"], r["pct_full"], r["cons_len"]))
        x = X0
        segs = r["_segs"]
        op5 = r.get("open5_unit", "-"); op3 = r.get("open3_unit", "-")
        if op5 not in ("-", ""):        # the element continues into a known unit beyond this end (stage 6b)
            s.append('<rect x="%.1f" y="%d" width="36" height="%d" rx="2" fill="none" stroke="#a23b3b" stroke-dasharray="4 3">'
                     '<title>continues into %s at the 5 end</title></rect><text x="%.1f" y="%d" font-size="10" fill="#a23b3b" '
                     'text-anchor="middle">%s?</text>' % (x - 40, y, BARH, html.escape(op5), x - 22, y + 14, html.escape(op5)))
        for fam, a, b in segs:
            if fam is None:
                w = a * PX
                s.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#d9d3c5"/>' % (x, y + 5, w, BARH - 10))
                s.append('<text x="%.1f" y="%d" font-size="10" fill="#666" text-anchor="middle">%d bp</text>' % (x + w / 2, y + BARH + 11, a))
                x += w
                continue
            w = (b - a + 1) * PX
            txt = "%s %d-%d" % (fam, a, b)
            plink = _plate(fam, species_code, aln_base, subfamilies)
            rect = '<rect x="%.1f" y="%d" width="%.1f" height="%d" rx="2" fill="%s"><title>%s</title></rect>' % (x, y, w, BARH, col[fam], html.escape(txt))
            if w > 50:
                rect += '<text x="%.1f" y="%d" font-size="10.5" font-weight="600" fill="#fff" text-anchor="middle">%s</text>' % (x + w / 2, y + 14, html.escape(txt))
            s.append('<a href="%s" target="_blank">%s</a>' % (html.escape(plink), rect) if plink else rect)
            full = _full(fam, names)
            if full:                       # the part as units: coloured units, grey = no unit
                s.append(_track(comp, full, a, b, x, y + BARH + 16, ucol))
            x += w
        if op3 not in ("-", ""):
            s.append('<rect x="%.1f" y="%d" width="36" height="%d" rx="2" fill="none" stroke="#a23b3b" stroke-dasharray="4 3">'
                     '<title>continues into %s at the 3 end</title></rect><text x="%.1f" y="%d" font-size="10" fill="#a23b3b" '
                     'text-anchor="middle">%s?</text>' % (x + 4, y, BARH, html.escape(op3), x + 22, y + 14, html.escape(op3)))
        y += ROWH
    s.append("</svg>")
    return ("<section class='card' id='hierarchy'><h2>Element hierarchy</h2>"
            "<p class='intro'>Elements built from two or more consensus units (a chain has three or more), found by flankscan from the junctions "
            "of the assigned copies: each part drawn to scale at its positions in its own consensus, grey = "
            "sequence neither unit covers (linker or middle). <b>accept</b> = at least 70 % of the layout's "
            "copies read as one full unit when the element is added to the bank; <b>open</b> = the element "
            "continues into a known unit beyond an end (dashed box), so it is the middle of a longer chain. The element name opens its "
            "own plate, each part the plate of that family. Source: <code>flankscan/hierarchy.tsv</code>."
            + (" The thin track under each part shows it as units (colours as below; <b>grey = no unit explains it</b>)."
               if comp else "") + "</p>"
            "<div style='overflow-x:auto'>" + "\n".join(s) + "</div>" + bank_block + "</section>")
