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


def section(run_root, species_code=None, aln_base=None, subfamilies=()):
    p = Path(run_root) / "flankscan" / "hierarchy.tsv"
    if not p.is_file():
        return ""
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
    if not rows:
        return ""
    rows.sort(key=lambda r: (r["verdict"] != "accept", -int(r["copies"] or 0)))
    fams = sorted({r["part1"] for r in rows} | {r["part2"] for r in rows})
    col = {f: PALETTE[i % len(PALETTE)] for i, f in enumerate(fams)}
    maxbp = max(int(r["part1_end"]) + int(r["gap"]) + int(r["part2_end"]) - int(r["part2_start"]) + 1 for r in rows)
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
        segs = [(r["part1"], int(r["part1_start"]), int(r["part1_end"]))]
        if int(r["gap"]) > 0:
            segs.append((None, int(r["gap"]), 0))
        segs.append((r["part2"], int(r["part2_start"]), int(r["part2_end"])))
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
            x += w
        if op3 not in ("-", ""):
            s.append('<rect x="%.1f" y="%d" width="36" height="%d" rx="2" fill="none" stroke="#a23b3b" stroke-dasharray="4 3">'
                     '<title>continues into %s at the 3 end</title></rect><text x="%.1f" y="%d" font-size="10" fill="#a23b3b" '
                     'text-anchor="middle">%s?</text>' % (x + 4, y, BARH, html.escape(op3), x + 22, y + 14, html.escape(op3)))
        y += ROWH
    s.append("</svg>")
    return ("<section class='card' id='hierarchy'><h2>Element hierarchy</h2>"
            "<p class='intro'>Elements built from two consensus units, found by flankscan from the junctions "
            "of the assigned copies: each part drawn to scale at its positions in its own consensus, grey = "
            "sequence neither unit covers (linker or middle). <b>accept</b> = at least 70 % of the layout's "
            "copies read as one full unit when the element is added to the bank; <b>open</b> = the element "
            "continues into a known unit beyond an end (dashed box), so it is the middle of a longer chain. The element name opens its "
            "own plate, each part the plate of that family. Source: <code>flankscan/hierarchy.tsv</code>.</p>"
            "<div style='overflow-x:auto'>" + "\n".join(s) + "</div></section>")
