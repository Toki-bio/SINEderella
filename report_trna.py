"""report_trna.py - the "Where the heads come from" section of the SINEderella report: tRNA / 5S / 7SL origin of every
consensus and of every stretch of its composites, from <step2_output>/trna/trna_heads.tsv (tools/trna_heads.sh).
Returns "" when the run has no such table."""
import csv
import html
from pathlib import Path


def _chip(label, kind, tip=""):
    col = {"ok": ("#e2efe9", "#1f6f5c"), "edge": ("#f5eed8", "#7a5c12"), "muted": ("#eceee8", "#68766f")}[kind]
    t = " title='%s'" % html.escape(tip, quote=True) if tip else ""
    return ("<span style='display:inline-block;padding:1px 7px;border-radius:9px;font-size:12px;background:%s;color:%s'%s>%s</span>"
            % (col[0], col[1], t, html.escape(label)))


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _calls(r):
    """every head of a query: [(label, kind, tip, span, tRNA-hit?)]. tRNA when the gene hit beats every shuffle (passing
    SINEBase's 60 % / 60 nt rule, else '?'); 5S / 7SL / 7SK from Rfam, with Rfam's span. A query can carry both
    (MEG-TR: tRNA-Val 1-55 + 5S 65-180)."""
    out = []
    bits, null = _num(r.get("bits")), _num(r.get("null_max_bits"))
    if r.get("best_tRNA", "-") != "-" and bits is not None and null is not None and bits > null:
        aa = r["best_tRNA"].split("-")[0]
        tip = "%s %% over %s nt (query %s) to tRNA-%s; %s bits vs at most %s for the shuffles; next amino acid %s" % (
            r.get("ident%"), r.get("aln_bp"), r.get("q_span"), r["best_tRNA"], r.get("bits"), r.get("null_max_bits"),
            r.get("2nd_isotype"))
        if r.get("SINEBase_60/60") == "yes":
            out.append(("tRNA-%s" % aa, "ok", tip, r.get("q_span", "-"), True))
        else:
            out.append(("tRNA-%s?" % aa, "edge", tip + " (below SINEBase's 60 % / 60 nt)", r.get("q_span", "-"), True))
    for hit in (r.get("Rfam", "-") or "-").split("; "):
        for rna, lab in (("5S_rRNA", "5S rRNA"), ("SRP", "7SL RNA"), ("7SK", "7SK RNA")):
            if hit.startswith(rna):
                f = hit.split()
                out.append((lab, "ok", "Rfam: " + hit, f[1] if len(f) > 1 else "-", False))
    if not out:
        out.append(("none detected", "muted", "no tRNA gene hit above the shuffled sequences, no 5S / 7SL in Rfam", "-", False))
    return out


def section(run_root):
    ss = sorted(Path(run_root).glob("step2/step2_output*/trna/trna_heads.tsv"), key=lambda q: q.stat().st_mtime)
    if not ss:
        return ""
    rows = list(csv.DictReader(open(ss[-1], encoding="utf-8"), delimiter="\t"))
    if not rows:
        return ""
    whole = [r for r in rows if ":" not in r["query"]]
    parts = {}
    for r in rows:
        if ":" in r["query"]:
            parts.setdefault(r["query"].split(":")[0], []).append(r)
    def cells(r, offset=0):
        """Head chips, their spans (in the consensus's coordinates), and identity / bits of the tRNA hit when there is one"""
        cs = _calls(r)
        spans = []
        for lab, kind, tip, sp, is_trna in cs:
            if sp not in ("-", "") and offset:
                try:
                    a, b = (int(x) for x in sp.split("-"))
                    sp = "%d-%d" % (offset + a - 1, offset + b - 1)
                except ValueError:
                    pass
            spans.append(sp)
        heads = "<br>".join(_chip(l, k, t) for l, k, t, _, _ in cs)
        trna = any(c[4] for c in cs)
        ident = r.get("ident%", "-") if trna else "-"
        bits = ("%s / %s" % (r.get("bits", "-"), r.get("null_max_bits", "-"))) if trna else ("- / %s" % r.get("null_max_bits", "-"))
        return heads, "<br>".join(html.escape(s) for s in spans), html.escape(ident), html.escape(bits)

    body = []
    for r in whole:
        heads, spans, ident, bits = cells(r)
        body.append("<tr><td><code>%s</code></td><td class='num'>%s</td><td>%s</td><td>%s</td><td class='num'>%s</td>"
                    "<td class='num'>%s</td><td class='small'>%s</td><td class='small'>%s</td></tr>"
                    % (html.escape(r["query"]), r["len"], heads, spans, ident, bits,
                       html.escape(r.get("Rfam", "-")), html.escape(r.get("tRNAscan", "-")[:70])))
        for p in parts.get(r["query"], []):
            seg, unit = p["query"].split(":", 1)[1].split("|")
            heads, spans, ident, bits = cells(p, int(seg.split("-")[0]))      # part coordinates -> composite coordinates
            body.append("<tr style='opacity:.92'><td style='padding-left:22px' class='small'>&#8627; %s <span class='muted'>(%s)</span></td>"
                        "<td class='num'>%s</td><td>%s</td><td>%s</td><td class='num'>%s</td><td class='num'>%s</td>"
                        "<td class='small'>%s</td><td class='small'>%s</td></tr>"
                        % (html.escape(seg), html.escape(unit if unit != "-" else "no unit"), p["len"], heads, spans, ident, bits,
                           html.escape(p.get("Rfam", "-")), html.escape(p.get("tRNAscan", "-")[:70])))
    return ("<section class='card' id='trna'><h2>Where the heads come from: tRNA, 5S, 7SL</h2>"
            "<p class='intro'>Every consensus, and under each composite every stretch of it as read by its units, compared with "
            "<b>tRNA genes</b> (ssearch36; the best tRNA, identity and its span), against <b>20 dinucleotide-preserving shuffles</b> "
            "of the same sequence (bits / best shuffle), plus Rfam (Infernal cmscan) and tRNAscan-SE. A head is called "
            "<b>tRNA-derived</b> when the gene hit passes SINEBase's rule (&ge; 60 % identity over &ge; 60 nt) and beats every "
            "shuffle; <b>?</b> = beats the shuffles but below 60 / 60; 5S / 7SL from Rfam. Rfam and tRNAscan-SE miss degraded "
            "heads; their anticodon calls on such heads are noise (hover a chip for the measurements, including the next amino "
            "acid). Source: <code>trna/trna_heads.tsv</code> (tools/trna_heads.sh).</p>"
            "<div style='overflow-x:auto'><table class='tbl'><thead><tr><th>Consensus / part</th><th>bp</th><th>Head</th>"
            "<th>Span</th><th>Identity %</th><th>Bits / best shuffle</th><th>Rfam</th><th>tRNAscan-SE</th></tr></thead><tbody>"
            + "".join(body) + "</tbody></table></div></section>")
