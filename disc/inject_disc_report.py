#!/usr/bin/env python3
"""Inject SINE-discriminator verdict columns into a SINEderella report."""
import argparse
import html
import os
import re
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verdict import verdict  # noqa: E402
from flank_uniqueness import scan as flank_scan  # noqa: E402
from overall import overall_call  # noqa: E402
from fix_alignments import is_seed  # noqa: E402

MSA = "https://toki-bio.github.io/MSA-viewer/"

VCHIP_CSS = """
.vchip { display: inline-block; font-size: .78rem; font-weight: 600;
         padding: 2px 7px; border-radius: 3px; white-space: nowrap; }
.vchip.v-ok { background: #e2efe9; color: #1f6f5c; }
.vchip.v-edge { background: #f5eed8; color: #7a5c12; }
.vchip.v-warn { background: #f6e8de; color: #a8501d; }
.vchip.v-muted { background: #eceee8; color: #68766f; }
"""

# One styled tooltip for the whole page. Native title= tooltips appear only after a delay, vanish
# after a few seconds and never on touch; this moves every title= into data-tip on load and shows it
# at once on hover (the texts themselves are unchanged, so every generator's titles keep working).
TIP_MARK = "sd-tip-js"
TIP_HTML = """<style id="sd-tip-css">
#sd-tip { position: fixed; z-index: 9999; max-width: 460px; background: #1f2933; color: #fff;
          font-size: 12.5px; line-height: 1.45; padding: 7px 10px; border-radius: 5px;
          box-shadow: 0 4px 14px rgba(0,0,0,.25); pointer-events: none; white-space: pre-line;
          display: none; }
[data-tip] { cursor: help; }
a[data-tip] { cursor: pointer; }
</style>
<script id="sd-tip-js">
(function () {
  function init() {
    var tip = document.createElement('div'); tip.id = 'sd-tip'; document.body.appendChild(tip);
    document.querySelectorAll('[title]').forEach(function (el) {
      if (el.closest('svg')) return;
      el.setAttribute('data-tip', el.getAttribute('title')); el.removeAttribute('title');
    });
    function place(e) {
      var x = e.clientX + 14, y = e.clientY + 16, w = tip.offsetWidth, h = tip.offsetHeight;
      if (x + w > window.innerWidth - 8) x = Math.max(8, window.innerWidth - w - 8);
      if (y + h > window.innerHeight - 8) y = Math.max(8, e.clientY - h - 10);
      tip.style.left = x + 'px'; tip.style.top = y + 'px';
    }
    document.addEventListener('mouseover', function (e) {
      var el = e.target.closest && e.target.closest('[data-tip]');
      if (!el) { tip.style.display = 'none'; return; }
      tip.textContent = el.getAttribute('data-tip'); tip.style.display = 'block'; place(e);
    });
    document.addEventListener('mousemove', function (e) { if (tip.style.display === 'block') place(e); });
    document.addEventListener('scroll', function () { tip.style.display = 'none'; }, true);
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
</script>
"""

TOP_RULE = ("Top 100 by bitscore: the firmly assigned copies (10/10 votes, above the bitscore bar) with "
            "the highest bitscore to the consensus. If fewer than 100 are firmly assigned, the plate is "
            "filled up to 100 with soft-assigned copies (found by this query but not unanimous), ranked "
            "by search score and marked [soft] in the row name.")
RAND_RULE = ("100 random copies of the firmly assigned set, drawn with a fixed random seed (%s), so the same "
             "copies come back on every rebuild. If fewer than 100 are firmly assigned, all of them are used "
             "and the plate is filled up to 100 with random soft-assigned copies, marked [soft] in the row "
             "name." % os.environ.get("RAND_SEED", "42"))
# With flank verdicts (SINEderella tools/flank_verdict.py; plates drawn by tools/plate_filter.py) the plates hold SINE loci
# only: no soft-assigned copies, no copies with shared, masked or untested flanks (SINEderella docs/AIM.md, 2026-10-18).
TOP_RULE_V = ("Top 100 by bitscore among the SINE loci of this subfamily: assigned copies whose flanks were tested unique "
              "(flankscan stage 9; not twins, arrays or the same flank at the same place in many copies). If fewer than "
              "100 are SINE loci, all of them are shown, topped up with its grey-zone copies (flanks masked or untestable, "
              "marked [grey]) ranked by bitscore.")
RAND_RULE_V = ("100 random SINE loci of this subfamily (assigned copies whose flanks were tested unique), drawn with a fixed "
               "random seed (%s), so the same copies come back on every rebuild. If fewer than 100, all of them are shown, "
               "topped up with random grey-zone copies (flanks masked or untestable, marked [grey])."
               % os.environ.get("RAND_SEED", "42"))
ROWS_NOTE = ("Row 1 (<subfamily>_extended) is the consensus rebuilt from these copies. Its lowercase letters "
             "are proposals, none applied: past the original's ends they are proposed additions; at an end "
             "inside the original they are original bases the copies do not carry (proposed trims), kept so "
             "nothing of the original is lost. Row 2 (<subfamily>) is the consensus the genome was searched "
             "with; the uppercase span of the copies is its span.")


def plate_counts(path: Path):
    """(copies, soft, grey) on a published plate: rows other than the consensus and the original."""
    n = soft = grey = 0
    first = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        name = line[1:].strip()
        if first is None:
            first = name
            continue
        if is_seed(name, first):
            continue
        n += 1
        soft += "[soft]" in line
        grey += "[grey]" in line
    return n, soft, grey


def load_proposals(aln_dir: Path) -> dict:
    """proposals.tsv (SINEderella tools/add_seed_row.py): plate -> row."""
    p = aln_dir / "proposals.tsv"
    if not p.is_file():
        return {}
    lines = p.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t") if lines else []
    return {r[0]: dict(zip(head, r)) for r in (l.split("\t") for l in lines[1:]) if r}


def proposal_note(row) -> str:
    """One line per proposed edit, with the copies' support (ungapped identity; unrelated DNA ~0.25)."""
    if not row:
        return ""
    out = []
    for end, lab in (("5", "5'"), ("3", "3'")):
        a = int(row.get("add%s_bp" % end) or 0)
        t = int(row.get("trim%s_bp" % end) or 0)
        if a:
            out.append("%s: %d bp proposed in lowercase (copies' identity to them %s; ~0.25 is unrelated DNA)"
                       % (lab, a, row.get("add%s_ungapped" % end) or "n/a"))
        if t:
            out.append("%s: the copies do not carry the original's %d bp at this end (identity %s) - "
                       "trim proposed, not applied" % (lab, t, row.get("trim%s_support" % end) or "n/a"))
    return ("\nProposed edits to the original:\n- " + "\n- ".join(out)) if out else "\nNo edits proposed."

# Re-use chip helpers from inject_oma_aln_section (same logic)
from inject_oma_aln_section import (  # noqa: E402
    status_element,
    status_flank_context,
    status_flanks,
    _chip,
)


def msa_href(raw_base: str, filename: str, title: str) -> str:
    url = raw_base.rstrip("/") + "/" + filename
    return (
        f"{MSA}?url={quote(url, safe='')}"
        f"&title={quote(title, safe='')}"
    )


def aln_name(species: str, sf: str, kind: str) -> str:
    """step8a naming: {species}_{subfam}_{kind}.aln.fa"""
    prefix = f"{species}_" if not sf.startswith(f"{species}_") else ""
    return f"{prefix}{sf}_{kind}.aln.fa"


def load_copy_counts(summary_tsv: Path) -> dict:
    counts = {}
    if not summary_tsv.is_file():
        return counts
    for line in summary_tsv.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sf = parts[0].strip()
        try:
            counts[sf] = int(parts[1]) + int(parts[2])
        except ValueError:
            try:
                counts[sf] = int(parts[3])
            except (ValueError, IndexError):
                pass
    return counts


def discover_subfams(aln_dir: Path, species: str) -> list:
    out = set()
    for p in aln_dir.glob("*_top100.aln.fa"):
        name = p.name.replace("_top100.aln.fa", "")
        if name.startswith(f"{species}_"):
            name = name[len(species) + 1:]
        out.add(name)
    return sorted(out)


def score_top100(aln_dir: Path, species: str, sf: str):
    for fn in (aln_name(species, sf, "top100"), f"{sf}_top100.aln.fa"):
        path = aln_dir / fn
        if path.is_file():
            try:
                return verdict(str(path))
            except Exception as exc:
                return {"error": str(exc)}
    return None


def scan_flanks(aln_dir: Path, species: str, sf: str, tier: str):
    kind = "top100" if tier == "top100" else "rand100"
    for fn in (aln_name(species, sf, kind), f"{sf}_{kind}.aln.fa"):
        path = aln_dir / fn
        if path.is_file():
            try:
                return flank_scan(str(path))
            except Exception as exc:
                return {"error": str(exc)}
    return None


def load_flank_verdicts(path: Path) -> dict:
    """SINEderella flank_verdict.tsv: subfamily -> row dict (verdict, assigned_copies, sine, grey, not_sine, reason)."""
    if not path or not path.is_file():
        return {}
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t") if lines else []
    return {r[0]: dict(zip(head, r)) for r in (l.split("\t") for l in lines[1:]) if r and r[0]}


def load_cons_lengths(bank: Path) -> dict:
    """consensus FASTA -> name: length in bp (letters only)."""
    out, name = {}, None
    if not bank or not bank.is_file():
        return out
    for line in bank.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = 0
        elif name:
            out[name] += sum(c.isalpha() for c in line)
    return out


def load_families(path: Path) -> dict:
    """SINEderella family table (families/families.clean.tsv: subfamily<TAB>family[...]) -> subfamily: family."""
    if not path or not path.is_file():
        return {}
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        p = line.split("\t")
        if len(p) >= 2 and p[0] and not line.startswith("#") and p[0] != "subfamily":
            out[p[0].strip()] = p[1].strip()
    return out


GREY_WORDS = {"untested": "the subfamily's flank uniqueness could not be tested",
              "not_tested": "the flank test was not run",
              "not_interspersed": "the subfamily's flanks are still shared at the expansion limit (not interspersed)"}


def _int(x) -> int:
    try:
        return int(x or 0)
    except ValueError:
        return 0


def build_section(species: str, subfams, aln_dir: Path,
                  copy_counts: dict, raw_base: str, fverdict: dict = None,
                  cons_len: dict = None, families: dict = None) -> str:
    dash = "<span class='muted small'>&mdash;</span>"
    props = load_proposals(aln_dir)
    fverdict = fverdict or {}
    cons_len = cons_len or {}
    families = families or {}
    top_rule, rand_rule = (TOP_RULE_V, RAND_RULE_V) if fverdict else (TOP_RULE, RAND_RULE)
    rows = []
    fam_of = lambda s: families.get(s) or families.get(f"{species}_{s}") or ""
    order = sorted(subfams, key=lambda s: (fam_of(s) == "", fam_of(s), s)) if families else sorted(subfams)
    ncols = 13 if fverdict else 10
    fam_seen = None
    for sf in order:
        t100 = aln_name(species, sf, "top100")
        r100 = aln_name(species, sf, "rand100")
        sub = aln_name(species, sf, "subfam")
        fv = fverdict.get(sf) or fverdict.get(f"{species}_{sf}")
        if families and fam_of(sf) != fam_seen:
            # a family header row with the family's totals (the run's family table, families/families.clean.tsv)
            fam_seen = fam_of(sf)
            members = [s for s in order if fam_of(s) == fam_seen]
            tot = {k: sum(_int((fverdict.get(s) or fverdict.get(f"{species}_{s}") or {}).get(k)) for s in members)
                   for k in ("assigned_copies", "sine", "grey", "not_sine")}
            name = html.escape(fam_seen) if fam_seen else "no family in the table"
            rows.append(f"<tr class='famrow'><td colspan='{ncols}'><strong>{name}</strong>"
                        f" &nbsp;<span class='small'>{len(members)} subfamil{'y' if len(members) == 1 else 'ies'}"
                        + (f"; {tot['assigned_copies']:,} assigned copies, {tot['sine']:,} SINE loci, "
                           f"{tot['grey']:,} grey, {tot['not_sine']:,} not SINE" if fverdict else "")
                        + "</span></td></tr>")
        cl = cons_len.get(sf) or cons_len.get(f"{species}_{sf}")
        len_cell = f"<td class='num'>{cl:,}</td>" if cl else f"<td>{dash}</td>"
        if fv:
            # SINE loci highlighted; assigned, grey and not SINE in their own columns
            tip = ("%s assigned copies: %s SINE loci (flanks tested unique), %s grey zone (flanks masked or untestable), "
                   "%s not SINE (flanks shared, or the same flank at the same place in many copies). Subfamily verdict: %s."
                   % (fv.get("assigned_copies", "?"), fv.get("sine", "?"), fv.get("grey", "?"), fv.get("not_sine", "?"),
                      fv.get("verdict", "?")))
            ns = _int(fv.get("sine"))
            n_cell = (len_cell + f"<td class='num'>{_int(fv.get('assigned_copies')):,}</td>"
                      f"<td class='hl num' title='{html.escape(tip, quote=True)}'>{ns:,}</td>"
                      f"<td class='num'>{_int(fv.get('grey')):,}</td><td class='num'>{_int(fv.get('not_sine')):,}</td>")
            if not (aln_dir / t100).is_file():
                # no SINE locus, no plate: the row stays, with the reason
                why = GREY_WORDS.get(fv.get("verdict"), "no copy with flanks tested unique")
                label = "Not SINE" if fv.get("verdict") == "not_interspersed" else "Grey zone"
                rows.append(
                    f"<tr><td><code>{html.escape(sf)}</code></td>{n_cell}"
                    f"<td>{dash}</td><td>{dash}</td><td>{dash}</td><td>{dash}</td><td>{dash}</td><td>{dash}</td>"
                    f"<td>{_chip(label, 'warn' if label == 'Not SINE' else 'muted', 'No plate: ' + why + '. ' + tip + ' ' + fv.get('reason', ''))}</td></tr>")
                continue
        else:
            n = copy_counts.get(sf) or copy_counts.get(f"{species}_{sf}")
            n_cell = len_cell + (
                f"<td class='hl'>{n:,}</td>" if n is not None
                else "<td class='muted small'>&mdash;</td>"
            )
        def link(fn, label, css="", rule=""):
            path = aln_dir / fn
            if not path.is_file():
                return dash
            tip = ""
            if rule:
                n_c, n_soft, n_grey = plate_counts(path)
                if n_c < 100:
                    label = label.replace("100", str(n_c))
                if n_soft:
                    label += f" ({n_soft} soft)"
                if n_grey:
                    label += f" ({n_grey} grey)"
                tip = ((f"{rule}\n\nThis plate: {n_c - n_grey} SINE loci + {n_grey} grey-zone copies (marked [grey]) = {n_c} copies."
                        if fverdict else
                        f"{rule}\n\nThis plate: {n_c - n_soft} firm + {n_soft} soft = {n_c} copies.")
                       + ("" if n_c >= 100 else " The run has no more copies of this subfamily.")
                       + "\n" + ROWS_NOTE + proposal_note(props.get(fn)))
            href = msa_href(raw_base, fn, f"{species} {sf} {label}") if raw_base else fn
            cls = "aln-link" + (f" {css}" if css else "")
            t = f" title='{html.escape(tip, quote=True)}'" if tip else ""
            return f"<a class='{cls}' href='{href}' target='_blank'{t}>{label}</a>"
        v = score_top100(aln_dir, species, sf)
        top_scan = scan_flanks(aln_dir, species, sf, 'top100')
        rand_scan = scan_flanks(aln_dir, species, sf, 'rand100')
        oc = overall_call(v, top_scan, rand_scan)
        if fv and fv.get("verdict") != "interspersed":
            why = GREY_WORDS.get(fv.get("verdict"), "flanks not tested")
            oc = {"label": "Not SINE" if fv.get("verdict") == "not_interspersed" else "Grey zone",
                  "kind": "warn" if fv.get("verdict") == "not_interspersed" else "muted",
                  "reasons": ["Subfamily flank verdict: " + why + ". The plates show its grey-zone copies [grey].",
                              "Discriminator call on the plate: " + oc["label"] + "."] + list(oc.get("reasons") or [])}
        rows.append(
            f"<tr><td><code>{html.escape(sf)}</code></td>{n_cell}"
            f"<td>{link(t100, 'top 100', '', top_rule)}</td>"
            f"<td>{link(r100, '100 random', 'orange', rand_rule)}</td>"
            f"<td>{link(sub, 'SubFam', 'green')}</td>"
            f"<td>{status_flanks(v)}</td>"
            f"<td>{status_flank_context(top_scan, rand_scan)}</td>"
            f"<td>{status_element(v)}</td>"
            f"<td>{_chip(oc['label'], oc['kind'], chr(10).join(oc['reasons']))}</td></tr>"
        )
    # consensus alignments written by SINEderella publish/align_for_publish.sh (2026-10-18)
    cl = []
    for fn, lab in ((f"{species}_consensuses.aln.fa", "All consensuses"),
                    (f"{species}_consensuses_sine.aln.fa", "Consensuses of SINE-confirmed subfamilies")):
        if (aln_dir / fn).is_file():
            k = (aln_dir / fn).read_text(encoding="utf-8", errors="replace").count(">")
            href = msa_href(raw_base, fn, f"{species} {lab}") if raw_base else fn
            cl.append(f"<a class='aln-link' href='{href}' target='_blank'>{lab} ({k})</a>")
    cons_links = ("<p><strong>Consensus alignments:</strong> " + " &nbsp; ".join(cl) + "</p>") if cl else ""
    intro = (
        "<p class='intro'>One row per subfamily. "
        "<strong>Links</strong> open 100-copy MAFFT alignments (50&nbsp;bp left, "
        "70&nbsp;bp right flanks) in the MSA viewer: the <strong>consensus rebuilt from the copies is "
        "row&nbsp;1</strong>, the consensus as searched is row&nbsp;2. "
        + ("The alignments hold <strong>SINE loci</strong> (assigned copies whose flanks were tested unique) and, when "
           "a subfamily has fewer than 100 of them, its <strong>grey-zone</strong> copies (flanks masked as repeats or "
           "untestable), marked <code>[grey]</code>; the button then shows their count. Copies with shared flanks, or the "
           "same flank at the same place in many copies, and soft-assigned copies are never shown. The <strong>SINE "
           "loci</strong> column counts SINE loci (hover for the split into SINE / grey zone / not SINE); a subfamily "
           "whose flanks could not be tested keeps <strong>Grey zone</strong> as its overall call. "
           if fverdict else
           "A subfamily with fewer than 100 firmly assigned copies is filled up with soft-assigned ones, marked "
           "<code>[soft]</code>; the button then shows the count (hover for the rule). ")
        + "Last four columns: automated <strong>SINE-discriminator</strong> checks "
        "on the top-100 view (<strong>Flanks</strong>, <strong>Flank context</strong>, "
        "<strong>Element</strong>, <strong>Overall</strong>).</p>"
        + cons_links +
        "<details class='legend'><summary>How these alignments were built</summary>"
        "<ul class='small muted'>"
        "<li><code>publish/align_for_publish.sh</code>: step7 → border loop → step8a.</li>"
        "<li><code>rebuild_consensus_row.py</code>, <code>boundary_justify.py</code>, "
        "<code>trim_display_flanks.py</code> (SINE-discriminator).</li>"
        "<li>Verdict chips: <code>verdict.py</code> on 50L/70R publish geometry.</li>"
        "</ul></details>"
    )
    thead = (
        "<th>Subfamily</th><th title='Length of the consensus as searched (bp)'>Consensus bp</th>"
        + ("<th title='Copies assigned to this subfamily by the vote'>Assigned</th>"
           "<th title='SINE loci: assigned copies whose flanks were tested unique'>SINE loci</th>"
           "<th title='Grey zone: flanks masked as repeats or untestable, or the subfamily could not be tested'>Grey</th>"
           "<th title='Not SINE: flanks shared with another copy, the same flank at the same place in many copies, "
           "or the subfamily flanks still shared at the expansion limit'>Not SINE</th>"
           if fverdict else "<th>Copies</th>")
        + f"<th title='{html.escape(top_rule, quote=True)}'>Top 100 by bitscore</th>"
        f"<th title='{html.escape(rand_rule, quote=True)}'>100 random copies "
        f"(seed {html.escape(os.environ.get('RAND_SEED', '42'))})</th>"
        "<th>SubFam (chunk consensuses)</th>"
        "<th>Flanks</th><th>Flank context</th><th>Element</th>"
        "<th title='Verdict score on the top-100 plate. Below 30 copies no score is given: the call is "
        "Doubtful if the copies show negative evidence, Cannot assess if they show none. Hover a chip "
        "for its reasons.'>Overall</th>"
    )
    return (
        "<section class='card' id='alignments'>"
        "<h2>Subfamily alignments</h2>"
        + intro
        + "<table class='tbl'><thead><tr>" + thead + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


FAMROW_CSS = "\n.famrow td { background: #e8eef5; border-top: 2px solid #9fb3c8; }\n"


def ensure_vchip_css(text: str) -> str:
    if ".famrow" not in text:
        text = text.replace("</style>", FAMROW_CSS + "</style>", 1)
    if ".vchip" in text:
        return text
    return text.replace("</style>", VCHIP_CSS + "</style>", 1)


def ensure_tooltip(text: str) -> str:
    if TIP_MARK in text:
        return text
    i = text.rfind("</body>")
    return text + TIP_HTML if i == -1 else text[:i] + TIP_HTML + text[i:]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("report", type=Path)
    ap.add_argument("aln_dir", type=Path)
    ap.add_argument("species", help="Species code prefix (ccr, oma, …)")
    ap.add_argument("--raw-base", default="",
                    help="Published raw URL prefix for MSA viewer links")
    ap.add_argument("--summary-tsv", type=Path, default=None,
                    help="summary.by_subfam.tsv for copy counts")
    ap.add_argument("--verdict-tsv", type=Path, default=None,
                    help="SINEderella flank_verdict.tsv: copies column = SINE loci, a row for every subfamily")
    ap.add_argument("--bank", type=Path, default=None, help="consensus bank as searched: the consensus length column")
    ap.add_argument("--families", type=Path, default=None,
                    help="SINEderella families/families.clean.tsv: rows grouped by family under a family row")
    args = ap.parse_args()

    fverdict = load_flank_verdicts(args.verdict_tsv) if args.verdict_tsv else {}
    subfams = sorted(set(discover_subfams(args.aln_dir, args.species))
                     | {sf[len(args.species) + 1:] if sf.startswith(args.species + "_") else sf
                        for sf, r in fverdict.items() if (r.get("assigned_copies") or "0") != "0"})
    if not subfams:
        sys.exit("no *_top100.aln.fa in %s" % args.aln_dir)

    counts = load_copy_counts(args.summary_tsv) if args.summary_tsv else {}
    text = args.report.read_text(encoding="utf-8")
    text = ensure_vchip_css(text)
    text = ensure_tooltip(text)
    new_sec = build_section(args.species, subfams, args.aln_dir,
                            counts, args.raw_base.rstrip("/"), fverdict,
                            load_cons_lengths(args.bank), load_families(args.families))
    pat = r'<section class=[\'"]card[\'"] id=[\'"]alignments[\'"].*?</section>'
    text, n = re.subn(pat, new_sec, text, count=1, flags=re.DOTALL)
    if n != 1:
        # insert before overview if no alignment section yet
        text, n = re.subn(
            r'(<section class="card" id="overview">)',
            new_sec + r"\1",
            text, count=1,
        )
    if n != 1:
        sys.exit("could not insert alignment section")
    args.report.write_text(text, encoding="utf-8")
    print("patched %d subfamilies → %s" % (len(subfams), args.report))


if __name__ == "__main__":
    main()
