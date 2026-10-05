"""report_blocks.py - the "Sequence similarity between consensuses" section of the SINEderella report.

The sequence-part counterpart of the Element hierarchy: where the hierarchy draws composites as their parts from
the copy junctions (flankscan), this shows which stretches of one consensus are similar to which stretches of
another, computed from the consensus sequences alone (tools/consensus_blocks.py: ungapped blocks of >= 30 bp and
>= 80 % identity, both strands, simple repeats masked). A matrix gives, for every row consensus, the share of its
length that lies in blocks shared with the column consensus (directional: r9 in r7 is not r7 in r9); clicking a
cell draws the two consensuses with the shared blocks as ribbons and lists their coordinates.

section(run_root) returns "" when the run has no consensus file, numpy is missing or there are fewer than two
consensuses, so the report is unchanged then.
"""
import html
import json
import os
import re
import sys
from pathlib import Path

MINLEN, MINID, MAXN = 30, 80.0, 40


def lab(n):
    m = re.match(r"(.+?)_\d+seqs$", n)
    return m.group(1) if m else n


def _bank(run_root):
    for rel in ("consensuses.clean.fa", "results/consensuses.fa", "consensuses.fa"):
        p = Path(run_root) / rel
        if p.exists() and p.stat().st_size:
            return p
    return None


WORDS = {"TWO_VERSIONS": "two separate SINEs", "UNLINKED_ENDS": "one SINE, variable end", "SINGLE_MODE": "one length",
         "UNRESOLVED": "not decided", "NOT_TESTED": "not tested (too few copies)"}


def _length_versions(run_root):
    """table of results/length_variants/summary.tsv (tools/length_variants_run.py), '' if the run has none"""
    p = Path(run_root) / "results" / "length_variants" / "summary.tsv"
    if not p.exists():
        return ""
    lines = [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8")]
    if len(lines) < 2:
        return ""
    h = lines[0]
    rows = []
    for l in lines[1:]:
        r = dict(zip(h, l))
        rows.append("<tr><td>%s / %s</td><td class='n'>%s</td><td><b>%s</b></td><td>%s</td></tr>" % (
            html.escape(lab(r["short"])), html.escape(lab(r["long"])), html.escape(r.get("consensus_identity", "")),
            WORDS.get(r["verdict"], r["verdict"]), html.escape(r.get("notes", ""))))
    return ("<h3 style='margin-top:16px'>Is the shorter one a separate SINE?</h3><p>For every consensus that is the 5&prime; part of a longer one, the copies of both are tested: "
            "do they end at two separate places, do the bases inside follow the length, does each end have its own TSD (<code>docs/LENGTH_VARIANTS.md</code>).</p>"
            "<table class='tbl'><thead><tr><th>Pair (short / long)</th><th class='n'>Identity %</th><th>Answer</th><th>Note</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


AUDIT_WORDS = {"MATCH": "matches the bank", "SHORTER": "rebuild is shorter", "LONGER": "rebuild is longer", "DIVERGED": "differs from the bank",
               "UNSTABLE": "the two rebuilds differ", "SKIPPED": "too few copies", "FAILED": "rebuild failed"}


def _consensus_audit(run_root):
    """table of results/consensus_audit/summary.tsv (tools/consensus_audit.py), '' if the run has none"""
    p = Path(run_root) / "results" / "consensus_audit" / "summary.tsv"
    if not p.exists():
        return ""
    lines = [l.rstrip(chr(10)).split(chr(9)) for l in open(p, encoding="utf-8")]
    if len(lines) < 2:
        return ""
    h = lines[0]
    rows = []
    for l in lines[1:]:
        r = dict(zip(h, l))
        s1 = [k for k in h if k.startswith("rebuilt_bp_s")]
        m1 = [k for k in h if k.startswith("mismatches_s")]
        bp = " / ".join(r.get(k, "") or "-" for k in s1)
        mm = " / ".join(r.get(k, "") or "-" for k in m1)
        rows.append("<tr><td>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td><b>%s</b></td></tr>" % (
            html.escape(lab(r["family"])), html.escape(r.get("copies", "")), html.escape(r.get("bank_bp", "")), html.escape(bp), html.escape(mm),
            AUDIT_WORDS.get(r.get("verdict", ""), html.escape(r.get("verdict", "")))))
    return ("<h3 style='margin-top:16px'>Is each consensus what its copies say?</h3><p>Every consensus is rebuilt from the copies assigned to it (random subsamples, "
            "two seeds) and compared with the bank sequence (<code>docs/CONSENSUS_AUDIT.md</code>). A shorter or longer rebuild is not an error by itself: "
            "it says that fewer, or more, than 30 % of the copies carry a stretch.</p>"
            "<table class='tbl'><thead><tr><th>Consensus</th><th class='n'>Copies</th><th class='n'>Bank bp</th><th class='n'>Rebuilt bp (2 seeds)</th>"
            "<th class='n'>Mismatches</th><th>Answer</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _arrays(run_root):
    """table of results/array_flag.tsv (tools/array_flag.py), '' if the run has none"""
    p = Path(run_root) / "results" / "array_flag.tsv"
    if not p.exists():
        return ""
    lines = [l.rstrip(chr(10)).split(chr(9)) for l in open(p, encoding="utf-8")]
    if len(lines) < 2:
        return ""
    h = lines[0]
    rows = []
    for l in lines[1:]:
        r = dict(zip(h, l))
        flag = "<b>tandem array</b>" if r.get("flag") == "ARRAY" else "&ndash;"
        rows.append("<tr><td>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td>%s</td></tr>" % (
            html.escape(lab(r["family"])), html.escape(r["copies"]), html.escape(r["copies_in_arrays"]), html.escape(r["pct_in_arrays"]),
            html.escape(r.get("null_pct", "")), html.escape(r.get("excess_pct", "")), html.escape(r["arrays"]), html.escape(r["median_spacing_bp"]), flag))
    return ("<h3 style='margin-top:16px'>Do the copies sit in tandem arrays?</h3><p>Copies in a run of at least five on one contig with regular spacing "
            "(gaps up to 6 kb, within a factor of 5 of the run's median) are units of an array, not independent insertions. The share is compared with "
            "the same copies spread at random over the contigs (chance): a dense family forms regular runs by chance (VES in <i>Taphozous</i>, one copy per "
            "3 kb, 56 % of its copies in runs, nearly all chance). A family whose share exceeds chance by 20 points or more is marked, and its plates take "
            "independent copies first (<code>tools/array_flag.py</code>).</p>"
            "<table class='tbl'><thead><tr><th>Family</th><th class='n'>Copies</th><th class='n'>In arrays</th><th class='n'>In arrays %</th><th class='n'>Chance %</th>"
            "<th class='n'>Excess</th><th class='n'>Arrays</th><th class='n'>Median spacing bp</th><th>Answer</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _twins(run_root):
    """table of results/flank_twins.tsv (flankscan/fs9_run.sh after assignment), '' if the run has none"""
    p = Path(run_root) / "results" / "flank_twins.tsv"
    if not p.exists():
        return ""
    lines = [l.rstrip(chr(10)).split(chr(9)) for l in open(p, encoding="utf-8")]
    if len(lines) < 2:
        return ""
    h = lines[0]
    rows = []
    for l in lines[1:]:
        r = dict(zip(h, l))
        try:
            pct = float(r.get("pct_twin", "0") or 0)
            nt = int(r.get("twin1", "0") or 0) + int(r.get("twin2", "0") or 0)
            nmask = int(r.get("masked", "0") or 0) + int(r.get("untestable", "0") or 0)
            ncop = int(r.get("copies", "0") or 0)
        except ValueError:
            pct, nt, nmask, ncop = 0.0, 0, 0, 0
        if nt == 0:
            ans = "flanks unique" if nmask < 0.5 * ncop else "mostly untestable (repeat-masked flanks)"
        elif pct >= 20:
            ans = "<b>%d copies (%.0f %%) share flanks: not independent insertions, see the [twin] rows of the plates</b>" % (nt, pct)
        else:
            ans = "%d copies share flanks ([twin] on the plates); the rest are independent" % nt
        rows.append("<tr><td>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td>%s</td></tr>" % (
            html.escape(lab(r["family"])), html.escape(r["copies"]), html.escape(r.get("twin1", "")), html.escape(r.get("twin2", "")),
            html.escape(r.get("masked", "")), html.escape(r.get("untestable", "")), html.escape(r.get("pct_twin", "")), ans))
    return ("<h3 style='margin-top:16px'>Are the flanks of the copies unique?</h3><p>Independent insertions have unrelated flanks. For every family "
            "the 100 bp next to each copy, read outward from the junction, are compared between all its firmly assigned copies (flank scan stage 9, "
            "<code>flankscan/fs9_twins.sh</code>); 20-mers that occur more than 20 times in the genome are masked first, so known repeats are never "
            "evidence. Two copies are twins when a flank is colinear and at least 85 % identical (tier 1) or 70 % (tier 2) over at least 50 bp from "
            "the junction on: segmental duplications, array units the satellite screen did not catch, copies carried inside another element. "
            "A copy with too few testable bases after masking is <i>masked</i>, one with a flank cut by a contig end <i>untestable</i>. Twin copies are "
            "marked [twin] on the plates and the top 100 takes one copy per twin group first.</p>"
            "<table class='tbl'><thead><tr><th>Family</th><th class='n'>Copies</th><th class='n'>Twins (&ge; 85 %)</th><th class='n'>Twins (70-85 %)</th>"
            "<th class='n'>Masked</th><th class='n'>Untestable</th><th class='n'>Twins %</th><th>Answer</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _satellites(run_root):
    """table of results/satellites/indication.tsv (tools/satellite_stage.py, run inside step 1), '' if the run has none"""
    p = Path(run_root) / "results" / "satellites" / "indication.tsv"
    if not p.exists():
        return ""
    lines = [l.rstrip(chr(10)).split(chr(9)) for l in open(p, encoding="utf-8")]
    if len(lines) < 2:
        return ""
    h = lines[0]
    rows = []
    for l in lines[1:]:
        r = dict(zip(h, l))
        fa, fb = r.get("flag_A", "-") == "SAT_A", r.get("flag_B", "-") == "SAT_B"
        nlong = int(r.get("kindB_long_runs", "0") or 0)
        nver = int(r.get("kindB_verified_arrays", "0") or 0)       # the runs whose units are near-identical: what `verified` excludes
        # the Answer follows the unit check, not the share flag: human chr21 Alu (one copy per 4 kb, clustered in GC-rich isochores)
        # reads 20.6 % regular spacing above a uniform null with 0 of 551 runs verified - clustered copies, not an array (2026-10-05)
        if fb and nver:
            kb = "<b>tandem array of a longer unit</b> (%d of %s regularly spaced runs verified by unit identity)" % (nver, r["kindB_runs"])
        elif fb:
            kb = ("regular spacing %s %% above a uniform null but <b>no run verified by unit identity</b>: clustered copies "
                  "(the null places copies uniformly along a contig; a family that favours GC-rich regions exceeds it), not an array" % r["kindB_excess_pct"])
        elif nver:
            kb = "<b>%d arrays verified by unit identity</b> inside a dispersed family" % nver
        elif nlong:
            kb = "%d long regularly spaced runs, none verified by unit identity" % nlong
        else:
            kb = ""
        ans = " and ".join(x for x in (("<b>SINE-derived satellite</b> (%s loci)" % r["kindA_loci"]) if fa else "", kb) if x) or "&ndash;"
        rows.append("<tr><td>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td><td>%s</td></tr>" % (
            html.escape(lab(r["consensus"])), html.escape(r["full_hits"]), html.escape(r["kindA_loci"]), html.escape(r["kindA_monomers"]),
            html.escape(r["kindA_largest"]), html.escape(r["kindB_runs"]), html.escape(r.get("kindB_verified_arrays", "")),
            html.escape(r.get("kindB_long_runs", "")), html.escape(r["kindB_excess_pct"]), ans))
    nex = 0
    ex = Path(run_root) / "results" / "satellites" / "excluded_hits.bed"
    if ex.exists():
        nex = sum(1 for _ in open(ex, encoding="utf-8"))
    return ("<h3 style='margin-top:16px'>Satellites: SINE sequence in tandem arrays</h3><p>Before extraction and assignment, the hits of each consensus were "
            "screened for tandem arrays (<code>docs/SATELLITES.md</code>): loci where the monomer is a part of the SINE (TRF on the hit windows, unit aligned to the "
            "consensus, at least 4 monomers) and arrays of regular spacing whose unit is longer than the SINE; a regularly spaced run counts as an array "
            "when its units are near-identical (median unit identity &ge; 85 %%, &ldquo;verified&rdquo;), whatever the family's share. %s "
            "Loci and monomer consensuses: <code>results/satellites/</code>.</p>"
            "<table class='tbl'><thead><tr><th>Consensus</th><th class='n'>Full-length hits</th><th class='n'>Satellite loci</th><th class='n'>Monomers</th>"
            "<th class='n'>Largest locus</th><th class='n'>Regular runs</th><th class='n'>verified arrays</th><th class='n'>long runs</th><th class='n'>Excess over chance %%</th><th>Answer</th></tr></thead><tbody>"
            % ("%d hits in satellite loci were removed from the SINE analysis (kept in <code>excluded_hits.bed</code>)." % nex if nex else
               "No hits were removed.") + "".join(rows) + "</tbody></table>")


def _import_blocks():
    """tools/consensus_blocks.py: next to this file (repo, or a run dir with tools/ copied), under SINEDERELLA_TOOLS / SINEDERELLA_BIN,
    or next to the running step6 script. Returns (module, reason-string when missing)."""
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [os.path.join(here, "tools"), os.environ.get("SINEDERELLA_TOOLS", ""),
             os.path.join(os.environ.get("SINEDERELLA_BIN", ""), "tools"),
             os.path.join(os.path.dirname(os.path.abspath(sys.argv[0] or ".")), "tools")]
    for d in cands:
        if d and os.path.exists(os.path.join(d, "consensus_blocks.py")) and d not in sys.path:
            sys.path.insert(0, d)
    try:
        import consensus_blocks as cb
        return cb, ""
    except Exception as e:
        return None, "consensus_blocks not importable (%s)" % e


def section(run_root):
    # the four decision tables (length versions, arrays, satellites, consensus audit) are shown whenever their files exist, with
    # or without the block matrix: before 2026-10-05 a missing consensus_blocks.py or a one-consensus bank dropped them all
    tables = _length_versions(run_root) + _arrays(run_root) + _satellites(run_root) + _twins(run_root) + _consensus_audit(run_root)
    bank = _bank(run_root)
    cb, why = _import_blocks()
    names = cb.rd(str(bank))[0] if (bank is not None and cb is not None) else []
    if bank is None or cb is None or len(names) < 2:
        if cb is None:
            sys.stderr.write("WARNING: similarity blocks skipped (%s); the decision tables are still shown\n" % why)
        if not tables:
            return ""
        reason = ("the bank has one consensus" if (bank is not None and cb is not None) else
                  "no consensus file" if bank is None else why)
        return ("<section class=\"card\" id=\"blocks\">\n  <h2>Sequence similarity between consensuses</h2>\n"
                "  <p style='opacity:.75'>No block matrix: %s.</p>\n  %s\n</section>\n" % (html.escape(reason), tables))
    return _blocks_section(run_root, bank, cb, tables)


def _blocks_section(run_root, bank, cb, tables):
    names, seqs = cb.rd(str(bank))
    note = ""
    if len(names) > MAXN:
        keep = sorted(range(len(names)), key=lambda i: -len(seqs[i]))[:MAXN]
        keep.sort()
        tmp = Path(run_root) / "results" / "_blocks_input.fa"
        tmp.parent.mkdir(exist_ok=True)
        tmp.write_text("".join(">%s\n%s\n" % (names[i], seqs[i]) for i in keep))
        d = cb.compute(str(tmp))
        tmp.unlink()
        note = " The bank has %d consensuses; the %d longest are compared." % (len(names), MAXN)
    else:
        d = cb.compute(str(bank))
    N, LEN = d["names"], d["lens"]
    B = [b for b in d["blocks"] if b[4] - b[3] + 1 >= MINLEN and b[7] >= MINID]
    out_dir = Path(run_root) / "results"
    try:
        out_dir.mkdir(exist_ok=True)
        with open(out_dir / "consensus_blocks.tsv", "w") as fh:
            fh.write("a\tb\tstrand\ta_start\ta_end\tb_start\tb_end\tidentity\n")
            for i, j, st, a0, a1, b0, b1, idn in d["blocks"]:
                fh.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n" % (N[i], N[j], st, a0, a1, b0, b1, idn))
    except OSError:
        pass
    cov = {}
    for i, j, st, a0, a1, b0, b1, idn in B:
        cov.setdefault((i, j), set()).update(range(a0, a1 + 1))
        cov.setdefault((j, i), set()).update(range(b0, b1 + 1))
    order = list(range(len(N)))
    rows = ['<table class="bm"><tr><th></th>' + "".join('<th class="bh">%s</th>' % html.escape(lab(N[c])) for c in order) + "</tr>"]
    for a in order:
        row = '<tr><th class="n" style="text-align:right">%s <span style="opacity:.6;font-weight:400">%d</span></th>' % (html.escape(lab(N[a])), LEN[a])
        for b in order:
            if a == b:
                row += '<td style="background:#8882"></td>'
                continue
            f = len(cov.get((a, b), ())) / float(LEN[a])
            if f == 0:
                row += '<td class="bc" data-a="%d" data-b="%d"></td>' % (a, b)
                continue
            row += ('<td class="bc" data-a="%d" data-b="%d" style="background:rgba(59,125,221,%.2f);color:%s" title="%s: %d %% of its %d bp lie in blocks shared with %s">%d</td>'
                    % (a, b, 0.12 + 0.8 * f, "#fff" if f > 0.55 else "inherit", html.escape(N[a]), round(100 * f), LEN[a], html.escape(N[b]), round(100 * f)))
        rows.append(row + "</tr>")
    rows.append("</table>")
    data = {"names": [lab(n) for n in N], "lens": LEN, "blocks": B}
    first = 0
    second = 1
    if B:   # open on the pair with the largest shared block
        top = max(B, key=lambda b: b[4] - b[3])
        first, second = top[0], top[1]
    js = """<div id="bpair" style="margin-top:14px"><b>Pair view.</b> <select id="bA"></select> against <select id="bB"></select>
<svg id="bsvg" viewBox="0 0 700 150" width="100%" style="max-width:700px;display:block;margin-top:8px"></svg><div id="btab" style="font-size:12px;margin-top:6px"></div></div>
<script>
(function(){var D=__DATA__;var A=document.getElementById('bA'),B=document.getElementById('bB'),S=document.getElementById('bsvg'),T=document.getElementById('btab');
D.names.forEach(function(n,i){A.add(new Option(n+' ('+D.lens[i]+' bp)',i));B.add(new Option(n+' ('+D.lens[i]+' bp)',i));});
function g(a,b){var r=[];D.blocks.forEach(function(x){if(x[0]==a&&x[1]==b)r.push({a0:x[3],a1:x[4],b0:x[5],b1:x[6],s:x[2],id:x[7]});else if(x[0]==b&&x[1]==a)r.push({a0:x[5],a1:x[6],b0:x[3],b1:x[4],s:x[2],id:x[7]});});return r}
function draw(){var a=+A.value,b=+B.value,sc=640/Math.max(D.lens[a],D.lens[b]),x0=30,h='';
h+='<rect x="'+x0+'" y="20" width="'+D.lens[a]*sc+'" height="12" fill="#8885" rx="2"/><rect x="'+x0+'" y="108" width="'+D.lens[b]*sc+'" height="12" fill="#8885" rx="2"/>';
h+='<text x="'+x0+'" y="14" font-size="12" fill="currentColor">'+D.names[a]+'</text><text x="'+x0+'" y="138" font-size="12" fill="currentColor">'+D.names[b]+'</text>';
var bl=g(a,b),rows=[];if(a==b)bl=[];
bl.forEach(function(k,i){var c=['#3b7ddd','#d9822b','#2f9e6f','#c0554d','#8e5fc7'][i%5];
var p=[[x0+(k.a0-1)*sc,32],[x0+k.a1*sc,32],[k.s=='+'?x0+k.b1*sc:x0+(k.b0-1)*sc,108],[k.s=='+'?x0+(k.b0-1)*sc:x0+k.b1*sc,108]];
h+='<polygon points="'+p.map(function(q){return q.join(',')}).join(' ')+'" fill="'+c+'" fill-opacity="'+(0.25+0.5*(k.id-78)/22)+'" stroke="'+c+'"/>';
rows.push('<tr><td>'+k.a0+'&ndash;'+k.a1+'</td><td>'+k.s+'</td><td>'+k.b0+'&ndash;'+k.b1+'</td><td>'+(k.a1-k.a0+1)+' bp</td><td>'+k.id+' %</td></tr>')});
S.innerHTML=h;T.innerHTML=a==b?'':(bl.length?'<table><tr><th>'+D.names[a]+'</th><th>strand</th><th>'+D.names[b]+'</th><th>length</th><th>identity</th></tr>'+rows.join('')+'</table>':'No block of 30 bp or more at 80 % or better.')}
A.onchange=B.onchange=draw;document.querySelectorAll('td.bc').forEach(function(t){t.onclick=function(){A.value=t.dataset.a;B.value=t.dataset.b;draw();document.getElementById('bpair').scrollIntoView({block:'nearest'})}});A.value=__A__;B.value=__B__;draw();})();
</script>""".replace("__DATA__", json.dumps(data)).replace("__A__", str(first)).replace("__B__", str(second))
    return """<style>.bm{border-collapse:collapse;font-size:11px}.bm th{text-transform:none;letter-spacing:0}.bm th.bh{writing-mode:vertical-rl;transform:rotate(180deg);font-weight:500;padding:2px;white-space:nowrap}.bm td.bc{width:30px;height:22px;text-align:center;border:1px solid #8883;cursor:pointer}.bm td.bc:hover{outline:2px solid #d9822b}</style>
<section class="card" id="blocks">
  <h2>Sequence similarity between consensuses</h2>
  <p>The sequence-part view of the hierarchy: which stretches of one consensus are similar to which stretches of another, from the sequences alone (ungapped blocks of at least 30 bp and 80 %% identity, both strands, simple repeats masked).
  Each cell is the share of the <b>row</b> consensus (its length is given after the name) that lies in blocks shared with the <b>column</b> consensus; empty = no block. It is directional: a short consensus inside a long one is 100 %% in that row and a small share in the other. Click a cell to see the two consensuses with their blocks.%s</p>
  <div style="overflow-x:auto">%s</div>
  %s
  <p style="font-size:12px;opacity:.75">All blocks (&ge; 20 bp, &ge; 78 %%) are in <code>results/consensus_blocks.tsv</code>.</p>
  %s
</section>
""" % (note, "".join(rows), js, tables)
