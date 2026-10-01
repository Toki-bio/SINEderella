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


def section(run_root):
    bank = _bank(run_root)
    if bank is None:
        return ""
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "tools"))
    try:
        import consensus_blocks as cb
    except Exception as e:
        sys.stderr.write("WARNING: similarity blocks skipped (%s)\n" % e)
        return ""
    names, seqs = cb.rd(str(bank))
    if len(names) < 2:
        return ""
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
""" % (note, "".join(rows), js, _length_versions(run_root))
