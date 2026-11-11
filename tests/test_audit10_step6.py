"""Regression test for audit 10, finding 1: step6 ran SINEplot.py and built `sineplot_section`, but the page
template had no place for it, so the SINEplot panel was never on the report. The section must be on the page,
after the built-in PCA card, with the SINEplot HTML intact inside the iframe srcdoc; nothing with --no-sineplot.

STEP6_REPORT_PY=<path> runs the same tests against another copy of step6_report.py (e.g. `git show main:step6_report.py`).
"""
import importlib.util
import os
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
STEP6 = Path(os.environ.get("STEP6_REPORT_PY") or ROOT / "step6_report.py")

# what SINEplot.py could write: & and " (the two characters an srcdoc="..." attribute cannot hold raw), entities, ', < >
SINEPLOT_HTML = ('<!doctype html><html><head><title>A & B</title></head><body class="x" data-q=\'s\'>'
                 '<p>a &amp; b &quot;c&quot; 1 < 2 > 0</p><script>var s = "a&b" + \'c\';</script></body></html>')


@pytest.fixture()
def R():
    spec = importlib.util.spec_from_file_location("step6_report_audit10", STEP6)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _Page(HTMLParser):
    """Order of the sections, the nav links and every iframe of a report page."""
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids, self.iframes, self.links, self.open_cards = [], [], [], []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in ("section", "details") and "card" in (a.get("class") or "").split():
            self.ids.append(a.get("id"))
        if tag == "section":
            self.open_cards.append(a.get("id"))
        if tag == "iframe":
            self.iframes.append((a, self.open_cards[-1] if self.open_cards else None))
        if tag == "a" and (a.get("href") or "").startswith("#"):
            self.links.append(a["href"])

    def handle_endtag(self, tag):
        if tag == "section" and self.open_cards:
            self.open_cards.pop()


def _build(R, tmp_path, monkeypatch, sineplot, stub):
    """build_html() on the smallest run folder it accepts, SINEplot replaced by a stub; returns (text, parsed page)."""
    run = tmp_path / "run"
    s2 = run / "step2" / "step2_output"
    s2.mkdir(parents=True, exist_ok=True)
    (run / "consensuses.clean.fa").write_text(">A\nACGTACGTACGTACGTACGT\n")
    (s2 / "assignment_full.tsv").write_text("c1\tA\t40\nc2\tA\t36\n")
    (s2 / "assigned.fasta").write_text(">c1|A\nACGTACGTACGTACGTACGT\n>c2|A\nACGTACGTACCTACGTACGT\n")
    calls = []

    def fake_iframe(*a, **k):
        calls.append(a)
        return stub

    monkeypatch.setattr(R, "build_sineplot_iframe", fake_iframe)
    monkeypatch.setattr(R.shutil, "which", lambda *a, **k: None)      # no mafft: the built-in PCA card shows its note
    out = tmp_path / "report.html"
    R.build_html(run, out, inline_plotly=False, max_table_rows=50, embed_images=False,
                 sineplot=sineplot, sineplot_max=5, threads=1)
    text = out.read_text(encoding="utf-8")
    page = _Page()
    page.feed(text)
    return text, page, calls


def test_sineplot_section_is_on_the_page(R, tmp_path, monkeypatch):
    text, page, calls = _build(R, tmp_path, monkeypatch, sineplot=True, stub=SINEPLOT_HTML)
    assert len(calls) == 1
    assert "sineplot" in page.ids, "the SINEplot result is built and then not put on the page"
    assert page.ids.index("sineplot") == page.ids.index("pca") + 1          # right after the built-in PCA card
    assert page.ids.index("sineplot") < page.ids.index("run")
    assert len(page.iframes) == 1
    attrs, card = page.iframes[0]
    assert card == "sineplot"
    assert attrs["srcdoc"] == SINEPLOT_HTML                                 # the browser gets back exactly what SINEplot wrote
    assert "embed" in attrs["class"].split()
    raw = text[text.index('srcdoc="') + 8:]
    raw = raw[:raw.index('"')]
    assert "&amp;amp; b &amp;quot;c&amp;quot;" in raw and "A &amp; B" in raw and 'class=&quot;x&quot;' in raw
    assert "#sineplot" in page.links and page.links.index("#sineplot") == page.links.index("#pca") + 1


def test_no_sineplot_leaves_nothing(R, tmp_path, monkeypatch):
    text, page, calls = _build(R, tmp_path, monkeypatch, sineplot=False, stub=SINEPLOT_HTML)
    assert calls == []
    assert page.iframes == [] and "<iframe" not in text
    assert "sineplot" not in page.ids and "#sineplot" not in page.links
    assert "pca" in page.ids


def test_sineplot_without_result_gives_a_note(R, tmp_path, monkeypatch):
    text, page, calls = _build(R, tmp_path, monkeypatch, sineplot=True, stub=None)
    assert len(calls) == 1
    assert page.ids.index("sineplot") == page.ids.index("pca") + 1
    assert page.iframes == []
    assert "SINEplot PCA skipped" in text and "#sineplot" in page.links


def test_options_are_accepted(R, tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(R, "build_html", lambda *a, **k: seen.update(k))
    assert R.main([str(tmp_path), "--no-sineplot", "--sineplot-max", "5"]) in (0, None)
    assert seen["sineplot"] is False and seen["sineplot_max"] == 5
    seen.clear()
    assert R.main([str(tmp_path)]) in (0, None)
    assert seen["sineplot"] is True
