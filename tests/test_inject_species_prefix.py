"""disc/inject_disc_report.py: one row per subfamily, linked to its plates, whatever the subfamily name looks like (pytest).
step8a always names plates {species}_{subfam}_{kind}.aln.fa. Brittle star seeds named abr_seedN (species abr) gave a
page with every subfamily twice (abr_seed1 and seed1), no links and no checks (2026-11-19)."""
import os
import re
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(__file__), "..", "disc", "inject_disc_report.py")
ALN = ">cons\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n>c1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n"


def _run(tmp_path, species, run_names, plates):
    aln = tmp_path / "aln"
    aln.mkdir()
    for p in plates:
        (aln / p).write_text(ALN)
    v = tmp_path / "flank_verdict.tsv"
    v.write_text("subfamily\tverdict\tassigned_copies\tsine\tgrey\tnot_sine\treason\n"
                 + "".join(f"{n}\tinterspersed\t10\t8\t1\t1\tok\n" for n in run_names))
    rep = tmp_path / "report.html"
    rep.write_text("<html><head><style></style></head><body>"
                   "<section class='card' id='alignments'></section></body></html>")
    subprocess.run([sys.executable, TOOL, str(rep), str(aln), species, "--raw-base", "https://x/aln/",
                    "--verdict-tsv", str(v)], check=True, capture_output=True, text=True,
                   env=dict(os.environ, PYTHONIOENCODING="utf-8"))
    return rep.read_text()


def _rows(page):
    return re.findall(r"<tr><td><code>([^<]+)</code>", page)


def test_names_beginning_with_species_code(tmp_path):
    page = _run(tmp_path, "abr", ["abr_seed1", "abr_seed2"],
                ["abr_abr_seed1_top100.aln.fa", "abr_abr_seed1_rand100.aln.fa", "abr_abr_seed1_subfam.aln.fa",
                 "abr_abr_seed2_top100.aln.fa"])
    assert sorted(_rows(page)) == ["abr_seed1", "abr_seed2"]
    for f in ("abr_abr_seed1_top100", "abr_abr_seed1_rand100", "abr_abr_seed1_subfam", "abr_abr_seed2_top100"):
        assert f in page, f


def test_plain_names(tmp_path):
    page = _run(tmp_path, "ccr", ["x1"], ["ccr_x1_top100.aln.fa"])
    assert _rows(page) == ["x1"]
    assert "ccr_x1_top100" in page


def test_old_single_prefix_plates_still_found(tmp_path):
    page = _run(tmp_path, "ccr", ["ccr_y"], ["ccr_y_top100.aln.fa"])
    assert _rows(page) == ["ccr_y"]
    assert "ccr_y_top100" in page
