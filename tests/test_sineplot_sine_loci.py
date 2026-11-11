"""SINEplot card and built-in PCA draw only SINE loci (2026-11-11, item 45 of docs/AUDIT_2026-11-08.md).

A run with loci_verdict.tsv: only copies with verdict SINE reach the SINEplot input (all.fa) and the PCA sample; grey and
not_SINE copies do not. A run without it: every assigned copy, as before. SINEplot gets --consensuses <consensus FASTA>,
so its subfamilies are the consensus names, not query names shorter than 15 characters. With no subfamilies/ folder the
copies come from assigned.fasta.

STEP6_REPORT_PY=<path> runs the tests against another copy of step6_report.py (e.g. `git show main:step6_report.py`).
"""
import importlib.util
import logging
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
STEP6 = Path(os.environ.get("STEP6_REPORT_PY") or ROOT / "step6_report.py")

CONS = ["MEG-RL", "Rhin_SINE_long1"]   # second name: 15 characters
COPIES = [("chr1:100-300(+)", "MEG-RL", "SINE"), ("chr1:900-1100(-)", "MEG-RL", "not_SINE"),
          ("chr2:100-300(+)", "Rhin_SINE_long1", "SINE"), ("chr2:700-900(+)", "Rhin_SINE_long1", "grey")]


@pytest.fixture()
def R():
    spec = importlib.util.spec_from_file_location("step6_report_sineloci", STEP6)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(tmp_path, verdicts=True):
    run = tmp_path / "run"
    s2 = run / "step2" / "step2_output"
    s2.mkdir(parents=True)
    (run / "consensuses.clean.fa").write_text("".join(">%s\nACGTACGTAC\n" % c for c in CONS))
    (s2 / "assigned.fasta").write_text("".join(">%s|%s|250\nACGTACGTAC\n" % (c, sf) for c, sf, _ in COPIES))
    if verdicts:
        (s2 / "loci_verdict.tsv").write_text("copy\tsubfamily\tassignment_status\tverdict\treason\n" + "".join(
            "%s\t%s\tassigned\t%s\tx\n" % (c, sf, v) for c, sf, v in COPIES))
    return run, s2


def _sineplot(R, monkeypatch, s2):
    """build_sineplot_iframe with SINEplot.py, ssearch36 and the subprocesses stubbed; returns (names in all.fa, argv)."""
    seen = {}
    monkeypatch.setattr(R, "find_sineplot", lambda d: Path("SINEplot.py"))
    monkeypatch.setattr(R.shutil, "which", lambda x: "/bin/" + x)

    def fake_run(cmd, **kw):
        if cmd[0] == "ssearch36":
            seen["names"] = [l[1:].strip() for l in open(cmd[-1]) if l.startswith(">")]
        else:
            seen["argv"] = cmd
            Path(cmd[cmd.index("-o") + 1]).write_text("<html>plot</html>")
        return subprocess.CompletedProcess(cmd, 0, stdout=b"", stderr=b"")

    monkeypatch.setattr(R.subprocess, "run", fake_run)
    kw = {"keep": R.read_sine_loci(s2)} if hasattr(R, "read_sine_loci") else {}
    out = R.build_sineplot_iframe(s2, None, 400, 1, **kw)
    assert out == "<html>plot</html>", "SINEplot was not run (no subfamilies/ folder must fall back to assigned.fasta)"
    return seen["names"], seen["argv"]


def test_sineplot_draws_only_sine_loci(R, tmp_path, monkeypatch):
    run, s2 = _run(tmp_path)
    names, argv = _sineplot(R, monkeypatch, s2)
    assert sorted(names) == sorted(CONS + ["chr1:100-300(+)", "chr2:100-300(+)"])
    assert argv[argv.index("--consensuses") + 1] == str(run / "consensuses.clean.fa")


def test_sineplot_without_verdicts_draws_all_assigned(R, tmp_path, monkeypatch):
    run, s2 = _run(tmp_path, verdicts=False)
    names, _ = _sineplot(R, monkeypatch, s2)
    assert sorted(names) == sorted(CONS + [c for c, _, _ in COPIES])


def test_pca_samples_only_sine_loci(R, tmp_path, monkeypatch, caplog):
    run, s2 = _run(tmp_path)
    monkeypatch.setattr(R.shutil, "which", lambda x: "/bin/" + x)
    monkeypatch.setattr(R.subprocess, "run", lambda cmd, **kw: (_ for _ in ()).throw(RuntimeError("stop after sampling")))
    kw = {"keep": R.read_sine_loci(s2)} if hasattr(R, "read_sine_loci") else {}
    with caplog.at_level(logging.INFO):
        try:
            R.build_pca_fig(run, s2, n_per_sf=200, threads=1, **kw)
        except RuntimeError:
            pass
    assert any("PCA: 2 copies" in r.getMessage() for r in caplog.records), [r.getMessage() for r in caplog.records]


def test_sine_loci_note(R):
    assert "no <code>loci_verdict.tsv</code>" in R.sine_loci_note(None)
    assert "(2 SINE loci)" in R.sine_loci_note({"a", "b"})
