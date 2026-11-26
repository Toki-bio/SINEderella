"""step6 must not write a page that silently lacks parts (2026-11-26: a page built with a Python lacking numpy lost its verdict
columns, diagrams and similarity blocks with only log warnings). preflight() lists missing modules/programs; apply_incomplete()
stops before writing (exit 1) or, with --allow-incomplete, marks the page. STEP6_REPORT_PY=<path> tests another copy."""
import builtins, importlib.util, os, shutil, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
STEP6 = Path(os.environ.get("STEP6_REPORT_PY") or ROOT / "step6_report.py")


def load():
    spec = importlib.util.spec_from_file_location("step6_report_loud", STEP6)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.INCOMPLETE.clear()
    return m


def test_missing_module_and_programs_are_reported(monkeypatch):
    R = load()
    monkeypatch.setitem(sys.modules, "report_trna", None)               # import of a None entry raises ImportError
    monkeypatch.setattr(shutil, "which", lambda name: None)
    got = "\n".join(R.preflight(sineplot=True, script_dir=None))
    assert "report_trna" in got and "mafft not in PATH" in got and "ssearch36 not in PATH" in got


def test_consensus_blocks_not_found_is_reported(monkeypatch):
    R = load()
    import report_blocks
    monkeypatch.setattr(report_blocks, "_import_blocks", lambda: (None, "consensus_blocks not importable (simulated)"))
    monkeypatch.setattr(shutil, "which", lambda name: "/bin/" + name)
    got = R.preflight(sineplot=False)
    assert any("consensus_blocks" in g for g in got)


def test_incomplete_build_stops_without_page():
    R = load()
    R.incomplete("similarity blocks: simulated")
    with pytest.raises(SystemExit) as e:
        R.apply_incomplete("<html><body><p>x</p></body></html>", allow=False)
    assert e.value.code == 1


def test_allow_incomplete_marks_the_page():
    R = load()
    R.incomplete("similarity blocks: simulated")
    out = R.apply_incomplete("<html><body class='a'><p>x</p></body></html>", allow=True)
    assert out.index("Incomplete report") > out.index("<body class='a'>") and "similarity blocks: simulated" in out


def test_control_complete_build_is_unchanged():
    R = load()
    doc = "<html><body><p>x</p></body></html>"
    assert R.apply_incomplete(doc, allow=False) == doc
