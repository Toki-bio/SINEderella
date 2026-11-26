"""step6 report: copies that step2 leaves unassigned are soft-assigned in step3 (soft_assigned column of summary.by_subfam.tsv);
the page must not call them unassigned. Found 2026-11-26 on the Sicista sbe58 report: headline '236,947 (10.9%) remained unassigned'
while the Soft column summed to exactly 236,947 (0 left). STEP6_REPORT_PY=<path> runs the tests against another copy."""
import importlib.util, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEP6 = Path(os.environ.get("STEP6_REPORT_PY") or ROOT / "step6_report.py")
spec = importlib.util.spec_from_file_location("step6_report_softcount", STEP6)
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
HDR = ["subfam", "firm_assigned", "soft_assigned", "total_assigned"]


def test_soft_and_left_counts():
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, HDR, [["A", "50", "20", "70"], ["B", "20", "6", "26"]])
    assert s["soft"] == "26" and s["left"] == "4"
    h = R.funnel_html(s)
    assert "soft-assigned in step3" in h and "left unassigned</td><td>4" in h and "search hit)</td><td>26" in h


def test_sbe58_numbers_leave_zero_unassigned():
    # the real run's counts: firm 1,935,999; step2 unassigned 236,947; Soft column sum 236,947
    s = {"total": "2172946", "assigned": "1935999", "unassigned": "236947"}
    R.add_soft_counts(s, HDR, [["X", "1935999", "236947", "2172946"]])
    assert s["left"] == "0"
    assert "&#9492; left unassigned</td><td>0" in R.funnel_html(s)


def test_control_without_soft_column_is_unchanged():
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, ["subfam", "firm_assigned"], [["A", "70"]])
    assert "soft" not in s and "left" not in s
    h = R.funnel_html(s)
    assert "soft-assigned" not in h and "Not assigned in step2" in h


# 2026-11-27: edge cases reported by a GLM read of fee9d5e (relayed by another session), reproduced here. When the Soft column
# exists, a bad cell, a short row or soft > unassigned must not give silently wrong numbers: it is recorded as incomplete
# (step6 then stops without a page unless --allow-incomplete, commit 250cc8a).
def _fresh():
    if hasattr(R, "INCOMPLETE"):
        R.INCOMPLETE.clear()


def test_bad_soft_cell_is_recorded():
    _fresh()
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, HDR, [["A", "50", "x", "70"], ["B", "20", "6", "26"]])
    assert "soft" not in s and "left" not in s
    assert any("soft_assigned" in m for m in R.INCOMPLETE)


def test_short_row_is_recorded_not_undercounted():
    _fresh()
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, HDR, [["A", "50"], ["B", "20", "6", "26"]])
    assert "soft" not in s
    assert any("soft_assigned" in m for m in R.INCOMPLETE)


def test_soft_above_unassigned_is_recorded():
    _fresh()
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, HDR, [["A", "50", "40", "90"]])
    assert "left" not in s
    assert any("exceeds" in m for m in R.INCOMPLETE)


def test_control_blank_row_is_ignored():
    _fresh()
    s = {"total": "100", "assigned": "70", "unassigned": "30"}
    R.add_soft_counts(s, HDR, [["A", "50", "20", "70"], [], ["B", "20", "6", "26"], [""]])
    assert s["soft"] == "26" and s["left"] == "4" and not R.INCOMPLETE
