"""Species name in the report title and heading (owner 2026-11-28): --species-name / $SPECIES_NAME, saved to
<run>/species_name.txt and reused by later rebuilds; no name = page exactly as before. STEP6_REPORT_PY=<path> tests another copy."""
import importlib.util, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEP6 = Path(os.environ.get("STEP6_REPORT_PY") or ROOT / "step6_report.py")


def load():
    spec = importlib.util.spec_from_file_location("step6_report_species", STEP6)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_no_name_no_file(tmp_path):
    assert load().resolve_species_name(tmp_path, None) is None
    assert not (tmp_path / "species_name.txt").exists()


def test_given_name_is_saved_and_reused(tmp_path):
    R = load()
    assert R.resolve_species_name(tmp_path, "  Atelerix albiventris ") == "Atelerix albiventris"
    assert (tmp_path / "species_name.txt").read_text(encoding="utf-8") == "Atelerix albiventris\n"
    assert R.resolve_species_name(tmp_path, None) == "Atelerix albiventris"
    assert R.resolve_species_name(tmp_path, "") == "Atelerix albiventris"


def test_new_name_replaces_saved(tmp_path):
    R = load()
    R.resolve_species_name(tmp_path, "Atelerix albiventris")
    assert R.resolve_species_name(tmp_path, "Erinaceus europaeus") == "Erinaceus europaeus"
    assert R.resolve_species_name(tmp_path, None) == "Erinaceus europaeus"


def test_template_lines():
    src = STEP6.read_text(encoding="utf-8")
    assert "<i>' + html.escape(species_name) + '</i> &mdash; ' if species_name else ''}SINEderella run report" in src
    assert "html.escape(species_name) + ' &mdash; ' if species_name else ''}SINEderella report &mdash;" in src
