"""Regression tests for the GLM audit fixes of 2026-10-19 (docs/AUDIT_FIXES_6.md)."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import consensus_bank_lib as cbl  # noqa: E402


class _Boom(str):
    """A sequence that fails part-way through being written."""
    def __getitem__(self, k):
        raise RuntimeError("disk full")


def test_write_fa_failure_keeps_old_bank(tmp_path):
    # finding 99: orient_consensus_bank.py without -o rewrote the bank in place; a failure part-way truncated it
    bank = tmp_path / "bank.fa"
    bank.write_text(">A\nACGT\n>B\nGGCC\n")
    with pytest.raises(RuntimeError):
        cbl.write_fa(bank, {"A": "ACGT", "B": _Boom("GGCC")})
    assert bank.read_text() == ">A\nACGT\n>B\nGGCC\n"


def _border_loop(tmp_path, monkeypatch, rebuilt):
    """border_loop_subfam.main() with every external tool replaced: both sides extend by 30 bp and come out independent."""
    sys.path.insert(0, str(ROOT / "publish"))
    import border_loop_subfam as B
    run = tmp_path / "run"
    (run / "step2/step2_output").mkdir(parents=True)
    (run / "consensuses.clean.fa").write_text(">SF\nACGTACGTAC\n")

    def side(name, bed, *a, **k):
        a[-1]["verdict_%s" % name] = "independent"
        a[-1]["stop_reason_%s" % name] = "independent (island_fraction=0.0)"
        return bed, 30
    monkeypatch.setattr(B.F, "load_top100", lambda *a: ["c1:100-200(+)|SF"])
    monkeypatch.setattr(B.F, "headers_to_bed", lambda h, w: str(tmp_path / "x.bed"))
    (tmp_path / "x.bed").write_text("c1\t70\t230\tc1:100-200(+)|SF\t0\t+\n")
    monkeypatch.setattr(B, "run_side", side)
    monkeypatch.setattr(B, "rebuild_consensus", lambda *a: (rebuilt, None if rebuilt else "sine_consensus failed"))
    monkeypatch.setattr(sys, "argv", ["border_loop_subfam.py", str(run), "SF"])
    B.main()
    import json
    return run, json.load(open(run / "border_loop/SF/border_loop.json"))


def test_border_loop_failed_rebuild_is_untested(tmp_path, monkeypatch):
    # finding 120: copies widened, consensus rebuild failed -> widened copies were published with the old consensus and
    # the sides stayed "independent"
    run, log = _border_loop(tmp_path, monkeypatch, None)
    assert not (run / "border_loop/SF/loci.adjusted.bed").exists()
    assert log["verdict_5prime"] == "untested" and log["verdict_3prime"] == "untested"


def test_border_loop_rebuilt_keeps_adjusted_bed(tmp_path, monkeypatch):
    run, log = _border_loop(tmp_path, monkeypatch, "ACGTACGTACGT")
    assert (run / "border_loop/SF/loci.adjusted.bed").exists()
    assert log["verdict_5prime"] == "independent"


def test_merge_rebuilt_without_record_stops(tmp_path):
    # finding 120: a rebuilt file without its record was skipped silently
    import subprocess
    (tmp_path / "rebuilt_consensus").mkdir()
    (tmp_path / "consensuses.clean.fa").write_text(">SF\nACGT\n")
    (tmp_path / "rebuilt_consensus/SF.fa").write_text(">other\nACGTAA\n")
    r = subprocess.run([sys.executable, str(ROOT / "publish/merge_rebuilt_consensus.py"), str(tmp_path)], capture_output=True, text=True)
    assert r.returncode != 0


def test_mask_from_run_empty_table_stops(tmp_path):
    # finding 133: an empty all_hits.labeled.bed gave a 0-byte mask and exit 0, so the "masked" run searched unmasked
    import subprocess
    d = tmp_path / "genome.clean_step1" / "searches"
    d.mkdir(parents=True)
    (d / "all_hits.labeled.bed").write_text("")
    r = subprocess.run([sys.executable, str(ROOT / "tools/mask_from_run.py"), str(tmp_path), "-o", str(tmp_path / "m.bed")],
                       capture_output=True, text=True)
    assert r.returncode != 0


def test_length_variants_failed_search_raises(tmp_path, monkeypatch):
    # finding 96: ssearch36 exits 0 with "***" on stderr when it cannot read a file; that was read as zero 5'-complete copies
    sys.path.insert(0, str(ROOT / "tools"))
    import subprocess
    import length_variants as lv
    (tmp_path / "g.fa").write_text(">c1\n" + "A" * 1000 + "\n")
    (tmp_path / "g.fa.fai").write_text("c1\t1000\t4\t1000\t1001\n")
    (tmp_path / "bank.fa").write_text(">L\nACGTACGTAC\n")
    (tmp_path / "hits.tsv").write_text("c1:100-200(+)\tL\t.\t.\tassigned\n")
    w = tmp_path / "w"

    def fake_check_call(cmd, *a, **k):
        (w / "copies.fa").write_text(">c0::c1:0-300(+)\n" + "A" * 300 + "\n")
    monkeypatch.setattr(lv.subprocess, "check_call", fake_check_call)
    monkeypatch.setattr(lv.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, "", " *** cannot open library\n"))
    with pytest.raises(RuntimeError):
        lv.collect(str(tmp_path / "g.fa"), str(tmp_path / "hits.tsv"), {"L"}, str(tmp_path / "bank.fa"), "L", str(w), 1)


def test_add_seed_row_failure_stops_publish():
    # finding 27: an add_seed_row failure was a WARNING and the plates went out without the consensus as searched
    for p, pat in (("publish_run.sh", 'die "add_seed_row.py failed'), ("publish/align_for_publish.sh", "add_seed_row failed: the plates")):
        txt = (ROOT / p).read_text(encoding="utf-8")
        line = next(l for l in txt.splitlines() if "tools/add_seed_row.py" in l and "python3" in l)
        assert pat in line and 'log "WARNING' not in line, p


def test_add_seed_row_no_plates_ok(tmp_path):
    import subprocess
    (tmp_path / "bank.fa").write_text(">A\nACGT\n")
    r = subprocess.run([sys.executable, str(ROOT / "tools/add_seed_row.py"), str(tmp_path / "bank.fa"), "X"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_run_sine_failed_rerun_not_stale(tmp_path):
    # findings 42 / 151: the second rebuild of a subfamily reuses work/sine; a failed run returned the first run's consensus
    sys.path.insert(0, str(ROOT / "publish"))
    import flank_border_iterate as F
    work = tmp_path / "w"
    (work / "sine").mkdir(parents=True)
    (work / "sine/final_consensus.fasta").write_text(">consensus_copies\nACGTACGTAAAA\n")
    (tmp_path / "c.fa").write_text(">a\nACGT\n")
    (tmp_path / "fail.sh").write_text("exit 1\n")
    seq, err = F.run_sine(str(tmp_path / "c.fa"), str(work), str(tmp_path / "fail.sh"))
    assert seq is None


def test_run_sine_not_converged_is_a_failed_rebuild(tmp_path):
    # finding 65 (owner decision 2026-10-19): run_sine asks sine_consensus.sh for strict mode; a script that writes a consensus
    # but exits 2 only when SINE_CONSENSUS_FAIL_NOTCONVERGED=1 (= did not converge) must give no consensus
    sys.path.insert(0, str(ROOT / "publish"))
    import flank_border_iterate as F
    (tmp_path / "c.fa").write_text(">a\nACGT\n")
    (tmp_path / "nc.sh").write_text('printf ">x\\nACGTACGT\\n" > final_consensus.fasta\n'
                                    '[ "${SINE_CONSENSUS_FAIL_NOTCONVERGED:-0}" = 1 ] && exit 2\nexit 0\n')
    seq, err = F.run_sine(str(tmp_path / "c.fa"), str(tmp_path / "w"), str(tmp_path / "nc.sh"))
    assert seq is None


def test_step4_diagnostic_picks_highest_run(tmp_path):
    # finding 75: a string sort picked step2_output_run9 over step2_output_run10
    import re
    src = (ROOT / "step4_diagnostic.py").read_text(encoding="utf-8")
    m = re.search(r"step2_dirs = sorted\(\[.*?\], (key=.*?)reverse=True\)", src, re.S)
    assert m, "sort call not found"
    key = eval("dict(%s)" % m.group(1).rstrip(", "), {"re": re})["key"]
    dirs = ["step2_output", "step2_output_run9", "step2_output_run10", "step2_output_run2"]
    assert sorted(dirs, key=key, reverse=True)[0] == "step2_output_run10"


def _plate(tmp_path):
    p = tmp_path / "SF_top100.aln.fa"
    body = "aaaaACGTACGTACGTcccc"
    p.write_text(">SF_extended\n%s\n>SF\n%s\n>c1\n%s\n>c2\n%s\n" % (body, body, body, body))
    return p


def test_cpg_div_skips_seed_row_and_refuses_no_element(tmp_path):
    # findings 33 / 94: row 2 (the consensus as searched) counted as a copy; a row 1 without element columns gave all-NA
    sys.path.insert(0, str(ROOT / "tools" / "cpg_div"))
    import cpg_div
    recs = list(cpg_div.parse_fasta(str(_plate(tmp_path))))
    assert len(cpg_div.compute_plate(recs)[1]) == 2
    with pytest.raises(ValueError):
        cpg_div.compute_plate([("SF", "acgtacgt"), ("c1", "acgtacgt")])


def test_tsd_curve_skips_seed_row(tmp_path):
    # finding 108: the seed-row skip was multiplied by 0
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "tools/tsd_curve.py"), str(_plate(tmp_path)), "5", "16"], capture_output=True, text=True)
    assert "copies 2 " in r.stdout, r.stdout + r.stderr


def test_sim_family_diag_underdelivered_stops(tmp_path):
    # finding 38: the shipped diag5 scenario placed 3 of 5 diagnostic positions in subfamily 4 without a word
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "tools/sim/sim_family.py"), str(tmp_path / "s"), "--seed", "1", "--nsub", "8",
                        "--sizes", "200", "--diag", "5", "--div-old", "0.12", "--div-young", "0.03"], capture_output=True, text=True)
    assert r.returncode != 0 and "diagnostic positions" in r.stderr


def test_prep_labeled_skips_empty_rows(tmp_path):
    # finding 129: an all-gap row was written and counted as a copy
    import subprocess
    d = tmp_path / "aln"
    d.mkdir()
    (d / "POS__SP__G1_3seqs.aln.fa").write_text(">CONSENSUS_G1\nACGT\n>c1\nAC-T\n>c2\n----\n>c3\nACGT\n")
    subprocess.run([sys.executable, str(ROOT / "tools/coseg_compare/prep_labeled.py"), str(d), "SP", str(tmp_path / "o")],
                   capture_output=True, text=True)
    labels = [l for l in open(next(tmp_path.glob("o*labels"))) if l.strip()]
    assert len(labels) == 2


def test_assembly_qc_paragraph(tmp_path):
    # finding 39 (owner decision 2026-10-19): a run-level line, never a verdict change; a skipped check does not read as SOLID
    sys.path.insert(0, str(ROOT))
    import step6_report as R
    q = tmp_path / "genome.clean_assembly_qc.tsv"
    q.write_text("metric\tvalue\tflag\nN50\t8000\tlow\ngap_pct\t0.1\tok\nverdict\tWARN\tWARN\n")
    p = R.assembly_qc_paragraph(tmp_path, {})
    assert "WARN" in p and "N50 8000 (low)" in p and "gap_pct" not in p
    q.unlink()
    assert "skipped" in R.assembly_qc_paragraph(tmp_path, {"ASSEMBLY_QC_VERDICT": "SKIPPED"})
    assert "not SOLID" in R.assembly_qc_paragraph(tmp_path, {})
