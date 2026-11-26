"""Regression tests for two small fixes of audit-fixes-10.

The scripts under test are taken from AUDIT10_ROOT when that is set (a copy of the tree with the files of main put back, to
show that the test fails on the old code), else from this repository. The extend_ends test needs bedtools (skipped without
it; since 2026-11-26 the tool is the shell script tools/extend_ends.sh); the other needs no external program.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("AUDIT10_ROOT", str(ROOT)))


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools")
def test_extend_ends_paths_with_space_and_shell_chars(tmp_path):
    """tools/extend_ends.sh (was extend_ends.py): every path is passed quoted, so a run directory or genome path with a
    space or a shell character does not break step 1. Real bedtools (the Python stand-in no longer applies)."""
    d = tmp_path / "run dir;$x"
    d.mkdir()
    merged = "chr1\t100\t200\t80\t100\t+\nchr1\t500\t600\t80\t100\t-\n"
    (d / "g.fa").write_text(">chr1\n" + "ACGT" * 250 + "\n")
    (d / "g.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    (d / "m.bed").write_text(merged)
    (d / "l.bed").write_text("chr1\t100\t200\tE\t100\t+\nchr1\t500\t600\tE\t100\t-\n")
    r = subprocess.run(["bash", str(SRC / "tools" / "extend_ends.sh"), str(d / "m.bed"), str(d / "l.bed"), str(d / "g.fa"),
                        str(d / "o.bed"), str(d / "x"), "--threads", "2"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (d / "o.bed").read_text() == merged                           # 2 windows per group: nothing is shared, nothing moves
    assert "EXTEND loci=2" in r.stdout
    assert (d / "x" / "ends_extension.tsv").exists() and (d / "x" / "loci_extended.tsv").exists()
    assert not (tmp_path / "x").exists() and not list(tmp_path.glob("o.bed*"))   # nothing written outside the run directory


def test_test2_props_prints_its_table(tmp_path):
    """disc/test2_props.py: main() reaches the table header ("%>10s" was not a format: ValueError)."""
    (tmp_path / "calls.tsv").write_text("set\tcorpus\tcall\ttext\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(SRC / "disc" / "test2_props.py")], cwd=str(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    head = [l for l in r.stdout.splitlines() if l.startswith("what you said")]
    assert len(head) == 1
    assert head[0] == "%-24s %3s " % ("what you said", "n") + " ".join(
        "%10s" % s for s in ("rowsplit", "rowsplit_frac", "rowsplit_dlen", "redge_shift", "redge_naive_shift",
                             "redge_tail_len", "dip_depth", "dip_width"))
