"""Regression tests for two small fixes of audit-fixes-10.

The scripts under test are taken from AUDIT10_ROOT when that is set (a copy of the tree with the files of main put back, to
show that the test fails on the old code), else from this repository. No external program is needed: bedtools is replaced by a
subprocess.run stand-in inside the test.
"""
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("AUDIT10_ROOT", str(ROOT)))


def test_extend_ends_paths_with_space_and_shell_chars(tmp_path, monkeypatch, capsys):
    """tools/extend_ends.py: bedtools gets every path as one argument and no shell is used, so a run directory or genome
    path with a space or a shell character does not break step 1."""
    spec = importlib.util.spec_from_file_location("extend_ends_audit10", str(SRC / "tools" / "extend_ends.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d = tmp_path / "run dir;$x"
    d.mkdir()
    merged = "chr1\t100\t200\t80\t100\t+\nchr1\t500\t600\t80\t100\t-\n"
    (d / "g.fa").write_text(">chr1\n" + "ACGT" * 250 + "\n")
    (d / "g.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    (d / "m.bed").write_text(merged)
    (d / "l.bed").write_text("chr1\t100\t200\tE\t100\t+\nchr1\t500\t600\tE\t100\t-\n")
    calls = []

    def fake_run(cmd, **kw):
        calls.append((cmd, kw))
        assert not kw.get("shell"), "a shell is used: %r" % (cmd,)
        assert isinstance(cmd, list), "the command is one string, not an argument list: %r" % (cmd,)
        if cmd[0] == "rm":
            shutil.rmtree(cmd[-1])
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    monkeypatch.setattr(sys, "argv", ["extend_ends.py", str(d / "m.bed"), str(d / "l.bed"), str(d / "g.fa"),
                                      str(d / "o.bed"), str(d / "x")])
    mod.main()
    bed = [c for c, _ in calls if c[0] == "bedtools"]
    assert [c[1] for c in bed] == ["intersect", "getfasta"]
    inter, getfa = bed
    assert inter[inter.index("-b") + 1] == str(d / "l.bed")              # each path is one intact argument
    assert getfa[getfa.index("-fi") + 1] == str(d / "g.fa")
    a_bed, win_bed = inter[inter.index("-a") + 1], getfa[getfa.index("-bed") + 1]
    assert a_bed.startswith(str(d / "x")) and a_bed.endswith("loci.bed")
    assert win_bed.startswith(str(d / "x")) and win_bed.endswith("win.bed")
    assert (d / "o.bed").read_text() == merged                           # no windows returned: nothing moves
    assert "EXTEND loci=2" in capsys.readouterr().out


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
