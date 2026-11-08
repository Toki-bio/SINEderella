"""Regression tests for the audit fixes of 2026-11-08 (docs/AUDIT_FIXES_7.md, findings in docs/AUDIT_2026-11-08.md).

Each test names its finding id. The scripts under test are taken from AUDIT7_ROOT when that is set (a copy of the tree with
the files of main put back, to show that the test fails on the old code), else from this repository.
External programs (minimap2, seqkit, ssearch36, ...) are replaced by small stubs, so the tests run without them.
"""
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("AUDIT7_ROOT", str(ROOT)))
BASH = shutil.which("bash")
need_bash = pytest.mark.skipif(BASH is None, reason="needs bash")


def _stub(bindir, name, body):
    """A small bash program named `name` in bindir."""
    bindir.mkdir(parents=True, exist_ok=True)
    p = bindir / name
    p.write_text("#!/usr/bin/env bash\n" + body + "\n", encoding="utf-8", newline="\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return p


def _bash(script, cwd, bindir=None, env_extra=None, args=()):
    env = dict(os.environ)
    if bindir is not None:
        env["PATH"] = str(bindir) + os.pathsep + env["PATH"]
    env.update(env_extra or {})
    return subprocess.run([BASH, Path(script).as_posix(), *args], cwd=str(cwd), env=env,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------------------------------------------------
# 41: sear, minimap2 path (queries >= 1000 bp): a crashed minimap2 read as "no raw alignments", empty BED, exit 0
# ---------------------------------------------------------------------------------------------------------------------
def _sear_minimap2(tmp_path, minimap2_body):
    b = tmp_path / "bin"
    for t in ("bedtools", "samtools"):
        _stub(b, t, "exit 0")
    _stub(b, "seqkit", "cat >/dev/null 2>&1; exit 0")
    _stub(b, "minimap2", minimap2_body)
    w = tmp_path / "w"
    w.mkdir()
    (w / "q.fa").write_text(">q\n" + "ACGT" * 300 + "\n", encoding="utf-8", newline="\n")   # 1200 bp: minimap2 mode
    (w / "gen.fa").write_text(">c1\n" + "A" * 5000 + "\n", encoding="utf-8", newline="\n")
    (w / "gen.fa.fai").write_text("c1\t5000\t4\t5000\t5001\n", encoding="utf-8", newline="\n")
    return _bash(SRC / "sear", w, b, args=("q.fa", "gen.fa", "0.8", "65", "50")), w


@need_bash
def test_41_sear_minimap2_crash_stops(tmp_path):
    r, w = _sear_minimap2(tmp_path, "echo 'minimap2: index failed' >&2; exit 1")
    assert r.returncode != 0, "a crashed minimap2 must not read as zero copies:\n" + r.stdout[-500:]


@need_bash
def test_41_sear_minimap2_no_hits_still_ok(tmp_path):
    r, w = _sear_minimap2(tmp_path, "cat >/dev/null; exit 0")
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    assert (w / "gen-q.bed").exists() and (w / "gen-q.bed").stat().st_size == 0


# ---------------------------------------------------------------------------------------------------------------------
# 72: step2 --add, tier-1 id -> BED: the contig was cut at the FIRST ':' ("ctg1:apple:100-200(+)" -> ctg1 0 200)
# ---------------------------------------------------------------------------------------------------------------------
def _awk_block(script, end_marker):
    """The awk program of step2's tier-1 id -> BED conversion, cut from the script as it is."""
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    end = next(i for i, l in enumerate(lines) if end_marker in l)
    start = max(i for i in range(end) if lines[i].strip() == "awk '{")
    prog = "\n".join(lines[start:end + 1])
    prog = prog[prog.index("'") + 1:]
    return prog[:prog.index("}' ") + 1]


@pytest.mark.skipif(shutil.which("gawk") is None and shutil.which("awk") is None, reason="needs awk")
def test_72_step2_tier1_id_with_colon_in_contig(tmp_path):
    prog = _awk_block(SRC / "step2_asSINEment.sh", '"$TMPDIR/tier1_present.ids" | sort -k1,1 -k2,2n')
    ids = tmp_path / "ids"
    ids.write_text("ctg1:apple:100-200(+)\nchr2:5-50(-)\n", encoding="utf-8", newline="\n")
    r = subprocess.run([shutil.which("gawk") or shutil.which("awk"), prog, str(ids)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in r.stdout.splitlines()]
    assert rows[0][:3] == ["ctg1:apple", "99", "200"], rows
    assert rows[1][:3] == ["chr2", "4", "50"], rows


# ---------------------------------------------------------------------------------------------------------------------
# 20: SINEderella canonicalize_consensus_bank(): a crash of canonicalize_consensus_bank.py was a WARNING, the run went on
# ---------------------------------------------------------------------------------------------------------------------
def _bash_function(script, name):
    """The text of one bash function `name(){ ... }` of a script, from its first line to the first line that is '}'."""
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(name + "(){") or l.startswith(name + "() {"))
    end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
    return "\n".join(lines[start:end + 1]) + "\n"


@need_bash
def test_20_canonicalizer_crash_stops_run(tmp_path):
    b = tmp_path / "bin"
    _stub(b, "python3", "echo 'Traceback: boom' >&2; exit 1")
    sd = tmp_path / "sd"
    sd.mkdir()
    (sd / "canonicalize_consensus_bank.py").write_text("", encoding="utf-8")
    (tmp_path / "bank.fa").write_text(">a\nACGT\n", encoding="utf-8", newline="\n")
    drv = tmp_path / "drv.sh"
    drv.write_text("set -euo pipefail\nlog(){ echo \"$*\" >&2; }\ndie(){ echo \"ERROR: $*\" >&2; exit 1; }\n"
                   "SCRIPT_DIR=%s\n" % sd.as_posix()
                   + _bash_function(SRC / "SINEderella", "canonicalize_consensus_bank")
                   + "canonicalize_consensus_bank %s\necho CONTINUED\n" % (tmp_path / "bank.fa").as_posix(),
                   encoding="utf-8", newline="\n")
    r = _bash(drv, tmp_path, b)
    assert "CONTINUED" not in r.stdout and r.returncode != 0, r.stderr
