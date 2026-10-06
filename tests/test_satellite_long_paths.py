"""ssearch36 cannot open a file whose path is longer than ~120 characters ("cannot open library"; TRF/ssearch36 calls of the satellite stage
all failed in the Sicista runs under /home/toki/sine_runs/Sicista/primary_hifiasm/, 2026-10-06). Both satellite verifiers must work in a
directory whose path is far longer than that. Needs ssearch36 on PATH."""
import os
import random
import shutil
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
sys.path.insert(0, HERE)
import satellite_kindB_verify as kb  # noqa: E402
import satellite_trf_verify as sv  # noqa: E402
from test_satellite_kindB_verify import toy  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("ssearch36") is None, reason="ssearch36 not on PATH")


def deep(tmp_path):
    d = tmp_path / ("run_" + "x" * 60) / ("genome.clean_step1_" + "y" * 40) / "satellites"
    d.mkdir(parents=True)
    assert len(str(d)) > 130
    return d


def test_kindB_in_a_long_path(tmp_path, monkeypatch):
    monkeypatch.delenv("SATELLITE_KINDB_SERIAL", raising=False)
    kb._PAIR_CACHE.clear()
    d = deep(tmp_path)
    g, runs, starts, expect = toy(d)
    assert len(str(d / "kindB_abcdefgh" / "tmpabcdefgh.l.fa")) > 125
    new = kb.verify_runs(runs, starts, g, str(d), threads=4)
    assert [new[i][3] for i in range(len(runs))] == expect
    kb._PAIR_CACHE.clear()
    ref = kb.verify_runs_serial(runs, starts, g, str(d), threads=2)
    assert [ref[i][3] for i in range(len(runs))] == expect


def test_kindA_verify_in_a_long_path(tmp_path, capsys):
    d = deep(tmp_path)
    r = random.Random(4)
    sine = "".join(r.choice("ACGT") for _ in range(250))
    cons = d / "toy.cons.fa"
    cons.write_text(">toy\n%s\n" % sine)
    unit = sine[130:250]                                    # a 120 bp monomer = the 3' part of the SINE
    rec = ("c1", 0, 1200, 120, 10.0, 95, 600, unit)          # the verifier reads the unit from field 7
    best = sv.verify([rec], str(cons), "ssearch36", str(d), 1)
    assert 0 in best and best[0][1] >= 100, best             # aligned over most of the monomer
    assert "WARNING" not in capsys.readouterr().err
