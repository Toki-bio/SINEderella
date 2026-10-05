"""tools/satellite_kindB_verify.py: the fast unit check (verify_runs) against the reference (verify_runs_serial) on a toy genome
with planted tandem arrays (single unit, dimeric, long period, near the 85 % threshold) and a run of ordinary copies.

The ssearch36 tests need ssearch36 on PATH (therioserver: conda env sinederella); the rest runs anywhere.
"""
import os
import random
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import satellite_kindB_verify as kb  # noqa: E402

HAVE_SSEARCH = shutil.which("ssearch36") is not None
need_ssearch = pytest.mark.skipif(not HAVE_SSEARCH, reason="ssearch36 not on PATH")


def rnd(r, n):
    return "".join(r.choice("ACGT") for _ in range(n))


def mutate(r, s, rate):
    return "".join(c if r.random() >= rate else r.choice([x for x in "ACGT" if x != c]) for c in s)


def toy(tmp_path):
    """one contig of planted blocks; returns genome path, runs, starts and the expected verdict per run"""
    r = random.Random(7)
    sine = rnd(r, 300)
    seq, hits, runs, expect = [rnd(r, 5000)], [], [], []
    pos = [5000]

    def add(piece, is_hit=False):
        if is_hit:
            hits.append((pos[0], pos[0] + len(sine)))
        seq.append(piece)
        pos[0] += len(piece)

    def block(units, verdict):
        first = len(hits)
        for u in units:
            add(mutate(r, sine, 0.03), True)
            add(u)
        add(mutate(r, sine, 0.03), True)            # the last hit closes the last unit
        s, e = hits[first][0], hits[-1][0]
        runs.append(("q", "chrA", s, e, int((e - s) / len(units)), len(units) + 1))
        expect.append(verdict)
        add(rnd(r, 20000))

    body = rnd(r, 2100)
    block([mutate(r, body, 0.01) for _ in range(8)], "ARRAY")
    x, y1, y2 = rnd(r, 800), rnd(r, 1055), rnd(r, 360)        # dimeric, like rsi MEG-RS (2 155 / 1 460 bp sharing ~800 bp)
    block([mutate(r, x + (y1 if i % 2 == 0 else y2), 0.01) for i in range(8)], "ARRAY")
    long_body = rnd(r, 6700)
    block([mutate(r, long_body, 0.02) for _ in range(10)], "ARRAY")
    near = rnd(r, 700)
    block([mutate(r, near, 0.03) for _ in range(8)], "ARRAY")   # ~94 % between units
    far = rnd(r, 700)
    block([mutate(r, far, 0.15) for _ in range(8)], "COPIES")  # ~72 % between units
    block([rnd(r, 2700) for _ in range(8)], "COPIES")                  # SINE + unrelated flank
    add(rnd(r, 5000))
    g = tmp_path / "toy.fa"
    g.write_text(">chrA\n" + "".join(seq) + "\n")
    starts = {"q": {"chrA": sorted(hits)}}
    return str(g), runs, starts, expect


def test_hits_in_is_the_linear_scan():
    r = random.Random(1)
    lst = sorted((r.randint(0, 10000), r.randint(0, 10000)) for _ in range(500))
    for _ in range(300):
        s = r.randint(-10, 10010)
        e = s + r.randint(0, 3000)
        assert kb._hits_in(lst, s, e) == [x for x in lst if s <= x[0] <= e]


@need_ssearch
def test_fast_equals_serial_on_planted_arrays(tmp_path, monkeypatch):
    monkeypatch.delenv("SATELLITE_KINDB_SERIAL", raising=False)
    kb._PAIR_CACHE.clear()
    g, runs, starts, expect = toy(tmp_path)
    ref = kb.verify_runs_serial(runs, starts, g, str(tmp_path), threads=2)
    new = kb.verify_runs(runs, starts, g, str(tmp_path), threads=4)
    for i, want in enumerate(expect):
        assert ref[i][3] == want, ("serial", i, ref[i])
        assert new[i][3] == want, ("fast", i, new[i])
        assert new[i][0] == ref[i][0]                       # units tested
        if want == "ARRAY":                                 # the high-identity pairs: the same best alignment
            assert abs(new[i][1] - ref[i][1]) < 0.5, (i, ref[i], new[i])


@need_ssearch
def test_each_pair_is_aligned_once(tmp_path, monkeypatch):
    """a narrow run inside a wide run, a second call for the same runs: no pair is aligned twice"""
    monkeypatch.delenv("SATELLITE_KINDB_SERIAL", raising=False)
    kb._PAIR_CACHE.clear()
    g, runs, starts, expect = toy(tmp_path)
    calls = []
    real = kb._align_query

    def counting(qseq, partners, tmpdir, ssearch="ssearch36"):
        calls.append(len(partners))
        return real(qseq, partners, tmpdir, ssearch)

    monkeypatch.setattr(kb, "_align_query", counting)
    q, c, s, e, unit, nh = runs[0]
    hs = kb._hits_in(starts["q"]["chrA"], s, e)
    narrow = (q, c, hs[2][0], hs[6][0], unit, 5)             # units 2-5 of the first array
    first = kb.verify_runs([runs[0], narrow], starts, g, str(tmp_path), threads=2)
    n_pairs = sum(calls)
    assert n_pairs == 7 + 6                                  # 8 units: 7 lag-1 + 6 lag-2 pairs, the narrow run's all among them
    again = kb.verify_runs([runs[0], narrow], starts, g, str(tmp_path), threads=2)
    assert sum(calls) == n_pairs                             # nothing new aligned
    assert first == again
    assert first[1][3] == "ARRAY"


def test_serial_switch(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(kb, "verify_runs_serial", lambda *a, **k: seen.append(1) or {})
    monkeypatch.setenv("SATELLITE_KINDB_SERIAL", "1")
    assert kb.verify_runs([], {}, "none.fa", str(tmp_path)) == {}
    assert seen == [1]
