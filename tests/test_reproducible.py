"""The same input must give the same output (docs/REPRODUCIBILITY.md).

ssearch36 -z 11 (and the default statistics on a library of fewer than 500 sequences) takes its statistics from shuffled
library sequences seeded by the clock, so bit scores and E-values of the same alignment differed between two identical
calls, and with them votes, thresholds, verdicts and plates. Every call of the workflow now uses fixed statistics (-z 3).
Tested here: the control (the old call differs between runs, the new one does not, with 1 or several threads), step 2 run
twice with different numbers of threads (assignment files identical byte for byte), the continuation check run twice, and
the continuation check in a directory whose path is longer than 120 characters (ssearch36 cannot open such a path).
The whole workflow twice: tests/test_reproducible_run.sh (server, minutes).

Needs ssearch36; step 2 needs seqkit, the continuation check bedtools.
"""
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
need = [t for t in ("ssearch36", "bash") if shutil.which(t) is None]
pytestmark = pytest.mark.skipif(bool(need), reason="needs %s" % need)

sys.path.insert(0, str(Path(__file__).resolve().parent))
B4 = "ACGT"


def toy_library(d, n=80, seed=5):
    """a query and a library of n diverged copies of it (10 or more sequences, under 500: ssearch36 shuffles them)"""
    r = random.Random(seed)
    cons = "".join(r.choice(B4) for _ in range(250))
    mut = lambda s, p: "".join(r.choice(B4.replace(c, "")) if r.random() < p else c for c in s)
    (d / "q.fa").write_text(">q\n%s\n" % cons)
    with open(d / "lib.fa", "w") as fh:
        for i in range(n):
            fh.write(">c%d\n%s\n" % (i, mut(cons, 0.05 + 0.3 * i / n)))


def search(d, z, threads):
    r = subprocess.run(["ssearch36", "-g", "-3", "-T", str(threads), "-Q", "-n", "-z", z, "-E", "10", "-m", "8", "q.fa", "lib.fa"],
                       capture_output=True, text=True, cwd=d)
    assert r.returncode == 0 and r.stdout, r.stderr
    return sorted(r.stdout.splitlines())


def test_fixed_statistics_give_the_same_scores_every_time(tmp_path):
    toy_library(tmp_path)
    runs = [search(tmp_path, "3", t) for t in (1, 1, 4, 4, 2)]
    assert all(x == runs[0] for x in runs)
    assert len(runs[0]) >= 60


def test_shuffle_statistics_are_the_problem(tmp_path):
    """the control: the old call gives other bit scores from run to run, for the same alignments"""
    toy_library(tmp_path)
    runs = [search(tmp_path, "11", 1) for _ in range(4)]
    assert len({"\n".join(x) for x in runs}) > 1
    coords = lambda x: {tuple(l.split("\t")[:10]) for l in x}
    assert coords(runs[0]) & coords(runs[1])                     # the alignments are the same, the scores are not
    scores = lambda x: {tuple(l.split("\t")[:10]): l.split("\t")[11] for l in x}
    a, b = scores(runs[0]), scores(runs[1])
    assert any(a[k] != b[k] for k in a if k in b) or any(scores(runs[0]) != scores(x) for x in runs[2:])


def step2_toy(d):
    """three consensuses, two of them 2 % apart (copies of one are close to the other), and 120 copies"""
    r = random.Random(11)
    mut = lambda s, p: "".join(r.choice(B4.replace(c, "")) if r.random() < p else c for c in s)
    a = "".join(r.choice(B4) for _ in range(200))
    cons = {"TA": a, "TB": mut(a, 0.02), "TC": "".join(r.choice(B4) for _ in range(200))}
    (d / "cons.fa").write_text("".join(">%s\n%s\n" % kv for kv in cons.items()))
    with open(d / "ext.fa", "w") as fh:
        i = 0
        for k, v in cons.items():
            for p in (0.05, 0.15, 0.25, 0.4):
                for _ in range(10):
                    i += 1
                    fh.write(">ctg%d:1000-1200(+)\n%s\n" % (i, mut(v, p)))


@pytest.mark.skipif(shutil.which("seqkit") is None, reason="needs seqkit")
def test_step2_twice_identical(tmp_path):
    step2_toy(tmp_path)
    outs = []
    for name, threads in (("o1", "2"), ("o2", "5")):       # the number of threads must not change the result either
        r = subprocess.run(["bash", str(ROOT / "step2_asSINEment.sh"), "cons.fa", "ext.fa", threads, name],
                           capture_output=True, text=True, cwd=tmp_path, env=dict(os.environ, UNIT_RULE="0"))
        assert r.returncode == 0, r.stderr[-2000:]
        outs.append(tmp_path / sorted(p.name for p in tmp_path.glob(name + "*"))[0])
    for f in ("assignment_full.tsv", "all_votes.tsv", "assignment_stats.tsv", "assigned.fasta", "unassigned.fasta"):
        assert (outs[0] / f).read_bytes() == (outs[1] / f).read_bytes(), f
    rows = [l.split("\t") for l in (outs[0] / "assignment_full.tsv").read_text().splitlines()[1:]]
    assert len(rows) == 120                                # every copy has a row
    status = {x[0]: x[4] for x in rows}
    assert "assigned" in status.values()
    # the winner rule: a clear winner has 10 votes, a copy whose two best consensuses score within 5 % has 5
    for x in rows:
        assert x[3] in ("10", "5"), x
        assert (x[4] == "no_unanimous") == (x[3] == "5"), x
    assert any(x[3] == "5" for x in rows)                  # TA and TB are 2 % apart: some copies have no clear winner
    tc = [x for x in rows if x[1] == "TC"]
    assert tc and all(x[3] == "10" for x in tc)            # the unrelated consensus always wins clearly


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="needs bedtools")
def test_continuation_check_twice_identical_and_in_a_long_path(tmp_path):
    from test_continuation_check import make_run, read_tsv
    tool = str(ROOT / "tools" / "continuation_check.py")
    tables = []
    short = tmp_path / "a"
    deep = tmp_path / ("run_" + "x" * 60) / ("sub_" + "y" * 50)        # ssearch36 cannot open paths over ~120 characters
    assert len(str(deep / "results" / "continuation" / "work" / "own_lib.fa")) > 150
    for d, threads in ((short, "1"), (short, "4"), (deep, "2")):
        if not d.exists():
            d.mkdir(parents=True)
            make_run(d)                                                # seeded: the same toy run in both directories
        r = subprocess.run([sys.executable, tool, str(d), "--threads", threads, "--no-build"], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        tables.append((d / "results" / "continuation" / "continuation.tsv").read_bytes())
    assert tables[0] == tables[1] == tables[2]
    flag = [x for x in read_tsv(deep / "results" / "continuation" / "continuation.tsv") if x["flag"] == "yes"]
    assert len(flag) == 1 and flag[0]["target"] == "P"                 # the search worked in the long path
