"""edge_profile.py: copy names when the CONSENSUS record is not first, and flanks shared over the whole window (GLM audit 2026-10-19)."""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import edge_profile as E  # noqa: E402


def _aln(path, shared):
    r = random.Random(2)
    rb = lambda n: "".join(r.choice("ACGT") for _ in range(n))  # noqa: E731
    mut = lambda x, p: "".join(c if r.random() > p else r.choice("ACGT") for c in x)  # noqa: E731
    core, fl = rb(200), rb(600)
    rows = [(">c%d" % i, (mut(fl, 0.03) if shared else rb(600)) + mut(core, 0.05)) for i in range(30)]
    rows.append((">CONSENSUS", "-" * 600 + core))           # consensus LAST
    with open(path, "w") as fh:
        fh.write("".join("%s\n%s\n" % x for x in rows))


def test_unique_flanks_sharp_and_names(tmp_path):
    p = str(tmp_path / "u.fa")
    _aln(p, False)
    c = E.classify(p, "L", 600)
    assert c["pattern"] == "SHARP"
    assert c["groups"][0]["members"][0] == "c0"
    assert not any("CONSENSUS" in m for g in c["groups"] for m in g["members"])


def test_shared_flanks_not_final(tmp_path):
    p = str(tmp_path / "s.fa")
    _aln(p, True)
    c = E.classify(p, "L", 600)
    assert c["n_resolved"] == 0 and c["pattern"] == "UNRESOLVED"
    assert E.recommend([c])[0] != "FINAL"
