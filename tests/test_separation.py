"""tools/separation_test.py on synthetic aligned copies (no mafft)."""
import os
import random
import sys

import pytest

pytest.importorskip("scipy")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import separation_test as st  # noqa: E402


def make(tmp_path, diag_sites, n=120, seed=3):
    rng = random.Random(seed)
    anc = [rng.choice("ACGT") for _ in range(250)]
    other = anc[:]
    for p in range(diag_sites):
        other[p * 7] = rng.choice([b for b in "ACGT" if b != anc[p * 7]])
    seqs, rows = {}, ["Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo"]
    for sub, base in (("A", anc), ("B", other)):
        for i in range(n):
            s = [c if rng.random() > 0.06 else rng.choice("ACGT") for c in base]
            name = f"{sub}{i}"
            seqs[name] = "".join(s)
            rows.append(f"{name}\t{sub}\t{1000 - i}\t10\tassigned\t100\tNA")
    fa = tmp_path / "c.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    asg = tmp_path / "a.tsv"
    asg.write_text("\n".join(rows) + "\n")
    return str(asg), str(fa)


def test_real_split_holds_and_control_does_not(tmp_path):
    asg, fa = make(tmp_path, diag_sites=6)
    n, na, nb, out = st.run(asg, fa, "A", "B", n=50, draws=4, aligned=True, null=200)
    assert out[0][2]["holds"]
    assert all(r["holds"] for _, _, r in out)
    assert all(r["control_p_purity"] > 0.01 for _, _, r in out)      # splitting one subfamily in two is not a split


def test_identical_subfamilies_do_not_hold(tmp_path):
    asg, fa = make(tmp_path, diag_sites=0)
    n, na, nb, out = st.run(asg, fa, "A", "B", n=50, draws=4, aligned=True, null=200)
    assert not any(r["holds"] for _, _, r in out)
