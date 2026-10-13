"""tools/families.py: the family table (parsing, checks against the bank, parent consensus supplied / derived / member).

Runs anywhere; the one test that aligns with mafft is skipped without it.
"""
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import families as fm  # noqa: E402


def write(p, text):
    p.write_text(text)
    return str(p)


BANK = ">t1_1\nACGTACGTAC\n>t1-2\nACGTACGAAC\n>t3-1\nGGGGCCCCAA\n>lone\nTTTTTTTTTT\n"


def test_table_roundtrip_and_singletons(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    tab = write(tmp_path / "f.tsv", "# comment\nsubfamily\tfamily\tparent_consensus\nt1_1\tt1\nt1-2\tt1\nt3-1\tt3\n")
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o"), aligned=True)
    d = {r[0]: r[1] for r in rows}
    assert d == {"t1_1": "t1", "t1-2": "t1", "t3-1": "t3", "lone": "lone"}  # not listed = its own family
    assert list(fams) == ["t1", "t3", "lone"]
    src = {s[0]: s[2] for s in summ}
    assert src == {"t1": "derived", "t3": "member", "lone": "member"}


def test_no_table_is_flat(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    rows, fams, summ = fm.build(bank, None, str(tmp_path / "o"), aligned=True)
    assert [r[1] for r in rows] == [r[0] for r in rows]  # every consensus its own family
    assert all(s[2] == "member" for s in summ)


def test_unknown_subfamily_is_an_error(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    tab = write(tmp_path / "f.tsv", "t1_1\tt1\nt1-3\tt1\n")
    with pytest.raises(fm.FamilyTableError, match="t1-3"):
        fm.build(bank, tab, str(tmp_path / "o"), aligned=True)


def test_duplicate_row_and_name_collision(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    with pytest.raises(fm.FamilyTableError, match="twice"):
        fm.build(bank, write(tmp_path / "a.tsv", "t1_1\tt1\nt1_1\tt2\n"), str(tmp_path / "o"), aligned=True)
    # family named like a subfamily of another family
    with pytest.raises(fm.FamilyTableError, match="also the name"):
        fm.build(bank, write(tmp_path / "b.tsv", "t3-1\tX\nt1_1\tt3-1\n"), str(tmp_path / "o"), aligned=True)
    with pytest.raises(fm.FamilyTableError, match="need at least"):
        fm.build(bank, write(tmp_path / "c.tsv", "t1_1\n"), str(tmp_path / "o"), aligned=True)


def test_supplied_parent_wins_and_is_checked(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    par = write(tmp_path / "par.fa", ">T1PARENT\nACGTACGTACGT\n")
    tab = write(tmp_path / "f.tsv", "t1_1\tt1\tT1PARENT\nt1-2\tt1\n")
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o"), parents_path=par, aligned=True)
    assert fams["t1"] == "ACGTACGTACGT" and [s for s in summ if s[0] == "t1"][0][2] == "supplied"
    with pytest.raises(fm.FamilyTableError, match="neither"):
        fm.build(bank, write(tmp_path / "g.tsv", "t1_1\tt1\tNOPE\n"), str(tmp_path / "o2"), aligned=True)
    with pytest.raises(fm.FamilyTableError, match="two different"):
        fm.build(bank, write(tmp_path / "h.tsv", "t1_1\tt1\tT1PARENT\nt1-2\tt1\tt3-1\n"), str(tmp_path / "o3"), parents_path=par, aligned=True)


def test_parent_can_be_a_bank_entry(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    tab = write(tmp_path / "f.tsv", "t1_1\tt1\tlone\nt1-2\tt1\n")  # 'lone' is a bank sequence used as the parent
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o"), aligned=True)
    assert fams["t1"] == "TTTTTTTTTT"


def test_weighted_consensus_follows_copy_numbers():
    al = ["ACGT", "AGGT", "AGGT"]
    assert fm.weighted_consensus(al, [1, 1, 1]) == "AGGT"           # two against one
    assert fm.weighted_consensus(al, [10, 1, 1]) == "ACGT"          # copy number tips the column
    assert fm.weighted_consensus(["ACGT", "AC-T", "AC-T"], [1, 1, 1]) == "ACT"      # mostly gap: column dropped
    assert fm.weighted_consensus(["ACGT", "AC-T", "AC-T"], [5, 1, 1]) == "ACGT"     # weight keeps it (5/7 non-gap)
    with pytest.raises(fm.FamilyTableError):
        fm.weighted_consensus(["ACG", "AC"], [1, 1])


def test_derived_parent_uses_weights_file(tmp_path):
    bank = write(tmp_path / "bank.fa", ">a\nACGT\n>b\nAGGT\n>c\nAGGT\n")
    tab = write(tmp_path / "f.tsv", "a\tF\nb\tF\nc\tF\n")
    w = write(tmp_path / "w.tsv", "subfamily\tn\na\t1000\nb\t10\nc\t10\n")
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o"), weights_path=w, aligned=True)
    assert fams["F"] == "ACGT"
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o2"), aligned=True)
    assert fams["F"] == "AGGT"


def test_outputs_written(tmp_path):
    bank = write(tmp_path / "bank.fa", BANK)
    tab = write(tmp_path / "f.tsv", "t1_1\tt1\nt1-2\tt1\n")
    fm.build(bank, tab, str(tmp_path / "o"), aligned=True)
    for f in ("families.clean.tsv", "family_consensuses.fa", "families.summary.tsv"):
        assert (tmp_path / "o" / f).stat().st_size > 0
    assert fm.main(["--bank", bank, "--table", tab, "--out", str(tmp_path / "o2"), "--aligned"]) == 0
    assert fm.main(["--bank", bank, "--table", write(tmp_path / "bad.tsv", "zzz\tF\n"), "--out", str(tmp_path / "o3")]) == 2


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not on PATH")
def test_derived_parent_with_mafft(tmp_path):
    bank = write(tmp_path / "bank.fa", ">a\nACGTACGTACGTAAGGCC\n>b\nACGTACGAACGTAAGGCC\n>c\nACGTACGAACGTAAGGCC\n")
    tab = write(tmp_path / "f.tsv", "a\tF\nb\tF\nc\tF\n")
    rows, fams, summ = fm.build(bank, tab, str(tmp_path / "o"))
    assert fams["F"] == "ACGTACGAACGTAAGGCC"
