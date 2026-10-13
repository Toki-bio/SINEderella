"""tools/rollup_families.py on a synthetic run directory (no ssearch needed)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import rollup_families as rf  # noqa: E402


def w(p, text):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write(text)
    return str(p)


def make_run(root):
    s2 = os.path.join(root, "step2", "step2_output")
    w(os.path.join(s2, "assignment_full.tsv"),
      "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n"
      "c1\tA1\t300\t10\tassigned\t100\trunner=A2;runner_bs=290;runner_ratio=0.9667\n"   # firm, runner in same family
      "c2\tB1\t250\t10\tassigned\t90\trunner=A1;runner_bs=230;runner_ratio=0.9200\n"     # firm, runner in another family
      "c3\tA2\t200\t6\tno_unanimous\t100\tNA\n"                                          # split between A1 and A2 -> family-firm
      "c4\tA1\t150\t5\tno_unanimous\t100\tNA\n"                                          # split between families -> stays
      )
    w(os.path.join(root, "genome.clean_step1", "extracted.fasta"), ">c1\nACGT\n")
    w(os.path.join(root, "genome.clean_step1", "searches", "all_hits.labeled.bed"),
      "chr1\t0\t100\tA1\t1\t+\nchr1\t90\t200\tA2\t1\t+\nchr1\t500\t600\tB1\t1\t-\n")
    w(os.path.join(root, "fam.tsv"), "subfamily\tfamily\tparent_source\nA1\tA\tderived\nA2\tA\tderived\nB1\tB1\tmember\n")
    w(os.path.join(root, "fam.fa"), ">A\nACGT\n>B1\nTTTT\n")
    w(os.path.join(root, "votes.tsv"),
      "c3\tA1\t4\t800\t0\tNA\nc3\tA2\t6\t1100\t0\tNA\n"
      "c4\tA1\t5\t700\t600\tB1\nc4\tB1\t5\t650\t600\tA1\n")


def table(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return {r[0]: r for r in rows[1:]}


def test_rollup_relabels_and_promotes(tmp_path):
    root = str(tmp_path / "run")
    make_run(root)
    out = str(tmp_path / "fv")
    rf.rollup(root, os.path.join(root, "fam.tsv"), os.path.join(root, "fam.fa"), out, os.path.join(root, "votes.tsv"))
    t = table(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"))
    assert t["c1"][1] == "A" and t["c1"][4] == "assigned" and t["c1"][6] == "runner=NA;runner_bs=0;runner_ratio=0.0000"
    assert t["c2"][1] == "B1" and "runner=A;" in t["c2"][6]            # runner across families stays, relabelled
    assert t["c3"][1] == "A" and t["c3"][3] == "10" and t["c3"][4] == "assigned" and "layer2=unresolved" in t["c3"][6]
    assert t["c3"][2] == "1900"
    assert t["c4"][4] == "no_unanimous"                                # 5/5 across two families is not firm
    bed = open(os.path.join(out, "genome.clean_step1", "searches", "regions.by_subfam.bed")).read().splitlines()
    assert bed == ["chr1\t0\t200\tA", "chr1\t500\t600\tB1"]            # A1+A2 merge into one region of one family: no CONFLICT
    hits = open(os.path.join(out, "genome.clean_step1", "searches", "all_hits.labeled.bed")).read().splitlines()
    assert [h.split("\t")[3] for h in hits] == ["A", "A", "B1"]
    s = table(os.path.join(out, "family_view.summary.tsv"))
    assert s["A"][3] == "1" and s["A"][4] == "1" and s["A"][5] == "2"


def test_without_votes_only_subfamily_firm_count(tmp_path):
    root = str(tmp_path / "run")
    make_run(root)
    out = str(tmp_path / "fv")
    rf.rollup(root, os.path.join(root, "fam.tsv"), os.path.join(root, "fam.fa"), out, None)
    t = table(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"))
    assert t["c3"][4] == "no_unanimous" and t["c3"][1] == "A"


def test_no_table_equals_flat(tmp_path):
    root = str(tmp_path / "run")
    make_run(root)
    w(os.path.join(root, "flat.tsv"), "subfamily\tfamily\tparent_source\nA1\tA1\tmember\nA2\tA2\tmember\nB1\tB1\tmember\n")
    out = str(tmp_path / "fv")
    rf.rollup(root, os.path.join(root, "flat.tsv"), os.path.join(root, "fam.fa"), out, None)
    a = open(os.path.join(out, "step2", "step2_output", "assignment_full.tsv")).read().splitlines()
    b = open(os.path.join(root, "step2", "step2_output", "assignment_full.tsv")).read().splitlines()
    assert [x.split("\t")[:6] for x in a] == [x.split("\t")[:6] for x in b]
