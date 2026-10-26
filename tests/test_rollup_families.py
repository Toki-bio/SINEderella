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
    assert s["A"][3] == "1" and s["A"][4] == "0" and s["A"][5] == "1" and s["A"][6] == "2"


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


def test_near_tie_between_subfamilies_of_one_family(tmp_path):
    root = str(tmp_path / "run")
    make_run(root)
    full = os.path.join(root, "step2", "step2_output", "assignment_full.tsv")
    with open(full, "a") as fh:
        fh.write("c5\tA1\t300\t10\tassigned\t100\trunner=A2;runner_bs=297;runner_ratio=0.9900\n")   # 1 % apart: a tie
    out = str(tmp_path / "fv")
    rf.rollup(root, os.path.join(root, "fam.tsv"), os.path.join(root, "fam.fa"), out, None)
    t = table(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"))
    assert "layer2=tie" in t["c5"][6] and "tie_with=A2" in t["c5"][6] and "runner=NA" in t["c5"][6]
    assert "layer2=tie" not in t["c1"][6]                              # 3.3 % apart: not a tie at the default 2 %
    s = table(os.path.join(out, "family_view.summary.tsv"))
    assert s["A"][4] == "1"
    out2 = str(tmp_path / "fv2")
    rf.rollup(root, os.path.join(root, "fam.tsv"), os.path.join(root, "fam.fa"), out2, None, near_tie=0.05)
    t2 = table(os.path.join(out2, "step2", "step2_output", "assignment_full.tsv"))
    assert "layer2=tie" in t2["c1"][6]                                  # at 5 % it is


def test_partial_copies(tmp_path):
    """family-firm copies matching only a stretch of their consensuses get layer2=partial:<stretch> when they pass the stretch threshold"""
    root = str(tmp_path / "run")
    make_run(root)
    with open(os.path.join(root, "step2", "step2_output", "assignment_full.tsv"), "a") as fh:
        fh.write("c9\tA2\t2000\t10\tassigned\t900\tNA\n"            # A2 threshold 900 (A1 100: the family's lowest)
                 "c6\tA2\t1900\t6\tno_unanimous\t900\tNA\n"         # split, first 80 bp of both: partial, passes 900 x 80 / 300
                 "c7\tA2\t1500\t6\tno_unanimous\t900\tNA\n"         # split, but covers most of A2: not partial
                 "c8\tA2\t150\t10\trejected_low_bitscore\t900\tNA\n")   # partial, family-firm (>= 100) but below 240
    with open(os.path.join(root, "votes.tsv"), "a") as fh:
        fh.write("c6\tA1\t4\t800\t0\tNA\t1\t80\t0\t0\nc6\tA2\t6\t1100\t0\tNA\t1\t80\t0\t0\n"
                 "c7\tA1\t4\t600\t0\tNA\t1\t150\t0\t0\nc7\tA2\t6\t900\t0\tNA\t1\t290\t0\t0\n"
                 "c8\tA2\t10\t150\t140\tA1\t1\t80\t1\t80\n")
    w(os.path.join(root, "bank.fa"), ">A1\n%s\n>A2\n%s\n>B1\n%s\n" % ("A" * 200, "C" * 300, "G" * 100))
    out = str(tmp_path / "fv")
    rf.rollup(root, os.path.join(root, "fam.tsv"), os.path.join(root, "fam.fa"), out, os.path.join(root, "votes.tsv"),
              bank=os.path.join(root, "bank.fa"), tie_rel=0)
    t = table(os.path.join(out, "step2", "step2_output", "assignment_full.tsv"))
    assert t["c6"][4] == "assigned" and t["c6"][6].endswith("layer2=partial:A_partial1")
    st = table(os.path.join(out, "partial_stretches.tsv"))
    assert st["A_partial1"][1:] == ["A", "1", "80", "A1:1-80,A2:1-80"]
    assert t["c7"][4] == "assigned" and "layer2=unresolved" in t["c7"][6] and "partial" not in t["c7"][6]
    assert t["c8"][4] == "assigned" and "layer2=unresolved" in t["c8"][6]          # runner A1 at 93 % competes; below the stretch bar
    assert "layer2=unresolved" in t["c3"][6]                                       # old vote table without spans: unchanged
    p = table(os.path.join(out, "partial_copies.tsv"))
    assert set(p) == {"c6", "c8"}
    assert p["c6"][2] == "A2" and p["c6"][4] == "0.27" and p["c6"][6] == "240" and p["c6"][7] == "partial"
    assert p["c8"][3] == "A2:1-80,A1:1-80" and p["c8"][7] == "below_stretch_threshold" and p["c8"][8] == "-"
    assert p["c6"][8] == "A_partial1"
    s = table(os.path.join(out, "family_view.summary.tsv"))
    assert s["A"][5] == "4"                                                         # c3, c6, c7, c8 family-firm only


def test_group_partials():
    """copies whose spans on a shared consensus start and end close together form one stretch; head-only vs head+module apart"""
    items = [("R", {"P2": (2, 89), "P1": (2, 89)}),        # head, best P2
             ("R", {"P2": (2, 83), "C11": (2, 83)}),       # head, linked to the first through P2
             ("R", {"P1": (2, 93), "C11": (7, 93)}),       # head, linked to the first through P1 (ends 4 bp apart)
             ("R", {"P1": (2, 161), "C11": (2, 161)}),     # head + module: ends 70 bp past the head
             ("R", {"C11": (7, 141), "P1": (7, 141)}),     # head + module, ends 20 bp from the fourth: more than 13.5 allowed
             ("M", {"P2": (2, 89)})]                       # another family: never linked
    gi, groups = rf.group_partials(items)
    assert gi[0] == gi[1] == gi[2] and gi[3] != gi[0] and gi[5] not in (gi[0], gi[3])
    assert gi[4] != gi[3]                                  # 20 bp > max(10, 10 % of 135): two stretches
    big = groups[gi[0]]
    assert big[0] == "R" and big[1] == [0, 1, 2]
    assert big[2]["P2"] == (2, 89) and big[2]["P1"] == (2, 93) and "C11" in big[2]
    assert rf.stretch_label({"P2": (2, 89), "P1": (2, 93)}) == "P1:2-93,P2:2-89"
