"""tools/flank_verdict.py on planted cases (pytest). The expansion's end state must decide the SINE loci (docs/AIM.md)."""
import json
import os
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(__file__), "..", "tools", "flank_verdict.py")


def _run(tmp_path, step7_rows, loops, sfs=()):
    s2 = tmp_path / "step2" / "step2_output"
    s2.mkdir(parents=True)
    rows = ["Sequence\tSubfamily\tBitscore\tVotes\tStatus"]
    for sf in sorted({r[0] for r in (step7_rows or [])} | set(loops) | set(sfs)):
        rows += ["c%s%d\t%s\t100\t10\tassigned" % (sf, i, sf) for i in (1, 2)]
    rows.append("cX\tA\t10\t3\tno_unanimous")
    (s2 / "assignment_full.tsv").write_text("\n".join(rows) + "\n")
    if step7_rows is not None:
        (s2 / "boundary_refinement.tsv").write_text(
            "subfamily\tpopulation\tside\tboundary_bp\tstatus\n"
            + "".join("%s\tgeneral\t%s\t0\t%s\n" % r for r in step7_rows))
    for sf, d in loops.items():
        p = tmp_path / "border_loop" / sf
        p.mkdir(parents=True)
        (p / "border_loop.json").write_text(json.dumps(dict(d, subfam=sf)))
    out = subprocess.run([sys.executable, TOOL, str(tmp_path)], capture_output=True, text=True, check=True).stdout
    fv = {l.split("\t")[0]: l.split("\t")[1] for l in (s2 / "flank_verdict.tsv").read_text().splitlines()[1:]}
    lv = {l.split("\t")[0]: l.split("\t")[3] for l in (s2 / "loci_verdict.tsv").read_text().splitlines()[1:]}
    return out, fv, lv


def test_end_state_decides(tmp_path):
    step7 = [("A", "upstream", "confirmed"), ("A", "downstream", "confirmed"),
             ("B", "upstream", "undetermined"), ("B", "downstream", "undetermined"),
             ("C", "upstream", "undetermined"), ("C", "downstream", "confirmed"),
             ("D", "upstream", "undetermined"), ("D", "downstream", "undetermined"),
             ("E", "upstream", "undetermined"), ("E", "downstream", "insufficient_data")]
    loops = {
        # rsi MEG-RS before the fix: 3' had nothing left to scan and was called independent
        "B": {"stop_reason_5prime": "independent (island_fraction=0.0888)",
              "stop_reason_3prime": "independent (island_fraction=None)"},
        "C": {"verdict_5prime": "independent", "verdict_3prime": "independent"},
        "D": {"verdict_5prime": "still_shared", "verdict_3prime": "independent"},
    }
    out, fv, lv = _run(tmp_path, step7, loops)
    assert fv == {"A": "interspersed", "B": "untested", "C": "interspersed", "D": "not_interspersed",
                  "E": "not_interspersed"}
    assert lv["cA1"] == "SINE" and lv["cC2"] == "SINE"
    assert lv["cB1"] == lv["cD1"] == lv["cE1"] == "not_SINE"
    assert lv["cX"] == "not_SINE"                     # not assigned
    assert len(lv) == 11                              # every copy gets a verdict
    assert "copies SINE=4 not_SINE=7" in out


def test_no_flank_test_is_not_a_pass(tmp_path):
    out, fv, lv = _run(tmp_path, None, {}, sfs=("A", "B"))
    assert set(fv.values()) == {"not_tested"}
    assert set(lv.values()) == {"not_SINE"}
