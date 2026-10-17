"""tools/flank_verdict.py on planted cases (pytest). The expansion's end state and the per-copy flank test must decide the
SINE loci (docs/AIM.md); untestable flanks go to the grey zone, the same flank at the same place in many copies is not SINE
(owner's rule 2026-10-18)."""
import json
import os
import subprocess
import sys

TOOL = os.path.join(os.path.dirname(__file__), "..", "tools", "flank_verdict.py")


def _run(tmp_path, step7_rows, loops, sfs=(), n=2, stage9=None, extra=()):
    s2 = tmp_path / "step2" / "step2_output"
    s2.mkdir(parents=True)
    rows = ["Sequence\tSubfamily\tBitscore\tVotes\tStatus"]
    for sf in sorted({r[0] for r in (step7_rows or [])} | set(loops) | set(sfs)):
        rows += ["c%s%d\t%s\t100\t10\tassigned" % (sf, i, sf) for i in range(1, n + 1)]
    rows += ["%s\t%s\t100\t10\tassigned" % (c, sf) for c, sf in extra]
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
    for sf, (status, pairs) in (stage9 or {}).items():
        p = tmp_path / "results" / "flank_twins" / sf
        p.mkdir(parents=True)
        (p / "copy_status.tsv").write_text("".join("%s\t%s\t0\t-\t-\n" % kv for kv in status.items()))
        (p / "flank_groups.tsv").write_text("copy_a\tcopy_b\tside\tidentity\tlength\ttier\tclass\n"
                                            + "".join("%s\t%s\t%s\t95\t100\tT1\t%s\n" % r for r in pairs))
    out = subprocess.run([sys.executable, TOOL, str(tmp_path)], capture_output=True, text=True, check=True).stdout
    fv = {l.split("\t")[0]: l.split("\t") for l in (s2 / "flank_verdict.tsv").read_text().splitlines()[1:]}
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
    st9 = {sf: ({"c%s1" % sf: "unique", "c%s2" % sf: "unique"}, []) for sf in "ABCDE"}
    out, fv, lv = _run(tmp_path, step7, loops, stage9=st9)
    assert {k: v[1] for k, v in fv.items()} == {"A": "interspersed", "B": "untested", "C": "interspersed",
                                                "D": "not_interspersed", "E": "not_interspersed"}
    assert lv["cA1"] == lv["cC2"] == "SINE"
    assert lv["cB1"] == "grey"                        # untested subfamily: grey zone, even if stage 9 says unique
    assert lv["cD1"] == lv["cE1"] == "not_SINE"       # still shared at the expansion limit
    assert lv["cX"] == "not_SINE"                     # not assigned
    assert len(lv) == 11                              # every copy gets a verdict
    assert "copies SINE=4 grey=2 not_SINE=5" in out


def test_per_copy_rule(tmp_path):
    step7 = [("A", "upstream", "confirmed"), ("A", "downstream", "confirmed")]
    hot = ["h%d" % i for i in range(10)]             # same flank at the same place in 10 copies (one side)
    small = ["s1", "s2", "s3"]                       # shared flank in only 3 copies, masked
    status = {"cA1": "unique", "cA2": "masked", "tw": "twin1", "ar": "array", "un": "untestable"}
    status.update({c: "masked" for c in hot})
    status.update({c: "masked" for c in small})
    pairs = [(hot[i], hot[i + 1], 5, "flank_only") for i in range(9)]
    pairs += [(small[0], small[1], 3, "flank_only"), (small[1], small[2], 3, "flank_only")]
    pairs += [("tw", "cA1", 3, "twin")]               # a twin pair does not make its partner's group a hot spot
    extra = [(c, "A") for c in hot + small + ["tw", "ar", "un", "nores"]]
    out, fv, lv = _run(tmp_path, step7, {}, stage9={"A": (status, pairs)}, extra=extra)
    assert lv["cA1"] == "SINE"
    assert lv["cA2"] == "grey"                        # masked flanks: grey zone
    assert all(lv[c] == "not_SINE" for c in hot)      # same repeat in the same place in many copies: not SINE
    assert all(lv[c] == "grey" for c in small)        # in 3 copies only: grey
    assert lv["tw"] == lv["ar"] == "not_SINE"         # flanks shared with another copy
    assert lv["un"] == "grey" and lv["nores"] == "grey"
    assert fv["A"][6:9] == ["1", "6", "12"]           # sine / grey / not_sine per subfamily


def test_no_flank_test_is_not_a_pass(tmp_path):
    out, fv, lv = _run(tmp_path, None, {}, sfs=("A", "B"))
    assert {v[1] for v in fv.values()} == {"not_tested"}
    assert lv["cA1"] == lv["cB2"] == "grey"
    assert lv["cX"] == "not_SINE"
