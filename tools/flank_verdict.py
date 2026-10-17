#!/usr/bin/env python3
"""flank_verdict.py RUN_ROOT

Turns the flank-expansion results into verdicts that the outputs must obey (docs/AIM.md: interspersed repeats are copies
with UNIQUE flanks; a test that could not be done is not a pass).

The expansion: step7_boundary_refine.sh widens each subfamily's flanks (50 bp steps, up to 1 kb) until they fall to
background ("confirmed") or the limit is reached ("undetermined"); for undetermined subfamilies the publish border loop
(publish/border_loop_subfam.py) extends the element itself and rescans each side until it is independent or a cap is hit.
Until 2026-10-17 nothing read the end state: a subfamily whose flanks were still shared at the limit kept all its copies
as SINE loci, and a side with nothing left to scan counted as "independent" (rsi MEG-RS, a tandem-array family).

Per subfamily and side, the result that decides:
  - border_loop/<sf>/border_loop.json when the border loop ran for it (it ran because step7 was not confirmed):
      verdict_<side> = independent | still_shared | untested
      (older files without verdict_<side>: stop_reason "independent (island_fraction=None)" -> untested,
       "independent ..." -> independent, "... STILL BAD" -> still_shared, anything else -> untested)
  - else step7 boundary_refinement.tsv, population "general":
      confirmed -> independent, undetermined -> still_shared, anything else -> untested
Subfamily verdict: interspersed (both sides independent) | not_interspersed (a side still shared at the limit) |
untested (a side could not be tested) | not_tested (no flank test ran in this run).
Copy verdict: SINE only if the copy is assigned AND its subfamily is interspersed; every other copy is not_SINE with the
reason. Every copy of assignment_full.tsv gets a verdict.

Out (in the step2 output folder, next to boundary_refinement.tsv):
  flank_verdict.tsv  subfamily verdict side5 side3 source assigned_copies reason
  loci_verdict.tsv   copy subfamily assignment_status verdict reason
Prints one line: FLANK_VERDICT subfamilies interspersed=.. not_interspersed=.. untested=.. not_tested=.. copies SINE=.. not_SINE=..
"""
import glob
import json
import os
import sys

run_root = sys.argv[1].rstrip("/")
s2s = sorted(glob.glob(os.path.join(run_root, "step2", "step2_output*")), key=os.path.getmtime)
if not s2s:
    sys.exit("flank_verdict: no step2/step2_output* in %s" % run_root)
s2 = s2s[-1]
assign = os.path.join(s2, "assignment_full.tsv")
if not os.path.isfile(assign):
    sys.exit("flank_verdict: missing %s" % assign)

STEP7 = {"confirmed": "independent", "undetermined": "still_shared"}


def from_stop_reason(r):
    r = r or ""
    if r.startswith("independent (island_fraction=None)"):
        return "untested"
    if r.startswith("independent"):
        return "independent"
    if "STILL BAD" in r:
        return "still_shared"
    return "untested"


step7 = {}                                   # sf -> {"5": v, "3": v}
bt = os.path.join(s2, "boundary_refinement.tsv")
if os.path.isfile(bt):
    with open(bt) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[1] != "general":
                continue
            side = "5" if p[2] == "upstream" else "3"
            step7.setdefault(p[0], {})[side] = STEP7.get(p[4], "untested")

loop = {}
for js in glob.glob(os.path.join(run_root, "border_loop", "*", "border_loop.json")):
    try:
        d = json.load(open(js))
    except Exception as e:                   # an unreadable result is not a pass
        sf = os.path.basename(os.path.dirname(js))
        loop[sf] = {"5": "untested", "3": "untested", "why": "unreadable border_loop.json (%s)" % e}
        continue
    sf = d.get("subfam") or os.path.basename(os.path.dirname(js))
    v = {}
    for side, key in (("5", "5prime"), ("3", "3prime")):
        v[side] = d.get("verdict_" + key) or from_stop_reason(d.get("stop_reason_" + key))
        v["why" + side] = d.get("stop_reason_" + key) or "no result"
    loop[sf] = v

copies, nassigned = [], {}
with open(assign) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        copies.append((p[0], p[1], p[4]))
        if p[4] == "assigned":
            nassigned[p[1]] = nassigned.get(p[1], 0) + 1

tested_any = bool(step7) or bool(loop)
sfv = {}
for sf in sorted(nassigned):
    if sf in loop:
        src, v = "border_loop", loop[sf]
        why = "; ".join("%s' %s" % (s, v.get("why" + s, v.get("why", ""))) for s in ("5", "3"))
    elif sf in step7:
        src, v = "step7", step7[sf]
        why = "step7 general sample: 5' %s, 3' %s" % (v.get("5", "untested"), v.get("3", "untested"))
    else:
        src, v = "-", {}
        why = "no flank test result for this subfamily" if tested_any else "flank test (step7) not run in this run"
    s5, s3 = v.get("5", "untested"), v.get("3", "untested")
    if not v:
        verdict = "untested" if tested_any else "not_tested"
    elif "still_shared" in (s5, s3):
        verdict = "not_interspersed"
    elif s5 == "independent" and s3 == "independent":
        verdict = "interspersed"
    else:
        verdict = "untested"
    sfv[sf] = (verdict, s5 if v else "-", s3 if v else "-", src, why)

with open(os.path.join(s2, "flank_verdict.tsv"), "w") as fo:
    fo.write("subfamily\tverdict\tside5\tside3\tsource\tassigned_copies\treason\n")
    for sf, (verdict, s5, s3, src, why) in sfv.items():
        fo.write("%s\t%s\t%s\t%s\t%s\t%d\t%s\n" % (sf, verdict, s5, s3, src, nassigned[sf], why))

REASON = {"not_interspersed": "flanks still shared at the expansion limit (not interspersed)",
          "untested": "flank uniqueness could not be tested", "not_tested": "flank test not run"}
nsine = 0
with open(os.path.join(s2, "loci_verdict.tsv"), "w") as fo:
    fo.write("copy\tsubfamily\tassignment_status\tverdict\treason\n")
    for c, sf, st in copies:
        if st != "assigned":
            verdict, why = "not_SINE", "not assigned (%s)" % st
        else:
            sv = sfv[sf][0]
            verdict, why = ("SINE", "assigned; subfamily flanks unique") if sv == "interspersed" else ("not_SINE", REASON[sv])
        nsine += verdict == "SINE"
        fo.write("%s\t%s\t%s\t%s\t%s\n" % (c, sf, st, verdict, why))

tally = {k: 0 for k in ("interspersed", "not_interspersed", "untested", "not_tested")}
for verdict, *_ in sfv.values():
    tally[verdict] += 1
print("FLANK_VERDICT subfamilies " + " ".join("%s=%d" % kv for kv in tally.items())
      + " copies SINE=%d not_SINE=%d" % (nsine, len(copies) - nsine))
