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
Copy verdict (owner's rule, 2026-10-18), for an assigned copy of an interspersed subfamily, from flankscan stage 9
(results/flank_twins/<sf>/copy_status.tsv and flank_groups.tsv):
  - not_SINE: a twin / array copy (flanks shared with another copy, bodies alike), or the same flank sequence sits at the same
    place in many copies: the copy is in a group of >= SAME_PLACE_MIN (10) copies linked on one side by stage-9 "flank_only"
    pairs (shared from the junction outward, independent bodies);
  - grey: flanks masked (repeat) or untestable and not in such a group, or no stage-9 result for the copy;
  - SINE: stage 9 says unique and the copy is in no such group.
A copy of a not_interspersed subfamily is not_SINE; of an untested / not_tested subfamily grey; an unassigned copy not_SINE.
Every copy of assignment_full.tsv gets a verdict.

Out (in the step2 output folder, next to boundary_refinement.tsv):
  flank_verdict.tsv  subfamily verdict side5 side3 source assigned_copies sine grey not_sine reason
  loci_verdict.tsv   copy subfamily assignment_status verdict reason          (verdict = SINE | grey | not_SINE)
Prints one line: FLANK_VERDICT subfamilies interspersed=.. not_interspersed=.. untested=.. not_tested=.. copies SINE=.. grey=.. not_SINE=..
Env: SAME_PLACE_MIN=10
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
assigned_fa = os.path.join(s2, "assigned.fasta")
if not os.path.isfile(assign) and not os.path.isfile(assigned_fa):
    sys.exit("flank_verdict: neither %s nor %s" % (assign, assigned_fa))

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
if os.path.isfile(assign):
    with open(assign) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                copies.append((p[0], p[1], p[4]))
else:                                        # a run dir with assigned.fasta only (headers: copy|subfamily|bitscore)
    with open(assigned_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].strip().split("|")
                if len(p) >= 2:
                    copies.append((p[0], p[1], "assigned"))
for c, sf, st in copies:
    if st == "assigned":
        nassigned[sf] = nassigned.get(sf, 0) + 1

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

SAME_PLACE_MIN = int(os.environ.get("SAME_PLACE_MIN", 10))
FT = os.path.join(run_root, "results", "flank_twins")


def stage9(sf):
    """copy -> stage-9 status, and copy -> size of its largest one-side flank_only group. (None, None) if stage 9 did not run."""
    cs, fg = os.path.join(FT, sf, "copy_status.tsv"), os.path.join(FT, sf, "flank_groups.tsv")
    if not os.path.isfile(cs):
        return None, None
    status = {}
    for line in open(cs):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "copy":
            status[p[0]] = p[1]
    up = {}

    def root(x):
        while up[x] != x:
            up[x] = up[up[x]]
            x = up[x]
        return x
    if os.path.isfile(fg):
        for line in open(fg):
            p = line.rstrip("\n").split("\t")
            if len(p) < 7 or p[6] != "flank_only":
                continue
            a, b = (p[2], p[0]), (p[2], p[1])          # groups are formed per side
            for k in (a, b):
                up.setdefault(k, k)
            ra, rb = root(a), root(b)
            if ra != rb:
                up[ra] = rb
    size = {}
    for k in up:
        r = root(k)
        size[r] = size.get(r, 0) + 1
    group = {}
    for k in up:
        group[k[1]] = max(group.get(k[1], 0), size[root(k)])
    return status, group


SUBFAM_REASON = {"not_interspersed": ("not_SINE", "subfamily flanks still shared at the expansion limit (not interspersed)"),
                 "untested": ("grey", "subfamily flank uniqueness could not be tested"),
                 "not_tested": ("grey", "flank test not run")}
s9 = {sf: stage9(sf) for sf in sfv}
per = {sf: {"SINE": 0, "grey": 0, "not_SINE": 0} for sf in sfv}
allc = {"SINE": 0, "grey": 0, "not_SINE": 0}
with open(os.path.join(s2, "loci_verdict.tsv"), "w") as fo:
    fo.write("copy\tsubfamily\tassignment_status\tverdict\treason\n")
    for c, sf, st in copies:
        if st != "assigned":
            verdict, why = "not_SINE", "not assigned (%s)" % st
        elif sfv[sf][0] != "interspersed":
            verdict, why = SUBFAM_REASON[sfv[sf][0]]
        else:
            status, group = s9[sf]
            cst = status.get(c) if status is not None else None
            g = group.get(c, 0) if group is not None else 0
            if cst in ("twin1", "twin2", "array"):
                verdict, why = "not_SINE", "flanks shared with another copy (%s)" % cst
            elif g >= SAME_PLACE_MIN:
                verdict, why = "not_SINE", "same flank sequence at the same place in %d copies" % g
            elif cst is None:
                verdict, why = "grey", "no stage-9 flank result for this copy"
            elif cst != "unique":
                verdict, why = "grey", "flanks %s (repeat), not testable" % cst
            else:
                verdict, why = "SINE", "assigned; flanks unique"
        allc[verdict] += 1
        if st == "assigned":
            per[sf][verdict] += 1
        fo.write("%s\t%s\t%s\t%s\t%s\n" % (c, sf, st, verdict, why))

with open(os.path.join(s2, "flank_verdict.tsv"), "w") as fo:
    fo.write("subfamily\tverdict\tside5\tside3\tsource\tassigned_copies\tsine\tgrey\tnot_sine\treason\n")
    for sf, (verdict, s5, s3, src, why) in sfv.items():
        q = per[sf]
        fo.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n"
                 % (sf, verdict, s5, s3, src, nassigned[sf], q["SINE"], q["grey"], q["not_SINE"], why))

tally = {k: 0 for k in ("interspersed", "not_interspersed", "untested", "not_tested")}
for verdict, *_ in sfv.values():
    tally[verdict] += 1
print("FLANK_VERDICT subfamilies " + " ".join("%s=%d" % kv for kv in tally.items())
      + " copies SINE=%d grey=%d not_SINE=%d" % (allc["SINE"], allc["grey"], allc["not_SINE"]))
