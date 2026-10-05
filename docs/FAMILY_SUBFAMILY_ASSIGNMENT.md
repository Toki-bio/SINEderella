# Family first, then subfamily? An evaluation of two-layer assignment (2026-10-05)

Status: **evaluation only, nothing in the code changed.** Asked by the owner after the MEG-RS mini-genome finding (the same copies
voted MEG-RL on a small run and MEG-RS on the whole genome): "currently we try to assign copies into subfamilies, but probably it
should be a 2-layers process - since families are supposed to have sharp edges, first step should assign to families, and then
more fine-graded assignment be done with subfamilies based on initial family grouping ... a tectonic change and needs super
thorough evaluation."

## 1. Where this is in SINEderella today

Step 2, `step2_asSINEment.sh` (asSINEment). The bank is a **flat list**: every consensus, whether a family, a subfamily, a length
version or a composite, is one equal competitor.

* Every extracted copy is searched with every consensus, 10 cycles (`ssearch36 -z 11`, library parts of 20 000 copies). Per cycle the
  consensus with the highest bitscore wins.
* Firm = the same consensus wins all 10 cycles, and the summed score passes 0.45 x the 10th-best sum of that consensus.
* Otherwise the copy is a soft call (`no_unanimous`, `rejected_low_bitscore`): kept, labelled with its best query, but not counted as
  a member, not used by the array flag, the twin check, the length-version test, the consensus audit, or (except top-up) the plates.
* Step 3 then flags LEAK (runner-up >= 0.90 of the best) and CONFLICT.

So one label per copy, chosen among all consensuses at once. There is no notion of family in the code.

## 2. What was measured

The step-2 vote was reproduced exactly (same `ssearch36` call, same 20 000-copy parts) keeping, for every copy and cycle, the score
of each of its top six consensuses. The same ten cycles were then re-tallied two ways: **flat** (today) and **two-layer** (layer 1:
the family whose best member scores highest must win 10 of 10; layer 2: inside that family, the best member). Scripts and data:
`tools/vote_eval/` (`vote_full.sh`, `retally.py`); data and log on therioserver `~/tmp/votetest2/`.

Family definitions used (placeholders, see section 5): human chr21 = one family Alu (AluJb, AluSx, AluY: textbook subfamilies of
one family); mouse chr19 = B1 and B2 (two different families); rsi = two families, MEG (four MEG consensuses) and R (the 17 r-type
consensuses including the seven composites). The rsi split is coarse and only stands in for the owner's definition.

| data set | copies | firm today | family firm (layer 1) | splits across families | splits inside a family |
|---|---|---|---|---|---|
| rsi, whole genome, run A | 66 568 | 97.9 % | 100.0 % | 1 | 1 416 |
| rsi, whole genome, run B (same input) | 66 568 | 97.7 % | 100.0 % | 1 | 1 528 |
| rsi mini genome (351 copies, one library part) | 351 | 49.9 % | 100.0 % | 0 | 176 |
| rsi mini, part padded to 20 000 with shuffled decoys | 348 | 72.4 % | 99.4 % | 0 | 94 |
| hs21, Alu | 10 521 | 96.6 % | 100.0 % | 0 | 354 |
| mm19, B1 and B2 | 16 474 | 100.0 % | 100.0 % | 6 | 0 |

Reproducibility (run A versus run B, identical input, only the shuffles differ): the flat outcome of a copy (firm member, or not
firm) is the same for 99.37 % of copies; the family outcome for 100.00 %.

### What the vote failures are

| data set | inside a family: composite and a unit it contains | inside a family: two monomers | across families |
|---|---|---|---|
| rsi run A (1 417) | 1 175 (83 %): r5_r6_P26 vs r6 435, C11 vs P1 207, P48 vs P26 142, ... | 241 (17 %): r5 vs r8 99, r5 vs r6 30, r6 vs r8 28, r7 vs r9 27, ... | 1 |
| rsi mini (176) | 44 (25 %) | 132 (75 %): MEG-RL vs MEG-RS 123 | 0 |
| hs21 (354) | 0 | 354: AluSx vs AluY 222, AluJb vs AluSx 132 | 0 |
| mm19 (6) | 0 | 0 | 6: B1 vs B2 |

("composite and a unit": one of the two competing consensuses is a composite, P/C names, and the other is or shares one of its units.
In 363 of the 1 175 rsi cases the copy covers less than 60 % of the composite: a lone unit that local alignment scores equally
against the unit and against the composite containing it.)

## 3. What the numbers say

1. **The owner's premise holds in every data set.** Families have sharp edges for the vote: of 1 994 failed votes across four data
   sets, 7 are between families (1 in rsi, 6 between B1 and B2), the rest inside one. A family call is firm for 99.4-100 % of copies,
   reproduces 100 % between two runs, and does not depend on the library size (the 351-copy mini: 100 % family-firm against 50 %
   firm today). The MEG-RS / MEG-RL instability is a subfamily problem, not a family one.

2. **Two layers with the same bitscore rule at layer 2 change no subfamily call.** By construction the flat winner is the best member
   of the winning family, and the measured counts agree to the copy ("subfamily firm" = "flat firm" in all six rows; no copy is flat-firm
   but family-unsure). What two layers add is a firm **family** label for the ~2-4 % of copies that are soft today (1 416 rsi copies,
   354 Alu copies, all 176 of the mini's split copies), with the subfamily question kept open and stated, instead of the copy dropping
   out of counts, arrays, twins, audits and plates. That is a real gain for the family-level questions (is it a SINE, how many copies,
   arrays, flanks), but it does not by itself make subfamilies better.

3. **The subfamily layer needs its own criterion, and that is where the work is.** Whole-length bitscore is the wrong instrument for
   subfamilies: it sums agreement over every column, so private decay and a few diagnostic columns weigh the same, and two subfamilies
   97 % identical differ by a handful of points that the shuffling noise covers (rsi near-ties within 2 %: 21 115 cycle-copies inside
   families, 10 across). The definitions in the manuscript already say what layer 2 should read: a shared pattern of diagnostic
   mutations, plus the wave seen in divergence. Candidates in the code or the record:
   * the diagnostic positions of `step4_diagnostic.py` (MI / KL / random forest; needs labels, so it would be trained on the
     subfamily-firm copies and applied to the split ones);
   * the owner's own test (50 best copies of each subfamily, `mafft --reorder`, do they separate; memory note
     `sine_his_separation_test_50best`), as a check that two members deserve to be separate subfamilies at all;
   * a "subfamily unresolved" outcome as a first-class result (family-firm, subfamily between A and B), which is honest and costs nothing.

4. **Composites are not subfamilies, and a family/subfamily tree does not fit them.** In rsi, 83 % of today's failures are a composite
   against a unit it contains. That is a containment question (does the copy contain both units or only one), answered by alignment
   coverage, not by a vote. Where a composite sits in the hierarchy is an owner decision (its own family, by the Kramerov-Vassetzky rule
   that a change of structural parts makes a family; or a member of the family of its head unit; or outside the tree, handled by the
   flank scan's element model). Until that is decided, a two-layer design for rsi is undefined for 7 of its 21 consensuses.

5. **The 0.45 threshold must stay per subfamily (or be made relative), not per family.** Applied per family, 0.45 x 10th-best rejects
   41 014 of 66 567 rsi copies (today 1 782), 116 of 10 521 Alu copies (today 20): the family's top scores come from its longest or
   youngest member (composites of 300-600 bp, AluY), and short or older members (r9 105 bp, AluJb) fall below. The threshold was
   designed against short or decayed fragments of one consensus, so it belongs to layer 2, or should be expressed as score / self-score
   of the member (the step-3 similarity ratio), which is length-free.

6. **Library size still matters inside a family.** Padding the mini's one small part to 20 000 with shuffled decoys raised subfamily
   firmness from 50 % to 72 %; on the whole genome the small last part (6 568 copies) voted 96.0 % unanimous against 97.2-98.4 % in the
   full parts, padded 98.2 %. Equal parts or padding changed 1.3-1.7 % of whole-genome outcomes against a 0.6 % run-to-run noise floor.
   Two layers remove this from the family call (point 1) but not from the subfamily call; padding is a separate, small fix for layer 2.

## 4. What a two-layer SINEderella would touch

| part | change |
|---|---|
| bank | each consensus carries a family label (header field or a side table); a bank without labels = every consensus its own family, which reproduces today's results exactly (a property to test) |
| step 2 | the ten cycles stay as they are (no new search); the tally gives family (unanimous) then subfamily (rule of section 3.3); thresholds per subfamily; new statuses: family-firm/subfamily-firm, family-firm/subfamily-unresolved (with the candidates), not family-firm |
| `--add` / `--exclude` | adding a subfamily re-votes only that family's copies; adding a family re-votes copies whose family can change; simpler than today's overlap rule |
| step 3 | LEAK inside a family becomes expected; LEAK across families becomes the alarm |
| satellite screen, array flag, twin check | naturally per family (one array is now found once per member consensus) |
| length versions, consensus audit | per subfamily, unchanged |
| plates, report | family plates (all family-firm copies) beside subfamily plates; counts per family and per subfamily |
| downstream readers of `assignment_full.tsv` column 2 | 27 files read it (`grep` count); a new column for the family keeps them working, then each is moved deliberately |
| manuscript | the Assignment section, Table 1, principle 1 and the criteria text (families differ by parts, subfamilies by lineage) would finally match the code |

## 5. Decisions only the owner can make, before any code

1. **What the families are in each test bank.** For rsi: is "R" one family (all heads tRNA(Ile)-derived) or several (by 3' part:
   r3-type, r5/r6/r7/r8-type, r9, r10...)? Are r9 / r7 / r8 (length versions, TWO_VERSIONS) subfamilies or families? MEG: one family?
2. **Where composites go** (section 3.4).
3. **The layer-2 rule**: bitscore vote as now (then only the firm family label is new), diagnostic positions, or "unresolved" as the
   default with diagnostics as the reader's aid.
4. **Whether a family-firm, subfamily-unresolved copy counts** as a member of the family (for plates, arrays, twins) and of neither
   subfamily.

## 6. Proposed next steps (each cheap, none touches the pipeline)

1. With the owner's family definitions for the rsi bank, re-run `retally.py` on the existing cycles (minutes; the vote data are kept).
2. External oracle for layer 2 on hs21: the 354 Alu copies split between AluSx/AluY and AluJb/AluSx have RepeatMasker subfamily calls
   (UCSC rmsk, already used by the release check). Compare bitscore vote, diagnostic positions and "unresolved" against them.
3. Composite rule test on rsi: for the 1 175 composite-vs-unit splits, decide by alignment coverage and check on plates.
4. Only after 1-3: a design for step 2, toy-tested, then the mini genome, then one whole-genome run.
