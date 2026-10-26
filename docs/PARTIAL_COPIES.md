# Partial copies (branch partial-copies, 2026-10-26)

## Problem
A copy that matches only a stretch of its consensuses (a lone unit of composites, a 5' or 3' piece) cannot pass step 2:
- the bar of a consensus is 0.45 x its 10th-best 10-cycle sum, set by full-length copies; a 155-bp copy against a
  300-600 bp composite stays below it (`rejected_low_bitscore`);
- a stretch shared by several composites scores almost the same against each of them, so `-z 11` noise decides the cycle
  winner and the vote is split (`no_unanimous`).
rsi: r1_9seqs (155 bp, 21 copies) is in the bank only for this; without it its copies end as family-firm Rhin with
`layer2=unresolved` (rerun therioserver `/home/toki/rhin/rsi_nor1_20261026/`).

## Rule
1. Step 2 keeps, per copy and consensus, where the copy's best HSP lies on the consensus (m8 columns 7-8) in its best cycle,
   and the same for the runner-up: `vote_summary.tsv.gz` columns 7-10 (columns 1-6 unchanged). No new search.
2. `tools/rollup_families.py`, for copies that become family-firm with the subfamily open: competitors = the family's
   consensuses that won a cycle, plus the winner's runner-up when it scores >= 90 % of it. Partial when the copy covers
   < 60 % of every competitor and >= 50 bp of the best one. Bar = the best consensus's threshold x stretch length /
   consensus length. Passing copies: `layer2=partial:<cons:start-end,...>`; all candidates in
   `family_view/partial_copies.tsv` (status `partial` or `below_stretch_threshold`).
3. Nothing else changes: subfamily-firm copies, ties, the subfamily view and `-z 11` are untouched.

## Test (therioserver `/home/toki/rhin/rsi_nor1_20261026/partial_test/`, `mini.sh`)
20,000 random copies of the no-r1 run plus the 19 copies overlapping the 21 r1 loci of rsi_main_20261022; step 2 (this
branch) + family roll-up, bank without r1.
- r1 loci: 18 of 19 labelled partial: 10 with r1_r3_P2 best (head, `r1_r3_P2:2-83` to `2-113`, 82-112 bp), 8 longer
  (head + module, best C11 or P1, 119-160 bp); 1 is a full r1_r3_P1 copy (subfamily-firm). Coverage 0.20-0.37, every one above its
  stretch bar (685-1401 against 355-671).
- All copies: 232 partial (1.2 %), 1 below its stretch bar, 941 still unresolved, of 1,173 subfamily-open; partial copies
  of C11 53, P1 37, P2 36, P18 36, r6 17, P26 14, others <= 8.
- Tests: pytest 131 passed on therioserver; step-2 shell tests pass.

## Open
- Labels are per copy (`r1_r3_P2:2-89`, `2-83`, `2-93`): the report needs them grouped into stretches (overlap on the best
  consensus) before they can be counted and plotted.
- The family bar is the lowest member threshold (rsi Rhin: 158, from r4_32seqs with 12 copies), lenient.
- Not in the subfamily list (loci_verdict) or the report yet.

## Stretches (commit 9ad52fa)
Partial copies are grouped: two copies of one family are linked when, on a consensus both match, their starts and their
ends each differ by at most max(10 bp, 10 % of the shorter span); a stretch is a connected group, named
`<family>_partial<n>` (largest first), with the median span on every consensus that at least half of its copies match.
`layer2=partial:<stretch id>`; `family_view/partial_stretches.tsv` (Stretch, Family, Copies, Median_length, Spans);
`partial_copies.tsv` gets the copy's Stretch_id. On the test data: 232 copies, 41 stretches (21 of one copy); r1 loci in
Rhin_partial2 (head, P1 / P2 2-89, 9 loci), Rhin_partial7 (head + module, C11 / P1 2-141, 5), Rhin_partial11 (P2 2-116, 2),
Rhin_partial1 (P1 2-200, 1) and one single-copy stretch. Assignments, votes and family calls identical to the ungrouped run.

## What the rule cannot reach: step 1 keeps only hits covering >= 80 % of a query
`sear` keeps a hit only when it covers >= Slf (0.8) x the query length. A copy that is a part of a composite and not 80 %
of any bank entry is never extracted, so no step-2 rule sees it (rsi without r1: 3 of the 21 r1 loci).

Test (therioserver `/home/toki/rhin/rsi_nor1_20261026/blocktest/`, `blocktest.sh`, `blockcmp.py`, `blockcmp.txt`): step 1
search (same code, cleaned genome and settings as the no-r1 run) with the 18 composite stretches >= 50 bp of
`composite_units.tsv` as queries of their own; hits merged into loci (bedtools merge) and compared with the no-r1 run's
loci before its satellite screen (67,299). The satellite screen of the test itself did not run (trf not on PATH outside
the orchestrator), so no satellite loci were removed from the test's side.
- 71,533 piece loci; 1,570 overlap no locus of the bank without r1 (2.3 %).
- All 21 r1 loci are among the piece loci, the 3 the bank without r1 missed included.
- New loci per query (a locus found by several queries counts for each): P18 290-426 (r8 module) 492, C11 472-589
  (r3-type end) 345, P18 152-289 (the 10-bp-unit middle piece) 334, the r1-type head of C11 / P1 / P2 283 / 281 / 247,
  every other stretch <= 70.
- Not yet known: whether the new loci are SINE loci. They need the vote and the flank verdict (a full run with the pieces
  in the search bank); the middle piece of P18 is a tandem of 10-bp units and may mostly find simple repeats.

## Full run on the whole genome (therioserver `/home/toki/rhin/rsi_nor1_partial_20261026/run_20261026_183325`, 9ad52fa)
Bank without r1, otherwise as rsi_nor1_20261026; finished 19:50, report built.
- `vote_summary.tsv.gz` carries the spans; 690 partial copies (4 below their stretch bar) in 40 stretches.
- Largest stretches: Rhin_partial1 196 copies (C11 2-180 = P1 2-168: the r1-type head and module), Rhin_partial2 145
  (C11 337-606), Rhin_partial3 139 (P26 2-145 = P48 1-141, the r5-type head), Rhin_partial4 125 (P18 283-418),
  Rhin_partial5 34 (r10 2-81 = P18del 2-81); the other 35 hold 8 copies or fewer.
- The 21 r1 loci of rsi_main_20261022: 16 in Rhin_partial1, 1 in Rhin_partial38, 1 a full r1_r3_P1 copy (tie), 3 not
  extracted (step 1, see above).
- Subfamily assignments against the run without the branch (copies matched by overlap, `moved.py`): every subfamily
  within 44 copies (r1_r3_P1 14,457 -> 14,413), the vote's run-to-run spread.
- Grouping at this scale: copies of the head alone (2-89) and of head + module (2-161) are joined through copies of
  intermediate length into one stretch (single linkage). Whether they should stay apart is open.
