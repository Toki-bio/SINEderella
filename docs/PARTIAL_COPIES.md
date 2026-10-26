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
