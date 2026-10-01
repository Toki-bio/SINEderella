# Task L - divergence difference next to the LEAK ratio (decision D3)

## 2026-10-01 - tools/leak_div_example.py created

Reads the step3 table `all_sines.bedlike.ALL.tsv` (step3_postprocess.sh: table built in
section 1, LEAK-annotated in section 2). Layout used: best bitscore col 7, LEAK flag col 8
(runner_ratio >= 0.90), runner-up subfam/bitscore/ratio in the col-10 note (`runner=`,
`runner_bs=`, `runner_ratio=`), single-ssearch bitscore and ratio cols 11-12; copy id
chr:start-end(strand) as in unassigned.tsv.

The table carries bitscores only - no alignment identity, no alignment length - so the
tool says so and falls back to the columns that exist (the scores): divergence proxy =
1 - bitscore / consensus self bitscore, each consensus's self bitscore recovered from its
copies as sim_bitscore / sim_ratio. Note keys `best_id=`/`runner_id=` would be used as
true 1 - identity if a table ever carried them.

Output: number of LEAK copies; top N by divergence difference (runner-up much more
diverged than the best) with copy id, families, both bitscores, ratio, both divergences,
difference; the share of LEAK copies a >5-points difference rule would unflag.

`python tools/leak_div_example.py --selftest` builds a tiny 5-row table in the real column
layout in a temp dir, runs the same logic, asserts the counts/ranking/share, prints
SELFTEST OK.

Still open for D3: run it on a real step3 table and take the top rows to him as the
concrete example the decision asks for.
