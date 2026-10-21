# Families + subfamilies mode

Optional. Without `--families` nothing changes. On main since 2026-10-13: robust mode (all subfamily consensuses searched, flat vote,
then a family roll-up) and fast mode (`--families-mode fast`, below); the two-round subfamily call is not made (Not done). With a
family table the page's "Subfamily alignments" table groups the subfamilies under one row per family with the family's totals
(assigned copies, SINE loci, grey, not SINE; 2026-10-21).

## Use

    SINEderella genome.fa consensuses.fa --families families.tsv [--parents parents.fa] [--family-weights counts.tsv]
    SINEderella --resume --run RUN --families families.tsv

`families.tsv` (tab separated, `#` comments, header optional): `subfamily  family  [parent_consensus]`. A consensus not listed is its own
single-member family. Parent consensus: supplied (third column, name in `--parents` or in the bank) > derived (MAFFT of the members,
weighted plurality, weights from `--family-weights`, else equal) > the member itself for single-member families. Checked by
`tools/families.py` (typos, duplicates, name collisions are errors).

## What it writes

* `<run>/families/` : `families.clean.tsv`, `family_consensuses.fa`, `families.summary.tsv`.
* `<run>/family_view/` : a run directory in which every label is the family (`tools/rollup_families.py`), post-processed by the normal
  step 3 (LEAK, CONFLICT, soft assignment, `summary.by_subfam.tsv` = per family). Subfamily results are untouched.
* step 2 now keeps `step2_output/vote_summary.tsv.gz` (per copy, per consensus: votes, bitscore sum, runner-up). The roll-up sums the votes
  of a copy over the members of a family: 10 of 10 in one family and a score above the family's lowest threshold = family-firm, tagged
  `layer2=unresolved` (subfamily left open). A runner-up from the copy's own family is not a LEAK at family level.

## Toy tests (therioserver `~/tmp/families_toy`, 2026-10-13)

| toy | families | subfamily-firm | family-firm | note |
|---|---|---|---|---|
| hs21 Alu, 10 537 copies | Alu = AluJb, AluSx, AluY | 10 147 | 10 522 (+375) | 15 copies fail the bitscore bar; LEAK 0 at family level |
| Sicista, 18 019 copies (<=600 per consensus and status) | B1fam = B1, B1-dID, pB1 | 2 232 | 2 746 (+514) | 217 copies stay no_unanimous (split across families); runner-up LEAKs 1 177 per subfamily -> 296 per family |

The Sicista group tests the composites question: B1-dID and pB1 grouped with B1 as members of one family. Their copies largely resolve to the
family; what the subfamily vote splits (B1 / B1-dID / pB1) is the `layer2=unresolved` population. Not yet done: single-member-family test of the
composites, sweep of the near-tie threshold, separation test, fast mode.

## Separation test (`tools/separation_test.py`, first push)

`separation_test.py --assignment assignment_full.tsv --copies extracted.fasta --a SUBA --b SUBB [--n 50 --draws 20]`. Tier 1: the 50 best (bitscore) firm
copies of each; tier 2: resampled random 50 (default 20 draws). Each set is aligned together (MAFFT); k=5 nearest-neighbour purity and the number of
label runs along the UPGMA leaf order are compared with 500 label shuffles; the control splits one subfamily at random in two and must look like noise.
"Holds" = both p <= 0.01. No threshold is calibrated; his calls decide.

First results on the toys (10 random draws each; control never looked like a real split):

| pair | best-50 purity / runs | random-50 holds |
|---|---|---|
| AluJb vs AluSx | 0.81 / 3 | 10/10 |
| AluSx vs AluY | 0.66 / 3 | 10/10 |
| AluJb vs AluY | 0.90 / 3 | 10/10 |
| B1 vs B1-dID | 0.62 / 9 | 10/10 |
| B1 vs pB1 | 0.91 / 5 | 10/10 |
| B1-dID vs pB1 | 0.64 / 28 | 7/10 |
| B1 vs B4 (different families, positive control) | 1.00 / 2 | 10/10 |

Reading: the weakest pair is B1-dID vs pB1 (28 runs, 7/10): the candidate for "not separable at this depth", to be judged by him. Purity is lowest where
one consensus is a composite containing the other. Purity alone is not enough for ranking pairs; the run count separates "blocks" from "interleaved".

## Layer 2: near-ties between subfamilies of one family

A subfamily call (10/10 for one subfamily) whose runner-up is another subfamily of the SAME family is a **near-tie** when its score advantage is
below `min(10 %, max(--near-tie 0.02, --tie-rel 0.5 x divergence of the two consensuses))`; divergence = share of differing columns among the
columns where both consensuses have a base (pairwise MAFFT). Such copies keep their subfamily label, are tagged `layer2=tie;tie_with=<other>` in
the family view and counted in `family_view.summary.tsv` (`subfamily_near_tie`). Env in the orchestrator: `NEAR_TIE`, `TIE_REL`. Sweep on the
toys (`--near-tie 0.02`): Alu tie_rel 0 / 0.5 / 1 / 2 = 1.3 / 1.9 / 9.2 / 45.9 % of the subfamily-firm copies; Sicista B1 group 0.5 / 4.6 / 4.6 / 4.6 %
(the 10 % cap). Divergence B1-dID vs B1 is 23 % over shared columns (55 % when gaps counted: the inserted element is not a score difference per column).
0.5 is a starting value; his calls on which copies are really between two subfamilies decide it.

## Fast mode (`--families-mode fast`, full runs)

Bank of the run = the family consensuses (one query per family); afterwards `tools/fast_split.sh` assigns each family's firm copies with the ordinary
step 2 against the consensuses of that family's subfamilies only (`<run>/subfamily_split/subfamily_assignment.tsv`). Toy (hs21 Alu, family "Alu" = AluJb, AluSx,
AluY, full run from the genome): 9 181 copies firm in both modes, all 9 181 with the same subfamily; fast mode finds 137 copies the robust run does not and
misses 197 it finds (the one-query-per-family recall difference), and leaves more copies without a unanimous subfamily (AluSx 852 vs 276). With a
single query per family the family-level counts depend on how well the parent consensus represents all its subfamilies.

## Not done

* The two-round vote as a separate mode (family votes first, then the subfamily vote among members) is what the roll-up with `vote_summary.tsv.gz` already
  computes for the family call; a subfamily call by plurality among members for copies that are family-firm but not subfamily-firm is NOT made (they stay `layer2=unresolved`).
* Composites (B1-dID, pB1) as single-member families (the other branch of the test): only the grouped case was run.
* Diagnostic weighting / edge downweight / COSEG (second push).
* Byte-identity of a run without `--families` with the old code cannot be tested (ssearch36 -z 11 is deliberately non-deterministic); the only change on that path
  is the extra `vote_summary.tsv.gz` file in the step 2 output, and `families_step` returns immediately without `--families`.
