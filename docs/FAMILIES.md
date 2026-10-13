# Families + subfamilies mode (branch `families-subfamilies`)

Optional. Without `--families` nothing changes. Status 2026-10-13: robust mode (all subfamily consensuses searched, flat vote, then a
family roll-up) is implemented and toy-tested. Fast mode (search family consensuses only) and the two-round vote are NOT yet implemented.

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
