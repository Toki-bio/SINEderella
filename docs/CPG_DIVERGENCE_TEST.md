# CpG-adjusted divergence: the cross-species test (Task C)

Related: docs/BORROWED_TRIAGE.md, "CpG-adjusted (Kimura) divergence" (section 5, item 2) and decision D2.
Tools: tools/cpg_div/cpg_div.py (one plate), tools/cpg_div/run_all.py (all species).
Data: tools/cpg_div/summary.tsv.

**STATUS: the numbers in this document are PENDING.** The method, the decision
rule and the limitations below are final. Every number must be taken from
tools/cpg_div/summary.tsv; nothing is invented. To produce the table:

    python tools/cpg_div/run_all.py

Per-copy detail for a single plate:

    python tools/cpg_div/cpg_div.py C:/work/Tal/SPECIES/alignments/FAMILY.top100.aln.fa

## Why this test

Mammalian SINE heads are CpG-rich, and methylated CpG deaminates fast (C->T on
one strand, G->A on the other). A CpG-rich family therefore accumulates
transitions at CpG sites much faster than its true substitution rate, so raw
identity overstates its divergence — and its apparent age. Two families of the
same real age can sit at different points of a divergence profile just because
one is CpG-richer. docs/BORROWED_TRIAGE.md (decision D2) asks for a thorough
test on multiple species and SINEs before anything in the pipeline changes.
This is that test. It changes nothing in SINEderella; it measures the size of
the effect on the published plates.

## What was measured

Plates: C:/work/Tal/*/alignments/*top100.aln.fa, one species per folder. Plate
format: row 1 is the consensus, the other rows are copies; lower case = flank,
upper case = element; '-' = gap.

For every copy, over the element columns (consensus upper-case and not a gap)
where the copy itself is not a gap:

- raw_p: mismatches / aligned positions — the mismatch proportion the current
  divergence profiles are built from.
- cpg_adj_p: the same ratio after RepeatMasker's CpG-adjusted rule. CpG sites
  are consensus positions that are the C or the G of a CG dinucleotide
  (consensus gaps removed). At those sites a transition C->T or G->A counts
  1/10 of a substitution; every other change counts 1.
- k2p_adj: Kimura 2-parameter distance, K = -0.5 ln(1-2P-Q) - 0.25 ln(1-2Q),
  with P (transitions) and Q (transversions) taken after the same CpG rule;
  NA when the logarithm saturates.
- n_aligned and n_cpg_sites: how much of the copy was usable, and how much
  of that is CpG.

Per plate, summary.tsv records: species, family, n_copies, the consensus CpG
share (CpG sites / element columns of the consensus), the median raw_p, the
median cpg_adj_p, the Spearman rank correlation between raw and adjusted
values across the copies, and the share of copies that move to a different
bin when the copies are split into 5 equal-size divergence bins by rank.

Two questions decide the outcome:

1. Size: how far does the median divergence drop per family, and does the
   drop track the family's consensus CpG share?
2. Order: does any copy change its divergence bin? If nothing moves, every
   profile keeps its shape and only the axis would compress.

## Results

PENDING — paste the rows of tools/cpg_div/summary.tsv into the table, then
write the three sections below from the table alone.

| species | family | n_copies | consensus_cpg_share | median_raw_p | median_cpg_adj_p | spearman_raw_vs_adj | bin_change_share |
|---|---|---|---|---|---|---|---|
| PENDING | PENDING | PENDING | PENDING | PENDING | PENDING | PENDING | PENDING |

### How large is the CpG effect?

PENDING — per species and family, in plain words: median raw_p -> median
cpg_adj_p, as a percentage of the raw value, next to that family's consensus
CpG share.

### Which families are affected most?

PENDING — rank the families by the drop in the median and by the bin-change
share; name the top ones and state what they have in common (expected: the
highest consensus CpG shares).

### Should the divergence profiles be CpG-adjusted?

PENDING — the answer follows the rule:

- If the bin-change share is 0 (or ~0) for every family: no copy changes its
  bin, so every profile keeps its shape; the adjustment would only compress
  the divergence axis. Then either keep the profiles raw and say so in the
  Methods, or switch profile and axis to cpg_adj_p — a presentation choice,
  not a change in results.
- If some families have moved copies: those profiles change shape; name
  them. For them the answer is yes, at least for the published profile.

## What this test does NOT show

- No independent age estimate. raw_p and cpg_adj_p are both relative measures
  inside a family. The test says which copies and families look less diverged
  once CpG deamination is discounted; it does not date any insertion and does
  not calibrate divergence to absolute time.
- The 1/10 weight is an assumption, not a measurement. It is RepeatMasker's
  convention; no methylation data was used, and the true C->T rate at
  methylated CpGs varies with lineage and time.
- Only the top-100 plates. Each plate carries the copies SINEderella publishes
  for that family, not the whole family; unassigned and filtered copies are
  absent, so a family with a biased plate gets a biased median.
- Only assigned, known families. Families the pipeline is unsure about are
  not on plates, so the survey covers only what SINEderella already trusts.
- k2p_adj is reported, not validated. It is the standard Kimura correction
  applied after the CpG rule; it was not compared with any other substitution
  model.
- Nothing about the pipeline. No SINEderella output changed for this test;
  whether the profiles change is a separate decision (above).
- One genome per species. Assembly quality and copy-number differences
  between species are not controlled for.
