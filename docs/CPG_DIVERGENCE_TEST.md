# CpG-adjusted divergence: the cross-species test (Task C)

Related: docs/BORROWED_TRIAGE.md, "CpG-adjusted (Kimura) divergence" (section 5, item 2) and decision D2.
Tools: tools/cpg_div/cpg_div.py (one plate), tools/cpg_div/run_all.py (all species).
Data: tools/cpg_div/summary.tsv.

**STATUS: numbers filled 2026-10-01 from tools/cpg_div/summary.tsv.** To reproduce the table:

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

Run 2026-10-01 on 414 plates (all `C:/work/Tal/*/alignments/*top100.aln.fa`; table: `tools/cpg_div/summary.tsv`).
The aggregates below were computed from that table (Claude, independent of the GLM text); 371 plates have values,
39 give NA (species saq 9, ccr 8, teu 6, toc 6, dmo 5, gpy 5: not yet explained, probably plates without
upper-case element columns or with a different row layout; to be checked before any claim about those species).

### How large is the CpG effect?

Across the 371 plates the median CpG-adjusted divergence is a median of **10 %** lower than raw (10th to 90th
percentile of plates: 2 % to 38 %). Per species (median over its plates, drop of median_cpg_adj_p relative to
median_raw_p, and median consensus CpG share):

| species | plates | drop | consensus CpG share |
|---|---|---|---|
| hum | 64 | 35 % | 0.151 |
| rsi_v3 / rsi_v2 / rsi | 17 / 17 / 14 | 32 / 31 / 29 % | 0.10 |
| rle | 6 | 27 % | 0.080 |
| rmi | 6 | 15 % | 0.150 |
| zeb | 15 | 14 % | 0.108 |
| tim (Timema) | 14 | 13 % | 0.096 |
| tbr | 6 | 10 % | 0.097 |
| eri | 8 | 2 % | 0.026 |
| timb | 55 | 5 % | 0.067 |
| sco | 25 | 5 % | 0.075 |

The relation to the consensus CpG share is weak (Pearson r = 0.33 over the 371 plates): the drop depends on how
much of the divergence is transitions at CpG sites, which the CpG share of the consensus alone does not predict.

### Which families are affected most?

The mammalian plates (human, the rsi and Rhinolophus sets) are affected most, 27-35 %; invertebrate and low-GC
sets least (about 5 %). The table per family is in `summary.tsv` (sort by the two divergence columns).

### Should the divergence profiles be CpG-adjusted?

The rule written above gives a clear answer for the mammalian sets: the rank correlation raw vs adjusted is high
(median 0.96), but the median share of copies that change one of five equal-size divergence bins is **0.26**, so
the profiles change shape and not only the axis. For the mammalian families the answer is yes, at least as a second
track beside the raw profile; for the invertebrate plates the effect is small. This is a recommendation for the
published profile, not applied: D2 stays open until the user decides, and the 39 NA plates need explaining first.

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
