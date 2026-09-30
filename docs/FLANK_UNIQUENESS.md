# Flank uniqueness over ALL copies of a family (design, 2026-09-30)

Status: **design, not implemented.** Agreed with Toki 2026-09-30: document first, then build a benchmark,
then implement and compare several approaches on real data before choosing one.

## Why the plates are not enough (his caveats, 2026-09-30)

Independent insertions have unrelated flanks; copies with shared flanks were multiplied some other way.
Today flank uniqueness is judged only on the plates (`report_flank_uniqueness.py`, the verdict's flank
check, `[array]` marks from `tools/array_order.py`). That has two blind spots:

1. **top100 selects for the problem.** The highest-scoring copies are often the youngest-looking ones, and
   "zombie" copies multiplied WITH their flanks - segmental duplications, satellitization (tandem
   arrays of SINE-containing units), carriage by another mobile element - look exactly like that. The
   top100 can be full of non-independent copies.
2. **rand100 cannot prove absence.** 100 random copies out of thousands are an insignificant sample: all
   100 flanks unique does not mean no group of copies in the family shares its flanks.

So flank uniqueness must be scanned on the level of **all genomic copies** of a family.

Existing partial measures and their limits:
- `[array]` (step8a, `tools/array_order.py`): >= 3 copies on one contig within 50 kb - proximity only,
  same contig only; pairs and cross-contig duplications escape it (PLATES.md: ttr MEG-RS pairs 7-11 kb
  apart; nle MEG-TR 8 copies sharing long flanks on different contigs).
- flankscan (SINEderella/flankscan): compares flanks with the consensus BANK only - a duplicated stretch
  of unique flank is invisible to it (GLM manuscript review round 2, task 08 #6; flankscan HANDOFF).
- step7 boundary walk: "undetermined" when copies never reach background - catches extension, not
  duplication.

## The central difficulty: repeats inside flanks

Flanks contain fragments of other transposable elements (LINEs, other SINE families, simple repeats).
An all-copies comparison will "find" thousands of flanks that share, e.g., a LINE fragment - that is not
the duplication we look for. Every approach must ignore genome-repetitive sequence: the assembly's soft
mask where present, TRF / dust (flankscan stage 2 already masks flank tandem repeats), and/or k-mers that
are frequent genome-wide.

## What counts as a shared flank (definition, 2026-09-30)

A copy duplicated together with its neighbourhood (segmental duplication, array unit, carried by another
element) has flanks identical to its twin's **from the junction with the SINE outward, colinear**. Two flanks
that merely contain the same LINE fragment match somewhere inside the flank, at unrelated offsets. So:

**shared flank = similarity that starts within a few bp of the junction in BOTH copies, runs outward in the
same order, over >= ~50 bp** (thresholds to be set on the benchmark). This separates duplication from the
repeat problem by definition, and it points to a cheap test: only the first ~100 bp next to each junction
need comparing (a junction-anchored k-mer or hash of the proximal flank).

## Ground truth: exhaustive alignment on a subset, not the `[array]` mark

`[array]` means >= 3 copies within 50 kb on one contig - proximity, not shared flanks - so it cannot be the
truth. Truth comes from exhaustive pairwise flank alignment (ssearch36) on a manageable subset (~2 000
copies: the array-marked and satellite-classed copies of the positive cases below plus random single copies
of clean monomers), scored with the definition above. Exact and slow, which is fine on a subset; the fast
approaches are then scored against it, and run on the full families for time and memory.

Marked copies available in the published plates (Tal `<sp>/alignments/<sp>_MEG-*_top100.aln.fa`, rows
with `[array]`; coordinates in the row names): rsi MEG-TR 73 of 102, tbr MEG-RS 17, ttr MEG-RS 11, nle
MEG-RS 7, nle MEG-TR 6, rsi MEG-RS 1; flankscan stage 2 satellite class: tbr MEG-RS 137 copies.

## Candidate approaches (CD-HIT excluded: global-identity greedy clustering, wrong question)

| | approach | idea | expected strength | expected weakness |
|---|---|---|---|---|
| A | low-copy k-mer index | index k-mers (k ~ 15-21) of all flanks; keep only k-mers rare in the whole genome (count with jellyfish/KMC); candidate pairs = copies sharing >= m rare k-mers; confirm by alignment | linear; the rare-k-mer rule solves the repeat problem directly | needs a genome k-mer count |
| B | map flanks back to the genome | minimap2 each copy's two flanks (~200 bp each, outside element + own tail); count strong secondary hits | answers "is this neighbourhood duplicated" directly, even when the duplicate has no SINE copy | repeat fragments multi-map - mask first |
| C | all-vs-all of flanks | minimap2 all-vs-all (or vsearch --usearch_global) among the flanks of one family | partial overlaps found | pair explosion on 600 k copies; same repeat problem |
| D | vsearch clustering | cluster flanks by identity | simple | global identity misses partial overlaps; slow at 600 k |

## Ground truth (real data, plus toy)

Positives (copies known to share flanks):
- tbr MEG-RS: tandem array, ~910 bp period (flankscan stage 2: 55 % of MEG-RS copies satellite units).
- rsi MEG-RS / MEG-TR / MEG-RL: one GC-rich array unit hit by three queries (PLATES.md).
- nle MEG-TR: 8 copies sharing long flanks across different contigs (PLATES.md).
- ttr MEG-RS: two tandem pairs 7-11 kb apart (PLATES.md).
- Toy genome: planted segmental duplications carrying copies (same contig and cross-contig), planted
  arrays - to be added to `flankscan/tests/make_toy.sh`.
Negatives (flanks should be unique): r9 single copies (rsi, clean monomer), other clean monomers; the toy's
single copies.

## Metrics per approach

- recall on the positives (per case), false-positive rate on the negatives;
- with and without repeat masking (to show the repeat problem is solved, not hidden);
- run time and memory on the worst case: tbr VES, 621 128 copies.

## Output (fixed, the same for every approach)

`flank_groups.tsv`: `copy side shares_with identity length` (one row per pair, side = 5 or 3), plus per
family: share of copies with a shared flank, and groups (connected copies). Plate rows get `[segdup]`
(flankscan HANDOFF item "Shared flank across contigs"; BORROWED_TRIAGE item 12).

## Plan and division of work

1. **Benchmark first (Claude):** ground-truth list with copy coordinates, toy segdups in make_toy.sh, a
   scoring script that reads `flank_groups.tsv` and prints recall / false positives / time.
2. **Implement A, B, C (GLM, via aider):** one narrow task each - a skeleton with fixed inputs, outputs and
   rules, gated by the benchmark check (the setting where GLM delivers; open design questions make it
   loop). Claude verifies every result independently.
3. **Run all on the real cases (Claude),** report numbers with the alignments of disagreeing copies.
4. **He chooses** the approach (and what counts as "duplicated" on the edge cases).

Related: `flankscan/fs8_singletons.sh` (stage 8) asks the same question for one family's single copies
with TSDs as end markers; its results are a first real test case.
