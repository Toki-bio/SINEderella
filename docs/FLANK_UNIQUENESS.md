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

---

## Implementation and test plan (Claude, 2026-10-01; GLM task U was stalling and is dropped in this form)

Decision: build it as a flankscan stage (bash/awk, readable, like stages 1-8), not three competing Python tools. The
benchmark below decides thresholds, and the approaches A/B stay as two modules of the same stage because they answer two
different questions.

### What exactly is asked (two classes, reported separately)

1. **Twin copies** (class T): copy A and copy B share flanks from the junction outward, colinear, >= 50 bp, identity >= 85 %.
   The two SINE copies are one event multiplied (segmental duplication, array unit, carried by another element). Mark `[twin]`.
2. **Copy inside a duplicated region** (class R): the flank of A matches somewhere else in the genome over >= 100 bp at >= 85 %, but
   the other place is not a copy of the same family at the same offset. The insertion may be independent; the region is
   not unique. Mark `[dup-region]`, do not count as non-independent.
The question "is this family made of independent insertions" is answered by class T only; class R is a warning.

### Traps the design has to survive (each is a benchmark case)

| trap | why it fools a naive test | rule that handles it |
|---|---|---|
| strand | an inverted duplicate looks unrelated | compare in element orientation: 5' flank with 5' flank, 3' with 3' (flank sequence as the copy's own orientation; handles both) |
| insertion hot spot in a repeat | independent insertions into the same LINE/simple repeat have identical proximal flanks | k-mers frequent genome-wide (> 20 copies) and soft-masked/TRF sequence are ignored; a pair needs >= 50 bp of UNmasked, low-copy flank |
| TSD / target-site preference | the first 5-15 bp next to the junction repeat by chance | window starts after the TSD; pair needs >= 50 bp, not 15 |
| indels | exact k-mers break at an indel | k = 12, key = (k-mer, offset from junction); pair = >= 8 shared keys on one diagonal band (band +-3 bp), then ungapped/gapped confirmation (ssearch36 on the pair only) |
| old duplicate | 80-90 % identity kills 20-mers | k = 12 and identity floor 85 %; below that = not called (report the floor) |
| divergent copies of the SINE itself | flanks fine, element differs | flanks only; element ignored |
| short flank | copy at a contig end has < 100 bp | pair allowed only if both flanks >= 50 bp; otherwise "untestable", counted separately (never counted unique) |
| N / low complexity | poly-N or (TA)n matches everything | masked before indexing |
| nested copies / adjacent copies | copy B lies inside A's flank (tandem or array) | exclude pairs closer than 2 kb on one contig from class T and report them as `[array]`, which exists already |
| huge families | 600 k copies, 1.2 M windows | index only the proximal 80 bp per side; bucket cap 20; sort-based, no all-vs-all |

### Pipeline (stage fs9, per family)

1. input: flankscan fs1 windows (strand-normalised, +-1000 bp) of the family's copies; keep proximal 80 bp each side, after the TSD.
2. mask: soft-mask/TRF (stage 2 masks already exist) plus genome-wide frequent k-mers (jellyfish count k = 12 on the genome, drop > 20).
3. keys: every 12-mer with its offset from the junction, `kmer TAB offset TAB copy TAB side`; `sort`; buckets > 20 dropped.
4. candidate pairs: from each bucket, all pairs; count shared keys per (pair, side, offset difference); keep >= 8 on a +-3 bp band.
5. confirm: ssearch36 of the two proximal 150 bp flanks; identity >= 85 % over >= 50 bp, starting <= 15 bp from the junction in both.
6. class R: minimap2 (or ssearch36) of each flank against the genome; second hit >= 100 bp, >= 85 %, not within 2 kb, not a copy of the family at the same offset.
7. output: `flank_groups.tsv` (copy_a copy_b side identity length class), union-find groups, per family: share of copies with a twin, number of
   groups, number untestable; plate rows get `[twin]` / `[dup-region]`; report line in the verdict ("x % of copies have a twin").

### Test plan

1. **Toy (planted truth), extend `tests/make_toy.sh`:** (a) 400 independent copies, unique flanks (negatives); (b) 6 segdups, 3 kb, copies inside:
   same contig, cross contig, inverted, 95 / 90 / 85 / 80 % identity; (c) tandem array of 10 units; (d) 60 insertions into one LINE-like repeat at the
   same offset (hot spot: must NOT be twins; must be suppressed by the frequency cap) and 60 at random offsets of it; (e) copies at contig ends
   (untestable); (f) copy inside a duplicated region where the SINE was inserted after duplication (class R, not T); (g) poly-N and (TA)n flank.
   Pass = recall 100 % at >= 90 % identity, every hot-spot copy unflagged, class R flagged R and not T, untestable counted, no unique copy flagged.
2. **Threshold calibration:** sweep k, band, min keys, identity floor on the toy; fix them where recall stays >= 95 % at 85 % identity and false pairs = 0;
   shuffled-flank control as in stage 8 (shuffle flanks between copies, expect 0 pairs).
3. **Exact reference:** on a subset of ~2 000 real copies (the array-marked ones plus random singles) run all-pairs ssearch36 of the proximal flanks and score
   the fast stage against it (recall and false positives), with and without masking.
4. **Real positives:** rsi MEG-TR (73 of 102 marked [array]), tbr MEG-RS satellite units, nle MEG-TR (8 copies on different contigs), ttr MEG-RS pairs;
   must be found as twins (tandem ones as `[array]`).
5. **Real negatives:** rsi r9 singles and other clean monomers: expect a low twin share; every flagged pair is inspected as an alignment (no score reported
   without the alignment, see the memory rule).
6. **Scale:** time and memory on tbr VES (621 128 copies); target: minutes on a 64-core node, memory dominated by the sort.
7. **Sensitivity of the verdict:** report the twin share per family for rsi r1-r10 with the thresholds' range, not a single number.

### Open design question for Toki

Where to put the line between "independent insertion" and "twin": identity floor 85 % over 50 bp is a first guess. Old duplicates (> 15 % divergent flanks) are
missed by design; do you want a second tier (70-85 %, class "possible old duplicate") reported separately?
