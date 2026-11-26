# Reproducibility: the same input gives the same output (2026-11-21)

Branch `reproducible`, not merged. Two steps: fixed statistics with fixed cuts (commit 34be237), then a cut per consensus
from a shuffle with a fixed seed (owner's choice of 2026-11-21, section 4b), which brought the loci back to the old
behaviour. What still differs from the old runs is in section 6; what is left to decide in section 8.
**Step 2 was changed again on 2026-11-22 (section 10): the best-scoring consensus wins, a close second is flagged and
no longer costs the copy its assignment.** Where sections 4, 6, 8 and 9 speak of `WIN_MARGIN` and of copies "without a
clear winner", they describe the version before that; the numbers that hold now are in section 10.

## 1. The problem

Two native runs on identical inputs disagreed. Test inputs (therioserver, `/home/toki/tmp/sd_container/e2e_ab3b29e/`):
`mini` = 11.4 Mb of the *R. sinicus* genome (68 sequences) with a bank of 21 consensuses; `toy` = the planted toy of
`tests/make_toy_satellite.py`. Three native runs of the old code on `mini` (natA, natB, natC; made by the container work,
commit ab3b29e = main plus path defaults), compared with `tests/compare_runs.py`:

| | natA | natB | natC |
|---|---|---|---|
| loci | 1,041 | 1,046 | 1,042 |
| assigned / no_unanimous / rejected_low_bitscore | 232 / 805 / 4 | 226 / 815 / 5 | 241 / 797 / 4 |
| verdict SINE / grey / not_SINE | 190 / 25 / 826 | 174 / 35 / 837 | 198 / 26 / 818 |
| plates | 39 | 39 | 41 |

| pair | loci in both | subfamily differs | status differs | verdict differs | plates: identical files |
|---|---|---|---|---|---|
| natA vs natB | 1,040 | 25 | 21 | 30 | 20 of 39 |
| natA vs natC | 1,040 | 56 | 11 | 13 | 24 of 39 |
| natB vs natC | 1,040 | 51 | 20 | 29 | 20 of 39 |

Over the three runs 112 of 1,048 loci (10.7 %) did not have the same subfamily, status and verdict in all three.
The `toy` runs agreed in every verdict and assignment but not in a single score line, and in 3 of 4 plates.

## 2. The cause

**ssearch36 takes its statistics from random shuffles, seeded by the clock.**

- With `-z 11` (used by `sear`, step 2, step 3, step 4, `length_variants.py`, `satellite_trf_verify.py`) every library
  sequence is also shuffled and scored, and the E-values and bit scores are fitted to those shuffled scores. Without
  `-z`, ssearch36 does the same by itself whenever fewer than 500 library sequences were scored
  (`comp_lib9.c:1079-1096`; `composite_units.py`, `continuation_check.py`, `fs9_bodycheck.py`, `trna_heads.sh`).
- The seed is the microsecond field of the clock at program start (`src/mrandom.c`, `my_srand`: `gettimeofday`,
  `tv_usec % 65535`). Each worker thread seeds its own generator the same way (`work_thr2.c:148`). There is no option
  for the seed.
- The alignments themselves do not change. Measured on `mini`, same command twice (therioserver, ssearch36 36.3.8g):

| call | threads | hits run a / b | same alignment (query, subject, coordinates) | identity differs | bit score differs | largest difference |
|---|---|---|---|---|---|---|
| `sear` call: 1 consensus vs 5,760 windows | 1 | 2,330 / 2,329 | 2,329 | 0 | 2,329 | 2.0 bits |
| same | 16 | 2,329 / 2,329 | 2,329 | 0 | 2,324 | 0.9 bits |
| step-2 call: 21 consensuses vs 1,042 copies | 1 | 8,942 / 8,965 | 8,926 | 0 | 8,626 | 19.1 bits |
| same | 16 | 8,967 / 8,955 | 8,933 | 0 | 8,704 | 40.5 bits |
| either call with `-z 3` | 1 and 16 | equal | all | 0 | 0 | 0 (same lines; with several threads in another order) |

So the hits near an E-value cut come and go, and every decision made on a bit score can come out differently:
which consensus wins a copy (step 2 votes), the summed score and the 0.45 threshold, the score in the FASTA headers that
ranks the top 100 copies of a plate and of the border loop, the continuation targets.

**A fixed seed does not help.** Tried first, because it would have kept every score on its old scale: ssearch36 run
under a stopped clock (libfaketime, and a 5-line preloaded `gettimeofday`), one thread. Twelve identical calls of the
step-2 search gave 8 different fits (one value 5 times, seven others once); the same with ssearch36 36.3.8j2 built from
the current source. Something besides the seed differs from run to run inside ssearch36, so the shuffles cannot be made
repeatable without changing its source.

## 3. Every decision that used a bit score or an E-value

Line numbers are those of main (0ebbaff). "cross" = scores of different queries are compared, which is where the shuffle
noise changes the outcome most. Calls of scripts that a full run does not reach are listed at the end.

| where | quantity | what is decided |
|---|---|---|
| `sear:451` | E <= 2 (`-z 11`) | whether an alignment is a hit at all; hits become loci (`sear:460-465`: merge, identity >= 65, length >= 0.8 x query) and raw hits for the satellite stage |
| `step1_search_extract.sh:240` | sear's maximum bit score | score column of `all_hits.labeled.bed` |
| `tools/extend_ends.py:119-121` (since 2026-11-26 `tools/extend_ends.sh`, function `best1`, same rule) | that score column (cross) | group of a locus for the end extension, which changes locus coordinates |
| `step3_postprocess.sh:401-422` | that score column (cross) | `Soft_Subfamily` of an unassigned locus |
| `tools/composite_units.py:59,66,77` (step 1, step 2, publish) | E <= 1e-5 (shuffle statistics: bank < 500 sequences) | which bank entries are composites: the pieces searched in step 1, the unit/composite rule of step 2 |
| `tools/satellite_trf_verify.py:225,243` | E <= 1e-2 (`-z 11`), best alignment by bit score | a TRF unit is SINE-derived (kind-A satellite locus); its region is excluded |
| `step2_asSINEment.sh:259` (`:689` incremental) | E <= 2 in each of 10 cycles | a copy without a hit in a cycle cannot reach 10/10 |
| `step2_asSINEment.sh:49-53` | bit score | best alignment of a copy on a consensus, and its span (unit rule, partial copies) |
| `step2_asSINEment.sh:67-75` | bit score (cross) | winner and runner-up of a cycle = the vote. Equal scores: whichever the awk array listed first |
| `step2_asSINEment.sh:332` | votes == 10 | unanimity |
| `step2_asSINEment.sh:357-360,378` | 10-cycle sum | TopN, threshold 0.45 x TopN, assigned / rejected_low_bitscore |
| `step2_asSINEment.sh:384-385`, `step3_postprocess.sh:302-308` | runner sum / winner sum >= 0.90 | LEAK flag |
| `step2_asSINEment.sh:411,461` | 10-cycle sum | the score in the headers of `assigned.fasta` and `subfamilies/*.fasta` |
| `step3_postprocess.sh:82,136` | E <= 100 (`-z 11`); the score read is the raw one | which copies get a similarity ratio |
| `step4_plots.sh:139,144-149` | E <= 100; best alignment by bit score | % identity histogram (display) |
| `step7_boundary_refine.sh:158,328` | header score | the top-100 population of the boundary test |
| `publish/flank_border_iterate.py:157-183` | header score | the 100 copies the border loop extends: rebuilt consensus, borders, flank verdict, `loci_verdict.tsv` |
| `step8a_extract_alignments.sh:447-459`, `tools/array_order.py:136-180` | header score | top-100 plate membership and order |
| `tools/continuation_check.py:113,135,212,319` | E <= 1e-5 (shuffle statistics), bit score (cross at 135) | own match, continuation target, flag, bank match of the candidate |
| `tools/length_variants.py:401,416` | E <= 1e-3 (`-z 11`), best alignment per copy | which copies enter the length-version test |
| `flankscan/fs9_bodycheck.py:101,109` | E <= 100 (shuffle statistics with 10 or more partners) | whether a pair has an identity; twin status of copies |
| `tools/trna_heads.sh:45,63`, `tools/trna_heads_summary.py`, `report_trna.py:29` | bit score against the null of shuffled queries | tRNA head call (only with `TRNA_REF`) |
| `tools/rollup_families.py:105-117,301-331` | votes, sums, runner ratio | family roll-up (only with `--families`) |

Already without shuffles: `tools/kindB_flank_check.sh` (`-z -1`), `step6_report.py:1105` (`-z 3`), `tools/trna_plates.py`
(one library sequence). Not reached by a full run and not changed: `sear_multi`, `asSINEment`, `step4_diagnostic.py`,
`tools/composite_scan.py`, `flankscan/fs3_partners.sh`, `fs8_singletons.sh`, `fs_lib.sh`, the serial path of
`satellite_kindB_verify.py`, `tools/vote_eval/`. They still use `-z 11`.

Other things that could differ between runs, found by reading the code and not fixed here: `sine_consensus.sh:121`
(`shuf` without a seed; the border loop uses it only when `conse` returns nothing), and `mafft --thread N` without
`--threadit 0` in `step4_plots.sh:178` and `tools/add_seed_row.py:130` (no iterative refinement, so probably the same
every time; not tested). After the fix two whole runs are identical on both test inputs, so none of these differed
there.

## 4. The fix

**Fixed statistics (`ssearch36 -z 3`) in every call of a full run, and in step 2 a clear-winner rule instead of ten
noisy votes.**

With `-z 3` ssearch36 uses fixed parameters (Lambda 0.318, K 0.130 for the +5/-4 matrix) and shuffles nothing. The bit
score is then a fixed linear function of the Smith-Waterman score, `bits = (0.318 x score + 2.04) / ln 2`: it depends on
the alignment only, not on the library, the run or the number of threads. With several threads ssearch36 writes the
same lines in another order, so the places that took "the first of equal scores" now sort the lines or name the tie-break.

Two things do not carry over and had to be set again:

**E-value cuts.** The fixed parameters give far smaller E-values for the same alignment (`-z 3 -E 2` reports 80 times
more hits than `-z 11 -E 2` on one window part), so each cut was moved to where the old cut lay. Measured on `mini`: the
old call three times, the `-z 3` call once, and for every alignment whether the old call reported it 3, 1-2 or 0 times.

| call | old | new | measured on `mini` |
|---|---|---|---|
| `sear` | `-z 11 -E 2` | `-z 3`, the consensus's own cut (section 4b); first version `-E 1e-10` | fixed cut, one part: of 865 hits always reported 855 kept, 364 never reported before come in; own cut: 864 kept, 59 come in |
| step 2 | `-z 11 -E 2`, 10 cycles | `-z 3`, the consensus's own cut, once; first version `-E 1e-7` | the weak alignments in question are runner-ups far below the winner |
| `composite_units.py` | default, `-E 1e-5` | `-z 3 -E 1e-17` | weakest always reported 5.8e-18, strongest never reported 6.4e-18 |
| `continuation_check.py` | default, `-E 1e-5` | `-z 3 -E 1e-16` | tails vs bank: 2,588 of 2,602 kept, 8 come in |
| `length_variants.py` | `-z 11 -E 1e-3` | `-z 3 -E 1e-17` | 229 of 229 kept, 1 comes in |
| `satellite_trf_verify.py` | `-z 11 -E 1e-2` | `-z 3 -E 1e-7` | one unit in the test run: not measured, set by analogy |
| step 3, step 4, `fs9_bodycheck.py`, `trna_heads.sh` | `-E 100` / `-E 100` / `-E 100` / `-E 10` | the same numbers with `-z 3` | these cuts are meant to pass everything |

The old cut was not one line on the fixed scale: under shuffles each consensus had its own threshold (bit score on the
fixed scale 57-70 for most consensuses, 110-125 for `r10_r8_P18`, above 73 for `MEG-T2`), because a consensus with
GC-rich or repeated stretches scores high against shuffled sequence too. A fixed cut cannot follow that, and the first version showed what it costs: `MEG-T2`, which had no hit under shuffles,
got 31 weak hits and 759 raw hits, the satellite stage flagged their arrays (satellite loci 25 -> 151) and excluded 336
intervals instead of 18, and 304 loci that were stable in the old runs disappeared.

### 4b. A cut per consensus, from a shuffle with a fixed seed

`tools/shuffle_cut.py` does outside ssearch36 what `-z 11` did inside it, the same way in every run:

1. **Shuffled library, once.** Up to 2,000 sequences (`SHUFFLE_N`) of the first part of the library, each with its
   residues shuffled. The generator is written out in the script (a 64-bit linear congruential generator and
   Fisher-Yates), so the bytes do not depend on the python version; seed 1. In `sear` the library is the first of the 20
   window parts (the parts are dealt round-robin, so it is a sample of the whole genome) and the file is kept with the
   parts for all queries; in step 2 it is the first part of the copies.
   The shuffle keeps the composition of each sequence, not its dinucleotides: that is what ssearch36 did, and the aim
   is the old behaviour. A dinucleotide shuffle would raise the cuts of CpG-poor or repeat-rich consensuses further; not
   tried.
2. **The same search** (`-z 3`, same matrix and gaps) of the consensus against the shuffled file.
3. **The cut.** The best bit score of each shuffled sequence follows an extreme-value distribution; its two parameters
   are taken from the 80 % and 95 % points of the scores, and the cut is the bit score that 2 shuffled sequences are
   expected to reach in a library of the size of one part: `cut = mu + beta * ln(N / 2)`. That is the meaning of
   "E <= 2" in the old call, and like it the cut rises with the size of the library.
   Only the top fifth of the scores is read because ssearch36 leaves out the lowest-scoring shuffled sequences (59-74 %
   were reported) and, with several threads, not always the same ones: a first version that used all reported scores
   wrote a step-2 cut file that differed between the 16- and the 10-thread run in one line (the assignments did not).
4. **Where it is applied.** `sear`: every alignment of the query must reach the cut (ssearch36 runs with `-E 1e-6`, which
   only limits its output). Step 2: every alignment of a consensus must reach that consensus's cut (`-E 1e-3`); the cuts
   are written to `consensus_cuts.tsv` in the step-2 output. A consensus with fewer than 30 scored shuffled sequences
   (a tiny library, as in the toy) keeps the fixed cut (`1e-10`, `1e-7`).
   Not applied in `composite_units.py`, `continuation_check.py`, `length_variants.py`, `satellite_trf_verify.py`: their
   libraries are a bank or a few hundred copies, their fixed cuts are in the table above.

Cuts on `mini` (288 windows per part), bit score on the fixed scale, with the old boundary measured in section 4:

| consensus | cut | under shuffles: weakest hit always reported / strongest never reported |
|---|---|---|
| MEG-RL | 53.1 | 56.8 / 50.9 |
| MEG-RS | 52.4 | 62.9 / 47.3 |
| MEG-TR | 56.6 | 60.5 / 55.0 |
| r1_r2_r3_r3_C11 | 67.5 | 67.3 / 63.6 |
| r5_r6_P26 | 59.9 | 64.1 / 63.2 |
| r10_r8_P18 | 93.3 | 110.4 / 124.6 |
| MEG-T2 | 76.7 | no hit / 78.4 |

Applied to the hits of all 21 consensuses on one part: 864 of the 865 always reported under shuffles are kept and 59
never reported come in (fixed cut: 855 and 364). On the whole window bank with the cut for that size: 16,610 of 16,619
kept, 631 come in (fixed cut: 16,591 and 4,525). So the cuts are a little lower than the old ones.

Cost: one more search per query against at most 2,000 windows (6 Mb), and one of the bank against at most 2,000 copies.
On `mini` it is not visible in the run time (section 6); on a large genome it is 6 Mb against 150 Mb per part times 20
parts: 0.2 % of the sequence searched (well under 1 % of the time; not measured on a large genome).

**The vote of step 2.** With fixed scores ten cycles are ten identical searches. What the ten cycles tested is asked
directly: the best consensus must be ahead of the runner-up by at least `WIN_MARGIN` = 0.05 of its bit score; equal
scores are broken by name. The margin is not asked where the unit/composite rule chose between the two (that rule gave
its winner all ten votes before). Calibration on `mini` (1,042 copies; ten old cycles against one fixed search):

| | n | margin (best - runner) / best |
|---|---|---|
| copies with 10/10 votes | 216 | median 0.168 (5 % at 0.021, 25 % at 0.112) |
| copies without | 826 | median 0.0217, half of them between 0.0208 and 0.0227, 95 % below 0.027 |

| margin | clear winners | agree with 10/10 |
|---|---|---|
| 0.02 | 963 | 26 % |
| 0.03 | 239 | 94.9 % |
| 0.04 | 223 | 95.7 % |
| **0.05** | 211 | **96.1 %** (193 both, 23 only 10/10, 18 only clear) |
| 0.07 | 188 | 96.2 % |
| 0.10 | 170 | 95.0 % |

The output files keep their layout. A clear winner is written as 10 votes (`Votes` 10, `Bitscore` = 10 x bit score,
so about 2.1 times the old sums); a copy without one as 5 votes for each of the two (`Votes` 5, status `no_unanimous`,
as before). The 0.45 x TopN threshold and the LEAK ratio are ratios of scores and are unchanged.

Not recalibrated, because they are ratios or rankings: 0.45 x TopN, LEAK >= 0.90, the top-100 rankings.

Changed files: `sear`, `step2_asSINEment.sh`, `step3_postprocess.sh`, `step4_plots.sh`, `tools/composite_units.py`,
`tools/continuation_check.py`, `tools/length_variants.py`, `tools/satellite_trf_verify.py`, `tools/trna_heads.sh`,
`flankscan/fs9_bodycheck.py`.

**Long paths (smaller bug, same branch).** ssearch36 cannot open a file whose path is longer than about 120 characters,
and `tools/continuation_check.py` passed full paths. It now copies the two files into one directory and gives bare
names, as the satellite tools do. Test: `tests/test_reproducible.py` runs it under a 160-character path.

## 5. The test

- `tests/test_reproducible.py` (pytest; needs ssearch36, seqkit, bedtools): the shuffled library is the same bytes in
  every run and keeps each sequence's composition, and a consensus with a repeated GC-rich motif gets a higher cut
  than a plain one; the fixed call gives the same lines at 1, 2 and 4 threads; the control, the old `-z 11` call, differs between runs; step 2 run twice with 2 and 5 threads writes
  byte-identical `assignment_full.tsv`, `all_votes.tsv`, `assignment_stats.tsv`, `assigned.fasta`, `unassigned.fasta`;
  the continuation check gives the same table at 1, 2 and 4 threads and in a long path.
- `tests/test_reproducible_run.sh GENOME BANK [WORKDIR] [SPECIES]`: the whole workflow twice with publish, with different
  numbers of threads, then `tests/compare_runs.py`: loci, hits, the cuts of step 1 and step 2, assignments, verdicts,
  published consensuses and every plate must be identical. Minutes; run on a server.
- `tests/compare_runs.py RUN_A RUN_B [--loci OUT.tsv]`: exit 0 only when the two runs are identical in those files.

pytest on therioserver (conda env sinederella): 172 passed.

## 6. Verification

Two native runs of the final code on each input, with 16 and 10 threads (therioserver, `/home/toki/tmp/repro/e2e_z3g/`):

| input | loci | hits | cuts (step 1, step 2) | assignment | verdict | published consensuses | plates |
|---|---|---|---|---|---|---|---|
| toy | identical (225) | identical | identical | identical | identical | identical | 4 of 4 identical files |
| mini | identical (1,071) | identical | identical (38 queries, 21 consensuses) | identical | identical | identical | 43 of 43 identical files |

`mini` takes 408-429 s (old code: 470-484 s; fixed cuts without the shuffle: 408-423 s). `toy`: 147-149 s (old 161-165 s).

**Final run against the old runs.** `toy`: every one of the 225 loci has the same subfamily, status and verdict as in
the three old runs (217 SINE, 2 grey, 6 not_SINE).

`mini`:

| | old natA / natB / natC | fixed cuts (34be237) | cut per consensus (final) |
|---|---|---|---|
| loci | 1,041 / 1,046 / 1,042 | 770 | 1,071 |
| intervals removed as satellite regions | 18 / 18 / not looked up | 336 | 18 |
| assigned | 232 / 226 / 241 | 242 | 245 |
| no_unanimous | 805 / 815 / 797 | 505 | 816 |
| rejected_low_bitscore | 4 / 5 / 4 | 23 | 10 |
| verdict SINE | 190 / 174 / 198 | 178 | 179 |
| verdict grey | 25 / 35 / 26 | 45 | 43 |
| verdict not_SINE | 826 / 837 / 818 | 547 | 849 |
| plates | 39 / 39 / 41 | 37 | 43 |

The final run differs from natA, natB, natC at 824, 845 and 793 loci (subfamily, status or verdict, or present in one
run only); 73, 94 and 42 of them had flipped among the old runs. **751 loci were the same in all three old runs and
are different now; in 700 of them only the name of the best consensus differs, and the verdict differs in 45.**
`MEG-T2` has no hit again, and the 304 MEG-RL loci that the fixed cuts had lost are back.

| loci | old (all three runs) | final | cause |
|---|---|---|---|
| 700 | MEG-RL, no_unanimous, not_SINE | MEG-RS, no_unanimous, not_SINE | Only the name written for a copy without a winner. These copies score about 2 % higher on MEG-RS (margin 0.021-0.023); under shuffles each consensus had its own score scale and most of the 10 votes went to MEG-RL (7 or 8 of 10 for most copies). Status and verdict unchanged. |
| 3 | no_unanimous, not_SINE | the same, other name | 2 MEG-RL -> MEG-TR, 1 r5_27seqs -> r5_r5_P48: as above. |
| 30 | absent | 18 not_SINE (5 r10_r8_P18, 1 r1_r3_P2 and 1 r9_15seqs rejected_low_bitscore; 3 r1_9seqs assigned; 8 without a winner), 8 grey, 4 SINE (3 r1_9seqs, 1 r1_r2_r3_r3_C11) | New loci. The cuts are a little lower than the old ones (section 4b: 59 more hits per part), so some weak hits now reach the length rule or join a neighbour. Not traced one by one. |
| 1 | r10_r8_P18, no_unanimous, not_SINE | absent | Not traced. |
| 7 | assigned, SINE | grey (6 with the same subfamily, 1 r8_83seqs -> r10_19seqs) | Not traced. The same number changed with the fixed cuts, so it is probably not the cut; the verdict comes from the top-100 copies and the border loop, which rank by score. |
| 1 | assigned, grey | assigned, SINE | Not traced. |
| 3 | assigned, SINE | no_unanimous, not_SINE (2 same subfamily, 1 r8_83seqs -> r1_r3_P1) | The margin rule (winner ahead by less than 0.05), by its definition; not looked at one by one. |
| 5 | no_unanimous, not_SINE | assigned (2 SINE, 1 grey, 2 not_SINE) | The margin rule (winner ahead by 0.05 or more although the 10 votes had split); not looked at one by one. |
| 1 | rejected_low_bitscore | no_unanimous | Not traced. |

Among the loci stable in the old runs, 10 lose the verdict SINE (7 to grey, 3 to not_SINE) and 7 gain it (4 new loci,
2 formerly without a winner, 1 from grey). The SINE total, 179, lies inside the old range (174-198). For comparison,
two old runs differed from each other at 62-90 loci, 13-37 of them in the verdict. The rows marked "not traced" are not
explained.

## 7. What is not restored

Under shuffles each consensus also had its own score scale, not only its own cut, and the winner of a copy was chosen
on those scales. Now the winner is the consensus with the highest Smith-Waterman score. For copies with a clear winner
this is the same consensus in all but a handful of loci (table above). For copies without one it changes the name that
is written (700 loci, MEG-RL -> MEG-RS), not their status or verdict. A score scaled by the consensus's shuffle fit
would bring the old names back; it is not written, because nothing reads the name of a copy without a winner as an
assignment.

## 8. For the owner to decide

1. **`WIN_MARGIN` 0.05.** It reproduces the 10/10 vote on 96 % of the `mini` copies. 0.03-0.07 do about the same; the
   group that the vote left unassigned sits at 0.021-0.023.
2. **The cuts are a little lower than the old ones** (30 new loci on `mini`, 4 of them SINE). They can be raised by
   asking for fewer expected shuffled hits (`--expect` in `tools/shuffle_cut.py`, now 2 as the old `-E 2`), or left.
3. **The name written for copies without a winner** (section 7): leave, or scale the scores per consensus.
4. **The status name `no_unanimous` and `Votes` 10 / 5** are kept so that every reader of `assignment_full.tsv` works
   unchanged. They now mean "no clear winner" and "clear winner / none".

Also: the scripts outside the full run (section 3, last paragraph) still use `-z 11`; 40 loci of section 6 (30 of them the new loci)
are not traced to a cause one by one.

## 9. Full genome: *R. sinicus* (2026-11-21)

Branch at d86b145, one run, 32 threads, no publish, the inputs and options of the latest complete run of main
(`/home/toki/rhin/rsi_final_20261102/run_20261101_235602`, main f2457a8: `rsi.fna`, `bank21_nor1.fa`, `--families
rsi_families_nor1.tsv`). New run: therioserver `/home/toki/tmp/repro/rsi_full/A/run_20261121_211005`. Only one old run
exists for this input, so how many loci would have differed between two old runs is not known here. A second copy for
byte identity was not run: the server was at a load of 80 of 128 before the start.

**Time.** 1 h 16 min 23 s, old 1 h 31 min 48 s. Step-1 searches 23 min 41 s (old 24 min 26 s). The shuffle searches
are a small part of that: 2 s for the shuffle of 2,000 windows and 15 s for 21 bank queries (measured apart, 32 threads),
so about half a minute for the 39 queries. Step 2: 1.6 min (old 9.3 min), of which shuffle and its search 4 s.

**Totals.**

| | old | new |
|---|---|---|
| loci | 68,580 | 68,191 |
| same locus name in both | 66,818 | |
| loci with no overlapping locus in the other run | 625 (lost) | 236 (new) |
| intervals removed as satellite regions | 196 | 181 |
| assigned | 64,697 | 61,081 |
| no_unanimous | 1,123 | 4,877 |
| rejected_low_bitscore | 2,760 | 2,233 |
| verdict SINE | 50,790 | 48,512 |
| verdict grey | 12,404 | 11,845 |
| verdict not_SINE | 5,386 | 7,834 |

8,065 loci differ in subfamily, status or verdict or are in one run only; the verdict differs in 6,173.

**Assigned copies and SINE verdicts per subfamily (old -> new).**

| subfamily | assigned | SINE | subfamily | assigned | SINE |
|---|---|---|---|---|---|
| r1_r3_P1 | 14,464 -> 13,830 | 11,183 -> 10,704 | r5_r3_P34 | 1,576 -> 1,531 | 1,207 -> 1,178 |
| r6_210seqs | 7,442 -> 7,388 | 5,934 -> 5,897 | MEG-RS | 729 -> 16 | 0 -> 13 |
| r1_r3_P2 | 6,940 -> 6,059 | 5,596 -> 4,899 | r10_r8_P18del | 505 -> 722 | 419 -> 590 |
| r9_15seqs | 6,937 -> 6,758 | 5,617 -> 5,467 | r5_r5_P48 | 464 -> 361 | 369 -> 290 |
| r7_133seqs | 6,862 -> 6,671 | 5,720 -> 5,556 | r10_19seqs | 351 -> 91 | 270 -> 65 |
| r5_r6_P26 | 5,677 -> 5,527 | 4,372 -> 4,258 | MEG-T2 | 103 -> 12 | 94 -> 10 |
| r10_r8_P18 | 4,308 -> 4,273 | 3,404 -> 3,378 | r3_58seqs | 40 -> 14 | 30 -> 9 |
| r8_83seqs | 3,385 -> 3,060 | 2,759 -> 2,495 | r2_3seqs | 29 -> 17 | 24 -> 16 |
| r1_r2_r3_r3_C11 | 2,866 -> 2,823 | 2,169 -> 2,136 | MEG-TR | 18 -> 11 | 16 -> 9 |
| r5_27seqs | 1,981 -> 1,897 | 1,589 -> 1,526 | r4_32seqs, MEG-RL | 12 -> 16, 8 -> 4 | 11 -> 16, 7 -> 0 |

**Causes, largest first.**

1. **The margin rule is stricter than the 10 votes were on a full genome: 3,415 copies.** Of the 63,861 copies that were
   assigned before and exist in both runs, 3,415 have no clear winner at `WIN_MARGIN` 0.05 (the 3,038 loci "assigned ->
   no_unanimous, same subfamily" of the comparison are among them: 1,835 SINE, 447 grey, 756 not_SINE before; all
   not_SINE now). Their margins: median 0.023, 10 % below 0.009, 90 % below 0.045; their old runner-up ratio had a median
   of 0.96. By subfamily: r1_r3_P2 748 (and 96 more now best on r1_r3_P1), MEG-RS 726, r1_r3_P1 625, r8_83seqs 190,
   r7_133seqs 175, r9_15seqs 154, r5_r5_P48 135, r5_r6_P26 131.
   The reason is the size of the library. The 10 votes failed when the shuffle noise exceeded the margin, and the noise
   of a fit to 20,000 copies is far smaller than of a fit to the 1,042 copies of `mini`, where the margin was set.
   On this genome no margin separates the two old groups (copies in both runs):

   | margin | old 10/10 and clear now | old 10/10, not clear now | old no_unanimous, clear now | neither | agreement |
   |---|---|---|---|---|---|
   | any difference (> 0) | 65,685 | 101 | 1,002 | 30 | 98.4 % |
   | 0.005 | 65,498 | 288 | 873 | 159 | 98.3 % |
   | 0.01 | 65,284 | 502 | 745 | 287 | 98.1 % |
   | 0.02 | 64,660 | 1,126 | 524 | 508 | 97.5 % |
   | 0.03 | 63,387 | 2,399 | 392 | 640 | 95.8 % |
   | 0.05 (now) | 62,021 | 3,765 | 266 | 766 | 94.0 % |

   Where both rules name a winner it is the same consensus in 60,050 of 60,446 copies.
2. **Winner by raw score: a few hundred copies change subfamily.** r10_19seqs -> r10_r8_P18del 143 (SINE in both),
   r10_r8_P18del -> r10_r8_P18 87 (assigned SINE -> rejected_low_bitscore), r1_r3_P2 -> r1_r3_P1 75, r8_83seqs ->
   r5_r5_P48 50. Not traced further.
3. **Loci lost: 625 with no overlapping locus in the new run.** 373 r10_r8_P18 and 25 r10_r8_P18del, all
   rejected_low_bitscore; 91 MEG-T2 (84 SINE, 4 grey, 3 not_SINE); about 70 rejected_low_bitscore loci and about 30 loci without a winner of other
   consensuses; about 30 assigned loci (8 r2_3seqs, 7 MEG-TR, 6 r1_r3_P1, 3 r1_r2_r3_r3_C11, 2 r8_83seqs, all SINE; 3 grey).
   18 of the 625 lie in satellite regions of the new run. The cause for the first two groups is the cut: `r10_r8_P18`
   has the highest cut of the bank (86.8 bits against 60-70 for most) and `MEG-T2` the third highest (75.8), and here they
   are higher than what the old shuffles asked for, the opposite of `mini`. `MEG-T2` had 115 hits and 1,074 raw hits
   in the old run and has 15 and 268 now.
4. **Loci new: 236.** about 120 rejected_low_bitscore (44 r10_r8_P18del), about 80 without a winner (46 MEG-RS), about 25 assigned
   (6 r4_32seqs, 4 MEG-RS, 4 r1_r3_P1, 2 r9_15seqs SINE; 7 grey). 64 lie in satellite regions of the old run.
5. **Other loci** (about 1,100 names in one run only that overlap a locus of the other run): the same copy with other
   ends, from a few more or fewer weak hits merged in. Hit counts per query agree within 0.1 % for most consensuses
   (r10_19seqs 65,812 -> 65,847, r1_r3_P1 22,476 -> 22,477); the short pieces gain more (`r1_r3_P1` piece 1-80:
   82,707 -> 86,943).

**Consensuses that score high against random sequence** exist on the full genome as on `mini`, and the cut finds the
same ones: `r10_r8_P18` 86.8, its piece 156-289 78.6, `MEG-T2` 75.8, `r1_r2_r3_r3_C11` 75.2, `r10_r8_P18del` 74.2,
`r5_r3_P34` 74.1, `r1_r3_P1` 73.5; the others 56-71 (the shortest piece, 50 bp, 48.4). Step 2 (2,000 shuffled copies): 50.7-69.9, the
same order.

**For the owner, from this run.**

1. `WIN_MARGIN`: 0.05 removes 5.3 % of the assigned copies of a full genome and nearly all of MEG-RS (729 -> 16). The 10
   votes behaved like a margin of about 0.02 on 1,042 copies and like "any difference" on 68,000. A smaller margin
   (0.01: 502 of the old 10/10 copies not clear, 745 formerly split copies clear) is closer to what full-genome runs did;
   no value reproduces both.
2. The cuts of `r10_r8_P18` and `MEG-T2` are higher than the old ones on a full genome (the fit through the 80 % and
   95 % points is extrapolated far into the tail when a part holds tens of thousands of windows). 84 MEG-T2 loci with the verdict SINE are gone. Whether these
   were real copies has to be looked at in their alignments.
3. Two runs of the full genome were not compared byte for byte.

## 10. Best score wins, close seconds flagged (2026-11-22/23; commits 6a90eb5, fbd15ff, f4e6811)

**Rule (owner, after section 9).** A copy goes to the consensus with the highest bit score; equal scores are broken by
name; the unit/composite rule is unchanged. There is no margin gate. When the second consensus is within `TIE_MARGIN`
(0.02) of the best and the unit/composite rule did not choose between the two, the copy is flagged: a row in
`ties.tsv` of the step-2 output (copy, best, second, margin, kind) and the tags `tie_with`, `tie_margin`, `tie_kind` in
the RunnerInfo column of `assignment_full.tsv`. `tie_kind` is `tie_family` when the two consensuses are in different
families of the `--families` table, otherwise `tie`. A flagged copy is assigned and gets its verdict like any other.
The 0.45 x TopN threshold is unchanged. Step 2 no longer writes the status `no_unanimous`; every copy with a hit has
Votes 10.

**What went wrong on the way, and the fixes.**

1. The flag was first written as an 8th column. `tools/rollup_families.py` unpacks 7 fields per row, so both full runs
   with `--families` stopped in the family roll-up. The small test runs had no `--families` and did not show it. The
   table is read by 20 scripts; it keeps 7 columns and the flag is in the tags, which are read by name.
   `tests/test_reproducible.py` now puts the step-2 rows through rollup_families' own unpacking and tag parser, and
   `tests/test_reproducible_run.sh` takes `FAMILIES=table` to run the workflow with `--families`.
2. With more copies assigned, the border loop reached `sine_consensus.sh`, whose subsampling used `shuf` without a seed:
   the rebuilt consensus of MEG-RS differed by one base between two `mini` runs. The seed is now the name of the
   subsample file.
3. Two full runs differed in three files of the family view (`sim_scores.tsv`: 58,336 and 58,225 rows). ssearch36 keeps
   a limited list of best library sequences and, on a larger library, drops low-scoring ones in an order that depends
   on its threads: one consensus against 68,345 copies gave 74,644, 73,620 and 75,923 alignment rows in three calls, in
   parts of 20,000 copies the same 242,211 rows at 32 and at 24 threads. Step 3 now searches in parts of 20,000 (step 2
   already did). `tests/compare_runs.py` also compares the step-3 tables and the family view.

**Plates.** Flagged copies are not kept out of the plates. A plate takes its copies from `assigned.publish.fasta`, whose
coordinates the border loop has changed, so a copy cannot be matched to `ties.tsv` by name there; and the old full-genome
runs had these copies in their plates too (the 10 votes passed nearly all of them, section 9). Keeping them out needs the
flag carried into that file; not done.

**Small inputs** (two runs each, 16 and 10 threads, `mini` with `--families`): identical in everything
`tests/compare_runs.py` reads. `toy`: the 225 loci as in the old runs. `mini`: 1,071 loci, 1,027 assigned, 44 rejected_low_bitscore, SINE 204, grey 40, not_SINE 827, 73 flagged copies, 41 plates (old runs: SINE 174-198, grey
25-35); the 700 MEG-RL/MEG-RS copies of section 6 are now assigned to MEG-RS with the verdict not_SINE.
pytest on therioserver: 172 passed (f4e6811).

**Full *R. sinicus* genome, f4e6811**, same inputs as section 9, two runs: 32 threads 1 h 10 min 37 s, 24 threads
1 h 15 min 14 s (old 1 h 31 min 48 s). Runs: therioserver `/home/toki/tmp/repro/rsi_full/B/run_20261123_072054` and
`/home/toki/tmp/repro/rsi_full/C/run_20261123_072054`.

*Byte identity.* `tests/compare_runs.py`: identical. All 2,038 files of each run directory, logs, manifest, report page
and timings left out: 50 differ. 42 are the PDF plots of step 4 (the PNG plots are identical). 6 hold the run directory
path, a date or a duration (`border_loop.json`, `consensus.path`, `run_info.txt`, `subfam_log.txt`, `summary.txt`,
`verdict_chain.done`). `vote_summary.tsv.gz` has the same content after decompression. One real difference is left:
`subfam_input/input.msf`, the alignment SubFam makes for itself (1,823 and 1,781 columns); the workflow does not read
it, and `input.clw` and the alignment built from it are identical. It comes from the vendored SubFam
(`mafft --thread` without `--threadit 0`) and was left alone.

*Totals.*

| | old (main f2457a8) | section 9 (margin 0.05) | now |
|---|---|---|---|
| loci | 68,580 | 68,191 | 68,191 |
| assigned | 64,697 | 61,081 | 64,900 |
| no_unanimous | 1,123 | 4,877 | 0 |
| rejected_low_bitscore | 2,760 | 2,233 | 3,291 |
| verdict SINE | 50,790 | 48,512 | 50,907 |
| verdict grey | 12,404 | 11,845 | 12,463 |
| verdict not_SINE | 5,386 | 7,834 | 4,821 |

*Assigned copies and SINE verdicts per subfamily (old -> now).*

| subfamily | assigned | SINE | subfamily | assigned | SINE |
|---|---|---|---|---|---|
| r1_r3_P1 | 14,464 -> 14,595 | 11,183 -> 11,287 | r5_r3_P34 | 1,576 -> 1,655 | 1,207 -> 1,265 |
| r6_210seqs | 7,442 -> 7,388 | 5,934 -> 5,897 | MEG-RS | 729 -> 770 | 0 -> 0 |
| r1_r3_P2 | 6,940 -> 6,813 | 5,596 -> 5,500 | r10_r8_P18del | 505 -> 722 | 419 -> 590 |
| r9_15seqs | 6,937 -> 6,915 | 5,617 -> 5,597 | r5_r5_P48 | 464 -> 641 | 369 -> 504 |
| r7_133seqs | 6,862 -> 6,883 | 5,720 -> 5,740 | r10_19seqs | 351 -> 92 | 270 -> 66 |
| r5_r6_P26 | 5,677 -> 5,730 | 4,372 -> 4,419 | MEG-T2 | 103 -> 12 | 94 -> 10 |
| r10_r8_P18 | 4,308 -> 4,274 | 3,404 -> 3,379 | r3_58seqs | 40 -> 14 | 30 -> 9 |
| r8_83seqs | 3,385 -> 3,282 | 2,759 -> 2,674 | r2_3seqs | 29 -> 17 | 24 -> 16 |
| r1_r2_r3_r3_C11 | 2,866 -> 3,021 | 2,169 -> 2,286 | MEG-TR | 18 -> 15 | 16 -> 12 |
| r5_27seqs | 1,981 -> 2,025 | 1,589 -> 1,626 | r4_32seqs, MEG-RL | 12 -> 17, 8 -> 19 | 11 -> 17, 7 -> 13 |

*Agreement with the old assignment* (66,818 copies with the same name in both runs):

| old | now | copies |
|---|---|---|
| assigned | assigned, same subfamily | 63,028 (98.7 % of the 63,861 old assigned) |
| assigned | assigned, other subfamily | 424 (r10_19seqs -> r10_r8_P18del 168, r8_83seqs -> r5_r5_P48 61, r1_r3_P2 -> r1_r3_P1 58) |
| assigned | rejected_low_bitscore | 409 (349 with another subfamily: 87 r10_r8_P18del -> r10_r8_P18) |
| no_unanimous | assigned | 656 (372 same subfamily, 284 other); verdicts SINE 516, grey 136, not_SINE 4 |
| no_unanimous | rejected_low_bitscore | 376 |
| rejected_low_bitscore | rejected_low_bitscore | 1,898 |
| rejected_low_bitscore | assigned | 27 |

*Close seconds.* 1,873 copies are flagged, all `tie`, none `tie_family`: no copy has its two best consensuses in
different families (MEG and Rhin). Of them 872 are assigned with the verdict SINE, 232 grey, 59 not_SINE, and 710 are
rejected_low_bitscore. By pair: r1_r3_P1 / r1_r3_P2 497, r1_r3_P1 / r1_r2_r3_r3_C11 433, r5_r5_P48 / r5_r6_P26 347,
r1_r3_P1 / r5_r3_P34 124, r8_83seqs / r5_27seqs 91; MEG-RS ahead of MEG-RL 43.
Margins of all 68,191 copies (from the runner-up score): 14.8 % below 0.005, 15.4 % below 0.01, 16.7 % below 0.02,
19.0 % below 0.03, 21.5 % below 0.05, 31.5 % below 0.10. So 11,408 copies have a second consensus within 0.02 and 1,873
are flagged; the others are pairs where the unit/composite rule chose (a unit and a composite that holds it align the
same way), which are left out of the flag on purpose. That reading of the 9,535 was not counted apart.

*What still differs from the old run.* 5,212 loci differ in subfamily, status or verdict or are in one run only; the
verdict differs in 4,219. The loci and their causes are those of section 9, points 2 to 5 (the searches are the same):
625 loci lost (398 rejected loci of r10_r8_P18 and r10_r8_P18del and 91 MEG-T2 loci, from their higher cuts), 236 new,
about 1,100 copies with other ends, a few hundred copies with another best consensus. Point 1 of section 9 is gone:
instead of 3,415 assigned copies lost, 656 copies whose 10 votes had split are assigned now (516 SINE).

**For the owner.**

1. `TIE_MARGIN` 0.02 flags 1,873 copies (2.7 %). It decides nothing; it can be moved freely.
2. Flagged copies are in the plates (above). Say if they should be kept out.
3. The cuts of `r10_r8_P18` and `MEG-T2` are unchanged, as agreed (84 MEG-T2 loci with the verdict SINE are still gone).
4. `input.msf` of SubFam differs between runs; nothing reads it. A fix belongs in SubFam.

## 11. The environment must not change the output either (2026-11-23)

Sections 1-10 made two runs on one machine identical. A native run and a run of the container image
(branch `container`) on the same input then still compared as different, and `tests/test_reproducible.py`
failed one test in both environments. Both were workflow bugs, fixed here; nothing was worked around in the
comparison.

### 11a. `sim_scores.tsv` was written in the host's collation

`tests/compare_runs.py` reported `sim_scores.tsv` of step 2 and of the family view as "same lines, other
order" (1,027 and 1,071 rows, none only in one run). The order of that file is the order of a sort of the
locus names, and the sort ran in the inherited locale: `en_US.UTF-8` on therioserver, `C.UTF-8` in the image.
glibc's `en_US.UTF-8` collation ignores `:`, `-` and `(`, the C collation does not, so the same rows came out
in two orders. Sorting the old native file under `LC_ALL=C` gave exactly the container's file (md5
`7f85a273…` both ways), which is what identified the cause.

The order is not cosmetic: `step6_report.stratified_sample_sim` reads this file in order and keeps a reservoir
sample of 3,000 copies per subfamily, so the divergence histogram of a subfamily with more copies than that
depends on it.

Fixed where the file is written: `step3_postprocess.sh` sorts `sim_scores.tsv` with `LC_ALL=C sort -k1,1`, and
every other sort in that script that writes a file (`self_bits.tsv`, `self_bits_real.tsv`, the soft-assignment
table, `summary.by_subfam.tsv`) and the cut-off table of `step1c_deplete.sh` run under `LC_ALL=C` too. Step 2
already exports `LC_ALL=C` for its whole run and `step8a_extract_alignments.sh` already sorted that way; step 3
was the gap.

### 11b. ssearch36 and paths longer than 120 characters, again

`tests/test_reproducible.py::test_step2_twice_identical` failed identically native and in the container with

```
***[ssearch36]  cannot open library …/test_step2_twice_identical0/.step2_assign_sOhi/consensx***
```

— the library name cut off mid-word at 120 characters. This is the limit already met in the Sicista runs of
2026-10-06 and fixed then in `tools/satellite_trf_verify.py`, `tools/satellite_kindB_verify.py` and later
`tools/continuation_check.py` (section 4, "Long paths"): ssearch36 truncates the file name, prints `***…***`
on stderr and **exits 0**. `step2_asSINEment.sh` makes its temp directory under `$PWD`
(`STEP2_TMP_BASE:-$PWD`) and handed ssearch36 the full path, so any run started from a deep directory died
there; under pytest it always did.

Every ssearch36 call of the repository was checked. Already safe: `sear` and step 1 (they `cd` into the search
directory and pass bare names), the flankscan stages (relative names), `tools/kindB_flank_check.sh`,
`tools/trna_heads.sh`, `tools/trna_plates.py` and the three tools fixed in 2026-10/11. Six calls still passed a
constructed path and were given bare names with the working directory set, the pattern of
`tools/continuation_check.py`:

| file | call |
|---|---|
| `step2_asSINEment.sh` | the shuffle search of `tools/shuffle_cut.py`, and the main search over `splits/part_*.fasta` |
| `step3_postprocess.sh` | the consensus self-alignment, and the per-copy search over `copies_parts/*` |
| `step4_plots.sh` | the %identity search per subfamily |
| `step6_report.py` | the SINEplot search (`cwd=work`) |
| `tools/length_variants.py` | the chunk search (`cwd=workdir`) |
| `flankscan/fs9_bodycheck.py` | the body-identity search per copy (`cwd=tmp`) |

Nothing about what is searched or scored changed. Not changed: the scripts a full run does not reach
(`asSINEment`, `sear_multi`, `tools/composite_scan.py`, `tools/vote_eval/`, `step4_diagnostic.py`), which also
still use `-z 11` (section 3).

`tests/test_reproducible.py::test_step2_in_a_long_path` runs step 2 from a directory whose path is over 150
characters and checks the assignment table and an empty `ssearch.log`. On the code before this section it fails
with the error above; after it, it passes. The stub in `tests/test_sineplot_sine_loci.py` now resolves the
ssearch36 arguments against the `cwd` the call sets, as ssearch36 does.

### 11c. Verification

therioserver, conda env `sinederella` and image `/home/toki/tmp/containers/sinederella_<container commit>.sif`
(image locale `LANG=C.UTF-8 LC_ALL=C.UTF-8`); runs under `/home/toki/tmp/sdc3/`.

- pytest: **173 passed** natively and in the image (172 before, plus the new test).
- `tests/test_reproducible_run.sh` (the whole workflow twice with publish, 8 and 12 threads) on `toy` and on
  `mini` with `--families`, in both environments: all four passed.
- `tests/compare_runs.py`, all four native/container pairs (A/A, B/B, A/B, B/A) on both inputs: exit 0.

| comparison | toy | mini |
|---|---|---|
| native A vs native B | identical | identical |
| container A vs container B | identical | identical |
| native vs container (4 pairs) | identical | identical |

`toy`: 225 loci, SINE 217 / grey 2 / not_SINE 6, 4 plates. `mini`: 1,071 loci, 1,027 in `sim_scores.tsv`,
73 flagged close seconds, SINE 204 / grey 40 / not_SINE 827, 41 plates — the numbers of section 10. Times:
`mini` 671 s native, 673 s in the image (both environments running at once).

### 11d. Left open

- Only `step2_asSINEment.sh`, `step3_postprocess.sh` and `step8a_extract_alignments.sh` force `LC_ALL=C`. The
  other sorts of the workflow (`sear`, `SINEderella`, the flankscan stages) were identical in both locales on
  these inputs, but nothing guarantees that. Exporting `LC_ALL=C` once in `SINEderella` would cover them all;
  it is not done here because it would change more than the two bugs this section is about.
- The toy and `mini` are small. Two runs of the full *R. sinicus* genome in the two environments were not made.

## 12. One collation for the workflow, the commit of every run, and the full genome in both environments (2026-11-23)

The two points left open in section 11d, and the full-genome check they were waiting for. Code: commit `0bac7e7`.

### 12a. `export LC_ALL=C` once, in `SINEderella`

Section 11a fixed the sorts of step 3 one by one. Three scripts forced the collation (`step2_asSINEment.sh`,
`step3_postprocess.sh`, `step8a_extract_alignments.sh`); every other sort of the workflow ran in whatever locale the
caller happened to have. `SINEderella` now exports `LC_ALL=C` before anything else (line 77). The per-script settings
are kept: they also protect a stage that is run on its own.

Every stage is an ordinary child process (`bash "$SCRIPT_DIR/step..."`, `python3 "$SCRIPT_DIR/tools/..."`), nothing in
the repository sets or unsets a locale variable afterwards, and nothing is started with `env -i`, so the export reaches
all of them. Checked at run time, not by reading the code: `/proc/<pid>/environ` of the 23 live processes of the two
full runs below, native and in the image.

| process | native run | container run |
|---|---|---|
| `rsi.sh`, `/usr/bin/time`, the Apptainer starter and `squashfuse_ll` (above `SINEderella`) | `LC_ALL` unset, `LANG=en_US.UTF-8` | unset, `LANG=C.UTF-8` |
| `SINEderella` itself | `LC_ALL` unset at exec, set by the script | the same |
| `step1_search_extract.sh`, its two `tee`, `sear`, `ssearch36`, `awk`, `sort -V` | `LC_ALL=C` | `LC_ALL=C` |

So the host locale no longer reaches a sort: in the image `LANG` stays `C.UTF-8` and natively `en_US.UTF-8`, and both
collate under `C`.

**What it changed.** Toy and `mini` (with `--families`) were run with the code before (`6fdb10a`) and after
(`0bac7e7`), on the same inputs, same 8 threads, publish on, natively: `/home/toki/tmp/sdc4/{toy,mini}/{before,after}/`.

`tests/compare_runs.py` reports the two runs identical on both inputs — every file it reads: loci, hits, the cuts of
step 1 and step 2, the assignment, `ties.tsv`, the step-3 tables, the family view, the verdicts, the published
consensuses and all 41 plates. `toy` 225 loci (SINE 217 / grey 2 / not_SINE 6, 4 plates); `mini` 1,071 loci, 1,027 in
`sim_scores.tsv`, 73 flagged close seconds, SINE 204 / grey 40 / not_SINE 827, 41 plates — the numbers of section 10.
**No result changed.**

Of all 977 files of the `mini` run directory, 72 differ, and each is accounted for:

| files | why |
|---|---|
| 38 PDF plots of step 2 | identical after the `CreationDate` / `ModDate` lines are dropped (the PNG plots are identical, as in section 10) |
| 12 `__pycache__/*.pyc` | bytecode of two source trees at two paths |
| 11 logs, `manifest.txt`, `results/report.html`, `summary.txt`, `verdict_chain.done`, `consensus.path`, two `border_loop.json` | identical once the commit, the run path, the Plotly `div` ids, the dates and the durations are masked |
| `run_info.txt` | one line: the PID |
| `subfam_input/subfam_log.txt` | two lines: SubFam's `mktemp` directory names |
| `vote_summary.tsv.gz` | the same content; gzip stores the mtime in its header |
| `publish/verdict_chain.sh` | the changed script itself, copied into the run directory |

One real change, cosmetic: `printf "DATE\t%s\n" "$(date)"` in the manifest no longer follows the locale.
`Mon Nov 23 12:27:16 PM MSK 2026` became `Mon Nov 23 12:27:16 MSK 2026`, and `results/report.html` shows that string.
It is a timestamp for a reader, nothing reads it, and the C form is the same on every host. Left as it is.

### 12b. A run always records the code that made it

The container recorded its commit and a native run from an unpacked tarball recorded `PIPELINE_COMMIT unknown`: there
is no git in such a tree, and the one-line `PIPELINE_COMMIT` file was written by hand for some runs and forgotten for
others (the two full runs of section 10 both say `unknown`).

`tools/pipeline_commit.sh` now resolves it in one place, for `SINEderella` (`manifest.txt`) and
`publish/verdict_chain.sh` (`verdict_chain.done`):

1. `$PIPELINE_COMMIT` from the environment, if set and not empty — it wins, so a wrapper can override.
2. `git describe --always --dirty`, but only when the tree **itself** is a checkout (`TREE/.git` exists). Without that
   test git walks up, and a tree unpacked inside somebody else's repository would be stamped with that repository's
   commit.
3. `TREE/PIPELINE_COMMIT`, first non-empty line, CR and blanks stripped.
4. `unknown` — and the run writes two `WARNING` lines in its log saying how to fix it, which it did not do before.

`tools/pack_tree.sh` makes the tarball for a machine without git: `git archive` of the commit plus a
`PIPELINE_COMMIT` file, so nobody has to remember the hand step. It refuses a dirty checkout, because `git archive`
would silently pack the committed tree instead of what is there.

Checked, all three ways of installing the code:

| install | what the run records | checked |
|---|---|---|
| plain `git clone` of `main` + run | `PIPELINE_COMMIT 0bac7e7`, and `code 0bac7e7` in `verdict_chain.done` | therioserver, clone of `main` at `0bac7e7`, toy run `/home/toki/tmp/sdc4/clonetest/run/` |
| tarball of `tools/pack_tree.sh`, unpacked, no `.git` | `PIPELINE_COMMIT 0bac7e7` | the toy, `mini` and full-genome native runs below |
| image `sinederella_0bac7e7.sif` | `PIPELINE_COMMIT 0bac7e7`, from `/opt/SINEderella/PIPELINE_COMMIT` written at build | `apptainer exec`, and the container full run below |
| `git archive` tree with neither `.git` nor the file (the old failing case) | `unknown`, and two `WARNING` lines in the log | therioserver, `/home/toki/tmp/sdc4/clonetest/run2/` |

`tests/test_reproducible.py::test_pipeline_commit_is_recorded_without_git` covers the five cases: no git and no file
(unknown), the file alone, the environment over the file, an empty file (unknown), a real checkout over a stale file,
and a tree unpacked inside another checkout.

### 12c. The full *R. sinicus* genome, native and in the image

Section 11c compared the two environments on the toy (0.5 Mb) and on `mini` (11.4 Mb) only. One run of the whole
*Rhinolophus sinicus* genome was now made natively and one in the image, with the code of 12a and 12b in place, both
at 24 threads, started together on an otherwise busy server (load 19 of 128 before the start, 56 at the peak).

Inputs and options of the section-10 runs: `rsi.fna`, `bank21_nor1.fa`, `--families rsi_families_nor1.tsv`, no
publish. Runs: therioserver `/home/toki/tmp/sdc4/rsi/{nat,cnt}/run_20261123_153329`; image
`/home/toki/tmp/containers/sinederella_0bac7e7.sif`. Native 1 h 13 min 11 s (12.3 GB peak), container 1 h 13 min 52 s
(section 10, one environment: 1 h 10 min 37 s at 32 threads and 1 h 15 min 14 s at 24).

**`tests/compare_runs.py`: the two runs are identical** — every file it reads, including the 1.6 million hit lines and
the family view:

| what | lines | native vs container |
|---|---|---|
| loci | 136,690 | identical |
| hits | 1,602,403 | identical |
| assignment, `assignment_full.tsv`, verdicts | 68,192 each | identical |
| `sim_scores.tsv` (step 3) / family view | 64,900 / 68,186 | identical |
| `all_sines.bedlike.ALL.tsv` / family view | 68,345 / 68,345 | identical |
| `ties.tsv` | 1,873 | identical |
| `unassigned.tsv` / family view | 3,446 / 155 | identical |
| `summary.by_subfam.tsv` / family view, `family_view.tsv`, `family_view.summary.tsv` | 22 / 3, 68,192, 3 | identical |
| cuts of step 1 and step 2 | 39 queries, 21 consensuses | identical |
| published consensuses | 42 | identical |
| loci, one by one | 68,191 in both, 0 in one only | subfamily differs 0, status differs 0, verdict differs 0 |

Plates are skipped (no publish, as in section 10).

**The totals are those of section 10, to the copy.**

| | section 10 (f4e6811, 2 runs, one environment) | now (0bac7e7, native) | now (0bac7e7, image) |
|---|---|---|---|
| loci | 68,191 | 68,191 | 68,191 |
| assigned | 64,900 | 64,900 | 64,900 |
| without a winner (`no_unanimous`) | 0 | 0 | 0 |
| `rejected_low_bitscore` | 3,291 | 3,291 | 3,291 |
| verdict SINE | 50,907 | 50,907 | 50,907 |
| verdict grey | 12,463 | 12,463 | 12,463 |
| verdict not_SINE | 4,821 | 4,821 | 4,821 |
| flagged close seconds (`ties.tsv`) | 1,873 | 1,873, all `tie` | 1,873, all `tie` |

So nothing in section 10 moved: not the locale change of 12a, not the commit change of 12b, and not the environment.
`manifest.txt` of both runs says `PIPELINE_COMMIT 0bac7e7` — the native run was started from an unpacked tarball, which
before 12b would have written `unknown`.

**pytest** (therioserver, 2026-11-23): 174 passed natively; in the image 173 passed and 1 skipped — the new commit test
skips its git half there, because the image has no git, which is the reason the `PIPELINE_COMMIT` file exists.

### 12d. Left open

- `SINEderella_multi` and the scripts that can be started on their own (`publish_run.sh`, `sear`, the flankscan stages)
  do not export `LC_ALL=C` themselves; three of them already did (12a). Started through `SINEderella` they all inherit
  it. A stage started by hand from a shell in another locale is on the user.
- The two full runs were compared with `tests/compare_runs.py`, not file by file over the whole run directory. Section
  10 did that for two native runs and found only PDF dates, paths, durations and SubFam's `input.msf`; the file-by-file
  comparison of 12a (977 files of a `mini` run) found nothing else either.
- The scripts a full run does not reach still use `-z 11` (section 3), and SubFam's `input.msf` still differs between
  runs (section 10, point 4).
