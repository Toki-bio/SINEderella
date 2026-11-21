# Reproducibility: the same input gives the same output (2026-11-21)

Branch `reproducible`. Not merged: the fix changes results beyond the old run-to-run noise (section 6), and three choices
in it are the owner's to make (section 8).

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
| `tools/extend_ends.py:119-121` | that score column (cross) | group of a locus for the end extension, which changes locus coordinates |
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
| `sear` | `-z 11 -E 2` | `-z 3 -E 1e-10` (`SEAR_EVALUE`) | one part: of 865 hits always reported 855 kept, 364 hits never reported before come in; whole window bank: 16,591 of 16,619 kept, 4,525 come in |
| step 2 | `-z 11 -E 2`, 10 cycles | `-z 3 -E 1e-7` (`STEP2_EVALUE`), once | the weak alignments in question are runner-ups far below the winner |
| `composite_units.py` | default, `-E 1e-5` | `-z 3 -E 1e-17` | weakest always reported 5.8e-18, strongest never reported 6.4e-18 |
| `continuation_check.py` | default, `-E 1e-5` | `-z 3 -E 1e-16` | tails vs bank: 2,588 of 2,602 kept, 8 come in |
| `length_variants.py` | `-z 11 -E 1e-3` | `-z 3 -E 1e-17` | 229 of 229 kept, 1 comes in |
| `satellite_trf_verify.py` | `-z 11 -E 1e-2` | `-z 3 -E 1e-7` | one unit in the test run: not measured, set by analogy |
| step 3, step 4, `fs9_bodycheck.py`, `trna_heads.sh` | `-E 100` / `-E 100` / `-E 100` / `-E 10` | the same numbers with `-z 3` | these cuts are meant to pass everything |

The old cut was not one line on the fixed scale: under shuffles each consensus had its own threshold (bit score on the
fixed scale 57-70 for most consensuses, 110-125 for `r10_r8_P18`, above 73 for `MEG-T2`), because a consensus with
GC-rich or repeated stretches scores high against shuffled sequence too. A fixed cut cannot follow that (section 6).

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

- `tests/test_reproducible.py` (pytest; needs ssearch36, seqkit, bedtools): the fixed call gives the same lines at 1, 2
  and 4 threads; the control, the old `-z 11` call, differs between runs; step 2 run twice with 2 and 5 threads writes
  byte-identical `assignment_full.tsv`, `all_votes.tsv`, `assignment_stats.tsv`, `assigned.fasta`, `unassigned.fasta`;
  the continuation check gives the same table at 1, 2 and 4 threads and in a long path.
- `tests/test_reproducible_run.sh GENOME BANK [WORKDIR] [SPECIES]`: the whole workflow twice with publish, with different
  numbers of threads, then `tests/compare_runs.py`: loci, hits, assignments, verdicts, published consensuses and every
  plate must be identical. Minutes; run on a server.
- `tests/compare_runs.py RUN_A RUN_B [--loci OUT.tsv]`: exit 0 only when the two runs are identical in those files.

pytest on therioserver (conda env sinederella): 171 passed.

## 6. Verification

Two native runs of the fixed code on each input, with 16 and 10 threads (therioserver, `/home/toki/tmp/repro/e2e_z3c/`):

| input | loci | hits | assignment | verdict | published consensuses | plates |
|---|---|---|---|---|---|---|
| toy | identical (225) | identical | identical | identical | identical | 4 of 4 identical files |
| mini | identical (770) | identical | identical | identical | identical | 37 of 37 identical files |

`mini` takes 408-423 s (old: 470-484 s): one search in step 2 instead of ten.

**Fixed run against the old runs.** `toy`: every one of the 225 loci has the same subfamily, status and verdict as in
the three old runs (217 SINE, 2 grey, 6 not_SINE).

`mini`:

| | old natA / natB / natC | fixed (both runs) |
|---|---|---|
| loci | 1,041 / 1,046 / 1,042 | 770 |
| assigned | 232 / 226 / 241 | 242 |
| no_unanimous | 805 / 815 / 797 | 505 |
| rejected_low_bitscore | 4 / 5 / 4 | 23 |
| verdict SINE | 190 / 174 / 198 | 178 |
| verdict grey | 25 / 35 / 26 | 45 |
| verdict not_SINE | 826 / 837 / 818 | 547 |
| plates | 39 / 39 / 41 | 37 |

The fixed run differs from natA, natB, natC at 866, 888 and 847 loci (subfamily, status or verdict, or present in one
run only). Of these, 72, 94 and 53 are loci that had flipped among the old runs. **794 loci were the same in all three
old runs and are different now.** This is not what was asked for (changes only at loci that flipped before), and each
group has a cause:

| loci | old (all three runs) | fixed | cause |
|---|---|---|---|
| 396 | MEG-RL, no_unanimous, not_SINE | MEG-RS, no_unanimous, not_SINE | Only the name of the best consensus of a copy without a winner. These copies score about 2 % higher on MEG-RS (margin 0.021-0.023); under shuffles most of the 10 votes went to MEG-RL (7 or 8 of 10 for most copies). Status and verdict unchanged. |
| 304 | MEG-RL, no_unanimous, not_SINE | absent | The satellite stage excludes their regions: with fixed statistics `MEG-T2` and a piece of `r10_r8_P18` have many weak raw hits in arrays there (section 7). All were not_SINE. |
| 62 | absent | 23 r10_r8_P18 rejected_low_bitscore (not_SINE); 21 MEG-T2 assigned (15 grey, 6 SINE); 15 MEG-RS no_unanimous; 2 r6_210seqs assigned (not_SINE); 1 MEG-RL | New loci from hits that the shuffle statistics had rejected (section 7). |
| 15 | 2 SINE, 5 grey, 8 not_SINE | absent | Loci in the newly excluded satellite regions, or made from weak hits now below the cut. Not traced one by one. |
| 7 | SINE | grey (6 with the same subfamily, 1 r8_83seqs -> r10_19seqs) | Not traced. The verdict comes from the top-100 copies and the border loop, which now see another set of copies. |
| 2 | assigned, SINE | no_unanimous, not_SINE | The margin rule (winner ahead by less than 0.05), by its definition; not looked at one by one. |
| 5 | no_unanimous, not_SINE | assigned (2 SINE, 1 grey, 2 not_SINE) | The margin rule (winner ahead by 0.05 or more although the 10 votes had split); not looked at one by one. |
| 3 | no_unanimous, not_SINE | the same, other name | 2 MEG-RL -> MEG-TR, 1 r5_27seqs -> r5_r5_P48: the name of the best consensus of a copy without a winner. |

Counted over the loci that were stable in the old runs: 11 lose the verdict SINE (7 to grey, 2 to not_SINE, 2 absent),
8 gain it (6 new MEG-T2 loci, 2 formerly without a winner). The SINE total, 178, lies inside the old range (174-198).
The rows marked "not traced" are not explained.

## 7. What the fixed scores cost

The shuffle statistics did one useful thing that fixed parameters do not: a consensus that scores high against shuffled
sequence of the same composition (GC-rich, with repeated motifs) got a higher threshold. On `mini`, `MEG-T2` had no hit
at all under shuffles and has 31 hits now, all weak (bit score 62-73 on the fixed scale; 65-71 % identity in the five
looked at), and 759 raw hits of any length; `r10_r8_P18` goes from 224 raw hits to 1,146. The satellite stage then flags
their arrays (satellite loci 25 -> 151) and excludes those regions with the loci in them. This is the largest part of the difference to the old runs. A
bank with such a consensus will gain weak hits on other genomes as well.

## 8. For the owner to decide

1. **`WIN_MARGIN` 0.05.** It reproduces the 10/10 vote on 96 % of the `mini` copies. 0.03-0.07 do about the same; the
   group that the vote left unassigned sits at 0.021-0.023.
2. **The `sear` cut `-E 1e-10`**, and what to do about consensuses that score high against random sequence (section 7).
   Options: leave it (the satellite stage and the verdicts deal with the extra hits); a stricter cut (1e-13 removes most
   `MEG-T2` hits and about 4 % of the weak hits of the other consensuses); or a threshold per consensus taken from a
   search against a shuffled copy of the library made with a fixed seed outside ssearch36 (keeps the old behaviour, costs
   a second search, not written).
3. **The status name `no_unanimous` and `Votes` 10 / 5** are kept so that every reader of `assignment_full.tsv` works
   unchanged. They now mean "no clear winner" and "clear winner / none".

Also: the scripts outside the full run (section 3, last paragraph) still use `-z 11`.
