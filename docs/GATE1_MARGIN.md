# Gate 1 keeps a margin above the shuffle noise floor

Owner, 2026-11-24: "require a margin in gate 1."

## What gate 1 is

A copy is verdicted `SINE` (docs/VERDICTS.md) only when three gates pass:

1. **assigned** - the copy's best local hit to the bank beats that consensus's own score cut;
2. the **subfamily's flanks** are independent on both sides (step 7, border loop);
3. **this copy's flanks** are unique (flankscan stage 9).

Gates 2 and 3 compare *flanks* between copies. Gate 1 compares a copy with a *consensus*. Nothing in the
chain compares the body of one copy with the body of another.

## Why a margin was needed

The cut of gate 1 comes from `tools/shuffle_cut.py`: the copies of the first search part are shuffled with a
fixed seed, scored against each consensus, and an extreme-value (Gumbel) fit gives

    cut = mu + beta * ln(nlib / expect)

`expect` was **2** - the bit score that 2 shuffled copies are expected to reach, which is what `E <= 2` meant
under ssearch36's own `-z 11` shuffles. That is the noise floor itself, so a copy that beat it by a hair was
assigned, and nothing downstream looked at the element again.

Measured on the full *Rhinolophus sinicus* run
`therioserver:/home/toki/tmp/repro/rsi_full/B/run_20261123_072054`
(genome `therioserver:/home/toki/rhin/genomes/rsi.fna`,
bank `therioserver:/home/toki/rhin/final_bank_20261026/bank21_nor1.fa`),
as `(copy bits - cut) / beta` for the weakest assigned copy of each subfamily:

| subfamily | copies | weakest copy, beta over cut | its E |
|---|---|---|---|
| r2_3seqs | 17 | **0.07** | 1.9 |
| r4_32seqs | 17 | **0.07** | 1.9 |
| MEG-TR | 15 | 2.75 | 1.3e-1 |
| MEG-T2 | 12 | 3.23 | 7.9e-2 |
| MEG-RL | 19 | 3.56 | 5.7e-2 |
| *(empty range 3.56 - 21.73)* | | | |
| r9_15seqs | 6915 | 21.73 | 7.3e-10 |
| MEG-RS | 770 | 24.29 | 5.7e-11 |
| r3_58seqs | 14 | 26.37 | 7.1e-12 |
| ... | | | |
| r1_r2_r3_r3_C11 | 3021 | 147.42 | 1.9e-64 |

Exactly those five subfamilies are the ones whose copies do not resemble each other. Pairwise identity of
every pair against a composition-matched null (each copy shuffled, same aligner, 3 replicates):

| subfamily | observed | null | excess |
|---|---|---|---|
| r4_32seqs | 74.5 | 73.6 | **+0.9** |
| r2_3seqs | 77.0 | 76.0 | **+1.0** |
| MEG-T2 | 76.9 | 75.7 | **+1.2** |
| MEG-TR | 76.2 | 75.0 | **+1.2** |
| MEG-RL | 78.0 | 75.9 | **+2.1** |
| r9_15seqs | 84.3 | 74.0 | +10.3 |
| r3_58seqs | 87.9 | 76.0 | +11.9 |
| r1_r3_P1 | 90.4 | 75.2 | +15.2 |
| MEG-RS | 94.3 | 76.9 | +17.4 |

A bare identity number says nothing here: 12 random genome windows score 70.1 on this scale (and 37.1 on a
MAFFT multiple alignment). Only the excess over the matched null separates the two groups, and it separates
them with nothing in between.

## The margin

`STEP2_EXPECT`, default **1e-7**, passed to `tools/shuffle_cut.py cut --expect`.

1e-7 is not a new number: it is the stringency the *other* branch of the same filter already demanded. A
consensus with no shuffle cut (fewer than 30 shuffled copies with a score) was always kept only at
`E <= 1e-7`. Both branches of gate 1 are now equally strict, and the single value is read from one variable
instead of two literals.

In beta units the margin is `ln(2 / 1e-7)` = **16.8**, inside the empty range 3.56-21.73. Because the margin
is expressed through the per-consensus Gumbel fit, it is automatically scaled for each consensus's length and
composition; a margin in bits would not be (the run's shortest consensus, r9_15seqs at 105 bp, clears its cut
by 142 bits on average, while MEG-T2 at 232 bp clears it by 16).

`STEP2_EXPECT=2` restores the old behaviour exactly.

## What this does and does not fix

- It stops a copy being **assigned** on a hit that is at the noise floor.
- It does **not** add a copy-against-copy test. Gates 2 and 3 still only compare flanks. A subfamily whose
  copies pass gate 1 with a wide margin but still do not resemble each other would not be caught; no such
  case is present in the rsi run, where the two measures agree.
- Step 1's own cut (which loci are found at all) is unchanged.

## Gate 1 is two gates, and the second one cannot see this problem at all

Step 2 applies two different thresholds:

- **1a, the shuffle cut** (`consensus_cuts.tsv`): absolute and properly calibrated per consensus, but it was set
  *at* the noise floor (`expect` = 2). This is what the margin above fixes.
- **1b, `thresholds.tsv`** (step2_asSINEment.sh:416 and again at :818): `MIN_REL_BITSCORE x TopN of that
  subfamily` - a threshold **relative to the subfamily's own best copies**. A copy below it is kept in the table
  with status `rejected_low_bitscore`.

1b is self-referential by construction: if every copy of a subfamily is noise, the top copy is noise too, and
the threshold moves down with it. MEG-T2's 12 copies score 76.2-86.3 bits and its 1b threshold is 34.5 bits, so
all 12 pass. **A relative threshold can never detect that a whole subfamily is noise.** Together with gates 2
and 3 looking only at flanks, that is the full reason 10 MEG-T2 copies came out verdicted SINE.

(Note the duplication: the same threshold logic exists twice in the script, at :416 and :818.)

## Open issue from the margin: a rejected copy loses its row

Measured on the step-2 toy fixture (3 consensuses x 4 divergence levels x 10 copies = 120), run both ways:

| STEP2_EXPECT | rows in `assignment_full.tsv` |
|---|---|
| 2 (old) | 120 |
| 1e-7 (new) | **115** |

The five copies are the most diverged ones (`TC_p0.4` x4, `TB_p0.4` x1). They are present in
`unassigned.fasta`, so nothing is lost from the run - but gate 1a filters hit *lines*, so a copy with no
surviving line never reaches `assignment_full.tsv`, and `publish/verdict_chain.sh` writes one verdict row per
row of that table. Such a copy therefore gets **no verdict**, which contradicts docs/AIM.md ("every initial
copy verdicted"). `tests/test_reproducible.py::test_step2_twice_identical` and `::test_step2_in_a_long_path`
both assert 120 rows and now fail - correctly: they are catching this, not complaining about the margin.

On the rsi run this did not arise, because step 1's own cut is stricter than the old `expect` = 2, so every one
of the 68,191 loci had a row. With the margin, some will not.

### Fixed (owner chose option 1, 2026-11-24)

`run_votes` keeps each such copy's best sub-cut hit in `below_cut.tsv`, and **both** assignment paths
(`run_assignment`, `run_incremental_assignment`) write it as a row with status `rejected_low_bitscore` and
`Threshold` `NA`, tagged `below_cut=1` in `RunnerInfo`. Purely additive: nothing above the cut changes, and these
copies never enter `unanimous.tsv`, so the relative thresholds of `thresholds.tsv` do not move.

Verified on the toy fixture, both ways:

| | rows | assigned | rejected_low_bitscore |
|---|---|---|---|
| `STEP2_EXPECT=2` | 120 | 91 | 29 |
| `STEP2_EXPECT=1e-7` | **120** | 91 | 29 |

All 120 input copies present exactly once, none missing, none duplicated. The five copies the margin pushed below
the cut (`TC_p0.4` x4, `TB_p0.4` x1) carry `below_cut=1`. Their final status is unchanged, because the relative
threshold 1b had already rejected them - the margin moved them from being rejected by 1b to being rejected by 1a.

The two tests that caught the problem (`test_step2_twice_identical`, `test_step2_in_a_long_path`) pass
**unchanged** - the assertion of 120 rows was not relaxed. Whole suite: 175 passed.

### Still not done

The effect of the margin on the full rsi genome has **not** been measured; a full run is needed to say how many
loci change status there, and the five noise subfamilies' loci are the point of the exercise.
