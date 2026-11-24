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
