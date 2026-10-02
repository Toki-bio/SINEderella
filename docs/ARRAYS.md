# Tandem arrays: copies that are units of an array, not independent insertions

Near-identical units of a tandem array rank first by bitscore and share their flanks. A plate of such copies shows one repeated unit
~100 times, its shared flanks read as a continuation of the element, and the report can still call the family a SINE (found on the
bat MEG-RS plates 2026-09-28 and again on *R. sinicus* MEG-RS 2026-10-02: 1 526 of 1 715 copies in 30 clusters, 77 of the top 100
from two arrays with spacings of 2.1 and 2.9 kb, flanks 90 % identical over 300 bp each side).

## Detection (`tools/array_order.py`, `regular_runs`)
A tandem run is at least `REG_MIN` = 5 consecutive copies on one contig, each gap at most `REG_GAP` = 10 kb and within a factor
`REG_RATIO` = 2 of the run's median gap. Applied to **all** loci of a family, not only the 300 best-scoring ones (the earlier rule,
neighbours within 50 kb among the first 300, is kept and the union is used: it catches arrays of large period, e.g. hla MEG-RL ~27 kb).
Calibration: simulated dispersed copies (exponential gaps, mean 20 kb, the densest realistic case) are marked in 1.5 % of copies,
at mean 50 kb in 0 %; on the 21 families of the *R. sinicus* run, 20 families have 0.0-0.2 % of copies in runs and MEG-RS 70.8 %.
Unit-length variants of one array (1.4 / 2.0 / 2.3 kb) stay in one run (tests/test_array_order.py).

## What it changes
* **Plates** (step8a): the top 100 takes every independent copy first, then the best copy of each array, then further array copies;
  array rows are marked `[array]` on every plate (top100 and rand100).
* **Flag** (`tools/array_flag.py`, run after assignment): `results/array_flag.tsv` per family: copies in arrays, share, number of arrays,
  largest array, median spacing; `ARRAY` when at least 20 % of the copies are in arrays.
* **Report:** the overall chip of a flagged family reads "Tandem array" instead of "SINE", with the numbers in its tooltip; the
  Similarity section has the table. Nothing is removed from the assignment.

## Limits
A family that is an array is not thereby "not a SINE": its units may contain a SINE-derived core (MEG-RS: 135 bp core in a ~2 kb unit).
The call needs the copies outside the arrays (few for MEG-RS in *R. sinicus*) or another species. The flank scan (TRF, window 1 kb)
cannot see periods above ~1.1 kb; this rule can, because it uses copy positions.
