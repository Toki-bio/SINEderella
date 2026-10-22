# Copies that continue past their ends (tools/continuation_check.py)

## Why

On rsi (*Rhinolophus sinicus*, 2026-10-22) the r10_19seqs plates showed a conserved G-rich, then A-rich stretch in the
"right flank". It was not flank: 535 of the 761 r10 copies were a composite the bank did not have, r10 joined directly to
the r8 module of r10_r8_P18 (P18 positions about 290-423, without P18's middle piece 151-291). The vote gave them the label
of the r10 unit; the plates cut them at r10's end, so the module looked like flank; the flank tests read SINE sequence as
flank; and about 216 more such copies matched neither r10 nor P18 well enough and were left unassigned. Nothing in the
workflow reported it. Adding the missing consensus (r10_r8_P18del, built from those copies) to the bank moved 729 copies to
it in the rerun (rsi_main_20261022): 347 from r10, 131 from P18, 216 that had been unassigned, 25 not in the old copy list.

## What the step does

After the verdicts and the consensus audit, for every consensus with at least 20 copies, separately for its assigned copies
and for the unassigned copies whose best (too weak) match it was (`rejected_low_bitscore`, `no_unanimous`):

1. up to 400 copies (seeded) are read from the genome with 150 bp on both sides, in the copy's own orientation;
2. the copy's own consensus is located in it (ssearch36, forward strand);
3. the sequence past each end of that match (the rest of the locus plus the margin) is searched against the bank, and the
   best hit that starts within 20 bp of the end of the own match and is at least 50 bp long is kept;
4. a (consensus, group, side, target) row is flagged when at least a fifth of the tested copies (those in which the own
   consensus was found; `tested` of `sampled` in the table), and at least 10, continue into the same consensus;
5. for a flag into another consensus, a candidate consensus "own match + continuation" is built from those copies (up to
   200) with `tools/vendor/sine_consensus.sh`, two seeds, and compared by a global alignment with every bank consensus
   that ssearch36 hits (not only the best local scorer): if one of them matches it end to end (mismatches + gap columns <= max(8, 5 % of its length)), the bank already has it and the
   copies were assigned elsewhere; otherwise a consensus is missing;
6. a flag into the copy's own consensus is a dimer or tandem when the continuation restarts at or before the end of the own
   match, and an internal deletion of the consensus (P without a middle piece) when it resumes more than 20 bp past it; a
   deletion gets a candidate like a flag into another consensus, its bank answer starting "deletion inside P".

`inside_locus` (median bp of the continuation inside the copy's coordinates, of its length) tells two cases apart: mostly
inside = the locus is long but carries the label of a shorter consensus (the rsi r10 case); mostly outside = the locus
coordinates stop early.

Nothing in the run is changed: copies, verdicts and the bank stay as they are. The table goes on the report ("Sequence
similarity between consensuses" section; `--resume` of a finished run rebuilds the report too). Adding a candidate to the bank is a new run (or `--add`), decided by the user.

Out: `results/continuation/continuation.tsv`, `candidates.fa`, `<candidate>.copies.fa`; log `RUN/continuation.log`.
`SKIP_CONTINUATION=1` skips the step; a failure is a warning.

## Checks

- Tests (`tests/test_continuation_check.py`): a planted genome with unit A, composite P = A + M + B in the bank and 60
  copies of A+B (not in the bank) labelled A. Flagged in both forms, locus stopping at A's end (continuation outside) and
  locus covering A+B (inside); target P 291-420; the candidate is A+B and "not in the bank"; with A+B in the bank the
  answer is "already in the bank as AB"; the 20 unassigned A+B copies with a long M (400 bp, not bridged by ssearch36) are
  reported as a deletion inside P, not a dimer; the report table.
- Real data (therioserver `/home/toki/tmp/cont_oracle`): on rsi_main_20261021b the step flags r10 3' into P18 270-415 in
  234 of 400 sampled assigned copies (59 %) and in 36 % of the unassigned ones, as found by hand. On rsi_main_20261022, with
  r10_r8_P18del in the bank, that flag is gone; a second case appears: 129 of r10's 315 copies (41 %) continue into P18's
  middle piece (positions 174-293) entirely outside their loci, so those loci stop at r10's end.
- After the GLM-audit fixes (therioserver `/home/toki/tmp/cont_oracle2`): on 21b the r10 flag stands (233 of 400, P18
  270-415), and the new deletion rule finds the same missing composite from the other side: 72 of 269 tested unassigned
  P18 copies are "deletion inside r10_r8_P18" (P18 1-77, then resuming at 253), i.e. P18 without its middle piece. On run
  22 (P18del in the bank) that row is gone; r10 into P18 173-291 outside the locus stays (128 of 315).
- The first version read past the locus end, not past the own match; on the real data it missed the r10 case because the
  long r10 loci already covered the r8 module. The test for the long-locus case was added for that.

## Limits

- A candidate is built only from sequence within 150 bp past the own match, and the builder writes a base only where at
  least half of the copies carry it; when the copies' continuations end at different points the candidate is shorter than
  the full element (rsi: 226 and 197 bp). It is a lead for review, not a finished consensus.
- Continuations shorter than 50 bp, or starting more than 20 bp from the end of the own match, are not counted.
- The own consensus is located by a local alignment: when the missing piece is short, ssearch36 aligns across it with a gap
  and the copy shows no continuation (a deletion variant of the own consensus is then seen only by step 7 and the plates).

## GLM audit 2026-10-22

Ten findings. Fixed: the report was not rebuilt on `--resume`; tails and candidates were searched with ssearch36 `-b 5` (a
junction hit to a sixth consensus was never seen); the candidate was compared only with the best local hit; `inside_locus`
was not clamped at the locus start (3') and end (5'); copies without an own match were left in the share's denominator and
not shown; the table header read "Copies" over the group; the unassigned branch had no test, and writing it showed internal
deletions reported as dimers. Not changed, by design: the step is report-only (adding a consensus is the user's decision);
the real-data check is on the author's server, the planted tests are the reproducible part.
