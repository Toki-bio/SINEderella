# Consensus audit: does each consensus say what its copies say?

A consensus of the bank is a hypothesis; the copies assigned to it are the evidence. `tools/consensus_audit.py` rebuilds every
consensus from its own assigned copies and compares the result with the bank sequence. It never changes the bank: the table
(`results/consensus_audit/summary.tsv`, shown in the report under "Sequence similarity between consensuses") is a decision input.

## Method

* Copies: `results/assigned.fasta`, split by the family label in the id (`locus|family|bits`). A family needs at least 20 assigned copies
  (`--min-copies`); fewer are `SKIPPED` (test it on a species with more copies).
* Builder: `tools/vendor/sine_consensus.sh`, a seeded copy of `sine_consensus.sh` of github.com/Toki-bio/SINE_consensus (pinned to
  commit 29719f5, 2026-03-13). Each round draws 100 random copies (`-n`), aligns them (MAFFT `--auto`), calls a mini-consensus, adds it to a
  master set, realigns the master and calls the next consensus; it stops when the change between rounds is below 1 % (at least 5 rounds),
  with an anchor phase (copies of the best consensus so far mixed into later draws) if it oscillates. A base is called where at least
  30 % of **all** copies carry it (gaps are in the denominator, minimum coverage 10 %): this trims ragged flanks and keeps a real conserved tail.
* Changes to the upstream script: `-r SEED` (the subsample of round N uses seed SEED×1000+N, so a run is reproducible); MAFFT with
  `--threadit 0` (see PUBLISH_WORKFLOW.md); an error unless awk is gawk (the column counting uses arrays of arrays).
* Each family is rebuilt with two seeds (`--seeds 1,2`), each rebuild is aligned globally to the bank sequence (match +1, mismatch −1,
  gap −1) and the table gives bp, mismatches and gap columns, and the mismatches between the two rebuilds.

## Verdicts

| Verdict | Meaning |
|---|---|
| `MATCH` | both rebuilds within 4 mismatches of the bank and within max(8 bp, 5 %) of its length |
| `SHORTER` | a rebuild is shorter: the copies are fragmentary, or the bank has a stretch that fewer than 30 % of the assigned copies carry |
| `LONGER` | a rebuild is longer: a stretch is carried by at least 30 % of the copies but is not in the bank (a longer version in the pool, or a tail the bank lacks) |
| `DIVERGED` | more than 4 mismatches against the bank |
| `UNSTABLE` | the two rebuilds differ from each other by more than 4 mismatches |
| `SKIPPED` / `FAILED` | too few copies / the builder did not finish |

## What it is and is not

It is an independent check and a way to see where a consensus rests on few copies. It is not a replacement for the bank sequence and
convergence of the bootstrap is not evidence of a family (chimeric or boundary-contaminated alignments also converge). A `SHORTER` or
`LONGER` verdict is a prompt to open the plates of that family, not an error. Example: in rsi, r7 rebuilds ~25 bp longer than its bank
sequence, because many r7 copies carry the extra 3′ stretch of r8 (length-version test: r7/r8 are two versions, linkage 0.56).

## rsi (rsi_fresh2, 21 consensuses, two seeds)

Result of the tool on the run: 11 MATCH (r5 with 0 mismatches including the 3′ tail `GTTCCCCAATATTCCCCAATAAAA`, r6, r8, r9, r10, P1, P2, P26, P34, P48, C11),
5 SHORTER (MEG-T2, MEG-TR, r1, r3, P18: few or fragmentary copies), 1 LONGER (r7, +25 bp, see above), 1 DIVERGED (MEG-RS: 9 mismatches, 114 bp
against 135; its assigned pool is mixed with MEG-RL-like copies), 3 SKIPPED (MEG-RL, r2, r4: fewer than 20 copies).
Two rebuilds of one family differ by 0-2 mismatches here (0-7 in the first, unseeded measurement, `rsi_fresh2/consensus_audit.tsv` in Tal):
the subsampling is random by design, only a fixed seed repeats it.

## Use

* In the pipeline: runs after the length-version test; `SKIP_CONSENSUS_AUDIT=1` skips it, `CONSENSUS_AUDIT_JOBS` (default 8) sets the parallel rebuilds.
* By hand: `python3 tools/consensus_audit.py RUN_DIR [--jobs 8] [--min-copies 20] [--seeds 1,2] [--subsample 100]`. Needs mafft and gawk.
* Tests: `tests/test_consensus_audit.py` (alignment counts, every verdict, copy splitting; an end-to-end toy when mafft and gawk exist,
  which must give the same rebuild twice for the same seed).
