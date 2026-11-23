# Consensus audit: does each consensus say what its copies say?

A consensus of the bank is a hypothesis; the copies assigned to it are the evidence. `tools/consensus_audit.py` rebuilds every
consensus from its own assigned copies and compares the result with the bank sequence. It never changes the bank: the table
(`results/consensus_audit/summary.tsv`, shown in the report under "Sequence similarity between consensuses") is a decision input.

The same copies answer a second question, added 2026-11-23: **do they resemble each other?** Every other gate of the workflow
compares a copy with a consensus - the per-consensus score cut of step 1, the flank-uniqueness test, the per-copy verdict - so a
consensus whose *composition* alone attracts hits can collect copies that each pass on their own and have nothing to do with one
another, and the report would still list them as a subfamily. The audit now measures that directly and prints it beside the rebuild.

## Method

* Copies: `results/assigned.fasta`, split by the family label in the id (`locus|family|bits`). A family needs at least 20 assigned copies
  (`--min-copies`); fewer are `SKIPPED` (test it on a species with more copies).
* Builder: `tools/vendor/sine_consensus.sh`, a seeded copy of `sine_consensus.sh` of github.com/Toki-bio/SINE_consensus (pinned to
  commit 29719f5, 2026-03-13). Each round draws 100 random copies (`-n`), aligns them (MAFFT `--auto`), calls a mini-consensus, adds it to a
  master set, realigns the master and calls the next consensus; it stops when the change between rounds is below 1 % (at least 5 rounds),
  with an anchor phase (copies of the best consensus so far mixed into later draws) if it oscillates. A base is called where at least
  50 % of **all** copies carry it and no other base ties it (gaps are in the denominator; a column with a base in fewer than 30 % of the
  copies is skipped): this trims ragged flanks and keeps a real conserved tail. (Until 2026-10-21 this page and the tool's docstring said
  30 % and 10 %; the script's defaults have been 50 % and 30 % since it was added, and the audit passes no other values.)
* Changes to the upstream script: `-r SEED` (the subsample of round N uses seed SEED×1000+N, so a run is reproducible); MAFFT with
  `--threadit 0` (see PUBLISH_WORKFLOW.md); an error unless awk is gawk (the column counting uses arrays of arrays).
* Each family is rebuilt with two seeds (`--seeds 1,2`), each rebuild is aligned globally to the bank sequence (match +1, mismatch −1,
  gap −1) and the table gives bp, mismatches and gap columns, and the mismatches between the two rebuilds.

## Identity between the copies (copy against copy, not copy against consensus)

* **Sample.** The assigned copies of a family are sorted by id and, if there are more than `--identity-sample` (60), a seeded draw
  (`--identity-seed`, fixed at 20261123) takes 60 of them; the drawn ids are sorted again. Same assigned set → same 60 copies →
  same numbers. 60 copies give 1770 pairs, which is plenty for a median; the cost is one alignment per family.
* **Alignment.** `mafft --threadit 0 --auto --quiet --preservecase` on those copies - the same binary and the same determinism
  setting as the rebuild. The alignment is kept, so the number can be read off a plate:
  `results/consensus_audit/copies/<family>.sample.aln.fa`.
* **Identity of one pair.** Over the columns where **both** rows carry A/C/G/T: matches ÷ shared columns. A gap or an N in either
  row removes the column from both the numerator and the denominator, so an insertion in one copy neither helps nor hurts the pair.
  A pair sharing fewer than 30 columns is dropped. The table gives the median and the 10th-90th percentile over all pairs, with
  `copies_compared` and `copy_pairs`.
* **Copy set.** The copies the *rebuild* uses, that is every copy in `results/assigned.fasta` with that label, whatever its per-copy
  verdict - so the identity column and the rebuild columns rest on the same evidence. A family needs at least 4 assigned copies;
  below that the cell is empty. This is deliberately **not** restricted to the flank-unique copies: it is a measurement table, not
  a plate to decide a locus from.
* **Composition.** `bank_gc`, `copy_gc_median` (median GC of the sampled copies) and `rebuilt_gc_s<seed>` are in the table because
  composition is the usual reason a consensus draws unrelated copies: a GC-rich consensus matches GC-rich genomic sequence.
* **It is computed for every family that has copies, including the ones with fewer than `--min-copies`** whose rebuild is `SKIPPED`
  (that is where MEG-T2 in rsi sits: 12 assigned copies, no rebuild, and an identity the rebuild could never have shown).
* **No threshold, and no verdict depends on it.** The `verdict` column is still only about the rebuild. Where the line between a
  real family and a composition artefact falls is a decision to be made from the distribution across real subfamilies, not a
  constant in this script.
* In the report the audit table is **sorted by this median, lowest first**, with a bar per row, so a consensus whose copies do not
  resemble one another is at the top of the table without anything being flagged.

## Verdicts

| Verdict | Meaning |
|---|---|
| `MATCH` | both rebuilds within 4 mismatches of the bank and within max(8 bp, 5 %) of its length |
| `SHORTER` | a rebuild is shorter: the copies are fragmentary, or the bank has a stretch that fewer than half of the assigned copies carry |
| `LONGER` | a rebuild is longer: a stretch is carried by at least half of the copies but is not in the bank (a longer version in the pool, or a tail the bank lacks) |
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
* By hand: `python3 tools/consensus_audit.py RUN_DIR [--jobs 8] [--min-copies 20] [--seeds 1,2] [--subsample 100]
  [--identity-sample 60] [--identity-seed 20261123]`. Needs mafft and gawk.
* Tests: `tests/test_consensus_audit.py` (alignment counts, every verdict, copy splitting, the pairwise-identity arithmetic on a
  hand-made alignment plus a random-sequence and a one-ancestor control; an end-to-end toy when mafft and gawk exist, which must give
  the same rebuild **and the same table** twice for the same seed). `tests/compare_runs.py` compares
  `results/consensus_audit/summary.tsv` between two runs of the same input, so both draws staying reproducible is checked by
  `tests/test_reproducible_run.sh` like the rest of the run.
