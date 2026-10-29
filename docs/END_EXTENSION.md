# End extension in step 1 (tools/extend_ends.py, 2026-10-29)

## The defect
A locus is the union of the hits of all bank consensuses, so its end is where the longest consensus that hit it stops.
When the genome's element goes on past every bank consensus that hit it, step 1 extracts a cut copy, and everything
downstream (SubFam plate, vote, flanks) sees the cut copy.

Found on hla (*Hipposideros larvatus*) with the rsi bank, peeling the step 1 SubFam plate: two groups of shorties
(r5-type head + linker + GGCAGGATGGCTCAGCT + A tail, 128-134 bp; 75 chunks, 3,750 copies) looked as if they varied by the
presence of that tail. In the genome, 94 % and 98 % of the copies carry it. Whether the extracted copy had it depended on
which consensus hit the copy: 1,583 of 1,642 and 393 of 636 copies with the tail had a hit of r5_r3_P34, whose second module
matches the tail; 27 of 708 and 39 of 764 cut copies had one. After the tail comes the A tail and then unique flank (Tal
`hla_rsibank_20261026/peel/README.md`).

The continuation check (docs/CONTINUATION.md) does not catch this: it runs after step 2, finds only continuations that
match another bank consensus, and changes nothing.

## What step 1 does now
After the satellite regions are removed and before extraction, `tools/extend_ends.py` moves the ends of `merged_hits.bed`
(the original is kept as `merged_hits.bed.before_extension`; `END_EXTENSION=0` skips the step):

1. loci are grouped by the label of their highest-scoring hit;
2. the 100 bp past each end are read in the element's orientation, never past the neighbouring locus or the contig end;
3. per group and side, every 10-mer is counted once per window (A/T-only and one-base-dominated 10-mers are not counted);
   a 10-mer is shared when it is in at least max(20, 2 % of the group's windows) windows; a group with fewer than 10
   shared 10-mers extends nothing;
4. each end is walked outward in 10-bp steps; a step is element when at least 20 % of its counted 10-mers are shared;
   up to two failed steps in a row are bridged (one mutation removes all 10-mers overlapping it); the new end is the far
   end of the last shared 10-mer, plus the A run that follows;
5. a walk still in element at 100 bp is not an end (array, dimer, longer element): the locus is left as it is and
   counted as `at_cap_not_extended`.

Out (`<step1>/end_extension/`): `ends_extension.tsv` (per group and side: windows, threshold, shared 10-mers, loci
extended, median and max bp, at the cap, `flag` when >= 20 % of the group grew by >= 10 bp: the bank consensus is shorter
than the element there), `loci_extended.tsv` (every moved locus). One line in the log: `EXTEND loci=.. extended5=..
extended3=.. flagged=..`.

## Tests
- `tests/test_extend_ends.py` (toy, 8 % divergence): 100 cut copies, 100 complete copies of the same element, 100 copies
  of another: >= 75 cut copies regain the tail and A tail, none goes more than 5 bp into the flank, no complete copy moves,
  only the cut element is flagged.
- hla, real loci of hla_rsibank_20261026 (233,376), against the shortie oracle (tail found by aligning the genome copies):
  cut copies reaching the tail end 566 of 764 and 621 of 708; more than 45 bp past it 2 and 3; copies that had the tail and
  moved 10 of 636 and 9 of 1,642. Flagged 3': r1_r3_P1, r1_r3_P2, r5_r5_P48, r9_15seqs (the consensuses that hit the
  shorties). 1 min 18 s, 1.4 GB.
- rsi control (its own bank; rsi_main_20261026, 68,607 loci): 89 loci moved (MEG-RS and r6 only), nothing flagged; MEG-RS
  loci mostly stop at the cap (698 of 762), as tandem arrays should.

Settings tried on hla (cut copies reaching the tail end, 764 / 708): k12 share 0.5: 287 / 269; k12 0.25: 356 / 336;
k10 0.3: 392 / 406; k10 0.2: 452 / 457; k8 0.3: 469 / 513 but 3.5 x the 5' extensions; k10 0.2 with two bridged steps:
566 / 621 (chosen).
