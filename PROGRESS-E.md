# Task E - end evidence for the open composite candidates

## Goal

`tools/cand_ends.py` reads candidate alignments (row 1 = candidate consensus,
dashes where only flank; other rows = copies, lower case = flank, upper case =
element) and reports end evidence for the 5' and 3' end of every candidate.
The evidence for the P9/P10/P21/P41 candidates and for chain_C5 then goes
into `docs/OPEN_CANDIDATES.md`, which says per candidate what the ends look
like and whether, and by how many bp, the consensus should be extended.

## Status

- [x] `tools/cand_ends.py` written (pure Python 3, standard library only).
- [ ] Check on the small built alignment prints PASS.
- [ ] Tool run on the P9/P10/P21/P41 candidates and on chain_C5.
- [ ] `docs/OPEN_CANDIDATES.md` written (one table plus one paragraph per
      candidate, evidence numbers from the tool output only).

## Log

- `tools/cand_ends.py`: for each end it reports the consensus end column,
  the first outside column whose top-base share falls under 60 %, the
  carrying shares 10/30/60 columns beyond the end, the commonest 12-mer
  outside the end with its share, A-run / oligo-G dominance (run >= 8 in
  >= 40 % of the copies) and one verdict: closed, open <share>, poly-A
  tail or oligo-G tail.
