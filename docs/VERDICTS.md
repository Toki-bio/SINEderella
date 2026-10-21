# Verdicts: every copy SINE, grey zone or not SINE

The aim of SINEderella (docs/AIM.md) is the list of SINE loci of a genome, with every initial copy verdicted. An
interspersed repeat is a set of copies with unique flanks, so the verdict is a flank verdict: a copy is a SINE locus only
when its flanks were tested and found unique. A test that could not be done is not a pass.

The verdicts are written by `tools/flank_verdict.py`, run by `publish/align_for_publish.sh` after step 7 and the border
loop, before the plates are drawn. Plates, page and counts read them; nothing downstream decides on its own.

## Outputs (step2 output folder)

| file | one row per | columns |
|---|---|---|
| `loci_verdict.tsv` | copy (every copy of `assignment_full.tsv`) | copy, subfamily, assignment_status, verdict (`SINE`, `grey`, `not_SINE`), reason |
| `flank_verdict.tsv` | subfamily | subfamily, verdict, side5, side3, source, assigned_copies, sine, grey, not_sine, reason |

Copy names encode `_` as `@U@` (the sanitized genome headers).

## Subfamily verdict: are the subfamily's flanks independent?

Each side (5′, 3′) is tested on the copies of the subfamily:

- **step 7** (`step7_boundary_refine.sh`) widens the flank window in 50 bp steps up to 1 kb until pairs of copies are no more
  alike there than random genome pairs (`confirmed`), or the limit is reached (`undetermined`). A side is shared when more
  than background + 5 % of copy pairs have elevated identity. It runs from **7 copies**: with 21 pairs one duplicated pair
  stays under that threshold, so a pair of twins (a copy-level matter, below) cannot condemn the subfamily; below 7 the
  subfamily is `insufficient_data`. Families of up to 17 copies (150 pairs or fewer) have every pair compared once. Until
  2026-10-21 the gate was 20 copies and smaller families were grey only for their size.
- **the border loop** (`publish/border_loop_subfam.py`) runs for subfamilies step 7 did not confirm: it extends the element
  and rescans each side until it is independent or a cap is hit.

The deciding result per side is the border loop's when it ran, otherwise step 7's: `independent`, `still_shared` or
`untested`. Subfamily verdict: `interspersed` (both sides independent), `not_interspersed` (a side still shared at the
limit), `untested` (a side could not be tested), `not_tested` (no flank test ran).

## Copy verdict

| copy | verdict | reason written |
|---|---|---|
| not assigned (no unanimous vote, low bitscore) | not_SINE | `not assigned (...)` |
| assigned, subfamily `not_interspersed` | not_SINE | subfamily flanks still shared at the expansion limit |
| assigned, subfamily `untested` / `not_tested` | grey | subfamily flank uniqueness could not be tested |
| assigned, interspersed subfamily, stage 9 says twin / array | not_SINE | flanks shared with another copy |
| ... in a group of >= 10 copies sharing the same flank at the same place | not_SINE | same flank sequence at the same place in N copies |
| ... stage 9 says masked or untestable | grey | flanks masked (repeat), not testable |
| ... no stage-9 result | grey | no stage-9 flank result for this copy |
| ... stage 9 says unique | **SINE** | assigned; flanks unique |

Stage 9 is flankscan's flank-twin test (`flankscan/fs9_run.sh`, `fs9_twins.sh`, docs/FLANK_UNIQUENESS.md): each copy's
flanks are compared with those of the other copies of its family, after a genome-wide 20-mer repeat mask. It runs on every
family with **2** or more assigned copies (20 until 2026-10-21). `SAME_PLACE_MIN` (default 10) sets the group size.

## Where the verdicts are used

- **Plates** (`step8a_extract_alignments.sh`, `tools/plate_filter.py`, docs/PLATES.md): top 100 and random 100 are drawn
  from SINE loci. A subfamily with fewer than 100 SINE loci is topped up with its grey copies, marked `[grey]` (owner
  2026-10-18: grey copies are still aligned and shown). not_SINE copies are never on a plate. A subfamily with no SINE
  locus and no grey copy gets no plate; its table row says why.
- **Page**, table "Subfamily alignments" (`disc/inject_disc_report.py`): per subfamily the consensus length, assigned
  copies, SINE loci, grey and not SINE; a plate button reading "top 26 (11 grey)" holds 15 SINE loci and 11 grey copies.
  When the run has a family table (`--families`, docs/FAMILIES.md) the rows are grouped under one row per family with the
  family's totals. The overview paragraph gives the run totals and names the subfamilies that are not interspersed or
  could not be tested.
- **Counts**: the page's copy columns count SINE loci; `summary.by_subfam.tsv` and the step 2 tables still count
  assigned copies.

## Without verdicts

`step8a` stops when `loci_verdict.tsv` is missing (a plate would show copies whose flanks were never tested), unless
`PLATE_VERDICT=off` is set; the plate toy test (`tests/toy/run_toy.sh`) uses that for its first part only.
