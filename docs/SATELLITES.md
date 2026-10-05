# SINE-derived satellites and SINE-containing tandem arrays: a screen before the SINE analysis

Design note of 2026-10-03, then the record of the implementation (sections 5b–5h). Sections 1–4 are the design as written before any
code; what runs today is `tools/satellite_stage.py` inside step 1 (5e), with the kind-A gate + windowed TRF (5d), the kind-B unit check
(5h) and the Tal viewer (5g). Not implemented: the characterisation of 4.3 (monomer as a SINE candidate, per-locus analysis).
Status 2026-10-05: the stage is wired in every mode (`--add` screens the added consensuses with `--only`); audit in `docs/AUDIT_2026-10-05.md`.

## 1. Why a separate stage

A satellite whose monomer contains a SINE sequence is not a SINE family, and SINEderella's logic breaks on it:

* **Resources.** Every monomer of every array is a hit, is extracted with flanks, is assigned, enters the plates and the
  flank scan. In *R. sinicus* the 135 bp MEG-RS consensus gives 1 715 assigned copies, 1 592 of them array units; the
  stage 8 border loop walked 1 038 bp into the array and gave up on both sides.
* **Logic.** The acceptance criteria (copies as independent insertions, flanks no more similar than random loci, a 3′ end
  at the genome junction) are built for interspersed elements. Array units share their flanks by definition, the "3′ end"
  is the next monomer, and the top-100 plate fills with the near-identical units, which rank first by bitscore. The
  report then shows one repeated unit a hundred times and the element column reads "Strong" (rsi MEG-RS before 2 Oct).
* **Biology.** Satellites arise by unequal recombination and are homogenised by concerted evolution, live in subtelomeric
  or pericentromeric heterochromatin, and amplify by mechanisms unrelated to retrotransposition. Their copy number, their
  age signal (homogeneity grows with array length) and their chromosomal distribution answer different questions than a
  SINE family's. Mixing the two in one count misstates both.

So the aim is to take satellite loci out of the SINE analysis **before** assignment, and to record and, when wanted,
characterise them separately. The screen must be cheap when a genome has no such satellites, because most do not.

## 2. Two kinds, from the two examples we have

| | A. SINE-derived satellite (sSat) | B. long tandem unit that contains a SINE piece |
|---|---|---|
| Example | sSat1–3 of squamates (Vassetzky, Kosushkin & Ryskov 2023, *Mobile DNA* 14:21) | MEG-RS array in *R. sinicus* (this project); tbr MEG-RS, ~910 bp period |
| Monomer | **shorter than the SINE**: 67–290 bp, a 5′, middle or 3′ portion of a 190–350 bp SINE, sometimes plus 14–62 bp of unknown origin | **longer than the SINE**: 2 155 bp in rsi, of which 107 bp (positions 2–108) are MEG-RS/RL/TR-like at 85 %, the rest unrelated to the bank; 67.6 % GC |
| Array | 4 to ~1 000 monomers per locus; a locus often has a longer *leading* or *trailing* monomer carrying the SINE up to its terminus | 5 to 307 units per array; units 99.9 % identical (u1 vs u2 over 2 155 bp); unit-length variants coexist (rsi: ~1.3, 1.45, 2.16, 4.7 kb) |
| Homogeneity | 77–97 % between monomers, higher in longer loci; locus-specific subvariants alternate | near-complete within an array (concerted evolution far advanced) |
| Where it sits | subtelomeric (snakes, leopard gecko); not centromeric (no CENP-B motif) | 17 contigs, large arrays over 150–560 kb; position on chromosomes not yet looked at |
| How the SINE search sees it | each monomer is a **partial** hit (< 80 % of the consensus) and is dropped by the length rule, or a few monomers merge into one odd locus | each unit gives a **full** hit of the short consensus: a clean, abundant, well-supported "family" |
| What finds it | hits of one consensus < 100 bp apart in runs (the 2023 method); TRF on a 1 kb window sees the period | regular spacing of full hits over the whole family (`tools/array_order.py`); TRF cannot (period above its window, and above its 2 000 bp maximum) |

Case B is the dangerous one for SINEderella, because it looks like a strong family. Case A mostly hides from the search
(the monomers are too short) and shows up as leakage, odd loci and a few composite-looking copies; its danger is to the
copy counts and to the peel, where merged monomer runs appear as long "copies".

## 3. What SINEderella does today

* `MANUAL.md` §16.9 argues **against** a TRF pre-filter: the search needs ≥ 80 % of the consensus length at ≥ 65 % identity,
  which a short-period repeat cannot produce, and a filter would also throw away satellite-adjacent real SINEs. Both points
  stand for short periods. Neither covers case B, where the hit is a complete SINE-length sequence embedded in a longer
  unit; and the fear of discarding real SINEs is answered by *recording* rather than deleting (section 4).
* The flank scan (`flankscan/fs2_trf.sh`) runs TRF (periods 2–2 000 bp) on 1 kb windows around every assigned copy and
  classes a copy `satellite` when a repeat spans it with ≥ 50 bp on both sides. It finds case A arrays (period near the
  SINE's length or shorter) and the A tails; it reported 0 % satellite for rsi MEG-RS, because a 2.1 kb period does not fit
  the window.
* `tools/array_order.py` + `tools/array_flag.py` (2 Oct 2026, `docs/ARRAYS.md`): runs of ≥ 5 copies with gaps ≤ 6 kb within
  a factor 5 of the run's median gap, over all copies of a family. Marks `[array]` rows, puts independent copies first on the
  top-100 plate, flags a family with ≥ 20 % of copies in arrays as "Tandem array" in the report. It finds case B after
  assignment, i.e. after the resources were spent, and it does not characterise the unit.
* The 2023 paper's own method (per family, all hits ≥ 20 % of the consensus length; merge hits < 100 bp apart; loci > 500 bp
  are candidates; ≥ 4 tandem SINE-derived units within 100 bp = sSat; monomers aligned and compared within and between
  loci) is in Perl/bash scripts outside SINEderella and finds case A only.
* `CLsat_workflow` (Toki-bio) is the satellite-side counterpart: TRF motifs (period 80–300, ≥ 10 copies) → simple-repeat
  filter → canonical rotation → vsearch clustering → family consensus → **seed expansion**: walk outward from a seed in
  3-monomer steps with a 3-copy tandem query until the ssearch signal is lost, then merge and rescan at monomer
  resolution → locus catalogue → flank orthology across assemblies. Its expansion-until-signal-lost is exactly what the
  SINE border loop lacks for arrays, and its monomer-resolution scan gives the per-locus monomer count and homogeneity.

## 4. Design of the screen ("step 0s", before assignment)

Principle: **indicate always, characterise on demand, exclude by default, never delete.** Satellite loci leave the SINE
analysis but are written to their own table and FASTA so that nothing is lost and the satellite can be studied with the
CLsat tools.

### 4.1 Indication (cheap, always on)

Input: the step 1 hit table (all `ssearch36` hits of every consensus, before the length rule, with coordinates and
strand) and the genome index. No TRF, no alignment. Two tests, one per kind:

* **A, monomer runs.** Per consensus and contig, sort hits; a run is ≥ 4 consecutive hits of the same consensus on the
  same strand with gaps ≤ 100 bp (the 2023 rule), spanning > 500 bp. Partial hits count (they are the monomers). Report
  the run's span, monomer count, median hit length, and which consensus positions the hits cover (the SINE part of the
  monomer).
* **B, regular spacing.** `regular_runs()` of `tools/array_order.py` on the full-length hits of each consensus (≥ 5 copies,
  gaps ≤ 6 kb, within a factor 5 of the median): the unit length is the median gap. Calibrated on rsi: 20 families 0 %,
  MEG-RS 92.8 %; simulated dispersed copies at 20 kb mean spacing 1.1 %.

Output `results/satellites/indication.tsv`: one row per family with copies, copies in runs (A and B), share, number of
loci, largest locus, unit length(s), and a flag `SAT_A` / `SAT_B` / `-` at ≥ 20 % of copies (the same threshold as the
array flag). Cost: seconds. A genome without SINE satellites costs nothing more than this table.

### 4.2 Exclusion (default when a flag is set)

Loci in flagged runs are written to `results/satellites/loci.bed` and `loci.fa`, and the corresponding hits are **removed
from the hit set that goes to extraction and assignment**. The family itself stays in the bank: its copies outside the
arrays are the dispersed SINE, if any (rsi MEG-RS: 123 copies outside arrays; whether they are real insertions is then a
normal SINEderella question with a normal-sized plate). `SKIP_SATELLITE_SCREEN=1` disables the stage, `SATELLITE_EXCLUDE=0`
keeps the loci in the analysis but still writes the table (for comparison runs). Every removed copy is listed with its
locus, so the count "assigned copies" can always be reconciled with "assigned + in satellite loci".

Why before assignment and not after: assignment is the costly step (ten cycles × all consensuses × all copies) and the
array units are the majority of the copies of such a family; the plates, border loop and flank scan then never see them.

### 4.3 Characterisation (on demand, `--satellites`)

For each flagged family, borrowed from CLsat_workflow and the 2023 paper, on the flagged loci only (not genome-wide):

1. **Monomer.** Case A: the hits themselves, cut at the run's spacing; case B: consecutive full hits of the same orientation
   cut at the median unit length (rsi: 12 units, 2 155 bp each, aligned with no gaps). Consensus by plurality; canonical
   rotation (lexicographically minimal, as CLsat Step 5) so that loci can be compared.
2. **Structure of the monomer.** `ssearch36` of the monomer against the bank (which SINE, which positions: rsi unit 2–108 =
   MEG head), against itself (internal repeats), dust/GC (rsi 67.6 % GC), simple-repeat content (CLsat Step 4 filter).
3. **Per locus.** Monomer count, mean monomer identity (within locus), leading/trailing monomers that carry the SINE
   terminus (the 2023 pattern), unit-length variants (rsi 1.3 / 1.45 / 2.16 / 4.7 kb), subvariants that alternate.
4. **Between loci.** Monomer identity between loci (homogeneity grows with locus length in the paper: 82 → 94 %), loci per
   contig, position relative to contig ends (subtelomeric?), and, when several assemblies of related species are available,
   flank orthology as in CLsat step 04.
5. **Boundaries.** CLsat's seed expansion (3-monomer steps, 3-copy tandem query, until the signal is lost) gives the array
   ends; this replaces the SINE border loop for these loci.

Output: `results/satellites/<family>.unit.fa`, `<family>.loci.tsv`, `<family>.monomers.aln.fa`, a short HTML section in the
report ("Satellites") with the table, the monomer structure drawn against the SINE consensus, and a histogram of monomers
per locus. TRF is used only here, on the loci, with the period raised to its 2 000 bp maximum; for units above that the
hit spacing is the period.

### 4.4 What stays in SINEderella proper

* The flank scan keeps its `satellite` class and TRF on 1 kb windows for the short-period cases around individual copies
  (A tails, dimers, small case-A nuclei with 2–3 monomers that the screen's ≥ 4 rule does not take).
* `array_flag.py` stays as the post-assignment check (it catches what the screen missed or a run with the screen off) and
  the report keeps the "Tandem array" chip.
* The user decides the family's status: dispersed SINE with a satellite derivative (snakes: Squam3 and sSat3 coexist),
  satellite only (rsi MEG-RS, if its 123 dispersed copies do not hold up), or a SINE that happens to sit in a segmental
  duplication (the `SQ2_NEXT_SESSION_PROPOSAL.md` table separates those by contig distance and recurrence; the screen
  reports distance so that the user can).

## 5. Test plan (separate, before any wiring)

Toy first (`tests/`): plant in a synthetic genome (a) a case-A array of 30 monomers of the 3′ 120 bp of a 250 bp toy SINE
with a leading full copy, (b) a case-B array of 20 units of 1 800 bp each containing one full toy SINE, (c) 200 dispersed
copies, (d) a segmental duplication of a region with 3 copies, (e) a dimer. The screen must flag a and b, report their
monomer length and SINE positions within 5 bp, leave c, d and e alone, and the excluded loci plus the remaining hits must
equal the original hit set.

Then the two real cases:

| Case | Genome | Family | Expected |
|---|---|---|---|
| B | *R. sinicus* (`~/rhin/genomes/rsi.fna`, run `~/rhin/rsi_fresh2/run_20261002_112935`) | MEG-RS (135 bp) | SAT_B, 21 arrays, unit 2 155 bp (variants 1.3–4.7 kb), 107 bp MEG at unit 2–108, 99.9 % unit identity; 123 copies stay for the SINE analysis |
| B | *Taphozous* (tbr, `~/chiro/`) | MEG-RS | SAT_B with ~910 bp unit (seen 2026-09-28 on the plates) |
| A | a snake, e.g. Indian cobra NN_10x_BNG, or Gekko japonicus V1.1 | Squam3 (consensus in `C:\work\glm-harness\sine-sequences\squam3_consensus.fasta`; subfamily Squam3C for snakes) | SAT_A: sSat3 monomer 67 bp = Squam3C positions 42–108 (cobra), thousands of loci, up to 320 monomers, 89 % mean monomer identity; Squam3 itself stays a dispersed family |
| A | Gekko japonicus | Squam2 | sSat2: 139 bp monomer (head + one third of the body), 79 % of loci |
| negative | *Scalopus*, *Timema* runs | all | no flag |

The squamate genomes are not on therioserver as far as my notes say (the sq2 work ran on "monsoon"); the accessions are in
the 2023 paper's Methods. Compare the screen's loci with the paper's counts per species before trusting it.

## 5b. First real test (2026-10-03, KIT, *Gekko japonicus*) — `tools/satellite_screen.py`

Indication only (hit coordinates, no alignment): 9 toy tests pass (planted A and B arrays, trimer/dimer and dispersed controls, dense
dispersed background, a real array inside a dense background). On the gecko genome (results in Tal `satellites/gja_test/`):

* Existing hit BEDs (80 % length rule applied): kind A 0.0-0.1 %, kind B excess over the chance null 7.5-10 %: nothing flagged, and
  nothing expected, since sSat monomers are shorter than the length threshold.
* Squam3A searched again down to 20 % of the consensus length (`sear ... 0.2 65 0`, as in the 2023 paper): 430 322 hits, **270 monomer
  runs** (>= 4 monomers < 100 bp apart; median monomer 141 bp, the size range the paper gives for gecko sSat3; longest runs 56 / 22 / 10
  monomers; 194 of the 270 have exactly 4) holding **1 289 hits = 0.3 % of the family's hits**. Kind B: 29.4 % in regular runs but
  28.3 % expected by chance (one hit per 6 kb): excess 1.1 %, not flagged.
* Consequences. (1) Kind A is found in seconds by coordinates alone. (2) A family-share threshold (>= 20 %) can never flag a family
  like Squam3A in gecko: the satellite loci are a fraction of a percent of its hits, yet they are real loci that should leave the SINE
  analysis. Exclusion for kind A therefore has to act on **loci**; the family share is only meaningful for kind B (rsi MEG-RS 92.8 %).
  (3) The null for kind B is essential in hit-rich families (without it 11-29 % of hits look "regularly spaced").
  (4) `sear` merges and filters its hits; a final version should read the raw `ssearch36` hits, and the monomer type (which part of the
  SINE the monomer covers; the user's `coordinates_by_consensus.sh` classes alpha-epsilon) needs the query coordinates of each hit.

## 5c. Second real test (2026-10-03, KIT, Indian cobra, Squam3C) and the revision it forces

Ground truth by Tandem Repeats Finder over the whole genome (44 parts, ~9 min on 44 cores): **2 580 Squam3C-derived tandem loci**, 21 799
monomers, 1.59 Mb; period 67 bp in 1 750 (the sSat3 unit of snakes; dimers of 133-134 bp in 113); 4-9 copies in 2 096 loci, 10-19 in 388,
20-49 in 79, 50-99 in 13, up to 559 monomers in one. The paper's "thousands of loci, up to ~320 monomers" is reproduced in order of
magnitude. (Squam3C used here is the reconstruction from the paper's Additional File 2; its positions 42-108 are the 67 nt unit with box B.)

The hit-spacing screen on `sear` hits found 156 runs (1.5 % of hits) and overlaps only 186 of the 2 580 loci (7 %), although 88 % of the
loci carry some `sear` hit. Cause: `sear` merges neighbouring monomers into one hit (hits 69-76 bp long at a spacing of 134 bp: two monomers
per hit) and most loci are short (4-9 monomers), so "4 consecutive hits" rarely holds. **Hit spacing alone is not a reliable kind-A detector.**

Revised design of the indication stage (cheap first, TRF only where it can pay off):
1. **Gate (seconds, hits only):** windows = hit clusters (hits within 300 bp, any hit length, both strands) of each consensus; the number and
   total size of windows is small (cobra: 54 269 hits = ~11 Mb of windows against 1.8 Gb). A consensus with no cluster of >= 2 hits has
   no kind-A candidates and costs nothing more. Kind B (regular spacing with the chance null) stays as it is: it worked on rsi MEG-RS.
2. **Verify (minutes, TRF on the windows only):** `trf windows.fa 2 5 7 80 10 40 300`; a locus = a TRF record with >= 4 copies whose unit aligns to
   the consensus (or to its doubled sequence). Gives the monomer length, copy number and TRF match per locus; no genome-wide TRF.
3. **Classify (characterisation):** align the unit to the SINE consensus: positions covered (cobra: 42-108 of Squam3C), full vs partial SINE
   in the monomer, leading/trailing monomers carrying the SINE terminus.
Cost check: the genome-wide TRF used as ground truth took ~9 min on 44 cores; the windowed version is two orders of magnitude smaller.
Open point: `sear` is not the right hit source for the final tool (merging); use raw `ssearch36` hits from step 1.

## 5d. The gate + windowed TRF stage works (`tools/satellite_trf_verify.py`, 2026-10-03)

Implemented as designed in 5c; defaults: a window (+-1 kb, merged) around every hit, TRF on the windows only (period >= 55, >= 4 copies),
unit aligned to the consensus (>= 45 bp), best record per locus; `--max-window-mb 600` falls back to clusters of >= 2 hits for hit-rich
families (and says so). Results on KIT (details and files: Tal `satellites/gate_trf_test/`):

* **Cobra / Squam3C:** 2 052 loci, 15 949 monomers, period 67 bp, SINE part 30-110, in 29 s over 103 Mb of windows (6 % of the genome);
  **89.7 % of the 2 580 loci found by genome-wide TRF** (90 % of the short loci with 4-9 monomers, 85 % of those with >= 20).
* **Gecko / Squam3A:** 537 loci, 2 865 monomers, periods 114-126 bp, in 21 s (fell back to >= 2-hit clusters, 112 Mb of windows).
* **Toy genome:** both planted satellites with the planted SINE part, no false loci from dispersed copies, dimers, trimers or a segmental duplication.
* Lessons: 100 bp of padding recovered only 60 % (partial windows cut arrays); 1 kb recovers 90 %. Misses are loci without any nearby hit
  (`sear` drops diverged monomers): the final tool should use raw `ssearch36` hits.
* Kind B (units longer than the SINE) stays with `satellite_screen.py` (regular spacing + chance null); a 20-unit toy array was 7 % of the hits
  and the screen reported it as a locus run although the family-share flag stayed off: **report loci, not only a family flag**.

Not done yet: monomer classification beyond the aligned SINE part, exclusion of the loci from the hit set, wiring into SINEderella (step
between search and assignment), peel input, the report section and the viewer, borderline cases on real data (dimers, composites, SINEs
without TSDs), raw-hit input, a snake and a gecko run against the paper's per-species counts.

## 5e. Wired into SINEderella (2026-10-03)

* `sear` now also writes `gen-<q>.rawhits.tsv` beside `gen-<q>.bed`: the merged hits at the homology cut with **no length rule**
  (partial copies, satellite monomers), in genome coordinates and with the original contig names. The filtered `gen-<q>.bed` is unchanged.
* `step1_search_extract.sh` runs `tools/satellite_stage.py` right after the searches, before `all_hits.labeled.bed`, `merged_hits.bed`,
  extraction, sampling and SubFam: kind A (`satellite_trf_verify`) on the raw hits, kind B (`satellite_screen`) on the full-length hits.
  Output in `<step1>/satellites/` (linked as `results/satellites/`): `indication.tsv` (per consensus), `loci.bed` (kind, consensus, locus,
  period or unit, monomers, SINE part), `units.fa`, `<q>.kindA.loci.tsv`, `excluded_hits.bed`.
* **Exclusion** (default on): hits of `gen-<q>.bed` overlapping a kind-A locus, or a kind-B run of a consensus flagged SAT_B, are removed
  (original kept as `gen-<q>.bed.before_satellites`), so extraction, the peel input (SubFam) and assignment never see them. The consensus
  stays in the bank. Env: `SKIP_SATELLITES=1` (stage off), `SATELLITE_EXCLUDE=0` (tables only), `SATELLITE_EXCLUDE_B=flagged|all|none`.
  The stage is skipped with a log line when `trf` is not on PATH; a failure inside it never stops step 1 (hits stay unfiltered, warning).
* The report's Similarity section gets a "Satellites" table (`report_blocks._satellites`).
* `SINEderella` copies the four tool files into the run dir (`tools/`) and exports `SINEDERELLA_TOOLS`, so a run is reproducible from its own copy.

## 5f. Bats, and the answer to "family share or every locus?" (2026-10-03)

Stage run on the existing tbr and rle runs (tables only; Tal `satellites/bats_test/`):
* **tbr MEG-RS**: 286 hits, two arrays of 81 and 109 units (926 / 912 bp): excess 66 %, SAT_B. The ~910 bp period seen on the plates.
* **rle MEG-RS**: 19 470 hits, excess 6.9 %, **not** SAT_B: the family is dispersed (it is the positive control of the length-version test),
  yet 35 regular runs, 26 of them with >= 20 units (9 with >= 50), one of 36 units with the **same 2 145 bp unit as the rsi array**; the chance
  null gives no run at all at this hit density. rle MEG-TR: 28 runs, 14 of 20-49 units, none expected. tbr VES (641 631 hits, one per 3 kb):
  chance alone gives 37 000 runs of 5-9 units and 180 of 20-49, observed 1 939 of 20-49 and 221 of >= 50 against 0 expected.

So the family-share flag is right for tbr and would leave the real rle arrays in the SINE analysis. Rule now in the stage (`--exclude-b long`;
**`flagged` stays the default**, see below): per consensus, the run-length distribution of the chance null is computed with the shares, and the **long minimum** is the smallest
run length L for which chance explains fewer than max(1, 5 %) of the observed runs with >= L units (rle MEG-RS: 5, tbr VES: 50; the same
calibration idea as the TSD minimum in flankscan stage 8). Runs at least that long are satellite loci on their own and are excluded in any
family; all runs of a SAT_B family are excluded as before. `flagged`, `all` and `none` remain as options; both modes stay testable.
Reported per consensus: `kindB_long_min`, `kindB_long_runs`; the report says "N long arrays inside a dispersed family".

*Correction after the rsi run (5h):* in rsi the calibrated minimum came out as 5 units for most r-families (45 000-83 000 hits each, one per
25-45 kb, so chance gives almost no regular run), and `long` would have excluded 65-167 five-to-six-unit clusters per family: local clusters of
ordinary copies, not satellites. Kind-B runs are geometric only; until the units are checked for sequence similarity (todo: align consecutive
units, accept >= 90 %) `long` needs a floor of 10 units (`--long-min-units`) and the default is `flagged`.

## 5g. Viewer (2026-10-03)

Tal `satellites/viewer.html` (https://toki-bio.github.io/Tal/satellites/viewer.html): loads a run's `results/satellites/loci.bed` (+ `indication.tsv`)
or a `*.loci.tsv` of `satellite_trf_verify.py`, and shows summary cards, the indication table, **where the monomers sit on the SINE** (loci per
consensus position, from the aligned SINE part), a monomers-per-locus histogram and a sortable, filterable locus table (kind, consensus,
min monomers, contig). Checked on the cobra (block at 34-110 of Squam3C) and rle (MEG-RS arrays of ~2 100, ~880 and ~1 500 bp units) data.
Not yet in it: per-locus monomer alignments and the per-locus analysis of the 2023 paper (needs the characterisation step, section 4.3).

## 5h. The two rules under test, settled by data (2026-10-04)

**Kind-B exclusion.** `tools/satellite_kindB_verify.py` cuts the units of every regularly spaced run from hit start to hit start and aligns
each unit with the next one and the one after (the rsi MEG-RS array is dimeric: 2 155 and 1 460 bp units alternate and share only ~800 bp,
so neighbours alone read as unrelated). Verdict ARRAY when the median unit identity is >= 85 %. Results: rsi MEG-RS 18 of 21 runs ARRAY at
92-99.8 % (the 3 others 88.7-89.3 %, arrays too under the 85 % rule), MEG-TR 1/1 (a 5.8 kb unit holding both heads), a 26-unit 2.6 kb array
at NC_142508.1:42.84-42.93 Mb carrying r9/r4/r7/r10 heads, and **1 850 of 1 862 five-to-nine-unit runs of the r-families are ordinary copies**
(25-55 % unit identity); rle MEG-RS 26 of 35 ARRAY (all 25 of the >= 20-unit runs), MEG-TR 28/28. So neither the family share (misses rle) nor
the run length (would take 1 850 clusters of copies in rsi) is the criterion: the unit sequence is. The stage now runs the unit check on every
run and excludes verified arrays (`--exclude-b verified`, default; the SAT_B share flag stays as information). Cost: a few minutes for rsi.

**Attribution of shared kind-A loci.** Per-consensus identity, monomer coverage and consensus coverage were tabulated for every shared locus
(`rsi_sat/satellites/attribution_test.tsv`). For the 140 bp x 58 locus: C11 95.2 % identity but only 45 % of C11 covered; r2 92.7 % with 93 %
of r2 covered; r1 92.1 % but only 56 % of the monomer; r4 84.9 %. The alignment-length rule picked C11 (the longest composite) for 16 of 17
shared loci. Rule now: among consensuses within 3 identity points of the best and covered by the monomer over >= 90 % of the alignment, the
consensus of which the monomer covers the largest share; ties -> the shorter consensus. The 140 bp satellite is then an **r2-derived satellite**
(the whole r2, 136 of 146 bp), the 136 bp x 4.6 locus likewise r2, the 59 bp x 4 locus r3. A first version of the rule read the wrong column
(score as identity) and changed nothing; found by recomputing the pick by hand on the per-consensus records.

**rsi_sat with both rules (tables v5, `rsi_sat/satellites/`):** 34 kind-A loci: 17 P18-derived (61-124 bp monomers of the r10-r8 junction
region), 4 r2 (the 136-143 bp whole-r2 satellites, 58 + 3 x 4-5 copies), 3 r3, 3 C11, 3 MEG-T2, 1 each r8, r7, P1, r10. Kind B: 42 verified
arrays (MEG-RS 21, MEG-TR 1, r9 15, r4 2, r10, r7, P26 1 each) against 1 862 runs that are copies; verified mode would remove 2 243 hits
(the run itself, in flagged mode, removed 2 066: the 177 more are the small r-family arrays, 92 of them r7, 176 r9, 44 r4).

## 5i. Speed of the kind-B unit check (2026-10-06)

**What was slow.** The focused rsi run of 2026-10-05 (`run_20261005_085343`, code 7f32adb) spent **13 h 24 min** in the satellite stage; the
run of the night before, without the 30 kb array tier, 23 min. All of it was the kind-B unit check; kind A (windowed TRF) took seconds per
consensus. TRF cannot take over kind B: TRF 4.10.0-rc.2 refuses a maximum period above 2 000 bp ("Value must be between 1 and 2000"), and
toy arrays of 2 500 and 7 000 bp units give no record, while the rsi MEG-RS unit is 2 155 bp and the wide tier allows 30 kb. Three costs
multiplied: (1) the 30 kb tier made 733 long runs (36 before) with units up to 30 kb, and Smith-Waterman grows with the product of the
lengths; (2) the check aligned every unit with every unit of its run (132 alignments for 12 units, of which the verdict uses the 21 at
lags 1 and 2), one run after the other; (3) with statistics on (`-z 11`, and equally `-z 1`, `-z 0` or the default) ssearch36 goes on
searching alternative, non-overlapping local alignments of each pair: one real 26 x 20 kb array pair gave 24 lines in 56 s, and with
`-z -1` (no statistics) 2 lines in 5.5 s with the same best alignment (99.92 % over 13 445 bp). Times per pair, -z 11 vs -z -1: 56 vs
5.5 s, 1.9 vs 0.13 s, 0.36 vs 0.03 s, 1.3 vs 0.09 s.

**What changed (`tools/satellite_kindB_verify.py`).** `verify_runs` aligns only the lag-1 and lag-2 pairs, each query unit against its
partners in its own `ssearch36 -z -1 -T 1` call, the calls in parallel (`--threads` processes), and keeps every pair it has aligned (a
narrow run inside the wide run that contains it, and repeated calls, align nothing new). The hits of a run are found by bisection instead
of a scan of the contig. The previous implementation stays as `verify_runs_serial` (`--serial`, or `SATELLITE_KINDB_SERIAL=1` for the
stage). One known difference: without the alternative alignments an unrelated pair can score lower (one rsi pair: 85.6 % over 215 bp
instead of 89.6 % over 303 bp), so COPIES identities can come out lower; the best alignment, which is what a real array's identity comes
from, is the same.

**Checks.** Toy (`tests/test_satellite_kindB_verify.py`, needs ssearch36): planted single-unit, dimeric (MEG-RS-like 2 155 / 1 460 bp),
7 kb-unit and near-threshold arrays, a diverged array and a run of ordinary copies give the expected verdicts with both implementations,
ARRAY identities within 0.5 point; no pair is aligned twice. Real data: all 2 739 kind-B runs of the 13-h rsi run (therioserver
`~/tmp/kbspeed/val/validate_rsi.tsv`): **9.5 min with 32 processes instead of 13 h 24 min; 0 verdicts changed (47 ARRAY, 2 692 COPIES);
the 47 ARRAY identities identical**; COPIES identities identical in 1 753 runs, the others -40 to +8 points (all below 85).

End to end (`~/tmp/kbspeed/e2e_7f32/`): the whole stage on copies of that run's step-1 inputs, with its own tool copies (7f32adb) and only
the verifier replaced, **718 s instead of 13 h 24 min**. Every kind-B row and verdict reproduces. The one row that differs is a kind-A
locus (r10_r8_P18, NC_142516.1:2 177 953-2 178 498, period 71): its 69 bp unit aligns to P18 over exactly 45 bp (`MIN_ALN` = 45) at 71 %,
and its E-value straddles `MAX_E` = 0.01 under the random shuffles of `-z 11` (10 repeats: 0.023-0.039; the 13-h run drew one <= 0.01).
Kind A was not changed; a locus that sits on both cut-offs is in or out by chance from run to run (with it come 1 excluded hit, 1 region,
the r10_r8_P18 indication row and 1 unit in `units.fa`).

The current code (2b4d35a: ed86b68's narrow runs verified as well) on the same inputs (`~/tmp/kbspeed/e2e_head/`): **717 s**; the narrow
runs cost almost nothing more, their unit pairs are those of the wide runs that contain them. It verifies 11 narrow arrays inside wide runs
that fail, among them the **2 168 bp MEG-RS array on NC_142509.1 (134 hits, 89 %)** that the 7f32adb run judged COPIES inside a wide run
(the mini-genome case of ed86b68, now on the whole genome); excluded hits 2 404 -> 2 543. The kind B of a run grows with `--threads`
(`THREADS`): every ssearch36 call is single-threaded and independent.

In a real orchestrated run (rsi, `d2d9e6a`, THREADS 32 on 32 pinned CPUs, therioserver `~/tmp/kbspeed/rsi_fast/run_20261006_000626`)
the stage ran **00:25:07-00:36:40, 11 min 33 s**; 2 543 hits excluded, the 2 168 bp MEG-RS array verified as a narrow run (ARRAY 89 %)
inside the wide run it was joined to (COPIES 82 %).

**Not done: a k-mer screen in front of ssearch36.** Skipping a pair is safe only when its identity is known to be below 70 (a median of 85
needs both middle values >= 2 x 85 - 100). On 8 005 real pairs (all 47 ARRAY runs and 700 random COPIES runs of rsi; therioserver
`~/tmp/kbspeed/calib.tsv`) the share of the shorter unit covered by shared k-mers does not give that: 31 of the 625 pairs scoring >= 70
have < 1 % of their positions in shared 16-mers, 7 share no 12-mer at all (scores 70.8-73.4, units 91-146 bp). At 70-85 % identity spread
along the unit the mismatches are 3-7 bp apart and no 12-20-mer survives, and an old, evenly diverged array with long units would fall
through the same way. A screen at 5 % cover would have skipped 68 % of the pairs (76 % of the Smith-Waterman cost), about 2x on top of
the above, without a guarantee; left out.

## 6. Decisions and requirements from the user (2026-10-03)

* **The SINE inside the satellite must still be detected and reported properly, and clearly separated from the
  satellite.** The screen removes satellite *loci* from the SINE analysis; it must not lose the SINE sequence they carry.
  For every flagged family the report states which SINE (or which part of it) the monomer contains, with coordinates on
  the SINE consensus, and the monomer consensus itself is run through SINEderella as a candidate ("can this monomer be a
  SINE in its own right?") — SINEderella handles shorter versions of an element, so a monomer that is a truncated SINE is
  a legitimate, testable hypothesis, not noise.
* **Borderline cases are the hard part**: dimers, composites and other tandems of two or three units that may be SINEs,
  and SINEs that lack the usual features (Squam2 has no proper TSDs). Rule: the screen excludes only runs of ≥ 4 monomers
  (the 2023 threshold); 2–3-unit tandems, dimers and composites stay in the SINE analysis and go through the flank scan
  (`docs/COMPOSITES.md`), marked, never dropped. Absence of a TSD is never a criterion for calling something a satellite.
* **Exclude by the family's share (≥ 20 %) or every flagged locus?** To be tested on both kinds (snake sSat3 next to a
  live Squam3 family; rsi MEG-RS). Implement both as options and compare the resulting SINE counts and plates.
* **Case-A monomer runs out of the peel input too?** Yes — but test first whether the monomer can be a SINE (previous
  point); the peel must still see a monomer that turns out to be a short SINE version.
* **Shorter vs longer monomer (A/B)?** The user's framing: the SINE is **fully or partially** inside the satellite
  monomer. A and B are then two ends of one range (B = full SINE inside a longer unit; A = part of a SINE as the monomer),
  and the detector pair (hit runs < 100 bp apart; regular spacing of full hits) covers the range. The report says, per
  family, how much of the SINE the monomer holds and how much of the monomer is SINE.
* **Characterisation depth:** first counts and classification with coordinates (indication table, loci BED, monomer
  consensus, SINE positions); then, finalising, the per-locus analysis of the 2023 paper and a satellite viewer like the
  CLsat viewer (`CLsat_workflow/viewer`).
* **Compute:** the squamate tests need the sq2 genomes on the monsoon server; see the session notes for the access
  question. The two test genomes (cobra NN_10x_BNG, Gekko_japonicus_V1.1) can also be downloaded to therioserver for the
  screen test alone.

## 7. Sub-task brief

Goal: implement sections 4.1–4.3 as `tools/satellite_screen.py` (indication + exclusion, standard library only) and
`tools/satellite_characterise.py` (needs mafft, ssearch36, TRF), wired as a stage between step 1 and step 2 with the env
vars above, with the toy test and the two real tests of section 5, and a `Satellites` section in the report. Read first:
this file, `docs/ARRAYS.md`, `flankscan/fs2_trf.sh`, `tools/array_order.py`, MANUAL §16.9, the 2023 paper (PMC10702118,
doi 10.1186/s13100-023-00309-2), CLsat_workflow `PIPELINE.md` and `02_seed_expansion_scan/clisat_expand_then_scan_merged.sh`.
The rsi unit consensus is in Tal `rsi_final/MEG-RS_array_unit_2155bp.fa`; the rsi arrays are in `rsi_final/array_flag.tsv`.
Do not change assignment, the flank scan or the plates; the screen only removes loci from the hit set and writes tables.

## References

Vassetzky NS, Kosushkin SA, Ryskov AP. SINE-derived satellites in scaled reptiles. Mob DNA. 2023;14:21.
doi:10.1186/s13100-023-00309-2 (PMC10702118). · Toki-bio/CLsat_workflow (CLsat/DarSat pipeline, 2026). ·
Gogolevsky KP, Vassetzky NS, Kramerov DA. 5S rRNA-derived and tRNA-derived SINEs in fruit bats. Genomics. 2009;93:494–500.
