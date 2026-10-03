# SINE-derived satellites and SINE-containing tandem arrays: a screen before the SINE analysis

Design note for a separate sub-task (2026-10-03). Nothing in this document is implemented yet, apart from the
tandem-array flag described under "What SINEderella does today". The sub-task brief is at the end.

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

## 6. Open questions for the user

* Threshold for exclusion: 20 % of a family's copies (as the array flag), or any flagged locus regardless of the family's
  share? (A snake genome has thousands of sSat3 loci next to a healthy Squam3 family: the loci should go, the family stays.)
* Should case-A monomer runs be removed from the **peel** input as well (they would otherwise appear as long chunk rows)?
* Is "monomer shorter than the SINE" (A) versus "unit longer" (B) the right split, or is the real distinction whether the
  SINE sequence is the whole monomer or a part of a larger unit?
* Characterisation depth: the per-locus analysis of the 2023 paper, or only the monomer and counts?

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
