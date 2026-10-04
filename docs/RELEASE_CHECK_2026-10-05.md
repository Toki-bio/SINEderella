# Release check 2026-10-05: end-to-end runs from a fresh clone

Purpose: before SINEderella is published, run it from a fresh `git clone` on cases with an outside answer, then on the Tal-page
species whose answers are the project's own earlier runs. Everything here ran on therioserver under `~/tmp/release_check/`
(`SD/` = the clone, `data/` = inputs, `hs21/`, `mm19/`, `tal/<species>/` = runs). Scoring script: `~/tmp/rc_compare.py`
(RepeatMasker comparison; copied below in section 5).

## 1. Install path

`git clone https://github.com/Toki-bio/SINEderella.git` at `d61d822` → every entry point executable (the clone of the morning had none,
fixed as D11 the same day), conda env `sinederella` (ssearch36, mafft 7.526, bedtools, samtools 1.24, seqkit, trf, gawk, EMBOSS,
numpy/matplotlib). `python -m unittest discover -s tests -p 'test_*.py'`: 62 tests pass.

## 2. Human chr21 (hg38, 46.7 Mb) against the UCSC RepeatMasker table

Bank: AluJb (312 bp), AluSx (312), AluY (311) from Dfam (DF000000007, DF000000047, DF000000002). RepeatMasker: 13,495 Alu on chr21,
of them AluJ* 3,208, AluS* 7,865, AluY* 1,884 (the rest AluYa5/Yb8/..., FLAM/FRAM, counted as Alu but not mapped to a bank member).
Overlap rule: a RepeatMasker locus is recalled when a SINEderella locus overlaps it by ≥ 50 bp.

| Run | Time | Firm loci | Recall (any extracted / firm) | Precision (firm loci on a target Alu) | Subfamily agreement (firm) |
|---|---|---|---|---|---|
| full, AluJb+AluSx | 20 min | 10,393 | AluJb 69 / 67 %, AluSx 90 / 89 % | 84.6 % (1,305 of the rest are AluY copies, not in the bank) | AluJb 98.8 %, AluSx 99.4 % |
| `--add` AluY | 10 min | 9,983 | AluJb 70 / 66 %, AluSx 90 / 84 %, AluY 90 / 90 % | 99.8 % (the rest: 15 SVA, 7SL, FLAM) | AluJb 98.9 %, AluSx 99.7 %, AluY 89.5 % |
| `--exclude` AluJb | 2 min | 10,173 | AluSx 90 / 87 %, AluY 90 / 89 % | 80.8 % | AluSx 100 %, AluY 92.7 % |
| `--resume` on the add run | 1 min | — | skipped steps 1, 2+3, 4 as complete; results rebuilt | — | — |
| `publish_run.sh` on the add run | 2.5 min | — | 9 plates, `proposals.tsv`, report with verdict columns; `DISC=~/SINE_discriminator/site` found; border loop: none flagged | — | — |

What the numbers say:

* **Recall is bounded by the length rule, not by divergence.** The missed RepeatMasker copies have the same divergence as the recalled
  ones (AluJb 16.9 vs 16.7 %, AluSx 11.5 vs 11.1 %) but a median length of 140–153 bp against 293–299: they are fragments below the
  80 % length requirement of `sear`. 97–99 % of the missed copies are shorter than 240 bp. A user who wants fragments must lower
  `0.8` in the sear call; the manual should say so (it does now, §6).
* **Assignment is "closest consensus in the bank".** With AluY absent, 1,305 AluY copies were extracted through AluSx and assigned to
  it (precision 84.6 %); adding AluY moved them (precision 99.8 %, AluSx firm recall 89 → 84 %). With AluJb excluded, ~2,000 AluJ
  copies were firmly assigned to AluSx (unanimous 10/10, above the 0.45 threshold: AluJ is 85–90 % identical to AluSx). A firm
  assignment states which bank member is closest, not that the copy belongs to that subfamily; a missing subfamily is absorbed by its
  nearest relative. The similarity ratio and the length-version test are the tools that show it.
* **Subfamily agreement** between SINEderella's firm label and RepeatMasker's name is 99–100 % for AluJb and AluSx; AluY 89.5 % (190 of
  1,836 firm AluY loci carry a RepeatMasker AluS* name, mostly AluSx/AluSc copies near the Y/S boundary).
* `--add`, `--exclude`, `--resume` and the publish flow completed without error; the only warning in all runs is the assembly-quality
  one (chr21 has N-rich edges). Two defects were found by looking at the outputs rather than the exit codes: the resume rebuilt
  `results/` without the satellites link (D25) and, worse, deleted the publish alignments and the report while doing so (D26); both
  fixed and re-verified on this run the same night (`--resume`, then a report-only publish: Satellites table present, plates intact).
* Satellite stage on Alu (D24): every Alu consensus was flagged `SAT_B` by the uniform-null share test (20.4–20.6 % excess) with 0 of
  ~550 regularly spaced runs verified by unit identity, 0 hits excluded. Alu is clustered in GC-rich isochores; the share flag is not
  evidence of arrays and the report now says so. Kind A found 2–6 small Alu-derived tandem loci per consensus (10–34 monomers).
* Consensus audit on the add run: AluSx and AluY rebuild within 4 mismatches of the Dfam consensus (MATCH); AluJb rebuilds DIVERGED
  (its chr21 copies are 17 % diverged and the pool holds AluJo/AluJr copies labelled AluJb), as the audit is meant to show.
* Length-version test: no candidate pair (the three Alus are the same length), so no table.

## 3. Mouse chr19 (mm39, 61.4 Mb), B1_Mus1 (148 bp) and B2_Mm2 (195 bp)

RepeatMasker: 14,822 B1 (repFamily Alu) and 9,412 B2 (repFamily B2) loci. Full run 15 min.

| | B1_Mus1 | B2_Mm2 |
|---|---|---|
| firm assigned | 6,402 | 3,049 |
| recall any / firm | 55.4 / 43.3 % | 72.0 / 33.4 % |
| precision (firm loci on a target) | 99.7 % overall | |
| subfamily agreement | 99.5 % | 100 % |
| missed copies: median divergence, length | 23.9 %, 100 bp (recalled: 13.8 %, 141) | 23.0 %, 114 bp (recalled: 19.1 %, 187) |

Recall is low for two reasons that the design states: short copies (the length rule again: 6,607 of 6,608 missed B1 are under 240 bp,
most under 120) and, for B2, the bank holds one of several mouse B2 subfamilies (B2_Mm1a, B2_Mm1t, B3, B3A are 75–90 % identical to
B2_Mm2): their copies are extracted (72 %) but fail the 0.45 bitscore threshold against B2_Mm2 (firm 33 %). The non-target extracted
loci are B4 (741), ID_B1 (741 as B4:ID_B1), RSINE1 and B4A: B4 is a B1–ID composite that carries a B1-like part, found by the B1 search.
The run is correct; the bank is incomplete for the question "all B1/B2 copies of the chromosome".

## 4. Tal-page species (rsi, tbr, rle)

Running at the time of writing (`~/tmp/release_check/tal/`): rsi with the 21-consensus bank (reference `~/rhin/rsi_fresh2/run_20261002_112935`),
then tbr and rle with `chiro_bank.fa` (references: the 2026-09-27 `run_add_*` runs). Results are appended here when they finish.

## 5. Scoring

`rc_compare.py RUN rmsk.tsv MAP` (MAP: `name=AluJ:AluJb,...` or `fam=B2:B2_Mm2`): RepeatMasker SINE loci of the mapped families as
bed; SINEderella loci from `results/all_sines.bedlike.ALL.tsv` (firm = `status=assigned`); recall by ≥ 50 bp overlap, precision by
overlap with any mapped RepeatMasker SINE, subfamily confusion by the RepeatMasker name with the largest overlap, and the divergence
(milliDiv) and length of recalled vs missed copies. The script is in the session scratchpad and on therioserver (`~/tmp/rc_compare.py`);
it is not part of the repository because it needs the UCSC tables.
