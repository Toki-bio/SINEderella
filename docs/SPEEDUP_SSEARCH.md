# Faster ssearch36-based steps, same results (2026-10-24)

Branch `speedup-ssearch` (from main 8d4961c). Not merged: the owner decides. Every number below was measured on therioserver
(128 cores, shared with the 22-bat rerun queue, which kept the load between 30 and 900 during the work; my jobs stayed at or
under 32 threads). Raw outputs: therioserver `/home/toki/ssspeed/` (section 9).

Rule used: a change is kept only if its outputs are identical to the old code's (byte for byte, or on every downstream
verdict and assignment). Where the old code is itself stochastic (step 1 E-values, the step 2 vote) the yardstick is the old
code's own run-to-run variation. Everything else is a proposal for the owner, off by default, or reverted.

## 1. Summary

| | old | new (default) | how |
|---|---|---|---|
| flank twins (fs9), tbr VES, 620,807 copies | 7 h 04 min, peak disk ~240 GB | ~25 min: 7.6 min to the body check (24 threads) + body check 17.2 min (4 workers, default -T); peak 2.0 GB; 9 min 53 s with the body-check proposal | seeds and seed pairs counted in memory instead of written as text and sorted; confirmation in ranges; union-find with path compression |
| flank twins (fs9_run), rsi, 22 families | 35.3 min | 9 min 19 s (4 body-check workers); 2 min 29 s with the body-check proposal | same |
| step 1 searches, rsi, 3 consensuses (2 runs each) | 259 / 238 s, peak 6.7 GB | 145 / 116 s, 3.3 GB | sear: each genome part was re-indexed after every search; window bank streamed into the parts |
| step 2 vote | unchanged | unchanged | parallel cycles changed the vote distribution: reverted (section 6) |

Outputs: identical in every comparison (section 4). The two big wins are not inside ssearch36: stage 9 wrote every pair of
copies that share a seed as a text row and sorted it (hits.p2.tsv 92 GB on tbr, 241 GB on ttr), and sear re-indexed every
155 MB genome part after every search although only its unused minus-bank mode reads that index.

Two findings about ssearch36 itself matter for reproducibility (section 5): its E-values come from clock-seeded shuffles (step
1 gives slightly different loci run to run), and for libraries under 500 sequences they depend on `-T`, i.e. on the number of
CPUs when `-T` is not given.

## 2. Where the time goes

Wall minutes per stage, from the timestamps in `pipeline.log` (therioserver `/home/toki/ssspeed/profile/stage_minutes.tsv`,
made by `/home/toki/ssspeed/stage_profile.py`). rsi THREADS=32; bats THREADS=24, five bats at a time.

| run | step1 | step2 | step3 + families | step4 | flank twins (fs9) | verdicts | publish | flankscan 1-7 | total |
|---|---|---|---|---|---|---|---|---|---|
| rsi (67,192 copies, 22 consensuses) | 29.6 | 7.6 | 4.2 | 6.4 | **35.3** | 7.7 | 10.6 | 16.7 | 128 |
| tbr (VES 620,807) | 21.6 | 7.3 | 8.4 | 0.8 | **426.1** | 0.7 | 3.2 | 8.8 | 482 |
| nle (VES 434,859) | 15.1 | 6.1 | 6.9 | 0.9 | **246.0** | 0.7 | 3.7 | 5.8 | 291 |
| mau (VES 425,267) | 18.9 | 4.7 | 6.3 | 1.0 | **188.2** | 0.2 | 2.6 | 5.5 | 232 |
| ntu (VES 344,170) | 17.6 | 4.4 | 7.1 | 0.8 | **170.9** | 0.5 | 2.9 | 4.1 | 213 |
| mev (VES 316,847) | 16.0 | 5.2 | 5.7 | 0.8 | **166.3** | 0.2 | 2.4 | 9.5 | 211 |
| cth (VES 91,807) | 11.1 | 0.8 | 7.0 | 0.7 | **45.2** | 0.1 | 5.9 | 0.7 | 76 |
| rle (MEG, 20k copies) | 6.5 | 0.3 | 0.7 | 1.0 | 15.2 | 0.1 | 2.2 | 1.0 | 31 |

Flank twins are 50-90 % of every bat run with a large VES family. Inside fs9 (cth VES, the old code alone, `bash -x`
timestamps, `/home/toki/ssspeed/bench/fs9_old_cth/`): repeat filter 228 s; pass 1 368 s (seed rows 72 s, pairs 68 s,
sort and count 96 s, confirmation 132 s); pass 2 526 s (confirmation 226 s); copy status 175 s; body check 43 s; 22 min 28 s
in all. tbr VES in the queue: both passes 3 h 57 min, copy status ~65 min, body check 2 h 01 min (serial calls; the parallel
body check 8d4961c came later that day).

Inside step 1 (rsi, `step1.stdout.log`): 22 searches 19.4 min, satellite screen 2.8 min, SubFam 2.3 min, MAFFT of 600
clusters 4.8 min. One search covers 20 genome parts. The ssearch36 call on one part takes 0.96 s at -T 32 (27 s CPU, 87 %
parallel efficiency), but after each part sear also ran `samtools faidx` + `sort -V` of that part: 1.4-2.0 s.

ssearch36 thread scaling (`/home/toki/ssspeed/bench/scale/`): a step-2 call (22 consensuses x 20,000 copies) 13.6 s at -T 32
(342 s CPU), 23.5 s at -T 16 (317 s), 40.1 s at -T 8 (303 s), 78.5 s at -T 4 (304 s); a step-1 call (1 consensus x 107,695
windows) 1.58 s at -T 32 (47 s CPU), 11.8 s at -T 4 (48 s).

### ssearch36 calls in the workflow

| where | call | library | -T | -z | order | share |
|---|---|---|---|---|---|---|
| sear (step 1) | 1 consensus x 1/20 of the genome in 3 kb windows, `-E 2 -m 8C` | ~55k windows | THREADS | 11 | 20 parts x consensuses in a row | ~1/3 of step 1 before the change |
| step2_asSINEment.sh | all consensuses x 20,000 copies, `-E 2`, 10 cycles | 20k | THREADS | 11 | cycles and parts in a row | 7.6 min rsi, 13.6 min ttr |
| step3_postprocess.sh | consensus x its copies, `-m 10` streamed into awk | per subfamily | THREADS | 11 | in a row | 9 s rsi |
| step4_plots.sh | consensus x <= 10,000 copies, then MAFFT --addfragments | per subfamily | THREADS | 11 | in a row | MAFFT and plots dominate |
| fs9_bodycheck.py | 1 copy x its 1-50 partners | tiny | none (= all CPUs) | default | parallel workers (8d4961c) | tbr 2 h serial |
| flankscan fs3, fs_lib (fold, ends) | consensuses x windows in chunks; candidate folds | chunks / tiny | T / none | 11 | in a row | fs3 193 s, fs5 491 s rsi |
| tools (kind-B, TRF verify, continuation, length variants, tRNA heads) | small | small | 1 or T | -1 / 11 / default | parallel where needed | minutes |

publish/ and disc/ (verdicts 7.7 min, publish 10.6 min on rsi) call no ssearch36.

## 3. Disk footprint

Largest files in the run directories (`find -size +200M`, `du`; read only).

| run | largest intermediate files | peak | at the end |
|---|---|---|---|
| ttr (VES 866,561), fs9 still running at 20:00 | hits.p2.tsv 241 GB, hits.p1.tsv 112 GB, uutils sort chunks 308 GB, keys.p1/p2 6.3 + 6.1 GB, flank20.tsv 2.5 GB | >= 659 GB so far | - |
| tbr (VES 620,807) | hits.p2 91.8 GB, hits.p1 45.3 GB, keys 3.7 + 3.5 GB, flank20 1.7 GB, sort chunks of the file being sorted | ~240 GB | 12 GB, of it flank_twins 3.2 GB |
| mau (VES 425,267) | flank20 1.3 GB (left in results/), flankscan windows.fa / windows.masked.fa / trf.raw 0.9 / 0.9 / 0.8 GB | (hits not kept) | 8.3 GB |
| cth (VES 91,807), old code alone | hits.p2 3.9 GB, hits.p1 2.8 GB, keys 0.6 + 0.6 GB, flank20 0.26 GB, sort chunks | ~12 GB | flank_twins 0.3 GB |
| rsi, step 1 (sear) | window bank genome.search.fa.2k.bnk 3.1 GB + its 20 parts 3.1 GB | 6.7 GB | 0 |

Why each big file existed:
* **hits.p1/p2.tsv** (fs9): one text row for every pair of copies that share a seed in the first 80 bp of a flank, then a full
  sort to bring one pair's seeds together. Every pair was written before anything was filtered.
* **keys.p1/p2.tsv** (fs9): every seed of every flank as a text row, sorted only to find equal seeds.
* **flank20.tsv** (fs9): every distinct flank 20-mer, sorted only to intersect it with the genome's k-mer count; never deleted.
* **cand / conf / pairs / seqs.raw** (fs9): intermediates no step reads, left in results/ (0.6 GB of cand files on tbr).
* **genome.search.fa.2k.bnk** (sear): the whole window bank written before splitting it into 20 parts, kept to the end.

After the changes, measured peak: tbr VES 2.0 GB, mau VES 0.87 GB, cth VES 0.33 GB, rsi fs9_run 31 MB (results/ with
assignment_full.tsv); results/flank_twins of rsi 238 MB -> 15 MB; rsi step-1 searches 3.3 GB. Not changed: flankscan windows.fa
and windows.masked.fa (read by stages 3, 5, 6, 8), trf.raw (5.3), genome.clean.fa (the run's genome).

## 4. Candidates

| # | change | where | measured | equivalence | verdict |
|---|---|---|---|---|---|
| 4.1 | seeds and seed pairs counted in memory per block of copies (new flankscan/fs9_seeds.py) instead of keys/hits files + sort | fs9_twins.sh steps 2-3 | cth: pass 1 236 s -> 14 s, pass 2 300 s -> 15 s (4 threads); mau: flanks (seqs.raw) -> cand.p1, repeat filter included, 50 min -> 3 min 46 s | cand.p1 / cand.p2 identical: cth 417,309 / 740,197 rows, mau 1,190,912 / 2,072,588, tbr 3,018,034 / 5,852,452 | accepted |
| 4.2 | repeat filter with the flank 20-mers held as 40-bit numbers (new flankscan/fs9_repmask.py), no flank20.tsv | fs9_twins.sh step 1 | cth 118 s -> 7.7 s (rep20), 111 s -> 7.5 s (mask) | rep20.tsv and the masked seqs.tsv identical (cth); every later file identical (cth, mau, tbr, rsi) | accepted |
| 4.3 | the confirmation awk over FS9_THREADS consecutive ranges of the candidates, outputs joined in order | fs9_twins.sh step 4 | tbr conf.p2 28 min -> 78 s (24 threads) | conf.p1 / conf.p2 identical (cth, mau, tbr: 544,636 / 1,674,446 rows) | accepted |
| 4.4 | union-find with path compression in the array and copy-status gawks | fs9_twins.sh | cth status 178 s -> 1.6 s | flank_groups.body and copy_status identical; the roots are unchanged by construction | accepted |
| 4.5 | delete the fs9 intermediates no step reads; FS9_KEEP=1 keeps them | fs9_run.sh | rsi flank_twins/ 238 MB -> 15 MB | grep: no reader outside fs9_twins.sh | accepted |
| 4.6 | sear: `samtools faidx` + `sort -V` of each part only when a minus bank is built (`--minus`; step 1 never asks for it) | sear | 1.4-2.0 s per part per consensus; rsi 3 consensuses 259 / 238 s -> 145 / 116 s | the removed commands wrote only `*_s.fai`, read only by the minus bank; loci differ old-new as much as old-old (below) | accepted |
| 4.7 | sear: window bank streamed into the 20 parts | sear | 3.1 GB less disk for rsi | the 20 parts byte-identical (cmp) | accepted |
| 4.8 | step 2: ten cycles at once (THREADS/10 threads per call; parts under 500 copies at -T THREADS) | step2_asSINEment.sh | rsi 497-574 s -> 404-435 s | the vote distribution changed (section 6) | **rejected, reverted** |
| 4.9 | body check: ssearch36 `-T 1` instead of no -T (= all CPUs) | fs9_bodycheck.py | 250 -> 35 ms CPU per call; tbr VES body check at 4 workers: default 1,031 s wall / 17,476 s CPU, -T 1 639 s / 2,493 s; -T 1 at 24 workers 2 min | verdicts identical in all 24 families tested; weak background pairs can differ (5.1) | owner: opt-in 2026-10-24, **default since 2026-10-26** (`BODYCHECK_SSEARCH_T=0` = no -T) |
| 4.10 | `-z 0` in the body check | fs9_bodycheck.py | no gain: it still shuffles (0.75 s per 20 calls, as -z 1) | - | rejected |
| 4.11 | `-z -1` (no statistics) in the body check | fs9_bodycheck.py | 0.14 vs 0.75 s per 20 calls | E-values become meaningless (1, 46) while `-E 100` decides which pairs count | rejected |
| 4.12 | prefilter (fasta36 ktup, blastn, minimap2) + ssearch36 only on candidates | sear, step 2 | not run | changes which copies are found (sear fast mode loses ~1 % of loci, docs/FAST_MODE.md); in step 2 -z 11 is the filter | rejected (rule 1) |
| 4.13 | other SW engines (parasail, SWIPE, opal) | anywhere | not run | the E-values come from ssearch36's own clock-seeded shuffles; no other engine reproduces them | rejected |
| 4.14 | another ssearch36 build | all | not run | FASTA36 has only an SSE2 kernel (no AVX2 path in the source), so a newer build gives no kernel gain; installed: bioconda fasta3 36.3.8 | not pursued |
| 4.15 | step 1: several consensuses at once in separate directories | step1_search_extract.sh | not run | after 4.6 a part call runs at 87 % efficiency at -T 32; sear's temporary names are not per query | not pursued |

Step-1 check for 4.6-4.7 (`/home/toki/ssspeed/bench/sear/`): two old and two new runs on rsi with MEG-RS, r1_r3_P1 and
r9_15seqs. The bit-score column differs in most rows between any two runs, old or new (clock-seeded -z 11). Loci
(contig, start, end, strand): MEG-RS 0 and r1_r3_P1 0 differences in every pair of runs; r9_15seqs (82,967 loci) 4 between
the two old runs, 0-6 between old and new runs.

## 5. Needs the owner's decision

### 5.1 ssearch36 `-T` changes E-values of small libraries

When a library gives fewer than `shuff_max` (500) scores, ssearch36 fits its statistics from shuffles of the library
(fasta36 comp_lib9.c line 1082), and the threaded build sizes the shuffle buffers `shuff_max / (2 x workers)` (compacc2.c
`buf_shuf_seq`). Measured on rsi MEG-T2 body-check background calls, 10 calls per setting (`/home/toki/ssspeed/bc_T.py`):

| pair (library size) | -T 1 | 2 | 4 | 8 | 16 | 32 | 44 | 64 | 128 = none |
|---|---|---|---|---|---|---|---|---|---|
| 72968061 vs 1807217 (14) | 10/10 E 14 | 10/10 | 10/10 | 4/10 | 3/10 | 0/10 | 0/10 | 0/10 | 0/10 |
| 1560800 vs 5106709 (12) | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 | 4/10 E 11 | 10/10 E 10-11 | 10/10 E 6 |
| 33962480 vs 91191978 (10) | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 | 1/10 | 8/10 | 10/10 E 4 |
| 76793976 vs 1807217 (7) | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 | 6/10 | 9/10 | 10/10 | 10/10 E 4 |

So which weak background pairs pass `-E 100` depends on `-T`, and with the current default (no -T) on the number of CPUs of
the machine. With `-T 1`: rsi MEG-T2 random_median 63.8 (old, 3 of 3 runs) -> 63.7 (2 of 3 runs, 63.8 in the third);
thresholds, pair classes and copy statuses identical in cth VES, tbr VES and all 22 rsi families
(`/home/toki/ssspeed/bench/bc_rep_MEG-T2/summary.txt`, `fs9run_rsi/compare.txt`, `fs9_tbr_VES/compare.txt`).
Options: (a) keep no -T: results depend on the machine, and 24 parallel workers x 129 threads gave load ~900; (b) `-T 1`:
fast, machine-independent, weak background rows change; (c) a fixed `-T 128`: the stored therioserver results, slow.
Recommendation (b), implemented as `BODYCHECK_SSEARCH_T=1`. The same mechanism applies to other small-library calls with no
`-T` or `-T THREADS` (fs_fold, fs_endcheck, step3/step4 on subfamilies under ~250 copies): their weak hits can change with
THREADS or with the machine today.

### 5.2 Step 1 is not reproducible run to run (existing behaviour)

Two runs of the old sear on rsi differ in the bit score of most rows and, for r9_15seqs, in 4 of 82,967 loci (4.6). The AIM
asks for reproducible output; ssearch36 has no seed option (mrandom.c seeds from gettimeofday microseconds).

### 5.3 trf.raw

flankscan fs2 leaves trf.raw (0.8 GB on mau) that only fs2 reads. Deleting it after fs2 would save that; not done here.

## 6. Step 2: parallel cycles failed the distribution test (rule 1)

Change tested: ten cycles at once, each call `-T 3` instead of `-T 30` (parts under 500 copies at -T 30). rsi, 22 consensuses,
67,192 copies in 4 parts, THREADS=30, 6 runs of each version interleaved (`/home/toki/ssspeed/bench/step2_rsi_all/`):

| metric | old (6 runs) | VOTE_PAR=10 (6 runs) |
|---|---|---|
| wall s | 497-574 | 404-435 |
| CPU s | 10,228-11,363 | 9,999-10,534 |
| assigned | 63,834-63,995 | 63,875-64,117 |
| best vote 10/10 | 65,380-65,559 | 65,430-65,687 |
| best vote 9/10 | 494-755 | 438-736 |
| **MEG-RS assigned** | **720-732** | **622-732** |
| **MEG-RS copies won in 6-9 cycles** | **4-16** | **2-114 (114, 112, 7, 2, 110, 74)** |
| **r10_r8_P18del assigned** | **721-735** | **709-723** (5 of 6 below 721) |
| r6_210seqs copies won in 6-9 cycles | 336-494 | 192-490 |
| per-copy agreement | old-old 98.26-98.82 % | new-new 98.42-98.75 %; old-new 98.10-98.93 % |

In 4 of 6 parallel runs one cycle sent 74-115 MEG-RS copies (those of the last part, 7,192 copies) to MEG-RL; in none of 6 old
runs. In isolated calls the settings did not differ (40 calls per setting on part 2: copies won per consensus with the same
medians within 1.5 % and the same spread; 60 per setting on part 4: MEG-RS 107-109 of the MEG copies in every call), so the
cause is not found. A step 2 shorter by ~1.3x does not justify a changed vote: reverted (f54bb09).

## 7. Whole-run estimate

From the stored stage times (section 2) with the measured stage changes; stages not changed are kept.

| | rsi | tbr |
|---|---|---|
| stored run | 128 min | 482 min |
| step 1: searches x 0.525 (sear test) | 29.6 -> ~20.4 | 21.6 -> ~18.4 (6 consensuses x 20 parts x ~1.6 s saved) |
| flank twins, default body check (4 workers) | 35.3 -> 9.3 | 426.1 -> ~26 (7.6 + 17.2 + small families) |
| flank twins with `BODYCHECK_SSEARCH_T=1` | 35.3 -> 2.5 | 426.1 -> ~11 (VES 9.9 + small families) |
| **total, default** | **~93 min** | **~78 min** |
| **total with `BODYCHECK_SSEARCH_T=1`** | **~86 min** | **~64 min** |

## 8. Commits on the branch (oldest first)

| commit | what | status |
|---|---|---|
| 183a286 | fs9: repeat filter and seed pairing in memory, confirmation in ranges; tests/test_fs9_engines.sh | accepted |
| fdc7ebf | fs9: union-find with path compression | accepted |
| 39f1375 | fs9_run: remove intermediates no step reads (FS9_KEEP=1 keeps them) | accepted |
| 56b3b41 | fs9_bodycheck: -T 1 | reverted by fa31519 (5.1) |
| 913de9a | sear: no per-part re-index, streamed window bank | accepted |
| fa31519 | revert of 56b3b41, with the -T measurement | - |
| 2954930 | fs9_bodycheck: BODYCHECK_SSEARCH_T (default unset = unchanged) | owner: opt-in; default 1 since 2026-10-26 (load ~1,400 from 32 parallel calls x 128 workers each) |
| 7b53a1c | fs9_seeds: overflow guard | accepted |
| 4210956 | step2: VOTE_PAR parallel cycles | reverted by f54bb09 (section 6) |
| f54bb09 | revert of 4210956, with the replicate result | - |

Tests: `bash tests/test_fs9_engines.sh` (both fs9 engines on a planted-case toy and a stress toy: every file identical);
existing `tests/test_step2_assigned.sh`, `tests/test_step2_dupnames.sh` pass. New dependency of flankscan stage 9: python3
with numpy (in the sinederella env); `FS9_AWK=1` runs the old awk path.

## 9. Raw outputs (therioserver)

* `/home/toki/ssspeed/SINEderella/`: my clone on speedup-ssearch; `/home/toki/ssspeed/SINEderella_8d4961c/`: the old code
* `/home/toki/ssspeed/profile/stage_minutes.tsv`: stage profile of rsi and the 22 bats
* `/home/toki/ssspeed/bench/fs9_old_cth/`: the old fs9 on cth VES with timestamps (README.txt)
* `/home/toki/ssspeed/bench/fs9_cth_VES_bc/`, `fs9_mau_VES/`, `fs9_tbr_VES/` (compare.txt), `fs9run_rsi/` (body check -T 1) and `fs9run_rsi_def/` (default body check) (compare.txt): new vs stored
* `/home/toki/ssspeed/bench/t_seeds/`, `t_repmask/`, `t_uf/`: single-step oracles on cth VES
* `/home/toki/ssspeed/bench/sear/` (summary.txt, cmp.txt), `/home/toki/ssspeed/bench/split/RESULT.txt`: step 1
* `/home/toki/ssspeed/bench/scale/`: thread scaling; `/home/toki/ssspeed/bench/bc_rep_MEG-T2/`, `/home/toki/ssspeed/bench/bc_tbr/`: body check
* `/home/toki/ssspeed/bench/step2_rsi_all/` (metrics.txt, megrs.txt; links to step2_rsi/ and step2_rsi_b/), `/home/toki/ssspeed/bench/s2win/`: step 2
* `/home/toki/ssspeed/src/fasta36/`: FASTA36 source (6d6718a), read for 5.1
