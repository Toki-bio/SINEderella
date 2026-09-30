# flankscan — handoff (2026-09-29, before a context compaction)

## His brief (verbatim essence)
Rewrite SINEderella's flank analysis and dimer / tandem / composite / satellite detection in
**bash/awk**, fully understandable to him, **thoroughly tested on real data**. His design:
extract all hits **once** with big flanks (~1000 bp) at the start; scan flanks for tandem repeats
with TRF tuned to catch everything from dimers up; after subfamily peeling, run flanks against all
candidate consensuses to detect composites, with a tolerance ("cancer allowed" = tolerance). A second
Claude session reviewed the design (pasted by him): keep core coords in headers, clamp at contig ends,
record inter-hit spacing, TRF periods 2–2000 (`trf x 2 5 7 80 10 20 2000 -h -ngs`), filter by copy
number, classify TRF hits by position (tail / spanning = satellite / flank background -> mask),
composites by **junction histograms in consensus coordinates** (partner, strand, core end, partner
start, gap; sharp peak with >= k copies vs a density null), linker conservation, mask consensus
poly-A tails / low complexity, handle homodimers, **A-tail insertion hotspot** (own class),
**nested insertion** (partner split around the core), misassigned copies (build composite
consensus, re-assign once), library should include tRNA / 7SL / 5S heads and LINE 3' ends; nhmmer
worth trying for old partners. He also said: delegate work to GLM.

## Plan — stages 0-6, one readable script each (SINEderella/flankscan/)
| stage | script | status |
|---|---|---|
| 1 | `fs1_extract.sh RUN OUT [1000]` — windows ±F once; loci.tsv with core_s/core_e, clamp5/3, gap5/3, nb5/3 | **written, toy 4/4 PASS** |
| 2 | `fs2_trf.sh OUT` — TRF, classes tail/head/satellite/core/partial/flank5/flank3, summary, masked windows | **toy 20/20 PASS** (2026-09-29, after tail-rule fix) |
| 3 | `fs3_partners.sh OUT CONS [T]` — masked consensuses (A tails + dust) vs masked windows (ssearch36 -z 11); units.tsv + junctions.tsv (per window and side: nearest partner, consensus coords at the junction, gap, gap seq, clamp) | **toy 11/11 PASS** |
| 4 | `fs4_junctions.sh OUT RUN [MINK=50]` — junction peaks per (family, side, partner, strand) vs density null; peaks.tsv (type, linker consensus + identity), copies.tsv, family_summary.tsv | **toy 19/19 PASS** (MINK 20 on toy) |
| 0 | `fs0_partnerlib.sh OUT.fa [TAXID=9397]` — partner library x.*: tRNAs, 7SL, 5S (Rfam cmemit), Dfam LINE 3' ends (clade + ancestors); fs3 4th arg | toy checks_6 PASS (2026-09-30) |
| 5 | `fs5_build.sh OUT RUN [NBEST=60] [MINK=50]` — candidate consensus per peak (mirrors folded, best 60 elements, L-INS-i majority), same-element fold | toy: planted elements 100 % (2026-09-30) |
| 6 | `fs6_reassign.sh OUT RUN` — re-assign ONCE with candidates added (fs3 rerun), per candidate % one full unit (accept >= 70 %), moves by old family | toy all accept, no false moves |

Tests: `tests/make_toy.sh` (gawk, seed 7) plants: 40 single TA; 40 dimers TA[1-130]+39 bp linker+TB
(two loci each, as SINEderella splits them); 20 chance TA-head+TC; 15 TA inserted in TC's A tail;
15 TA nested in TC; a 6-unit satellite (TA+300 bp); 20 TB with a (TA)15 tail; 1 TA near a contig
start (clamp). Truth in `truth.tsv`. `tests/run_tests.sh [stage...]` builds it, runs, PASS/FAIL;
stage checks go in `tests/checks_<stage>.sh` (sourced). Server wrapper: `~/tmp/up/fstest.sh`
(copies ~/tmp/up/fs/* to ~/tmp/flankscan, runs tests).

Stage-2 checks to write: tatail TB copies -> class tail, motif TA/AT (>= 18/20); satellite -> class
satellite (>= 4/6); single TA -> no satellite, no non-A tail; dimer/nested -> no satellite.
Then real data: rsi (run `~/rhin/rsi_comp/run_add_20260928_222535` and the original
`~/rhin/rsi/run_add_20260927_180847`), expected: tbr MEG-RS array (~910 bp period) for satellite;
rsi composites r1+39bp+r3, r5head+r6, r10+105bp+groupB (Tal rsi/REFINEMENT.md §10–13).

## What the Python prototype established (keep as ground truth)
- `tools/composite_scan.py`, `tools/build_composite.py` (SINEderella 014e2ae, 48c06c3); design doc
  `docs/COMPOSITES.md` (516a46d). rsi results: Tal `rsi/REFINEMENT.md` §9–14, `rsi/analysis/`.
- h/mid test: a downstream unit starting at its consensus 5' end = new copy; starting mid-consensus =
  piecewise match of one element (consensus defect).
- ssearch36 needs `-z 11` on an all-homolog library (default stats return nothing); -z 11 is
  non-deterministic (~1 % count noise).
- rsi: r1_r3 (396 bp) 18 329 copies but only 50 % full (r3 internal-repeat variant needs its own
  consensus); r5h_r6 85 % full; r10_groupB 83 % full. All rsi heads tRNA(Ile)-derived (Dfam), no
  7SL/5S; linker, r10_groupB middle, group-B part unknown. Test page: Tal `rsi_v2/`.

## Environment
therioserver; conda `sinederella` (trf, ssearch36, nhmmer, bedtools, samtools, seqkit, gawk 5.4,
dustmasker) and `rnatools` (Infernal, tRNAscan-SE); `~/refs/smallrna/Rfam.cm` (cmpressed); Dfam via
web API (urlencoded POST https://dfam.org/api/searches, organism "Homo sapiens" or "Myotis
lucifugus"; "Mammalia"/"Chiroptera" fail).
Rules: remote commands via uploaded script files only (PowerShell strips quotes: a `>` truncated a
plate once); temp under /home; toy-test before real runs; commit + push after fixes.

## GLM (delegated 2026-09-29)
Task `C:\work\glm-harness\tasks\flankscan_audit\flankscan_audit.md` (3 questions: TRF field parsing,
class exclusivity, off-by-one in fs1 core coords + fs2 masking), run with `GLM_READ_ROOTS="C:/work/glm-harness/tasks/flankscan_audit" node glm.js tasks/flankscan_audit/flankscan_audit.md`
(glm.js only reads inside GLM_READ_ROOTS - the first run failed on that);
output `C:\work\glm-harness\out\flankscan_audit.json` (log `out/flankscan_audit.log`).
**Verify every claim independently before acting** (memory: feedback_use_glm_narrow_tasks).
Next GLM candidates: write stage-2 checks (checks_2.sh) from the planted truth; review fs3/fs4
once written (<= 3 claims per task).
**Result (2026-09-29, verified by me):** no defects in fs1/fs2. Q1 TRF fields f[3] period, f[4]
copies, f[6] %match, f[8] score, f[14] motif - correct for -ngs (GLM first claimed an off-by-one,
then retracted it itself). Q2 classes mutually exclusive - correct; GLM's example (400,980) with
core 1001-1100 is NOT head as it said (980 < 986 = s-15) but flank5 - the code is right, its example
wrong. Q3 coordinates - no defect; independently confirmed by the toy (core = planted core, 197/197).

## Stage 2 result (2026-09-29, after compaction)
- GLM audit (out/flankscan_audit.json): Q1 TRF fields, Q2 class exclusivity, Q3 off-by-one -> all
  "no defect" (it retracted its own first Q1 claim). Verified independently: TRF rows show
  pct_match 89-100 / score 24..3614 in f[6]/f[8]; toy core + masking checks pass.
- GLM MISSED the real defect, found by looking at class x case counts: the core already contains
  the A tail, so `tail` required `b > e` and 24/40 single copies had their own (A)n classed
  `core`. Fixed: tail = a in [e-30, e+15] and b >= e-15; head mirrored (a <= s+15). Also the
  trf_summary.tsv header was sorted to the bottom - fixed.
- Toy after fix: own A tail = tail 30/40 single; (TA)n tails 20/20 (motif AT/TA); satellite 5/6
  (first unit = partial, nothing upstream); every A-tail insertion has an (A)n head 15/15 (the
  stage-4 atail signal); 0 satellite in dimer/nested/chance/atail; 602 N masked, all in flank repeats.
- Note: flank3 of dimer_left = the partner TB's A tail, which gets masked. Harmless for stage 3
  (consensus A tails are masked too) but stage 3 must not read N runs as gaps.
- tests/checks_2.sh has 20 checks incl. anti-vacuity ones (N > 0, >= 25/40 A tails found).

## Stage 3 result (2026-09-29)
- ssearch36 -m 8 DOES report several alignments of one consensus in one window (nested TC twice).
- Two junction problems found on the toy and fixed by rule, not by tolerance:
  1. Local alignments overrun a junction by 10-30 bp when the sequence across it scores positive
     (TC ran 31 bp into an inserted TA: TA 1-10 resembles TC 81-90). The prototype's rule "drop
     a hit overlapping a kept one by > 20 bp" lost the partner (nested found 6/15). Now: best hit
     first, weaker hits keep only their parts outside kept units (>= 40 bp, consensus coords
     linear along the alignment). Same weakness exists in tools/composite_scan.py.
  2. Consensus A tails are masked, so the main unit stopped where its tail starts and the
     neighbour's alignment ran back over the copy's own A tail (TC start 67 instead of 81). Now the
     main unit owns its own tail when it reaches its consensus tail (tag e): extended over the A run
     (>= 80 % A walker), neighbours trimmed; gap is counted after the copy's own tail (main_tail col).
- Toy: dimers linker 39+-3 >= 36/40 both halves (tested as main overrun + gap + partner start, so it
  holds however the aligners split the linker); chance TA head within 30 bp; A-tail insertion: TC
  tag e + gap >= 70 % A; nested TC facing ends 80-82|82 (truth 80|81); satellite gap ~300; single
  and (TA)n-tailed copies no partner; clamp reported.

## Stage 4 result (2026-09-29)
- Toy extended: + 30 homodimers (TB + fixed 20 bp + TB) and 30 piecewise pairs (TA[1-100] + TB[80-]);
  now 317 loci. Two peaks share one group (TB side 5 TA: dimer and piecewise) - found in turn.
- Peaks = gap mode, then main_j, then p_j, each +-10; kept if >= MINK and >= 10x chance
  (N_family x density_partner x 21 x 1/2). Type by the DOWNSTREAM unit's start: <= 15 composite
  (homodimer if same family), mid = piecewise; opposite strands = inverted.
- Fixes on the way: mode() returned the centre of the first window covering a tight cluster (10 bp
  off) -> now the commonest value inside the densest window. Tail ownership made symmetric in fs3:
  a partner whose 3' end faces the copy owns its A tail too (p_tail col 25); before, the homodimer
  gap read 34 on one side (partner tail + 20 bp) and 20 on the other. atail rule now = partner
  reaches its tail and owns >= 5 bp A run against the copy (a TA head cut before its tail next to a
  TC was called atail 3/20 before).
- Toy peaks exact: dimer TA..138 +31+ TB 1.. (= 39 bp linker, 100 % id), piecewise TA 100 | TB 80,
  homodimer gap 20 both sides with the planted linker; chance TC 20/20 chance, no peak.
- Open design point for him: an A-rich structural linker (Alu-like) and an A-tail insertion
  hotspot look alike; peaks.tsv reports linker_A, the call is his.
NEXT: real data - rsi (~/rhin/rsi/run_add_20260927_180847), expect r1+39bp+r3, r5head+r6,
r10+105bp+groupB; then tbr MEG-RS ~910 bp satellite.

## Real data 1: rsi (2026-09-29) - `fs_all.sh ~/rhin/rsi/run_add_20260927_180847 ~/tmp/fs_rsi 16 50`
67 418 copies, 2.5 min total (fs2 18 s with 16 parallel TRF parts, fs3 130 s). 47 peaks.
Reproduces the hand / prototype results (Tal rsi/REFINEMENT.md §5-11):
- r3 <- r1 composite 9 164 copies, gap 39, r1 ends 154 (full), linker consensus 96 % id  [r1 + 39 bp + r3]
- r3 <- r2 composite gap 0 (5 804); r1 -> r2 gap 31 (231)                                  [r1 ~ r2 + r3]
- r3 -> r3 piecewise 4 357: r3 ends 161, next part starts 114                            [r3 1-163 + 115-201]
- r6 <- r5 composite 3 513 (r5 ends 127, gap 10) + 1 109 (139, gap 0); r5 -> r6 1 154     [r5 head ~130 + r6]
- r3 <- r5 composite 678+393+110 [r5h + r3, 4.6 %]; r5 homodimer 73+52 [r5 + r5, 4 %]
- single shares: r3 2.4 % (proto 2.3), r6 54.6 (51), r5 49.6 (48), r9 93.6 (91)
Open, NOT verified (look at alignments before any claim):
- r10 + 105 bp + group B: only ~210 copies in peaks (r10 -> r8, gaps 100-121, linker id 0.5-0.65) vs
  ~4 100 compound copies per §10. Likely the known limit - group B is barely in this bank (r8 covers
  it partly), a part absent from the library shows as flank, not partner. r10 has 13.3 % "near".
- r3 <- r1 piecewise 3 457 (r1 1-79, then r3 from 29); prototype called a similar layout
  r10[h] + r3[e] (9.5 % of r3). r1/r10 heads both tRNA-derived; which wins may depend on masking.
- r7 -> r3 piecewise 411: 57 bp matching r3 145-201 right after r7's end; prototype said r7 99 %
  single. Possibly a 3' part missing from r7's consensus - check the alignment.
- MEG families mostly nomain (MEG-T2, MEG-TR 100 %): no core hit at E 1e-5 (weak-match limit), 162 copies.

## Real data 2: tbr (2026-09-29) and a search fix it forced
`fs_all.sh ~/chiro/tbr/run_add_20260927_143901 ~/tmp/fs_tbr2 16 50` - 621 939 copies, ~18 min.
- **Satellite: tbr MEG-RS array found** - 55 % of MEG-RS copies (137/249) classed satellite, TRF
  period 900-930 (62 at 920, 50 at 910) = the ~910 bp array. Unchanged across the fix below.
- **Bug found on tbr, fixed:** one ssearch36 over all 622 k windows returned 81 k hits; 91 % of VES
  copies had no core hit (sample: 257/3000), Rhin-1 got 2 hits. The same windows in a 3 000 or
  20 000 library: 100 % covered (VES assignment bits median 1066 - not weak copies). ssearch36
  loses hits when the library is huge and one family dominates. Fix in fs3: search in chunks of
  20 000 windows with a fixed -Z 20000 (E means the same in every genome). VES nomain 91.2 % -> 0 %.
  rsi re-run with the fix (~/tmp/fs_rsi2): same peaks within the -z 11 noise (r1+39+r3 9 194 vs
  9 164), 49 vs 47 peaks. Toy 52/52 PASS.
- Rhin-1 in tbr stays 92 % nomain: genuinely weak (best ssearch hit median 44 bits / 76 bp) - the
  documented weak-match limit, like MEG-T2/TR.
- VES: 0 peaks; 76 % single, 17.8 % near, 5.5 % chance. Observation (not a finding): 56 k VES have
  another VES <= 200 bp downstream, same strand, starting at its 5' end, gaps spread 0-50 bp;
  ~5x the density expectation in the 0-10 bp bin (rough), below the 10x peak rule. Head-to-tail
  VES neighbours - insertion preference? Needs his look.
Run times: rsi 67 k copies 3.5 min; tbr 622 k copies 18 min (fs2 TRF 8.5 min, fs3 7.7 min).
Results: therioserver ~/tmp/fs_rsi2, ~/tmp/fs_tbr2 (peaks.tsv, copies.tsv, family_summary.tsv).


## Session 2026-09-30: stages 5-6, partner library, direct tests
Full toy suite `run_tests.sh 2 3 4 5 6` = **70/70 PASS**. Commits a0862bc, then the library commit.
**Stage 5+6 on rsi (~/tmp/fs_rsi3, re-run after the fold fix) reproduce the hand rebuilds (Tal
rsi/REFINEMENT.md §10-12) without hand work** - checked by aligning candidate vs hand consensus
(~/rhin/rsi_comp/run_add_20260928_222535/consensuses.clean.fa):
- r1__r3_P1 (398 bp) = r1_r3: 99.75 % over 396/396. accept, 72 % of its elements one full unit.
- r5__r6_P26 (370 bp) = r5h_r6: 100 % over 360/360. accept 89.6 %; ~7 400 elements move in (hand 7 294).
- r10__r8_P18 (393 bp) ~ r10_groupB: 97.2 %, NOT identical - ours 1-393 aligns to hand 8-368
  (~30 bp indel, likely the middle-length variant). Built from a peak of only 87 copies, yet
  re-assignment moves **4 374 r10 elements** to it (3 775 one full unit), 469 stay r10 (hand: r10_groupB
  4 181, r10 5 188 -> 1 072). = the "misassigned copies" item of his brief, working.
- r1-head+r3 piecewise P2 (255 bp) 93 % accept - the prototype's r10[h]+r3[e] layout.
- The r3 internal-repeat variant (P9, r3 1-161 + 114-201) is only 2 % full: its copies have r1
  upstream and go to r1_r3. The 3-part element r1 + 39 + r3(with repeat) needs a SECOND round
  (stage 3-4 again with the kept candidates in the bank) - not done yet.
- Fold rule bug found and fixed on the way: "same element" first = 90 % of the SHORTER, which folded
  contained elements (r3 variant into r1_r3, r10__r8 403 bp into r10__r6 194 bp). Now 90 % of BOTH.
- Coordinate bug found on the toy first: window -> genome used the unit strand instead of the window
  strand (half the toy elements were minus-strand: 87 % consensuses, 62 elements instead of 40).
**r7 -> r3 piecewise (496) = an artifact, now removed** (looked at the P20 alignment): r3 145-201
shares r7's last ~40 bp (CTTGACTTGG...CCTGGAAAAACACACT); its alignment ran past r7's end into the
copies' A-rich tails (TAAATAAATAAAAGTT + A runs), and the trimmed >= 40 bp remainder was kept as a
partner. Rule added in fs3: a trimmed remainder >= 60 % A (or T) is a tail overrun and dropped.
**nhmmer rescue - tested directly, NOT implemented:** profile HMMs from the run's top100 alignments
(hmmbuild --hand, RF = the consensus row -> HMM length = consensus length exactly) against the 145
rsi nomain copies: 9 rescued (MEG-RS 4, TR 2, RL 3), MEG-T2 0/109. The MEG-T2 copies are 61-68 %
identical to the consensus (ssearch vs ONE copy E ~1e-3): old copies, not a search-engine limit.
SINEderella assigned them at bits 300-370 by its own votes; flankscan cannot place their junctions.
**Partner library:** therioserver `~/refs/partnerlib/partnerlib_9397.fa` (101 seqs, 35 kb: 88 Dfam
LINE 3' ends of Chiroptera + ancestors, 7SL, 5S, tRNAs). The local ~/refs/smallrna/hg38-tRNAs.fa is
TRUNCATED (123 records, Ala..Cys only - no Ile/Leu/...); a 120 s curl from GtRNAdb stops at the
same 123, so the old file was probably a cut-off download too (GtRNAdb drops the connection at
~23 kB, curl 56, even with resume). **He gave the source: local `C:\work\hg19-tRNAs.fa`** (GtRNAdb
hg19, 419 genes, 49 anticodons incl. Ile-AAT/GAT/TAT) -> copied to therioserver
~/refs/smallrna/hg19-tRNAs.fa, now fs0's default; library rebuilt (49 tRNAs + 7SL + 5S + 88 LINE ends).
Recorded in SINE_discriminator/DATA_LOCATIONS.md.
**Tail rule on rsi (~/tmp/fs_rsi4 vs fs_rsi3):** 49 vs 50 peaks, every other peak within the -z 11
noise; r7 -> r3 fell 496 -> 329 but did not vanish. The 329 left are the same thing, looked at: the
piece is r7's real 3' end TAAATAA(A)TAAAAGTT + A run, then 10-30 bp of unrelated flank (A share
0.44-0.57, below 0.6). No further heuristic added: stage 6 already rejects the candidate (0 % of its
copies re-assign to it - masked, it equals r7), and its alignment (cand/r7_133seqs__r3_58seqs_P20.aln.fa)
shows the actual finding: **the r7 consensus stops ~17 bp before the copies' structured tail
(ACACT|TAAATAAATAAAAGTT(A)n)** - a boundary note for r7, his call.

## 2026-09-30: stage 5 consensus built from flanked copies (systematic end correction)
His call on the rsi P1 plate: "consensus should be corrected". The element cut (A-tail rule) stopped
short: 90-100 % of copies continued 28 bp past r1_r3's 3' end (r3's simple-repeat tail
GTCCTGTTCCCCTTCCCCAATAAAATCT) and 5 bp before its 5' end (GGGCC - every r1-headed candidate; the r1
consensus itself likely lacks its first 5 bp). fs5 now cuts every copy with FL=100 bp flank (lowercase),
aligns with --preservecase, takes the majority over the uppercase span and extends each end while
>= SHARE=0.6 of ALL copies carry the top base (columns < SKIP=0.3 occupied are passed over).
Tests: toy 71/71 PASS; rsi fs5 rerun (~/tmp/fs_rsi5): P1 431 bp = the hand correction (Tal
rsi_v2/composites/corrected, same rule in correct_consensus.py); P6/P9/P26/P39/P40/P43 same lengths,
others within 1-5 bp. Side effect: the fold step now keeps 17 candidates (P34, P10, P12, P44 new; P13
folded) - stage 6 NOT yet rerun on ~/tmp/fs_rsi5. Known limit: an extension that reaches the flank end
(P13 3', +100) is not flagged yet - needs a longer flank or an "unresolved" mark.
- 2026-09-30 fold-rule flaw (rsi r8): stage 5 folds a candidate into the first kept one it matches at
  >= 90 % of both, in peak-size order. P42 (r8+r8) matched P26 (r5h_r6) at 90.6 % and was folded there,
  although it is 99.0 % identical to P43 (r8+r8, kept). Fix to do: fold into the MOST similar kept
  candidate, and fold P43-like twins together first.
- 2026-09-30 FIXED (his agreement): stage 5 merges by BEST hit - each candidate links to its best >= 90 %/90 %
  match among ALL candidates, linked groups keep their largest peak. rsi (~/tmp/fs_rsi6): r8 dimer group
  P42 kept + P43, P38; r5h_r6 group P26 + P37, P27-29, P48; 16 kept. Toy PASS.

## 2026-09-30: stage 7 (hierarchy) and stage 8 (singletons) - status
- fs7_hierarchy.sh: hierarchy.tsv + accepted.fa; report_hierarchy.py draws it in step6 (commit 5823986).
  rsi final report pipeline (therioserver ~/rhin/rsi_v3, log pipeline.log) adds the accepted candidates
  with --add and publishes; result goes to Tal rsi_v3/.
- fs8_singletons.sh FAMILY [CONTROL]: are a family's "single" copies standalone or missed composites?
  Groups S (singles), B (family copies in composites), A (control singles), <= 1000 each (seed 42).
  (1) relaxed partner search E <= 1 per flank, >= 20 bp, 300 bp each side; (2) full-length in consensus;
  (3) TSD with ViewAlign's detector, PORTED EXACTLY (MSA-viewer script.js _findBestTsdInFlanks;
  gawk port = JS on 300 test pairs, 300/300 identical). His notes: ViewAlign's defaults (4-20 bp, 20 %
  mismatch) are too relaxed - short motifs count by chance -> the minimum length is calibrated on
  shuffled pairs (<= 5 % chance), table in tsd_calibration.tsv, TSD_MIN overrides. If a family has TSDs,
  they mark the element's real ends: the boundary scan moves each end outward (<= 150 bp) - its best of
  ~300 positions needs its own chance calibration (toy without TSDs gave 75 % "moved out" before that).
  Toy: B = truncated TA heads, partner 3' 100 %; S = full, no partner, TSD at chance. rsi r3 (control
  r9) run in progress (~/tmp/fs_rsi6/singletons/r3_58seqs). NOT committed until the rsi run passes.
- Next agreed task: flank uniqueness over ALL copies - design in docs/FLANK_UNIQUENESS.md.
