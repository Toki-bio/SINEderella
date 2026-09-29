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

## Plan — four stages, one readable script each (SINEderella/flankscan/)
| stage | script | status |
|---|---|---|
| 1 | `fs1_extract.sh RUN OUT [1000]` — windows ±F once; loci.tsv with core_s/core_e, clamp5/3, gap5/3, nb5/3 | **written, toy 4/4 PASS** |
| 2 | `fs2_trf.sh OUT` — TRF, classes tail/head/satellite/core/partial/flank5/flank3, summary, masked windows | **toy 20/20 PASS** (2026-09-29, after tail-rule fix) |
| 3 | `fs3_partners.sh OUT CONS [T]` — masked consensuses (A tails + dust) vs masked windows (ssearch36 -z 11); units.tsv + junctions.tsv (per window and side: nearest partner, consensus coords at the junction, gap, gap seq, clamp) | **toy 11/11 PASS** |
| 4 | `fs4_junctions.sh OUT RUN [MINK=50]` — junction peaks per (family, side, partner, strand) vs density null; peaks.tsv (type, linker consensus + identity), copies.tsv, family_summary.tsv | **toy 19/19 PASS** (MINK 20 on toy) |

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
