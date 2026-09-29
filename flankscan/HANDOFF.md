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
| 4 | `fs4_junctions.sh OUT` — per family junction histograms, null from density, classes composite / chance / atail / nested / homodimer / piecewise / satellite; linker agreement | to write |

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
