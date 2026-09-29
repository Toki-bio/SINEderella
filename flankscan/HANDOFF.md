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
| 2 | `fs2_trf.sh OUT` — TRF, classes tail/head/satellite/core/partial/flank5/flank3, summary, masked windows | **written, NOT yet run** |
| 3 | `fs3_partners.sh OUT CONS` — masked consensuses (A tails + dust) vs masked windows (ssearch36 -z 11, both strands); per partner: side, family, strand, core cons end, partner cons start, gap, gap sequence | to write |
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
class exclusivity, off-by-one in fs1 core coords + fs2 masking), run with `node glm.js ...`;
output `C:\work\glm-harness\out\flankscan_audit.json` (log `out/flankscan_audit.log`).
**Verify every claim independently before acting** (memory: feedback_use_glm_narrow_tasks).
Next GLM candidates: write stage-2 checks (checks_2.sh) from the planted truth; review fs3/fs4
once written (<= 3 claims per task).
