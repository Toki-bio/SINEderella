# Audit fixes 7 (GLM + Fable audit 2026-11-08), branch audit-fixes-7

Scope: the 12 findings of the first table of docs/AUDIT_2026-11-08.md ("Real, in the pipeline path"), nothing else. Each was
re-read in the code before it was changed; all 12 were real. One commit per fix. Tests: tests/test_audit7.py (pytest; the
external programs are small stubs, so they run without ssearch36 / minimap2 / mafft). "Old" = the same test run against the
files of main (`git show main:<file>` into a copy of the tree, pointed to by AUDIT7_ROOT); "new" = this branch.
Old: 14 failed, 2 passed (the 2 are controls that must pass on both). New: 16 passed.

| id | file | fix | test (old -> new) |
|---|---|---|---|
| 41 | sear (minimap2 path, queries >= 1 kb) | exit codes of seqkit / minimap2 / awk read from PIPESTATUS; a failure stops sear (exit 1, last lines of their stderr); stderr kept in `<bank>.minimap_err.txt`, removed with the other per-query files | test_41_sear_minimap2_crash_stops: stub minimap2 exit 1, old rc 0 with an empty BED, new rc 1. Control test_41_..._no_hits_still_ok: minimap2 rc 0 with no output still gives an empty BED, rc 0 (both) |
| 72 | step2_asSINEment.sh (--add, tier-1 id -> BED) | the contig is everything before the last ':' ahead of `start-end(strand)`; an id that does not match that form stops step 2 | test_72_...: `ctg1:apple:100-200(+)`, old `ctg1 0 200`, new `ctg1:apple 99 200` |
| 20 | SINEderella canonicalize_consensus_bank() | a crash of canonicalize_consensus_bank.py is `die` (was WARNING + unmerged bank); SKIP_CANONICALIZE=1 stays the way to run without the merge. The "script missing" WARNING is unchanged (not part of the finding) | test_20_...: the function with a failing python3, old went on (CONTINUED, rc 0), new rc 1 |
| 53 | tools/rollup_families.py (--votes) | a family none of whose subfamilies has a threshold is not score-tested: its copies keep their status (were assigned with threshold 0); docstring says so | test_53_...: family C without a threshold, 10 family votes, old `assigned` threshold 0, new `no_unanimous`; family A (threshold 100) still promoted |
| 9 | SINEderella mode_resume | the "step 2+3 already complete" branch calls families_step after step 3, before step 4, as run_step2_step3 does | test_9_...: mode_resume with every step stubbed, old never calls families_step, new calls it before report_step |
| 16 | tools/satellite_trf_verify.py verify() | failure = exit != 0 or `***` on stderr; after the retry a failure raises SystemExit (was a WARNING and an empty result, read by satellite_stage.py as "no tandem loci"); satellite_stage.py already stops on a non-zero exit | test_16_... (rc 1; rc 0 with `***`): old returns {}, new SystemExit |
| 50 | tools/kindB_flank_check.sh (piece median) | the loop over stretches runs to hits-2 (hits per run from nhits.tsv, the same MAX_UNITS cap as the piece cutting), not to the count of flank-eligible hits | test_50_... on Fable's toy t50 (4 hits, 3 eligible, one 95 % piece at stretch 2): old `r 0 0` (BORDER -> COPIES), new `r 95 1` |
| D16 | flankscan/fs9_run.sh | `results/flank_twins/` and `flank_twins.tsv` of an earlier run are removed before stage 9 runs; a family that fails now has no stage-9 result, and flank_verdict.py makes its copies grey (it already does that when copy_status.tsv is missing). Also covers a stage skipped for want of jellyfish | test_D16_...: stub fs9_twins.sh fails, old keeps the earlier copy_status.tsv, new has none |
| 224 | tools/length_variants.py (+ report_blocks.py, length_variants_run.py docstring, docs/LENGTH_VARIANTS.md) | linkage counts the columns it measured; none measured -> index None (not 0.0) and verdict UNTESTED, note "linkage not tested: no column ... >= 50 copies of both modes". summary.tsv: verdict UNTESTED, linkage_index empty (the existing None handling); the report table words it "not tested (too few copies per mode for the linkage test)" | test_224_...: 40 short-mode copies + 160 long, no base differences: old UNLINKED_ENDS (index 0.00, 0 informative), new UNTESTED. Control with 60 short copies: UNLINKED_ENDS, index 0.0 on both. test_224_report_words_untested: old no word, new "not tested ..." |
| L94 | step8a_extract_alignments.sh (top-100 and soft top-up sorts) | `LC_ALL=C sort -t$'\t' -k2,2nr -k3,3 -k4,4n -k5,5n` (bitscore, then contig, start, end) | test_L94_...: old sort lines without LC_ALL (and at a bitscore tie the whole-line fallback put c2:100 before c2:99); new LC_ALL=C, ties c2:99 before c2:100. This machine's sort has no locale collation, so the locale difference itself could not be shown here |
| 38 | publish_run.sh | step-2 output = the newest step2_output* by modification time (`ls -dt | head -n1`), the same rule as SINEderella, step4_plots.sh, align_for_publish.sh, verdict_chain.sh, flank_verdict.py and rollup_families.py (first fix used the run number; changed in review so that every step uses one rule) | test_38_publish_picks_newest_step2_output: six dirs made in run order; old: name sort gives _run9 (Fable toy t38); new: _run11 |
| D17 | publish/align_for_publish.sh | `<SP>_consensuses.aln.fa` and `<SP>_consensuses_sine.aln.fa` are removed from the alignment dir and the HTML dir before they are written | test_D17_...: stub mafft, one interspersed subfamily, stale SINE file from an earlier publish: old published again, new absent |

## For the owner

- 38, "the newest run directory": the brief says the rule is the run number as an integer, as the other steps do. In the code,
  only step4_diagnostic.py (no pipeline caller) uses the run number. The pipeline steps (SINEderella families_step and --add,
  step4_plots.sh, align_for_publish.sh, verdict_chain.sh, tools/flank_verdict.py, tools/rollup_families.py) take the newest by
  modification time (`ls -dt`, `getmtime`). publish_run.sh now uses the run number. In a normal run the highest number is also
  the newest, so they agree; they disagree only if an older step-2 dir was touched later. Not changed here (outside the 12).
- 20: when canonicalize_consensus_bank.py itself is missing from the install, the run still goes on with a WARNING (unchanged).
- 224: the length-version verdict list now has UNTESTED next to NOT_TESTED (pair not run: too few assigned copies or beyond
  --max-pairs) and FAILED (search failed).

## Test suite

`python -m pytest -q` from the repo root (this Windows machine; no ssearch36 / mafft / bedtools, so those tests skip):
main 1 failed, 118 passed, 15 skipped; this branch 1 failed, 134 passed, 15 skipped (the 16 new tests). The one failure on both is tests/test_extend_ends.py
(it fails in a bedtools call on this machine), not touched here.
