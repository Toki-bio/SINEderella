# Audit fixes 6 (GLM audit 2026-10-19), branch audit-fixes-6

Findings counted as confirmed only when both checks (sv-NNN-a/b) say CONFIRMED. "No caller": nothing in the pipeline
(SINEderella, step*.sh, publish_run.sh, publish/*.sh, flankscan/fs9_run.sh) runs that code; it is a manual or research tool.
"Old/new" outcomes are measured by running the file from main (git stash / git archive main) and from this branch.

## Main pipeline: fixed

| id | file | fix | test (old -> new) | rsi real data |
|---|---|---|---|---|
| 45 (+148) | flankscan/fs9_bodycheck.py (fs9_run.sh -> fs9_twins.sh) | random body background drawn only from copies in no flank pair; fewer than 10 such copies: no check, flank result stays. Docstring `--quantile` -> `--alpha` (148) | tests/test_fs9_bodycheck_background.sh: 20 duplicated copies with one body + 12 diverged unpaired: old 0/20 stay twin1 (threshold 100.0), new 20/20 (threshold 63.8) | 77 of ~63,800 copies change, in 7 of 17 families (all twin -> flank_only/unique); within seed noise (below) |
| 113 | step2_asSINEment.sh (asSINEment itself is unused) | step 2 stops on two copies with one name (votes are kept per name) in both the full and the incremental path | tests/test_step2_dupnames.sh: old exit 0, new stops | 0 duplicate names among 67,118 extracted copies: no change |
| 120 | publish/border_loop_subfam.py, publish/merge_rebuilt_consensus.py | copies extended but no consensus rebuilt: no loci.adjusted.bed (copies keep their borders), an independent side becomes untested (flank_verdict.py -> grey); merge stops on a rebuilt file without its record | test_audit6.py::test_border_loop_failed_rebuild_is_untested, ::test_merge_rebuilt_without_record_stops (fail on old, pass on new) | rsi has one border loop (MEG-RS, 638 bp 5' extension) and it was rebuilt: no change |
| 42 / 151 (pattern in the publish path) | publish/flank_border_iterate.py run_sine (used by the border loop) | old outputs removed before sine_consensus.sh runs, exit code read; the second rebuild of a subfamily reuses the work dir and returned the first rebuild's consensus when it failed | test_audit6.py::test_run_sine_failed_rerun_not_stale | no failed rebuild in rsi: no change |
| 99 | consensus_bank_lib.write_fa (orient_consensus_bank.py, canonicalize, rebuild) | write to `<path>.tmp`, then rename | test_audit6.py::test_write_fa_failure_keeps_old_bank: old truncated bank, new old bank intact | output bytes identical |
| 133 | tools/mask_from_run.py (`--mask-run`) | no hit interval -> error | test_audit6.py::test_mask_from_run_empty_table_stops: old exit 0 + 0-byte mask | rsi not masked: no change |
| 96 | tools/length_variants.py, length_variants_run.py | ssearch36 exit != 0 or `***` on stderr raises; the pair gets verdict FAILED with the reason (was UNRESOLVED "fewer than 100 5'-complete copies") | test_audit6.py::test_length_variants_failed_search_raises | length_variants.log has no failure; summary (3 TWO_VERSIONS, 1 SINGLE_MODE, 1 UNLINKED_ENDS, 2 NOT_TESTED) unchanged |
| 27 | publish/align_for_publish.sh (main path), publish_run.sh (LEGACY_SEED_ROW=1) | add_seed_row failure stops the publish (was a WARNING); add_seed_row.py accepts zero plates so a run without plates does not stop | test_audit6.py::test_add_seed_row_failure_stops_publish, ::test_add_seed_row_no_plates_ok | no add_seed_row warning in the rsi publish logs: no change |
| 4 / 217 | sear (also used by extract_alignments.sh tier D, no caller) | output names use the query FILE name | tests/test_sear_query_path.sh: old "gen-/home/.../consensus.rawhits.tsv: No such file or directory", no BED; new gen-consensus.bed, 7 hits for 6 planted copies | step 1 passes bare file names: unchanged (toy rc 0) |

## No caller: small fixes

| id | file | fix | test |
|---|---|---|---|
| 44 | flankscan/fs_lib.sh fs_endcheck (fs_all.sh, stage 6b) | 5' flank cut by a contig end: hit coordinates shifted by EFL - length, so near test, link gap and chain extension use the true distance | tests/test_flankscan_fixes.sh: old '-', new open '5' |
| 58 | flankscan/fs7_hierarchy.sh | stops without ends.tsv (old: exit 0, every end "closed") | same test: old exit 0, new exit 1 |
| 116 | flankscan/fs1_extract.sh | rebuilds a .fai older than the genome | same test: old clamp3 400 (stale index), new 0 |
| 118 | flankscan/fs8_singletons.sh | stops unless FAM is exactly one record. The substring part of the finding does not happen: `seqkit grep -p` is an exact ID match (r1 does not pick r1b, in the test); the missing-record part did | same test |
| 26, 79 | step5_align_subfamilies.sh, step5_direct.sh | a failed mafft's partial .aln.fa is removed; only non-empty files counted | tests/test_step5_failed_mafft.sh: old 4 FAIL, new 4 PASS |
| 62, 63 | run_step5_wrapper.sh | success needs a non-empty .aln.fa; a run dir not named run_* is linked as run_<name> | same test |
| 41 | extract_subfam_only.sh | skip only non-empty outputs; failed output removed | read only (same pattern as 26/79) |
| 15 | sine_consensus.sh (the border loop calls it with numbers) | refuses a non-numeric subsample size / max_iters; MANUAL usage corrected | read only |
| 65 | sine_consensus.sh | last trace record FINAL_not_converged when it did not converge | tests/test_sine_consensus_label.sh: old FINAL_converged, new PASS. Whether the border loop should treat non-convergence as a failed rebuild: see "Needs a decision" |
| 17 | sine_consensus_smart.sh | keeps its log when it holds LOW_CONFIDENCE / not converged | read only |
| 69 | SINEderella_multi | warns when generate_sinekb_report.py is missing (I could not find that file in this repo; if it lives elsewhere, the warning names the path looked at) | read only |
| 71 | step1b_cluster_subfamilies_assist.sh | `set -eo pipefail` | read only |
| 75, 76 | step4_diagnostic.py | newest step2_output_runN by number; a failed ssearch36 leaves no sineplot_input.txt | test_audit6.py::test_step4_diagnostic_picks_highest_run |
| 125 | sine_pairwise_consensus.sh | -r must hold exactly one sequence | read only |

## Research tools: fixed

| id | file | fix | test |
|---|---|---|---|
| 28, 29, 91 | tools/coseg_compare/peel_sweep.sh | failed truth step / peel -> row FAILED_truth / FAILED_peel instead of zeros | tests/test_research_failures.sh (old 3 FAIL, new PASS) |
| 30 | tools/coseg_compare/rescore.sh | writes through scores.txt.tmp | same |
| 31 | tools/coseg_compare/run_core.sh | failed COSEG / SubFam / peel arms listed in scores.txt, exit 1 | read only |
| 32 | tools/coseg_compare/run_core.sh | SUBFAM may be the wrapper file (the default) or a directory | read only |
| 129 | tools/coseg_compare/prep_labeled.py | records with no sequence left out (were counted as copies) | test_audit6.py::test_prep_labeled_skips_empty_rows |
| 35 | tools/sim/order_sweep.sh | failed SubFam: FAILED row and status 1 | read only |
| 36, 37, 107 | tools/sim/run_sim.sh | failed sim_family.py / run_core.sh: FAILED marker, exit 1 | tests/test_research_failures.sh |
| 182 | tools/sim/run_sim.sh | OUT made absolute before the cd | read only |
| 105 | tools/sim/aggregate.py | runs without a result are named on stderr | read only |
| 38 | tools/sim/sim_family.py | stops when a subfamily gets fewer diagnostic positions than --diag. Seed 1 of the shipped diag5 scenario: old rc 0 with subfamily 4 at 3 of 5; new stops. So diag5 now fails: see "Needs a decision" | test_audit6.py::test_sim_family_diag_underdelivered_stops |
| 33, 94 | tools/cpg_div/cpg_div.py | row 2 = the consensus as searched is not a copy; a row 1 without element columns is an error (was an all-NA plate) | test_audit6.py::test_cpg_div_skips_seed_row_and_refuses_no_element |
| 108 | tools/tsd_curve.py | the seed row is skipped (the condition was multiplied by 0) | test_audit6.py::test_tsd_curve_skips_seed_row |

## Not a bug

- 89 tools/composite_scan.py, 109 SubFam annotate script: ssearch36 without `-n` on DNA. Measured on therioserver: a 200 bp
  DNA query against a library with the forward and the reverse-complement copy, without `-n`: both hits, 100 %, the
  reverse one with qstart 200 / qend 1; the output reports "200 nt". ssearch36 detects nucleotide sequence itself.
- 124 sear_multi:334-337: the offset is the window start, not start + end. `gsub(/-/," -")` turns "2001-5000" into
  "2001 -5000" and only positive fields are summed (checked: chr1_sliding:2001-5000 -> 2001, NC_0001.1_sliding:2001-5000
  -> 2001). sear:502-505 does the same and the main pipeline's coordinates are right (test_sear_query_path.sh).
- 78 step5_direct.sh: the repo's SubFam wrapper writes `<prefix>.clw` (SubFam:43-50); the finding read another SubFam.
  No caller either.

## No caller, not changed

- 5 extract_subfam_only.sh pools firm copies only, extract_alignments.sh firm + soft: both are manual tools.
- 6, 42 flank_border_consensus_test.py (manual comparison page). The stale-consensus pattern of 42 was fixed where the
  pipeline has it (run_sine above).
- 49 publish/flank_border_iterate.py main(): the prototype; the pipeline only imports its functions.
- 88 tools/compare_cons.py.
- 110 SubFam/annotate_aligned_by_ssearch_seqkit.sh: I could not find this file in this repository (it may be in the
  separate SubFam repository); not changed here.

## Needs a decision

- 112 `--add` from a masked source run: the new consensus is searched on the unmasked genome while the old families were
  searched masked, and the new manifest has no MASK_* lines. docs/MASKING.md:12,86,117 says --add deliberately refuses a
  mask (the overlap of old and new families is what the re-vote and CONFLICT need), but does not say what to do when the
  SOURCE run was masked. Options: refuse --add on a masked source run; search the new consensus on the source's
  genome.masked.fa; or keep the current behaviour and copy the MASK_* lines with a note. Evidence: SINEderella:1055-1113
  (mode_add links only genome.clean.fa), build_mask's only caller is mode_full.
- 39 Assembly QC verdict (SINEderella:915) is a log line only. Should a WARN assembly change anything (verdicts, a banner
  on the page)?
- 46 flankscan/fs9_dupregion.sh is run only by its toy test; docs/FLANK_UNIQUENESS.md says it is not wired. Wire it into
  fs9_run.sh / the plates, or retire it?
- 65 sine_consensus.sh not converging: the border loop takes its consensus as rebuilt. Should non-convergence count as a
  failed rebuild (the side then becomes untested, as for 120)?
- 38 the shipped diag5 scenario (tools/sim/scenarios.tsv: --nsub 8 --diag 5, length 250, spacing 10) cannot be delivered
  and now stops run_sim.sh. Lower --diag / --spacing, or lengthen the element?

- fs9 body check depends on its random seed: in r1_r3_P1 12-36 of ~18,000 copies change twin/unique verdict with
  the seed alone (the threshold is a far-tail quantile of only --random 2000 pairs). Raising --random (e.g. 20000)
  narrows it at the cost of more body alignments, and moves those borderline copies once. (Found in the review of 45.)

## Tests and runs

- pytest: 111 passed, 9 skipped (branch at f9b7570, rerun by Claude Code after the agent)
- Shell tests on therioserver (/home/toki/tmp/auditfix6/), old (main) vs new: test_fs9_bodycheck_background.sh,
  test_flankscan_fixes.sh (old 4 FAIL / new PASS), test_sear_query_path.sh, test_step2_dupnames.sh,
  test_sine_consensus_label.sh, test_step2_assigned.sh (new: passed). Locally: test_step5_failed_mafft.sh,
  test_research_failures.sh.
- flankscan toy suite (flankscan/tests/run_tests.sh 2 3 4 5 6 and test_chain.sh): all checks PASS on main and on the branch.
- Small full-pipeline toy (/home/toki/tmp/composite_toy/small, THREADS=16 SKIP_LENGTH_VARIANTS=1 SKIP_CONSENSUS_AUDIT=1
  SKIP_REBUILD_CONS=1 SKIP_ASSEMBLY_QC=1), main vs branch at efe81ef: both rc 0; assigned 153 vs 199 of 323 (the -z 11
  shuffle spread, earlier measured 134-196); stage 9 every copy unique in both (no twins in the toy).
- rsi fs9 body check (45), old vs new on every family's flank-only tables of run_20261017_234630 (read only, rerun under
  /home/toki/tmp/auditfix6/bc45): the agent's rerun stopped at r1_9seqs (no flank_groups.flank.tsv: the body check never ran on it); the rest was
  run under /home/toki/tmp/auditfix6_review/bc. Changed copies: r1_r3_P1 28, r10_r8_P18 12, r6_210seqs 9, r7_133seqs 8,
  r1_r3_P2 7, r9_15seqs 7, r5_r6_P26 6; the other 10 families 0. All in one direction (twin -> flank_only / unique),
  thresholds +0.1 to +0.3. Seed control on r1_r3_P1: the OLD code with seeds 2/3/4 changes 31/36/12 copies against
  its own seed 1 (threshold 93.2-93.5); the new code with seeds 2/3/4: 21/6/35 (93.2-93.4). The fix's 28 is within
  the seed spread: no real-data effect on rsi beyond noise; the toy shows the case it fixes.

## Commits on audit-fixes-6

- 74ff643 audit fixes 6: flankscan 44, 45, 58, 116, 118
- 79d5b3c audit fixes 6: legacy step 5 failed-mafft outputs (26, 41, 62, 63, 79), sear query path (4, 217)
- c9d464d audit fixes 6: step2 duplicate names, bank write, border-loop rebuild, mask, length versions, seed row (113, 99, 120, 133, 96, 27, 42/151 in run_sine)
- efe81ef audit fixes 6: manual tools (15, 17, 65, 69, 71, 75, 76)
- 88b411f audit fixes 6: research tools (28, 29, 30, 31, 32, 33, 35, 36, 37, 91, 94, 105, 107, 108, 182)
- f9b7570 audit fixes 6: sim_family --diag, prep_labeled empty rows, sine_pairwise -r (38, 129, 125)
- (this report)
