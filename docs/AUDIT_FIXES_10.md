# Audit fixes 10 (Qwen audit of main 8ca35e8, 2026-11-11), branch audit-fixes-10

Scope: the four findings of the outside audit, each read in the code first and found true. One commit per fix.
Tests: tests/test_audit10.py, tests/test_audit10_step6.py. "Old" = the same test on the files of main 8ca35e8; "new" = this
branch. pytest tests on this branch (Windows, no mafft/ssearch36/bedtools/gawk): 140 passed, 16 skipped.

| id | file | fix | test (old -> new) |
|---|---|---|---|
| 1 | step6_report.py | the SINEplot result was built (ssearch36 all-vs-all + SINEplot.py, up to 1800 s) and never put on the page: the template had no `{sineplot_section}` since step 6 was added. It is now a card after the PCA card, with a nav entry. With `--no-sineplot`: no card, no nav entry | SINEplot stubbed, page built from a minimal run folder: no card -> card after `#pca`, one iframe, srcdoc equal to the stub. Stub returns nothing: no card -> card with the skip note |
| 2 | tools/extend_ends.py | the two bedtools calls are argument lists, no shell (were `%s` strings with shell=True) | run directory `run dir;$x`, subprocess replaced: "a shell is used" -> pass. With bedtools 2.31.1 on therioserver, tests/test_extend_ends.py in a plain directory: pass on both; in `run dir;$x`: fail -> pass |
| 3 | disc/test2_props.py | header format `%10s` (was `%>10s`, not a format) | the script run on a header-only calls.tsv: ValueError -> exit 0, header printed |
| 4 | tests/test_extend_ends.py | skipped when bedtools is not on PATH, as the other tool tests | no bedtools: FAIL -> SKIPPED |

## For the owner

- 1 changes the page of a plain run (no `--publish`): it gains the SINEplot card. Published pages do not change, because
  publish_run.sh passes `--no-sineplot`.
- SINEplot.py is not in this repository. It is looked for at `$SINEPLOT_PY`, `$SCRIPT_DIR/SINEplot.py`,
  `$SCRIPT_DIR/SINEplot/SINEplot.py`, then PATH. Where it is not found the card holds the skip note.
- Not run: a real SINEplot.py, so the look of its page inside the 800 px iframe and the size of the report are not seen.
- The skip note names only "SINEplot.py not found / ssearch36 unavailable"; a timeout or a SINEplot.py error gives the same
  note. Not changed.
- The SINEplot copies come from the subfamily fastas with no verdict filter, so not_SINE copies are plotted
  (item 45 of docs/AUDIT_2026-11-08.md). Not changed.

## Follow-up, same day: SINEplot on the published page (branch sineplot-publish)

Owner's decisions: SINE loci only; fix the subfamily rule in SINEplot itself.

| change | file | check |
|---|---|---|
| SINEplot and the built-in PCA draw only copies with verdict SINE in loci_verdict.tsv (item 45 of AUDIT_2026-11-08); a line under each plot says so; a run without the file draws all assigned copies and says that | step6_report.py | tests/test_sineplot_sine_loci.py: on main the PCA samples 4 copies of the toy run, on the branch the 2 SINE loci |
| copies from assigned.fasta when step2_output has no subfamilies/ folder (was: SINEplot skipped) | step6_report.py | same test file |
| SINEplot gets `--consensuses <consensus FASTA>`: the subfamilies are those names. Before, SINEplot took query names shorter than 15 characters; rsi has a 15-character consensus (r1_r2_r3_r3_C11), which became a copy | step6_report.py; Toki-bio/SINEplot 23234c4 | toy scores: 2 of 3 subfamilies -> 3 of 3 |
| ssearch36 consensuses against all sequences, `-z 3 -b =N` (was all-vs-all with default statistics) | step6_report.py | rsi run_20261101_235602: all-vs-all wrote 2.3 GB in 23 min, unfinished; SINEplot reads only consensus rows. Default statistics: 3 of 21 consensuses with no hit (r5_27seqs, r5_r5_P48, r8_83seqs). -z 3: all 21, self score 2.30-2.33 bits/nt, two runs equal but for one 29-bit hit at the cutoff. -z 11 was random (71 of 5,642 copies changed best consensus) |
| figure at SINEplot's full 1500 x 1000 px in a sideways-scrolling box (legend and mode buttons were cut off) | step6_report.py | screenshot of the real page |
| publish_run.sh no longer passes --no-sineplot | publish_run.sh | |

Real run, therioserver, rsi run_20261101_235602 (report built to a scratch file, the run folder untouched): step6 109 s
(SINEplot part 9 s, built-in PCA 90 s), page 5.57 MB, SINEplot card 1.44 MB, 50,790 SINE loci, all 21 consensuses in the plot.

On therioserver: SINEplot is cloned at /home/toki/SINEplot (23234c4) and linked as
/home/toki/miniforge3/envs/sinederella/bin/SINEplot.py, so step6 finds it on PATH in the sinederella env; plotly 7.1.0 was
added to that env (pip; nothing else changed). Other servers need the same before their published pages get the card;
without it the card shows the skip note.
