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
