# Audit fixes 9 (Fable audit 2026-11-10, flankscan, low findings), branch audit-fixes-9

Scope: the nine findings under "Low" in docs/AUDIT_2026-11-10.md, nothing else. One commit per fix.
Tests: tests/test_audit9.sh (bash, real tools; minimap2 is a stub in block c). "Old" = the same test with
`FS=<main 8b49d3a>/flankscan`; "new" = this branch. Old: 10 FAIL, 1 PASS (a, a control: same rows on both). New: 11 PASS.
Suites on this branch: flankscan/tests/run_tests.sh 2 3 4 5 6: 71 PASS, 0 FAIL; test_chain.sh 9 PASS; test_fs9.sh (now with a
real minimap2), tests/test_flankscan_fixes.sh, test_fs9_engines.sh, test_fs9_bodycheck_background.sh, test_audit8.sh pass;
pytest tests: 150 passed.

| id | file | fix | test (old -> new) |
|---|---|---|---|
| a | flankscan/fs5_build.sh (distinct elements of a peak) | a span is compared with the kept spans of the 1 kb bins it touches, not with every kept span; same rows | 4,000 random spans: the rows of an every-pair comparison in python on both (630 elements). Time for 20,000 elements of one peak: 35.9 s -> 0.05 s |
| b | flankscan/fs8_singletons.sh | `$G.cal.tmp` is created before the awk, so a group with no conserved window no longer stops stage 8 at `cat`; a line on stderr names the group; its ends.tsv row has "-" | 30 unrelated "single" copies: exit 1 (`cat: S.cal.tmp: No such file`) -> exit 0, summary.tsv written. This finding was "read, not run" in the audit; it is run now |
| f | flankscan/fs8_singletons.sh (reference row) | the consensus from lib.in.fa up to its A tail (cons.tsv column 3), not the masked consensus with every N deleted; stops when lib.in.fa is missing | consensus with a 20 bp dust-masked stretch inside: row 1 of S.aln.fa 160 bp -> 180 bp |
| c | flankscan/fs9_dupregion.sh | no minimap2 on PATH or a non-zero exit stops the script (was `|| true`: empty PAF, "0 copies duplicated") | stub minimap2 with exit 1: exit 0 -> exit 1 |
| d | flankscan/fs3_partners.sh | dustmasker, ssearch36 and seqkit are looked for at the start | PATH without dustmasker: exit 127 with no message -> "fs3: dustmasker not found on PATH" |
| e | flankscan/fs9_twins.sh, fs9_dupregion.sh, fs9_run.sh | an index older than the genome is made again, as fs1_extract.sh does | index of a 500 bp genome dated 2020 beside the toy genome: kept -> rebuilt, in all three |
| i | flankscan/fs9_run.sh | a locus with strand `(+,-)` is taken as + (as fs1_extract.sh), so it gets a stage-9 status | the first toy locus written as `(+,-)`: not in copy_status.tsv -> in it. On old code this run also stops on the stale index of block e |
| g | flankscan/fs6_reassign.sh | the no-candidate header of reassign.tsv has `from_full`, `from_part` | 7 -> 8 columns |
| h | flankscan/fs6c_chain.sh | stops when candidates.tsv already holds chain rows or `extended:` parents | second run: exit 0 and chains.tsv cut to its header -> exit 1, chains.tsv unchanged |

## For the owner

- f changes the first row of the stage-8 alignments (it is aligned with the copies, not scored) for every consensus that has a
  dust-masked stretch inside. Not run on a real family.
- e and i are in the pipeline path (fs9_run.sh, fs9_twins.sh). i matters only for runs made with MERGE_STRAND_AWARE=0.
- h: a second run of stage 6c alone is now refused, not repaired. `fs_all.sh` is not affected.
- flankscan/tests/test_fs9.sh needs minimap2 now; before, it passed without one because of the `|| true` of c.
