# Audit fixes 8 (Fable audit 2026-11-10, flankscan), branch audit-fixes-8

Scope: findings 1-4 of docs/AUDIT_2026-11-10.md, nothing else; the low findings listed there are open. One commit per fix.
Tests: tests/test_audit8.sh (bash, real ssearch36 / bedtools / jellyfish; trf is a stub). "Old" = the same test with
`FS=<main ca54a61>/flankscan`; "new" = this branch. Old: 6 FAIL, 2 PASS (controls that pass on both). New: 8 PASS.
Suites on this branch, with the real TRF 4.09.1: flankscan/tests/run_tests.sh 2 3 4 5 6: 71 PASS, 0 FAIL (main: stage 2
stops; with a trf wrapper 66 PASS, 5 FAIL); test_chain.sh 9 PASS; test_fs9.sh, tests/test_flankscan_fixes.sh,
test_fs9_engines.sh, test_fs9_bodycheck_background.sh pass.

| id | file | fix | test (old -> new) |
|---|---|---|---|
| 1 | flankscan/fs3_partners.sh (new step 2b) | ssearch36 `-m 8CC` (CIGAR in column 13 of hits.m8); a hit with an I or D run of >= SPLIT (40) bp is cut there into pieces with their own coordinates on both sequences, the hit's identity and E-value and a share of its bitscore (aligned columns of the piece / of the hit) -> hits.split.m8, read by step 3. A hit without such a gap is passed on unchanged. The gap cost stays `-g -3` (decision of 2026-11-07) | 30 TI copies (15 per strand) nested in a host: nested 3 -> 30 of 30; host junction coordinates within 5 bp of 150 / 151: 6 -> 60 of 60 |
| 2 | flankscan/fs_lib.sh fs_fold, rule 1 | the aligned columns (`$4`) must reach 90 % of both candidates, not only the span | a 233 bp candidate that lacks the middle 148 bp of a 381 bp one: `same_as:long` -> kept |
| 3 | flankscan/fs2_trf.sh | a trf run is an error when its exit code is 126, 127 or 129-192 (killed by a signal); 253, which TRF 4.09.1 returns on success, is not | stub trf with exit 253: stage 2 exit 124 -> 0; stub trf killed by a signal: stops on both |
| 4 | flankscan/fs9_run.sh | the 20-mer count (`genome.clean.fa.k20gt20.tsv`) is made once with jellyfish for every genome size and passed to every family; hash size = genome size up to 300 Mb, 4G above as before; made again when older than the genome; an empty count gets one line with count 0. Without jellyfish: genomes up to 300 Mb use the awk count of fs9_twins.sh with a warning, larger ones skip stage 9 as before | small toy genome: no count file -> count file made; copy_status.tsv and flank_groups.tsv of fs9_twins.sh identical with the jellyfish count and with the awk count (386 copies, control) |

## For the owner

- Fix 1 changes stage 3 on every hit with a gap of 40 bp or more, also those that 12 + 4k bridged before 2026-11-07.
  Tested on toys only: no real genome was run here. Stages 1-8 are run by hand; the pipeline calls only stage 9.
- Fix 4 is in the pipeline path. On the toy the result is identical; on a real genome under 300 Mb it was not run.
  The run now needs jellyfish for the fast path (it already did for genomes over 300 Mb).
- Fix 3: check the exit code of trf on the servers (`trf windows.fa 2 5 7 80 10 20 2000 -h -ngs > /dev/null; echo $?`).
  If it is 129-192 on success there, the range in fs2_trf.sh is wrong for that build.
- Identical insertion sites: the nested toy still gives two junction peaks (TI + H composite / piecewise), on this branch
  and with the gap cost of before 2026-11-07 alike, because every toy copy sits at the same host position.
