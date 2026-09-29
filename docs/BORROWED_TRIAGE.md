# Borrowed ideas and defensive tricks: triage against the current code (2026-09-29)

Sources: the two internal catalogs of 2026-08-22/23 -
"Borrowed Ideas for SINEderella" (claude.ai artifact 60b1a3bf..., ~25 algorithmic proposals from AnnoSINE_v2,
HiTE, EarlGrey, RepeatModeler2, RepeatMasker) and "Defensive Tricks Catalog" (artifact ad4d2c14..., 85 small
defensive tricks from the same five codebases; source notes DRAGEN /staging/tmp/sine_tool_comparison/analysis/).
Both were written against SINEderella as of 2026-08-22 ("no consensus refinement, no structural signal
detection, no redundancy consolidation"); much of that changed since. Every verdict below was checked
against the code at commit 2fd73b8 (SINEderella, publish/, tools/, flankscan/).

Verdicts: **DONE** (already in the code - where), **DROP** (does not apply, or conflicts with the design -
why), **IMPLEMENT** (real gap - where and how), **DECIDE** (changes results or framing - his call).

## 1. DONE - already in the code

| idea / trick (source) | where in SINEderella |
|---|---|
| CR / whitespace stripping of input FASTA (RM #1, part of #2) | `SINEderella` `sanitize_fasta`: `gsub(/[ \t\r]/,"")` on every sequence line, CR stripped from headers |
| Guaranteed trailing newline, fixed wrap (RM #35) | `sanitize_fasta` rewraps at 60 |
| Records with empty header or sequence not passed on (RM #4, #8) | `sanitize_fasta` keeps only `hdr!="" && seq!=""` (silently - see IMPLEMENT 7) |
| Input-file and tool existence checked with a clear message (RM #34) | `SINEderella` l.75-99 (`die "Missing ..."`), `step2_asSINEment.sh` `need`/"Missing tools" |
| Non-zero-size checks, die on empty intermediate (RM #11) | `[[ -s ... ]] || die` throughout; `genome.clean.fa is empty (sanitize failed)` |
| Backup before overwrite (RM #13, #31) | `consensuses.clean.fa.pre_publish.bak`, `rebuild_consensus_bank.py --inplace` -> `.pre_rebuild.bak` |
| Fixed random seeds (general; RM2 #9 spirit) | `step4_plots.sh` `seqkit sample -s 42`; `step7` `random.seed(42)`; `step8a` `RAND_SEED` (default 42, printed) |
| Deterministic ordering of ties (RM2 #9) | top100 `sort -t$'\t' -k2,2nr` (step8a l.407): GNU sort without `-s` breaks ties by the whole line, so the order is reproducible; flankscan fs4 `mode()` tie fixed 2026-09-29 |
| Flank coordinates clamped to contig ends (EarlGrey #6) | step8a `bedtools slop -g genome.fai` + start `>= 0`; flankscan fs1 records `clamp5/clamp3` |
| Strand-aware flanks and coordinates (RM2 #11; Borrowed "orientation-aware boundary") | step8a `slop -s`, `getfasta -s`; flankscan fs1 core coordinates tested on minus-strand copies (toy 197/197) |
| Final partial chunk flushed (HiTE #3) | flankscan fs2/fs3 `seqkit split2` |
| External tool that writes to cwd (AnnoSINE #6) | flankscan fs2 runs TRF with `-ngs` to stdout |
| Temp-dir name collisions (RM2 #1) | `mktemp` everywhere (but see IMPLEMENT 5 for WHERE) |
| Silent-dropout detection (EarlGrey TEstrainer; Borrowed step 3) | every failed copy kept with a reason: `no_ssearch_data`, `no_unanimous_votes`, `rejected_low_bitscore`, soft label |
| Nesting labelled, nothing deleted (EarlGrey `filteringOverlappingRepeats`) | CONFLICT flag (step3); flankscan class `nested` (host continues across the copy) |
| Redundant overlapping hits (Borrowed RepeatMasker masklevel) | step1 merges all consensuses' hits into one locus set (`bedtools merge`) before extraction |
| Reverse-complement duplicate consensuses (Borrowed CD-HIT, in part) | `canonicalize_consensus_bank.py` (RC merge, `CANON_MIN_ID` 80) - deliberately NOT forward clustering, see DROP |
| Boundaries from the copies, not from where merge cut (Borrowed AnnoSINE per-position profile) | step7 boundary walk vs genome background; step8a continuation (+150 bp up to +600); rebuilt row with proposals |
| Tandem / simple-repeat detection (Borrowed "pre-filter", 4x) | as FLAGS, not a filter: step8a `[array]` + one copy per cluster; verdict `TANDEM_ARRAY`; flankscan TRF position classes, satellite class |
| Progressive masking (Borrowed RM2 main) | `SINEderella --mask-run`, `tools/mask_genome.py`, `sear -m` depletion |
| Runner-up as a number, not only a flag (Borrowed continuous confidence) | step3 runner-up ratio per copy; LEAK = ratio >= 0.90 is a threshold on that number |
| Segmental-duplication guard via extension cap (RM2 #14) | step7: a walk that never reaches background is reported `undetermined`, not a boundary at the cap |
| Complexity masking before comparison (Borrowed RM2 Refiner) | flankscan fs3 masks consensus A tails + `dustmasker` (assignment does not - see DECIDE 4) |

## 2. DROP - does not apply, or conflicts with the design

| idea / trick (source) | why dropped |
|---|---|
| Parser internals of the five tools: RM #3, 6, 7, 9, 10, 12, 14, 15, 16, 17, 26; RM2 #4, 5, 6, 8, 13, 15, 17; AnnoSINE #4, 5, 7, 8, 10, 11, 13; EarlGrey #1, 3, 8 | guard bugs in code SINEderella does not have (FastaDB in-place compaction, EMBL/Dfam library building, BLAST text parsing, pandas/Biopython/JS quirks) |
| Installer / build tricks: HiTE #6-10 | SINEderella has no installer or compile step |
| Windows / Cygwin paths and colons: RM #19, #20 | Linux-only pipeline (therioserver, KIT, DRAGEN) |
| Fork retry with backoff, batch retry: RM #23, #24 | a failed step should stop with its reason, not be retried silently; no fork-per-sequence design here |
| Progressive search-parameter relaxation on failure: RM #27 | silently changes what counts as a hit - against reproducible, stated criteria |
| Auto-compress large outputs: RM #33 | not a problem at current output sizes |
| Poly-A TSD false-positive filter: AnnoSINE #2 | SINEderella does not call TSDs (see next row) |
| Borrowed: TSD verification as a quality flag | design decision (manuscript, Plate construction): TSDs are for inspection only; using them in flags makes old families and families that never made TSDs look like failures |
| Borrowed: CD-HIT 80 % clustering of the consensus bank (3x) | would merge real subfamilies (rsi MEG-RS vs MEG-RL share 134 bp at 97.8 %; peeled subfamilies differ by a few diagnostic columns). Keep the RC merge only |
| Borrowed: iterative consensus refinement to convergence (3x) | rejected in the manuscript (Consensus construction): iterating the builder extends the edge outward, and chimeric or boundary-contaminated alignments also converge. Replaced by one rebuild from the plate copies with every change shown as a proposal |
| Borrowed: A/B-box regex pre-filter on hits | tRNA/7SL only (5S heads have A box / IE / C box), boxes decay in old copies; not a criterion |
| Borrowed: staged excision (fast strict pass skipping the vote) | bypasses the 10/10 vote, which IS the filter (memory: -z 11 is the filter) |
| Borrowed: gap-tolerant merge (`bedtools merge -d 50-100`) | would fuse the composite and tandem units that flankscan must keep apart |
| Borrowed: block sampling with used-count tracking | sampling is single-pass with a fixed seed; nothing iterative to track |
| Borrowed: local realignment of low-scoring blocks before sim_ratio | sim_ratio is a score ratio, not an alignment used downstream |
| Borrowed: backward-looking best-match check | the vote already scores every copy against every consensus |
| Borrowed: divergence-dependent similarity threshold | assignment already uses a relative per-subfamily threshold; the fixed step-1 floor (0.8 x 65 %) is documented as relaxable |
| Borrowed: RAM-aware thread capping | no OOM seen on therioserver (499 GB); revisit only if it happens |
| Borrowed: automatic overlap resolution (ProcessRepeats) | already rejected in the Borrowed page's own dissent: SINEderella flags and defers |

## 3. IMPLEMENT - real gaps, cheap, no change to results

1. **Duplicate sequence names -> die** (RM #5, RM2 #3). `sanitize_fasta` never checks. Duplicate contig
   names break `samtools faidx`/`bedtools getfasta` silently (first one wins); duplicate consensus names
   make votes ambiguous. Where: `sanitize_fasta`, count names, `die` with the duplicates listed.
2. **Uppercase in `sanitize_fasta`** (AnnoSINE #9, HiTE #1, EarlGrey #7). NCBI assemblies are usually
   soft-masked, so lowercase reaches every extracted copy. Case carries meaning downstream (MAFFT
   `--preservecase`; lowercase = flank / proposal on plates). The plate correction recases copies, so the
   published plates are safe, and I have not found a place that breaks today - this is a cheap guard.
3. **Non-sequence characters** (RM #2). Anything other than CR/space/tab passes into `genome.clean.fa`.
   Count non-IUPAC letters/bytes per file; die above zero with the first offenders.
4. **Tool versions in `manifest.txt`** (RM2 #6, #8). The manifest records paths, threads and
   date but no versions. Record ssearch36, MAFFT, bedtools, samtools, seqkit, SubFam, sear, TRF,
   dustmasker, python. Also needed for the manuscript (Methods "Computing environment" is PENDING).
5. **Temp files under the run dir** (RM2 #12, RM #32; memory: therioserver /tmp is RAM, /var/tmp the SSD).
   `step2_asSINEment.sh` uses `mktemp -d -t` -> `/tmp` unless TMPDIR is set; `SINEderella --add` uses
   `${TMPDIR:-/tmp}`. step8a already uses `$RUN_ROOT/.step8a_XXXXXX` and step1 works inside its own step
   directory - make step2 and `--add` do the same.
6. **Count invariants between steps** (AnnoSINE #12). step1: sequences in `extracted.fasta` == loci in
   `merged_hits.bed`; step2: assigned + unassigned == extracted. Die with both numbers when they differ.
7. **Say what `sanitize_fasta` dropped** (RM #4, #8). Log the number of records skipped for an empty
   header or sequence (today silent).
8. **`.gz` input** (RM #21, #22). A gzipped genome fails as "genome.clean.fa is empty". Read through
   `zcat` when the name ends in `.gz`, or die with "gzipped input - decompress first".
9. **Quote paths in Python `subprocess(..., shell=True)`** (EarlGrey #5, RM2 #7). 13 call sites in six
   files (`tools/build_composite.py`, `tools/composite_scan.py`, `publish/flank_border_iterate.py`,
   `publish/needs_border_loop.py`, `flank_border_consensus_test.py`, `step4_diagnostic.py`) build shell
   strings with `%` formatting; a path with a space breaks them. Use `shlex.quote` or argument lists.

10b. **Report consensuses contained in another** (Borrowed RM2 main "satellite/contained-consensus
    filtering"). Not as a filter - a family whose head is another family's head (rsi MEG-RS inside
    MEG-RL, 134 bp at 97.8 %) is real and must stay - but as a line in `tools/audit_consensus_bank.py`,
    which today reports RC duplicates, N content and divergence only. Explains split votes and soft
    calls before anyone reads the plates.

## 4. IMPLEMENT - flankscan into the pipeline (from today's work)

10. **Run flankscan after assignment in `--publish`** and write its classes onto the plates as marks
    (`[dimer-L]` / `[dimer-R]` / `[tandem]` / `[nested]` / `[satellite]`), next to `[array]` and `[soft]`,
    and leave them out of the continuation decision like `[array]`. This is the Borrowed "tandem
    pre-filter" done the SINEderella way: flags, not a filter. Design: `docs/COMPOSITES.md` section 4.
11. **Shared flank across contigs / segmental duplications** (RM2 #14 idea; GLM round-2 task 08 #6;
    PLATES.md "nle MEG-TR, 8 copies"). flankscan searches flanks only against consensuses, so duplicated
    unique flank is invisible. Add a stage: flank ends (e.g. 200 bp past each copy) of all copies against
    each other; copies whose flank matches another locus (other contig or > 50 kb away) get `[segdup]`.

## 5. DECIDE - changes results or framing (his call)

1. **Bitscore filter for small subfamilies** (found 2026-09-29, not from the catalogs):
   `step2_asSINEment.sh` l.296-314, N = min(10, unanimous count); with <= 10 unanimous copies the
   threshold is 0.45 x the lowest score and can reject nothing. Options: take N from all hits of the
   consensus, or a floor relative to the self-score (sim_ratio).
2. **CpG-adjusted (Kimura) divergence** for the divergence profiles (RM2 Refiner). Mammalian SINE heads are
   CpG-rich; raw identity overstates age. Matters for "waves inferred from divergence" (GLM round-2 note).
3. **Divergence difference as LEAK support** (RepeatMasker `isTooDiverged`): a close runner-up bitscore
   with a much more diverged runner-up alignment is not a real leak. A second number next to the ratio.
4. **Mask consensus A tails / low complexity inside assignment** (RM2 complexity adjustment). flankscan does
   it; asSINEment does not, so poly-A can lift bitscores toward any A-tailed consensus. Changes
   assignments - toy test and one real rerun first.
5. **blastx (TE protein) negative check for bank entries from other programs** (RM2 RepeatClassifier):
   a SINE has no protein-coding hit; a strong hit marks a LINE/LTR/DNA fragment. Useful for the
   AnnoSINE_v2 candidates (Timema, scorpion). Needs a protein database (Dfam/RepeatPeps).

## 6. Open code checks raised by the GLM manuscript review (round 2)

- Are `[array]` copies excluded from step7's best-copy boundary walk? (If not, arrays can push the edge out.)
- Does `sear` merge same-consensus hits inside one insertion, so a dimer or array unit is counted once?
- Would a composite family trip step4's 5'/3' chimera flag on every copy?
- How are the up to 10 000 copies of the subfamily plate sampled - is there a seed?

## 7. The remaining catalog tricks (every one of the 85 is placed)

| trick | verdict |
|---|---|
| RM #18 divergence division-by-zero, RM #29 zero non-ambiguous bases | DROP - no such divisions here; the only ratio (sim_ratio) divides by a consensus self-score, which is never 0 |
| RM #25 touch empty output, RM2 #16 remove empty masked output, EarlGrey #2, #9, #10 empty guards | DONE - empty intermediates stop the run (`[[ -s ]] \|\| die`) |
| RM #28 writability probe | DONE - run dir created with `mkdir -p ... \|\| die` at the start |
| RM #30 signal-safe `system()` wrapper | DROP - Bash steps with `trap cleanup EXIT`; runs are stopped by PID tree (memory rule) |
| RM2 #2 tool exits 0 while failing | covered by IMPLEMENT 6 (count invariants catch silently lost records, e.g. `getfasta` with stderr to /dev/null in step8a) |
| RM2 #10 delete-while-iterating | DROP - Perl-specific |
| RM2 #12 temp-dir writability | IMPLEMENT 5 |
| AnnoSINE #1 IUPAC -> N | DROP - ssearch36, bedtools and MAFFT accept IUPAC; IMPLEMENT 3 only rejects non-sequence bytes |
| AnnoSINE #3 truncate before append | DONE - every run writes into a new timestamped run dir |
| AnnoSINE #12 parsed vs header count | IMPLEMENT 6 |
| HiTE #2, EarlGrey #4 single-line FASTA | DONE where offsets are computed (flankscan `seqkit seq -w 0`); elsewhere samtools/bedtools index the 60-col file |
| HiTE #4 clean rebuild of output dir | DONE - new timestamped run dir per run |
| HiTE #5 absolute paths | DONE - manifest `readlink -f`; flankscan fs3 `readlink -f` on the consensus file |

Count, one verdict per trick: DONE 26, DROP 43, IMPLEMENT 16 (= 85; RepeatMasker 10/20/5, RepeatModeler2
5/6/6, AnnoSINE_v2 2/9/2, HiTE 4/5/1, EarlGrey 5/3/2). Borrowed ideas (27 incl. the dissent): DONE 10,
DROP 12, IMPLEMENT 1 (contained-consensus report), DECIDE 4.

## 8. His decisions (2026-09-30)

- D1 bitscore filter for <= 10 unanimous copies: **NO** - leave as is.
- D2 CpG-adjusted divergence: **very important question - needs a thorough test on multiple species
  and SINEs** before any change (not a code change yet; design the test).
- D3 divergence difference next to LEAK: **undecided - needs a concrete example** (find a LEAK copy
  where the two numbers disagree, show it to him).
- D4 mask consensus A tails in assignment: **NO** - many SINEs have other types of repeat ends.
- D5 TE-protein negative check: **NO**.
- Implement items 1-9 (hygiene / reproducibility): **YES**, delegated to GLM (aider), checked by Claude.
