# Published plates: what each row means, and what the pipeline may change

A *plate* is one published alignment of a subfamily: `<code>_<subfamily>_top100.aln.fa`,
`_rand100.aln.fa` (and `_subfam.aln.fa`, SubFam chunk consensuses). Plates are built by
`publish/align_for_publish.sh`: step7 boundary refinement → border loop → **step8a** extraction and
MAFFT → the SINE-discriminator chain (`boundary_justify` → `orient_publish_aln` →
`rebuild_consensus_row` → `trim_display_flanks` → `correct_published_aln`) → **`tools/add_seed_row.py`**.

This file records the conventions introduced on 2026-09-27/28 after Toki's review of the bat plates
(rsi MEG-RS, MEG-RL, r1_9seqs), and the measurements behind them.

## Rows

| row | name | what it is |
|---|---|---|
| 1 | `<subfamily>_extended` | the consensus **rebuilt from the copies on this plate** |
| 2 | `<subfamily>` | the consensus **exactly as the genome was searched with** (the original) |
| 3… | `ctg:start-end(strand)` | copies; `[soft]` after the name marks a soft-assigned copy |

Plates published before 2026-09-28 name row 1 `<subfamily>` and row 2
`<subfamily>_seed_as_searched`; `add_seed_row.py` converts them in place (no realignment).

**How row 2 gets there.** step8a aligns the copies together with the consensus, so right after
step8a's MAFFT the consensus row sits exactly over the copies. `tools/dup_original_row.py` then
renames it `<subfamily>_extended` and inserts an identical row `<subfamily>` as row 2 — only if that
row is the original as searched (`consensuses.clean.fa.pre_publish.bak`), i.e. not widened by the
border loop. Every chain step skips row 2 as a copy and carries it unchanged (never packed or
recased). The earlier way — adding the original at the very end with `mafft --add` — misplaced it
on rsi r1_9seqs: its 3′ part landed ~100 columns away, in columns no copy occupies, although in the
step8a alignment it sat on the copies at support 0.93. `add_seed_row.py` still uses `mafft --add`
when row 2 is missing (plates made before this, or a border-loop-widened consensus).

### Case in row 1: every automatic change is a *proposal*

The rebuild adds bases past the ends of the original and drops original bases the copies do not
carry. Neither is applied automatically any more (his decision, 2026-09-28: "trims also be
proposals only since I don't know how to do it safely"; extensions "should stay until user removes
them"). Row 1 therefore reads, against row 2 directly below it:

- **UPPERCASE** — what the copies support within the original's span.
- **lowercase outside the original's span** — a proposed **addition**.
- **lowercase at an end inside the original's span** — a proposed **trim**: an original base the
  copies do not carry, put back so nothing of the original is lost.

Copies are uppercase over row 1's uppercase span and lowercase elsewhere; only the case changes —
no copy is realigned, and no base is added or removed (checked on 312 plates, below).

### `proposals.tsv` (next to the plates)

One row per plate, rewritten on re-run (trim fields survive a re-run, since restored trims look
like row-1 letters the second time):

| column | meaning |
|---|---|
| `add5_bp`, `add3_bp` | proposed additions at each end |
| `add5_ungapped`, `add3_ungapped` | median per-copy identity of each copy's own **ungapped** flank, read outward from the original's edge, to the proposed bases. Unrelated DNA gives ~0.25; this is the number to read |
| `add5_support`/`_null` (and 3') | the same with gaps allowed, and with the proposal shuffled. Gapped alignment lifts random matches to ~0.5, and A-rich proposals score the same shuffled; kept for reference, not as a decision number |
| `trim5_bp`, `trim3_bp` | proposed trims |
| `trim5_support`, `trim3_support` | median per-copy identity to the trimmed original bases, in the alignment columns |

Measured on the bat corpus (312 raw top100/rand100 plates, re-extracted from the published ones):
- 132 additions of ≥ 5 bp sit at background (`ungapped` < 0.35) — e.g. rsi MEG-RL 5′ `tgggggaaata`,
  0.27; the column walk wrote these into the consensus silently before.
- 89 additions are clearly carried (≥ 0.6), up to +495 bp (cse MEG-RS, MEG-T2).
- trims proposed on 71 (5′) and 122 (3′) plates, median 14 bp, support ~0.25 — mostly original ends
  the copies really do not carry (SINEbase consensuses longer than the bat copies).

### What the verdict judges

`verdict.py` and `flank_uniqueness.py` judge the element **as rebuilt**: row 1 without its
lowercase *trim* proposals, *with* its lowercase additions (`fix_alignments.judged_span`). That is
what they judged before the marking, so their calibration stands. Judging only the uppercase span
was tried and is wrong: the A-tail and target-site duplications the rebuild placed at the ends fell
into the "flank", every copy's flank began with the same A-run, the shared-flank test collapsed the
core (rsi r3_58seqs 82 → 1 of 100 copies), the TSD score went to 0, and six real SINEs were called
Not SINE.

Row 2 is not a copy. Every per-copy measurement skips it (`fix_alignments.is_seed(name, cons_name)`:
the old `_seed_as_searched` tag, or the consensus name without `_extended`).

## Copies on the plate: 100 means 100

`top100` / `rand100` promise 100 copies. A subfamily with fewer than 100 **firmly** assigned copies
(10/10 votes, above the bitscore bar) is filled up to 100 with its **soft**-assigned copies
(`step2_output/unassigned.tsv`, `Soft_Subfamily` = the query that found the locus): top100 ranks
them by search score, rand100 draws them at random; firm copies always come first. Soft rows carry
` [soft]`. A subfamily with only soft copies still gets plates. `SOFT_TOPUP=0` restores firm-only
plates. Found on rsi: MEG-RS had 16 firm and 1,570 soft copies (they split their votes with MEG-RL,
whose first 134 bp are MEG-RS at 97.8 %), so its "top 100" held 16 and the verdict could not assess it.
The report's buttons show the real count ("top 16", "top 100 (84 soft)") and state the rule on hover.

## Sequence shared past the ends: kept aligned, and extracted far enough

His review of rsi r1_9seqs: past the right end the copies continue with the same sequence (an A-run
of variable length, then `gtcctggaagtacacactgttccccaataaagtcctgttcccc…`) for ~190 bp. That region
must be shown **aligned until similarity is lost**, and where copies run out of sequence before that,
the end (or its absence) cannot be established.

`continuation.py` (SINE-discriminator) measures it **per copy**, because MAFFT splits such a stretch
into blocks different copies occupy (column occupancy 0.1–0.9) and a column walk cannot follow it:
each copy's own bases are read outward from the original's edge, a base scores 1 when it equals its
column's majority (columns reached by ≥ 3 copies and ≥ 5 %), a 20-base window slides along the copy,
and the median over copies is the similarity profile. Unrelated flank sits at 0.40–0.47, shared
sequence at 0.9–1.0; similarity is lost below 0.60.

| status | meaning | what happens |
|---|---|---|
| `none` | below 0.60 at the edge | nothing |
| `ends` | falls below 0.60 while ≥ 50 % of copies still have sequence | `boundary_justify` keeps that stretch aligned; flanks are packed only beyond it; `correct_published_aln` neither extends nor repacks on that side |
| `unresolved` | copies run out (< 50 %) while still similar | as `ends`, and **step8a re-extracts** the plate with +150 bp on that side and realigns, up to +600 bp (`CONT_MAX_EXTRA`); `CONTINUATION=0` switches it off |

Recorded per plate and side in `continuation.tsv` next to the plates.

**How the re-extraction loop runs** (`align_until_resolved` in step8a, after the fixes of 2026-09-28):

1. Round 0 is the normal plate: 50/70 bp flanks, MAFFT L-INS-i (`--localpair --maxiterate 1000`).
2. `continuation.py --need` asks for +150 bp on each `unresolved` side. It decides on the
   **independent copies only**: rows marked ` [array]` are left out, because tandem-array units share
   their flank far past the element and would drive the loop to the cap for nothing. With fewer than
   3 independent copies it asks for nothing.
3. Each further round **re-extracts only the new 150 bp per side** (`extract_flank_align … extract`,
   no alignment) and `tools/extend_plate.py` aligns just those segments among themselves (L-INS-i,
   `--maxiterate 2`, ~100 × 150 bp, seconds) and appends the block to the plate on that side. Rows are
   matched to the new extraction **by coordinate containment** (the header coordinates are the extracted
   region, so they change every round; matching by name matched nothing and the loop ran to the cap on
   the toy). Rows MAFFT reversed (`_R_…`) take the other side's segment, reverse-complemented. Rows with
   no new segment (contig end) and the consensus rows get gaps.
4. When `continuation.py --need` asks for nothing more, the delivered plate is built **once** with full
   L-INS-i at the final flank size (only if the flanks were extended; otherwise round 0 is the plate).

Speed-up: the **subfamilies are aligned in parallel** (`PAR_SF`, default `THREADS`), each in its own temp dir; every MAFFT stays **single-threaded** (`MAFFT_THREADS=1`), so the delivered plate is the reference alignment. Threading MAFFT was tried and rejected — see the benchmark.

Why (his review, 2026-09-28): L-INS-i aligns every pair of sequences locally — 101 rows = **5,050
pairwise alignments** — plus up to 1,000 refinement iterations, and its cost grows ~L². Before any fix,
every round re-aligned the whole plate: 5 rounds from ~320 to ~1,520 bp cost ~50× one plate, twice per
subfamily, single-threaded, plates in sequence. cse MEG-RS rand100 alone ran > 30 min. A first fix
(7afbf65) used FFT-NS-2 for the decision rounds; the benchmark below showed FFT modes change the
continuation calls, so it was replaced by the incremental append. `--adjustdirection` is kept although
a GLM audit suggested dropping it: step8a extracts `+,-` loci as `+`, and MAFFT turns those copies round.

### Speed benchmark (2026-09-28, therioserver, `tests/bench/`)

Real copies of published cse plates (and rsi r1), re-extracted at base flanks and at +600 bp.
Q = fraction of the 1-thread L-INS-i reference's aligned residue pairs that the test reproduces;
cont = `continuation.py` status on the result.

| case | flanks | L-INS-i 1 thr | L-INS-i 16 thr | L-INS-i it2 16 thr | FFT-NS-i | FFT-NS-2 |
|---|---|---|---|---|---|---|
| MEG-RL | base | 0.3 s | 0.2 s, Q .95 | 0.2 s, Q .96 | Q .23, 5′ ends:2 ✗ | Q .18, 3′ ends:12 ✗ |
| MEG-RL | +600 | 3.5 s, 3′ ends:2 | 1.5 s, Q .98, 3′ none ✗ | 1.1 s, Q .96, 3′ ends:2 | Q .07 ✗ | Q .07 ✗ |
| MEG-RS | base | 8.7 s | 1.2 s, Q .96 | 1.0 s, Q .95 | Q .43, 3′ unresolved:50 ✗ | Q .43 |
| MEG-RS | +600 | **260.5 s** | 81.2 s, Q .47 | 21.5 s, Q .68 | 31.6 s, Q .15 | 13.7 s, Q .14 |
| MEG-T2 | base | **254.4 s** | 46.1 s, Q .07 | 15.7 s, Q .05 | 15.3 s, Q .01 | 4.6 s, Q .01 |

✗ = continuation call differs from the reference.

Pairwise Q alone was misleading, so each alignment was also scored **reference-free**
(`tests/bench/colq.py`): over the element columns (the CONSENSUS row's first..last letter), the fraction
of residues equal to their column majority, next to how many columns the element occupies.

| case | L-INS-i 1 thr | 16 thr | it2 16 thr | FFT-NS-i | FFT-NS-2 |
|---|---|---|---|---|---|
| MEG-RL base: element cols / agreement | 279 / .575 | 282 / .581 | 280 / .578 | 289 / .655 | 294 / .651 |
| MEG-RL +600 | 948 / .564 | 950 / .564 | 945 / .563 | 995 / .615 | 987 / .612 |
| MEG-RS base | 195 / .805 | 200 / .804 | 200 / .805 | 404 / .780 | 351 / .787 |
| MEG-RS +600 | 284 / .778 | 290 / .786 | 270 / .779 | 1227 / .704 | 1123 / .700 |
| MEG-T2 base | **671 / .477** | **1580 / .443** | 742 / .473 | 2037 / .484 | 2683 / .470 |

Element-only Q (`qscore.py --elem`): MEG-RS +600 rises from .47/.68 to **.86/.91** — the low whole-row
Q was the unrelated flank, which has no true alignment. MEG-RL / MEG-RS: the L-INS-i variants give
equally good, slightly different alignments. **MEG-T2 (copies ~1 kb): 16-thread L-INS-i is genuinely
worse** — the consensus row spreads over 1,580 columns instead of 671; `--maxiterate 2` stays close.
FFT modes spread the element over 2–4× the columns and change continuation calls.

**Decision:** the delivered plate keeps the reference setting (L-INS-i, 1,000 iterations, 1 thread);
speed comes from parallel subfamilies plus the incremental continuation (the appended 150 bp blocks use
`--maxiterate 2`; they only feed the decision, the final plate is re-aligned in full).

Examples (bat corpus): rsi r1_9seqs 3′ `ends` at +189 bp (79 % of copies still covered);
cse MEG-T2 5′ `unresolved` — identical at +371 bp when coverage drops; rda Rhin-1 and rsi MEG-RL:
`none`.

## Why not a better column walk

Before settling on proposals, a stricter walk was drafted and measured. At 20 copies the top-base
fraction of **random** DNA passes the walk's 0.45 cutoff in 17–25 % of columns (95th percentile of
the top-base fraction: 0.60 at 10 copies, 0.50 at 20, 0.43 at 30, 0.35 at 100), and MAFFT pulls
similar bases into the same flank column, so aligned flank passes more often still (rsi MEG-RL:
17 of 26 flank columns). A per-plate cutoff with a run rule was too strict on old low-copy families
(his 9-copy MEG-RL lost its real 3′ end), and a poly-A end rule fails because many SINEs do not end
in poly-A. No cutoff removes flank letters without also losing real ones; hence proposals, with the
support number next to each.

## Stray original end blocks

MAFFT sometimes aligns a few of the original's end bases on their own, far from the rest: SINEbase
Rhin-1's first `GGGG` sat ~40 columns left of its body on cth, cse, fho, mtu, msc and rna; VES's first
bases did the same on cth and rmi. Taking the original's span as first-to-last letter then stretched it
over those empty columns, and whatever the rebuild had put there (flank letters at support 3–4)
became UPPERCASE and was judged as element (rna Rhin-1: `gT--A--TAA--A--T` in front of the head,
score 95.5). `add_seed_row.mark` now takes the span from the original's **main block**: an end block
of ≤ 12 letters (`STRAY_MAX`) separated by > 5 empty columns (`STRAY_GAP`) from the next original
letter is shown as a trim proposal instead — **unless the copies carry it**: ≥ 50 % of them have
letters in its columns (`STRAY_OCC`) *and* match its letters at median identity ≥ 0.40 (`STRAY_ID`).

Both conditions were learned the hard way. Without the "carried" test, cth Rhin-1's original head
`GGGGCGGCCGGT` — 25 columns from its body, but so are the copies' heads — became a trim, every copy's
head became "shared flank" (85 of 100), and the call went SINE 96.8 → Not SINE 45 (rsi r4 the same).
Occupancy alone does not work on a published plate, where flanks are packed against the element and
fill a stray block's columns (mtu Rhin-1 `GGGG`: occupied, identity 0.25); identity does (cth head 0.58,
rsi r4 0.83, stray `GGGG` on rna/vmu/tni/cse/fho 0.00–0.25). A 15-column gap missed mtu's `GGGG`, 14
columns from the body.

The same main block is used by `continuation.py` (cth Rhin-1's last 6 original letters sat ~50 columns
past its body where 8 % of copies reach; measured from there, cover was 0.11 and the ~80 bp the copies
share were missed) and by the verdict's `judged_span`.

## Tandem arrays

A family whose best-scoring copies sit in one tandem array looks like a perfect SINE on a top100
plate: the array's units are near-identical *including their flanks*, so they rank first by bitscore,
fill the plate, agree at every column, and their shared flank reads as a "continuation" of the element
with ungapped support ~1.00. Bat corpus: the MEG-RS top100 plates of vmu, tbr, fho, tni and cse held
99, 98, 95, 87 and 82 copies from tandem clusters (spacing ~1–30 kb); rsi MEG-RS/MEG-TR/MEG-RL are one
GC-rich array unit (`…ccctgccgccccttgcccct…`) hit by three queries; rmi MEG-RS is one array on
CM093732.1 84.25–84.41 Mb.

Two measures:

- **step8a** (`tools/array_order.py`) marks loci in a tandem cluster (≥ 3 copies on one contig with
  neighbours ≤ 50 kb apart) and, for top100, **selects one copy per cluster** before any second one:
  the other members go after every independent locus, so they reach the 100 only if the family has
  too few independent copies. Plate rows from a cluster carry ` [array]`. (Display order on the
  plate is MAFFT's `--reorder` guide-tree order, not this order.)
  Clusters are looked for **only among the candidates that can reach the plate**: the first 300 ranked
  loci for top100, the 100 drawn rows for rand100 (`--limit`). The first version clustered all loci of
  a subfamily; an abundant family (100,000 copies in 2 Gb, one per ~20 kb) then had nearly every copy
  marked — 98 of 100 dispersed top100 copies in a test at genome density
  (`tests/toy/test_array_order.py`). The toy run caught it before the bat republish.
- **step8a continuation**: ` [array]` rows do not count when deciding to re-extract (above).
- **verdict** (`verdict.py`, from the plate row names) reports `TANDEM_ARRAY` when ≥ 10 % of the copies
  are in such clusters, and caps the call when ≥ 50 % are (`overall.py` counts it as negative
  evidence: "Doubtful").

The array share is not a verdict on the family: in ntu, mme and mtu the same MEG-RS consensus finds
real, dispersed copies with the head `GTCTACGGCCATACCAC` and tail `TGTAGGCTTT(A)n` mixed with array
units. Putting the independent loci first is what lets those be seen.

## `rebuilt_vs_orig`

`proposals.tsv` also records the identity of row 1 (as rebuilt from the copies) to row 2 over the
original's span. A low value means the copies are a *different* element that the query merely hits —
rmi's Rhin-1, VES and MEG-T2 plates are all the local tRNA SINE with the head
`GGGGATGCCGGGATAGCGCAGTGG`; lly and mev "Rhin-1" copies are other tRNA families.

## Reading 159 plates (2026-09-28)

Every top100 plate in the corpus views (159, all 25 bat species plus rsi's r-subfamilies: both edges, 10–25 copies each,
support and occupancy per column) was read by eye against the verdict and the proposals. The log is
`SINE_discriminator/plate_reading_2026-09-28.tsv` (plate, my 5′ and 3′ reading, whether the proposals
are right, whether the verdict agrees, note).

- **Verdict agrees** on 142 of 159, disagrees on 14, unclear on 3 (mau MEG-RS, msc MEG-RS, ttr MEG-RS).
  After the fixes of the day 11 of the 14 are resolved on the 312-plate corpus rerun (rsi r2/r10, cth
  Rhin-1, msc MEG-RS; cse MEG-RS/T2, hla MEG-TR, rmi MEG-RS, rsi MEG-TR, vmu MEG-RS -> Doubtful 45; rsi
  MEG-RS -> Not SINE 45). Open: hla MEG-RS (SINE 100 on 36 copies, about half junk, not bimodal enough
  for CONTAMINATED), mev MEG-RS (SINE 100 with array copies on top), ntu MEG-RS (Not SINE 45; a real
  MEG-RS diluted by arrays - needs step8a's array ordering, which the corpus cannot test). Rhin-1 and VES are real wherever they have ≥ 100 copies; MEG-T2
  is junk in every microbat; MEG families are all real SINEs in rle (the positive control, SINE 100 on
  all four).
- **Disagreements** are almost all tandem arrays (rsi MEG-RS/TR, rmi MEG-RS: called SINE/Cannot
  assess on array units) or real MEG-RS diluted by arrays (ntu MEG-RS: Not SINE 45 though ~half the
  copies carry the real head and tail). The views were made before `TANDEM_ARRAY` and `array_order`.
- **Proposals** are right on nearly every plate: A-tails and shared tails come out at 0.5–0.8, flank
  junk at 0.2–0.35. Wrong only where the copies are array units (support ~1.00 from identical flanks)
  and where a stray original block widened the span (fixed, above).
- Seen, not yet acted on:
  - additions longer than the continuation extent carry junk at their outer end (rle MEG-TR: +30 bp
    proposed, continuation `ends` at +19, the extra is `aaaaa` at support 4–5);
  - an A-tail addition at ungapped support 0.45 (rle MEG-T2) is dropped from the judged span; harmless
    there, but A-runs could always be kept;
  - support measured from `[array]` copies should not count toward `add*_ungapped`;
  - copies sharing a long flank *across different contigs* (nle MEG-TR, 8 copies) are a larger repeat,
    which the tandem test (same contig) does not see;
  - tandem PAIRS with shared flanks escape the >= 3 rule (ttr MEG-RS: two pairs, 7-11 kb apart);
  - several ttr MEG-RS copies continue past the A-run into a tRNA-like head - possibly a MEG-RS + tRNA
    SINE dimer; not examined further.

## Testing: a toy run before any real publish

Changed step8a / publish code is run end to end on a toy first (his rule, 2026-09-28: "always test long
term code in advance on toy example"). `tests/toy/`:

| file | what it does |
|---|---|
| `make_toy.py DIR` | builds a SINEderella run dir in seconds: a random 3 Mb + 2 × 100 kb genome, consensuses, `step2_output/assigned.fasta`, `unassigned.tsv` |
| `run_toy.sh` | runs step8a and `publish/align_for_publish.sh` on it, with PASS/FAIL checks; `SD=` / `DISCD=` point it at test copies of SINEderella / SINE-discriminator |
| `test_array_order.py tools/array_order.py` | tandem-array selection at genome density (100,000 random loci over 2 Gb + a 20-unit array) |

Toy families, one per branch: **TOYS** — dispersed + 5-unit array + soft copies (array marks, soft
top-up); **TOYA** — 8 array units sharing 300 bp + 4 independent copies (must NOT re-extract);
**TOYB** — shared 40 bp tail (continuation `ends`); **TOYL** — 200 bp shared past the 70 bp flank
(must re-extract once and end at ~195 bp); **TOYC** — soft copies only. Dispersed copies sit 60 kb
apart (a 12 kb spacing made every copy look like an array). `publish_run.sh` stops at step6 on the
toy (no step3 files); `inject_disc_report.py` was tested on a copy of the real rsi report instead.

Each check was shown to fail on the code it guards against: the old array_order marked TOYB/TOYC
copies; the old continuation re-extracted TOYA twice.

The final plates are checked too, not only step8a's output: the cse republish showed every
` [soft]` / ` [array]` mark gone from the published plates although step8a wrote them.
`correct_published_aln.py` read the plate with a FASTA reader that keeps only the first word of a
header and rewrote it. Fixed (it now keeps full names); `run_toy.sh` checks the marks after
`align_for_publish`. The verdict's TANDEM_ARRAY was not affected (it uses the coordinates in the
names), but the report's soft-copy counts were.

## Rebuilt consensi across species

Row 1 of every top100 plate, from all 25 bat species, aligned with the queries: Tal
`chiroptera/recreated/` (scripts in `scripts/`, page card on `chiroptera.html`, LOG in
`chiroptera/LOG.md`). The "Rhin-1" rebuilt in 15 non-rhinolophoid bats is one other element
(≥ 0.90 to each other, 0.58–0.66 to Rhin-1); real Rhin-1 is only in the Rhinolophoidea.
