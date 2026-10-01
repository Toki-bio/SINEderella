# Length versions: separate SINEs, or one SINE with a decayed 3′ end?

Two consensuses where the shorter is the 5′ part of the longer one (rsi r9 105 bp / r7 154 bp, MEG-RS 135 bp / MEG-RL 207 bp,
Gogolevsky et al. 2009) are either two monomer SINEs, each inserted as such, or one element whose copies end at different places
(decay, or a variable simple-repeat tail). Sequence similarity of the consensuses cannot tell these apart, and the bank cleaning
used to merge them (fixed 2026-10-02, see PUBLISH_WORKFLOW.md). `tools/length_variants.py` decides it from the **copies**.
It never merges or renames anything: the verdict is a decision input.

## The four measurements

Every copy of both families is cut from the genome with 100 bp of flank, aligned to the longer consensus (`ssearch36 -m 8CB`,
`-z 11`), and gets its consensus start/end and the bases that differ from the consensus (BTOP).

| | measurement | real versions | decay / variable end |
|---|---|---|---|
| A | histogram of the 3′ end of the 5′-complete copies; modes (peak ≥ 5 % of copies, ≥ 4× the surrounding median, ≥ 20 bp apart); `valley_ratio` = mean density between two modes / the lower peak | separate modes, empty valley (≤ 0.25) | one mode and a smooth tail |
| B | columns where a base is ≥ 0.25 more frequent in the short mode than in the long mode; each copy typed by its bases there; `linkage` = P(short type given short mode) − P(short type given long mode) | the internal bases follow the length (≥ 0.4) | no column differs, linkage 0 |
| C | TSD (fs8 rule: 5′ copy ≤ 4 bp before the element, 3′ copy ≤ 45 bp after its end, ≥ 7 bp, ≤ 20 % mismatch) after each mode's own end; real flank pairs against shuffled pairs (5′ of copy i with 3′ of copy i+1) | excess ≥ 10 points in every mode | the short copies are not insertions of their own |
| D | similarity of the flank past the short mode's end to the long consensus' extension, real against shuffled | none (identity at chance, 0 of n ≥ 0.75) | a decayed extension leaves residual similarity |

Verdicts: `TWO_VERSIONS` (A, B and C all hold), `UNLINKED_ENDS` (modes but no linkage: one monomer with a variable end),
`SINGLE_MODE` (no second mode), `UNRESOLVED` (too few copies, or the evidence disagrees; the report names the failed test).
D is reported, not part of the rule. A short mode can be a mixture (decayed long copies inside it): that lowers the linkage index but a
clear difference still counts, which is why B uses frequency differences and not fixed ones.

## Where it runs

* `SINEderella` runs `tools/length_variants_run.py` after the results directory is built (all three flows). Candidate pairs come from
  `consensus_bank_lib.find_length_variant_pairs`: the shorter is the 5′ part of the longer at ≥ 85 % identity over its whole length
  (edit-distance alignment, start offset ≤ 10), the longer is 10–120 bp longer, both families have ≥ 100 assigned copies, at most 12 pairs.
* Output: `results/length_variants/summary.tsv`, one `<short>__<long>.report.txt` and `.ends.tsv` (end histogram) per pair.
* `SKIP_LENGTH_VARIANTS=1` skips it; `LENGTH_VARIANTS_MAX_PAIRS` sets the cap.
* By hand: `tools/length_variants.py --genome G --hits assignment_full.tsv --families A,B --cons bank.fa --long B --out prefix`.

## Calibration on real data (2026-10-02)

Positive control, the case in Gogolevsky 2009: *Rousettus leschenaultii* (rle, Pteropodidae), MEG-RL 9,706 + MEG-RS 6,717 firm copies.
Pooled, blind to the family labels: modes at 135 and 204, valley ratio 0.072 (1,509 of 15,519 copies between), 6 diagnostic columns
(15, 25, 26, 35, 102, 121; consensus 25–26 is AA in MEG-RL, GT in MEG-RS), linkage 0.59, TSD excess 42 and 44 points,
0 of 3,536 short-mode copies with flank similarity ≥ 0.75 to the MEG-RL extension → `TWO_VERSIONS`.
About 41 % of the short mode carry the long type at the diagnostic columns (decayed long copies, or a mixed population): MEG-RS is the
majority of its mode, not all of it.

rsi (v7 bank, original r7 154 bp consensus):

| test | verdict | modes | valley | linkage | TSD excess |
|---|---|---|---|---|---|
| r9 + r7, long r7 | TWO_VERSIONS | 114, 151 | 0.086 | 0.97 (12 columns) | +46, +48 |
| r7 + r8, long r8 | TWO_VERSIONS | 152, 178 | 0.047 | 0.56 (8 columns; 42 % of the r8 mode carry the r7 type) | +45, +48 |
| r9 + r7 + r8, long r8 | TWO_VERSIONS (3 modes) | 113, 152, 178 | 0.09, 0.05 | 0.95 | +46, +47, +52 |
| r6, r8, r5, r3 alone | SINGLE_MODE | 222 / 178 / 175 / 198 | – | – | – |

Single families give one mode: the method does not invent versions where the end is only scattered.

## Limits

* The mode finder needs the shorter and the longer consensus in the same bank and copies assigned to them; a family that exists only as
  truncated copies of a longer one produces a short mode inside a single family (the same pooled test applies).
* `linkage` compares the first and last mode only; for three or more modes the report lists every mode's TSD and valley.
* Thresholds (`VALLEY_MAX`, `LINK_MIN`, `TSD_MIN`, `DIAG_DIFF`) were set on the rle positive control and the toy tests
  (`tests/test_length_variants.py`); they are not a statistical test. A verdict other than `TWO_VERSIONS` on a family with few copies
  means "not shown", not "absent".
* Hits with mixed strand in the assignment ID (`(+,-)`) are skipped; about 0.03 % of copies.
