# Parked work (reviewed 2026-11-01)

Every branch that is not in main, what it does, and why it is not merged. Inventory by GLM (glm-harness
`tasks/parked/`, `out/parked-*.json`), every claim below checked against the code and, where marked, by a test run.

| branch | what it does | status |
|---|---|---|
| `fs9-5prime-orientation` | fs9 repeat filter: 5' flanks are stored reversed, so their 20-mers were never the genome 20-mers and no 5' flank was ever masked. | **Merged 2026-11-01.** Fixed in the rep step, the awk mask and `fs9_repmask.py`; the body check's partner count k also comes from the flanks before the filter (`partners_unmasked.tsv`), so hot-spot copies keep their partners. Tests pass (test_fs9: HOT 0/40, S7 78 % found 2/2). rsi (64,133 copies): unique -> masked 3,007, twin -> masked 291, twin -> unique 35, masked -> unique 6. |
| `coverage-tiebreak` | step 2: unit vs composite by coverage, with a hand-declared COMPOSITES table | **Superseded and merged 2026-11-02 as the `unit-rule` branch**: same idea, intervals taken automatically from `composite_units.py --bank`, and the composite must also cover the unit's own stretches (the first version moved 1,522 of 1,927 rsi r5 copies into r5_r5_P48 because r5 is 97 % like r6). rsi: r6 6,212 -> 7,430 copies, r5 and P48 stable. Open: r10_r8_P18del vs r10_r8_P18 (180 P18del copies go to P18, 149 unassigned); check with alignments. |
| `junction-flanks` | border loop: core aligned once per side, flanks from the genome packed at the junction (owner's design, 2026-10-18); step8a delivers the joined plate | **Border-loop part merged 2026-11-01** (rsi 6 subfamilies: same verdicts except r5 3', where main's MAFFT left no flank column and said untested; extension no longer double-counted; 5-8x faster). The unused row rename was dropped (rows keep locus names). **Still open:** step8a's joined-plate delivery (instead of a final full L-INS-i) - needs a plate comparison first. Also fixed on main the same day: the border loop emptied its own bed from round 2 on (MEG-RS crashed). |
| `fast-mode` (rebased as `fast-mode-2`) | sear `--fast`: fasta36 ktup 6 first round + ssearch36 rescue round on the minus bank | **Tested on rsi 2026-11-02 and rejected as default, not merged.** Regular vs fast, same bank and code: step 1 found 68,603 vs 68,224 loci; fast missed 977 regular loci (419 weak P18 copies, but also 101 MEG-T2 and 34 r9 loci regular assigned) and found 598 others; SINE-verdict copies 49,340 vs 48,595 (-1.5 %); step-1 search and extraction about 41 vs 93 min (server load 300-450, not a clean timing). The documented gain was measured on a 1.65 Gb genome with one long query; the rsi bank has 21 short queries. A fair timing needs an idle server (KIT). The merged-span length rule of the original branch admitted 6 % loci regular rejects and was already dropped. |
| `fs9b-partners` | one-sided unknown partners of flank groups (`fs9b_partners.sh`, `fs9b_review.py`, toy) | **Not merged**: not wired; known defect (2026-10-17): "no stage-9 pair on the other side" is taken as "other side unique", while in rsi 6 of 7 candidates had that side untested (masked). |
| `glm-uniq-cpg-leak` (was local only; pushed 2026-11-01) | CpG divergence tool and survey (already in main), `tools/cand_ends.py`, `tools/leak_div_example.py` | **Obsolete**: the CpG parts are in main; the two research tools are not wired and not needed by the workflow. |
| `untested-minimap2-adaptive-thresholds` | looser minlen / merge distance for sear's minimap2 mode | **Obsolete**: main rewrote the minimap2 mode (2026-10-05); only queries >= 1 kb reach it. |

## One gap cost in ssearch36: tested 2026-11-07, decision open

The steps mix two gap costs: sear, step 2 and step4_plots pass `-g -3` (12 + 3k), the others pass nothing (12 + 4k). Full oma runs on
current main, one `ssearch36` wrapper per arm, no code changed, in DRAGEN `/staging/tmp/zoo_gap/arms`:

| compared | firm set Jaccard | firm copies | pctid median change, mean / max |
|---|---|---|---|
| same code twice (noise) | 0.9985 | 77,624 / 77,588 | 0.003 / 0.02 |
| current against all 12 + 3k | 0.9984 | 77,624 / 77,601 | 0.024 / 0.23 |
| current against all 12 + 4k | 0.979 | 77,624 / 72,597 | 0.46 / 4.1 |

All 12 + 3k is the same as now within noise, apart from a few pctid medians. All 12 + 4k loses 11,456 of 97,174 loci (161 go the other way). The 3,582 firm
copies among them are older and shorter (median identity 85.4 against 88.9, length 202 against 222) but no worse by verdict (SINE 44 %,
not_SINE 6 %, against 40 % and 7 % of all firm copies). One genome and one bank.
