# SubFam / SINEderella route / COSEG: comparison log

Companion to `COSEG_COMPARISON_DESIGN.md` (design and decision rules, fixed before any run). This file records
what was run, where, and what came out. Newest entries at the bottom. Everything runs on KIT under
`/data/V/toki/coseg_cmp/` (hand-curated sets: `<species>/`, simulations: `sim*/`, code clones: `code/`, `code2/`).

## How to reproduce

```
# hand-curated sets (needs the owner's aln_c on KIT)
bash tools/coseg_compare/run_compare.sh saq        # saq | ccr | teu | dmo
# simulated sets with a known tree
bash tools/sim/run_sim.sh OUTDIR "1 2 3" 4          # scenarios in tools/sim/scenarios.tsv
```

Arms (see `tools/coseg_compare/run_core.sh`): COSEG (`-k -d -m M`, Price's rule "drop" and "keep"), SubFam chunks at
several chunk sizes, SubFam chunks plus the peel (`SINE-discriminator/peel_features.py`) on the chunk consensuses.
Scoring: `score.py` (purity, homogeneity, completeness, V-measure, adjusted Rand index; on the copies each method
placed, and on the copies placed by all methods). Reference for COSEG: plurality consensus of all copies
(`colcons.awk`), so that no subfamily is favoured.

## Entries

### saq, 900 hand-curated copies, 9 groups (2026-10-08)

100 copies per group from the owner's `POS__saq__*` alignments. COSEG: 883 placed (17 dropped as truncated), 10 groups,
V 0.60, ARI 0.45. SubFam chunks (n 50 / 20 / 10): ARI 0.20 / 0.11 / 0.07, over-split by design (completeness 0.31–0.36).
SubFam + peel (n 20): 5 groups, 780 placed, V 0.36.

Reading: this set is too small for the peel, which was built for about 600 chunk consensuses from 30,000 copies
(900 copies give 45–90 chunks, and the peel minimum group is 5 chunks). The result is not a fair test of the route and is
not evidence against it. The fair regime is the simulation at 16,000 copies and more.

### Tools checks (2026-10-08)

* SubFam 1.2.1 awk consensus vs EMBOSS `cons`: 3,600 runs, 0 different; MSF writer vs EMBOSS `seqret`: 80 files, 0 different.
* `kmer_order.c` vs the Python/numpy ordering: 292 synthetic cases and real saq/toc sets, byte-identical; vs ViewAlign `kmer-tree.js`: identical on 10 comparisons.
* Robustness suite (`SubFam/tests/test_robust.sh`): found a silent exit on protein input, an MSF header that embeds the output directory, and a final alignment that depends on thread timing; fixed on the SubFam branch `tests-robust`.

### SubFam 1.2.2 reproducibility checks (2026-10-08, KIT)

* `tests/test_robust.sh`: 24 of 24 pass after the fixes (1 record, 2 records, fewer records than the chunk size, `-n 2`, duplicate ids, ids with spaces and pipes, lower case / N / IUPAC, an empty record, CRLF, wrapped lines with blank lines and no final newline, 120 identical sequences, gapped input, sequences shorter than k, mixed strands with and without `-r`, protein input (now refused with a message), `-c`, an output directory with a space, and determinism).
* Determinism on 3,000 real saq copies: the consensuses, the final alignment, the chunk table and the MSF file are byte-identical across runs with 1, 8, 8 and 16 threads (before `--threadit 0`, the final alignment differed between two 8-thread runs).
* The sinederella wrapper (`SubFam input.fasta 50`) on the same 3,000 copies: 60 chunk files, 60 consensuses, `input.clw`, `input.msf`, `input.chunks.tsv`, as the pipeline expects.

### Four hand-curated sets (2026-10-08, KIT, `run_compare.sh`)

100 copies per group from the owner's `POS__<sp>__*` alignments; ground truth = the owner's group of each copy. ARI over all
copies (unplaced copies count as singletons), best setting of each arm in brackets; full tables in `<sp>/scores.txt`.

| set | copies | groups | COSEG (best of m 20/10/5, drop/keep) | SubFam chunks (best n) | SubFam + peel (placed copies) |
| --- | --- | --- | --- | --- | --- |
| saq | 900 | 9 | 0.447 (keep, m 20) | 0.198 (n 50) | 0.195 (n 20; 780 of 900) |
| ccr | 800 | 8 | 0.316 (keep, m 20) | 0.191 (n 50) | 0.087 (n 10; 760 of 800) |
| teu | 600 | 6 | 0.429 (drop, m 20) | 0.247 (n 50) | 0.430 (n 20; 560 of 600) |
| dmo | 500 | 5 | 0.424 (keep, m 20) | 0.267 (n 50) | 0.416 (n 10; 430 of 500) |

Reading, without going beyond the data: (1) COSEG is the stronger or equal arm on all four sets; on teu and dmo the peel equals it
on ARI over all copies, on dmo with the peel placing 360\-430 of 500 copies (ARI on placed copies only: 0.56 against COSEG's 0.42);
on saq and ccr the peel is clearly behind. (2) SubFam chunks alone are over-split by design and are not a partition; their best ARI is
the coarsest chunk size. (3) These sets hold 500\-900 copies, so SubFam makes 25\-90 chunks and the peel, built for about 600
chunk consensuses, has little to work with; the numbers are a lower bound for the route, not a fair measure of it. (4) The truth is the
owner's chunk-level curation carried onto 100 sampled copies per group, and chunk-level labels of this kind had a purity ceiling of 0.88 on Timema
(SUBFAMILY_METHOD.md), so no method is expected to reach 1.0.

### Simulation with a known tree (2026-10-08, KIT, `tools/sim`, output `simA/`)

Eight subfamilies on a tree (parent = the previous subfamily with probability 0.5, otherwise a random earlier one), each adding
`--diag` new diagnostic substitutions to a 250 bp consensus; copies carry private substitutions at a per-copy divergence that falls from
`--div-old` to `--div-young` along the subfamily index. Ten small scenarios (1,600 copies, 3 seeds) and three large ones (16,000-18,450
copies, 1 seed). ARI over all copies, best setting of each arm (mean over seeds):

| scenario | what changes | COSEG | SubFam chunks | SubFam + peel (defaults) |
| --- | --- | --- | --- | --- |
| base | 2 diagnostic changes per subfamily, divergence 12 % to 3 % | 0.625 | 0.090 | 0.148 |
| diag1 | 1 change per subfamily | 0.330 | 0.045 | 0 (no group) |
| diag5 | 5 changes per subfamily | 0.313 | 0.092 | 0.241 |
| indel | 60 % of subfamilies use a 3 bp deletion | 0.676 | 0.076 | 0.223 |
| trunc | 30 % of copies cut at the 5' end | 0.556 | 0.075 | 0.162 |
| rare | three subfamilies of 40, 20 and 10 copies | 0.545 | 0.029 | 0 |
| conv | 10 % of copies carry a tract of another subfamily | 0.625 | 0.085 | 0.162 |
| old | divergence 22 % to 10 % | 0.515 | 0.021 | 0.071 |
| young | divergence 5 % to 1 % | 0.772 | 0.169 | 0.285 |
| revcomp | half the copies reverse-complemented | 0.278 (keep) | 0.072 | 0.015 |
| large (16,000) | as base | 0.688 | 0.025 | 0.145 |
| large_diag1 | as diag1 | 0.442 | 0.013 | 0 |
| large_rare | 4 x 4,000, 2,000, 300, 100, 50 | 0.742 | 0.006 | 0 |

The three `large` rows are single runs; the others are means over three seeds (`simA/summary.tsv`). In `revcomp`, COSEG "drop" places almost nothing because the
aligner is forward-only; "keep" is shown. More diagnostic changes (`diag5`) did not raise the COSEG score over `diag1` here (0.313 against 0.330), which I have not explained.

Why the peel finds so little here, from its own logs:
1. The chunks are not pure. SubFam chunks of the simulated sets have purity 0.45-0.58 at the default scenario (0.28 for old families, 0.62-0.77 for young ones), so the
   consensus of a chunk already mixes subfamilies; the k-mer ordering cannot separate subfamilies that differ by two substitutions in 250 bases under 3-12 % private divergence.
2. The peel skips columns where one base holds more than 80 % of the alignment (`PEEL_GLOBAL_CONS`); a diagnostic base held by one subfamily out of eight is below 20 %,
   so it never becomes a feature ("4 features for 320 chunks").
3. A block needs 3 co-occurring features (`PEEL_MIN_BLOCK`), more than a subfamily with two new changes can supply.

Peel sweep (`peel_sweep.sh`, 16 data sets, chunk sizes 10-100, 396 runs eligible with at least 40 chunks): mean ARI over all copies is 0.10 at the defaults
(`GLOBAL_CONS` 0.80, `MIN_SET` 5, `MIN_BLOCK` 3, `FEAT_JACCARD` 0.45); 0.18 with `GLOBAL_CONS` 0.90 and `MIN_BLOCK` 2 (either `FEAT_JACCARD` 0.30 or 0.45), the best single
setting over all data sets. No setting reached the COSEG scores. Best per data set is not a result (it picks the setting after seeing the answer); it is given in `peel_sweep*.tsv`.

### Does the ordering matter for chunk purity? (2026-10-08, KIT, `tools/sim/order_sweep.sh`, chunk size 20)

Chunk purity (share of copies whose chunk majority is their own group) for k-mer trees with k = 4, 5, 6 (default), 8, 10 and for the MAFFT guide tree, on seven simulated sets (seed 1) and the four hand-curated sets:

| set | k4 | k5 | k6 | k8 | k10 | MAFFT |
| --- | --- | --- | --- | --- | --- | --- |
| base | 0.515 | 0.501 | 0.521 | 0.499 | 0.524 | 0.566 |
| diag1 | 0.397 | 0.380 | 0.391 | 0.408 | 0.414 | 0.423 |
| diag5 | 0.474 | 0.472 | 0.496 | 0.510 | 0.491 | 0.535 |
| indel | 0.475 | 0.500 | 0.501 | 0.512 | 0.501 | 0.543 |
| trunc | 0.475 | 0.486 | 0.484 | 0.502 | 0.499 | 0.536 |
| old | 0.352 | 0.347 | 0.342 | 0.364 | 0.351 | 0.372 |
| young | 0.714 | 0.727 | 0.693 | 0.700 | 0.691 | 0.759 |
| saq | 0.487 | 0.550 | 0.574 | 0.631 | 0.622 | 0.507 |
| ccr | 0.514 | 0.595 | 0.591 | 0.596 | 0.580 | 0.565 |
| teu | 0.642 | 0.705 | 0.755 | 0.698 | 0.727 | 0.687 |
| dmo | 0.616 | 0.758 | 0.772 | 0.784 | 0.786 | 0.716 |
| mean | 0.515 | 0.547 | 0.556 | 0.564 | 0.562 | 0.564 |

The default (k = 6) is within 0.01 of the best ordering on average; the MAFFT tree is better on the simulations (by 0.01-0.05 in six of seven sets) and worse on three of the four
real sets (saq, teu, dmo), k = 8-10 better on saq and dmo. The differences are small and not consistent, so the choice of ordering does not explain the low chunk purity: with
two diagnostic changes in 250 bases under 3-12 % private divergence, no ordering tested puts more than about half of a chunk in one subfamily (eight equal subfamilies would give
about 0.2 by chance).
