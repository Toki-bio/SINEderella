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
