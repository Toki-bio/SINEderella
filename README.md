# SINEderella

A reproducible Bash pipeline for genome-wide identification, classification, and characterization of SINE transposable elements.

## Pipeline Steps

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `step1_search_extract.sh` | Search genome for SINE hits (≥0.9 length, ≥0.65 similarity), merge overlaps, extract per-family FASTA/BED |
| 2 | `step2_asSINEment.sh` | Assign sequences to known subfamily consensuses via repeated ssearch36 cycles with stability rule |
| 3 | `step3_postprocess.sh` | Postprocess assignments: compute per-subfamily statistics, generate summary tables |
| 4 | `step4_plots.sh` | Generate divergence plots and conservation statistics per subfamily |
| 5 | `step5_align_subfamilies.sh` | Cross-species subfamily alignment against consensus bank |
| 6 | `step6_report.sh` / `step6_report.py` | HTML report (`results/report.html`): subfamily table, divergence, plates, element hierarchy, similarity blocks with the length-version, array, satellite and consensus-audit tables |
| 7 | `step7_boundary_refine.sh` | Standalone/modular: per-subfamily boundary refinement — stepwise flank extension until a fraction-of-pairs-above-threshold test confirms background-level identity (or hits a 1000bp cap), writes `boundary_refinement.tsv` |
| 8a | `step8a_extract_alignments.sh` | Standalone/modular: builds real top100/rand100/subfam alignments per subfamily, using `boundary_refinement.tsv` (if present) to size each subfamily's flanks |
| 1c | `step1c_deplete.sh` | Standalone/modular, after step2+3: keeps the copies the bank does not explain (unassigned + each subfamily's copies below the P-th percentile of step3 `sim_ratio`, both only if ≥ 0.8 × their consensus length; shorter ones are set aside), runs SubFam with bin size 20 on that residual, writes `step1c/subfam_input/input.clw.al` for the §6.1 review; for subfamilies too rare to fill a chunk in step1's 30k sample. Close the loop with `SINEderella --add` |
| 8b | `step8b_publish_report.sh` | Standalone/modular: wires step8a's alignments into an existing `report.html` as MSA-viewer links |
| — | `publish_run.sh` | **Full publish:** step4 (if needed) → `publish/align_for_publish.sh` → step6 → SINE-discriminator inject. See [PUBLISH_WORKFLOW.md](PUBLISH_WORKFLOW.md). |

Steps 7/8a/8b are runnable independently against any completed run
(`RUN_ROOT`) — not yet wired into the `SINEderella`/`SINEderella_multi`
main-orchestrator auto-run sequence. See
[RESEARCH_DIRECTIONS.md](RESEARCH_DIRECTIONS.md) for the modular-step
design philosophy behind this.

## Main Wrappers

- **`SINEderella`** — Single-species 5-step pipeline
- **`SINEderella_multi`** — Multi-genome wrapper (processes species list from TSV)

## Search Engines

- **`sear`** — Single-query SINE search with fragment scanning (ssearch36 + bedtools)
- **`sear_multi`** — Multi-query variant: all consensuses in one ssearch36 pass per genome fragment

## Supporting Tools

| Script | Purpose |
|--------|---------|
| `asSINEment` | Standalone subfamily assignment engine (earlier version of step2 logic) |
| `SubFam` | Compression of many copies to chunk consensuses (k-mer guide-tree order, chunks, plurality consensus, final alignment); wrapper around SubFam 1.2.0 in `tools/vendor/SubFam.sh` with the old interface; the subfamily call itself is made on its output by a person (MANUAL §6.1) or the peel; the previous script is `tools/vendor/SubFam.old.sh` |
| `sine_consensus.sh` | Bootstrap consensus builder (gaps excluded from denominator) |
| `sine_consensus_smart.sh` | Enhanced consensus with variance detection and early stopping |
| `analyze_convergence.sh` | Post-run convergence quality analysis |
| `plot_subfamily.py` | Divergence/conservation plots (matplotlib) |
| `extract_alignments.sh` | Generate tiered representative alignments (core, best50, SubFam, evidence) |
| `extract_subfam_only.sh` | Extract subfamily-only sequences |
| `benchmark_sear.sh` | Benchmark sear vs sear_multi performance |
| `run_step5_wrapper.sh` | Step 5 batch runner |
| `step5_direct.sh` | Direct single-run subfamily alignment (simplified step 5) |
| `tools/length_variants_run.py` | After assignment: is a consensus that is the 5′ part of a longer one a separate SINE or the same element with a worn 3′ end (end modes, linkage, TSD on the copies); [docs/LENGTH_VARIANTS.md](docs/LENGTH_VARIANTS.md) |
| `tools/array_flag.py`, `tools/array_order.py` | Tandem arrays: copies in runs of regular spacing, a per-family flag and plate selection that takes independent copies first; [docs/ARRAYS.md](docs/ARRAYS.md) |
| `tools/consensus_audit.py` | After assignment: every consensus rebuilt from its own assigned copies and compared with the bank; uses `tools/vendor/sine_consensus.sh` (seeded copy of the current SINE_consensus bootstrap builder); [docs/CONSENSUS_AUDIT.md](docs/CONSENSUS_AUDIT.md) |
| `tools/satellite_stage.py` (+ `satellite_screen.py`, `satellite_trf_verify.py`, `satellite_kindB_verify.py`) | Inside step 1, before extraction: SINE-derived satellites (hit windows → TRF → unit vs consensus) and tandem arrays of a longer unit (regular spacing, unit identity) leave the hit set as loci, recorded in `results/satellites/`; [docs/SATELLITES.md](docs/SATELLITES.md) |
| `flankscan/fs9_run.sh` (+ `fs9_twins.sh`) | After assignment: the flanks of all firmly assigned copies of every family compared with one another; copies with shared flanks (segmental duplications, missed array units) are listed in `results/flank_twins.tsv`, marked `[twin]` on the plates and put after independent copies in the top 100; [docs/FLANK_UNIQUENESS.md](docs/FLANK_UNIQUENESS.md) |
| `tools/consensus_blocks.py`, `report_blocks.py` | Similarity blocks between consensuses (matrix + pair view in the report) |
| `tools/tsd_curve.py` | TSD share of a plate's copies against shuffled pairs, by minimum length |
| `canonicalize_consensus_bank.py`, `consensus_bank_lib.py` | Bank cleaning before step 1: RC duplicates merged (same-orientation pairs only at 98 %, lengths within 90 %), simple-repeat tail oriented 3′; length-version candidate pairs |
| `flankscan/` | Flank analysis in bash/awk: tandem repeats around copies, composites from junction peaks, candidate elements, chains, hierarchy, singletons/TSD, flank twins; `flankscan/HANDOFF.md` |
| `import_squamata_run.py` | Import SINEderella run results into SINEdb data format (requires sine-kb models) |

## Standalone and historical scripts

Not called by `SINEderella`, `SINEderella_multi` or `publish_run.sh`; kept for the analyses that used them and runnable by hand:
`asSINEment` (the earlier assignment engine), `sear_multi` and `benchmark_sear.sh`, `step4_diagnostic.sh/.py` (needs pandas, scipy,
scikit-learn), `step5_align_subfamilies.sh`, `step5_direct.sh`, `run_step5_wrapper.sh`, `run_subfam_per_sf.sh`, `extract_alignments.sh`
and `extract_subfam_only.sh` (the multi-run tier extractor that preceded step 8a), `sine_consensus.sh`, `sine_consensus_smart.sh`,
`sine_pairwise_consensus.sh`, `analyze_convergence.sh`, `step1b_cluster_subfamilies_assist.sh` with `cluster_assist.js` /
`subfam_cluster_lib.js` (MANUAL §6.1), `flank_border_consensus_test.sh/.py`, `publish/flank_border_iterate.py` as a script (its
functions are used by the border loop), `import_squamata_run.py`, `tools/composite_scan.py` and `tools/build_composite.py` (the Python
prototype of flankscan, docs/COMPOSITES.md), `tools/compare_cons.py`, `workflow.html`. The audit of 2026-10-05 (docs/AUDIT_2026-10-05.md)
did not re-test these.

## Dependencies

- `ssearch36` (FASTA36 package)
- `mafft`
- `bedtools`
- `samtools`
- `seqkit`
- `cons` and `seqret` (EMBOSS; only the old `tools/vendor/SubFam.old.sh` needs them, SubFam 1.2.0 does not)
- `trf` (Tandem Repeats Finder; the satellite stage is skipped with a log line without it), `gawk` (consensus audit, flankscan), `dustmasker` (flankscan)
- `jellyfish` (conda `kmer-jellyfish`): 20-mer counts for the flank-twin check on genomes over 300 Mb; without it that check is skipped with a log line on large genomes
- Python 3 with numpy and matplotlib (plots, report panels, consensus blocks)
- Tests: `python -m pytest` from the repo root (`pytest.ini`; the end-to-end tests need mafft, gawk, samtools and are skipped without them)

## Documentation

- [MANUAL.md](MANUAL.md) — Full user manual
- [docs/LENGTH_VARIANTS.md](docs/LENGTH_VARIANTS.md) — separate SINE or decayed 3′ end; [docs/CONSENSUS_AUDIT.md](docs/CONSENSUS_AUDIT.md) — consensus vs its own copies; [PUBLISH_WORKFLOW.md](PUBLISH_WORKFLOW.md) — publish flow, env vars, reproducible alignments (`mafft --threadit 0`)
- [QUALITY_FLAGGING_README.md](QUALITY_FLAGGING_README.md) — Consensus convergence QC system
- [PLAN_alignment_viewer.md](PLAN_alignment_viewer.md) — Alignment tier design rationale
- [ALIGNMENT_DEPLOYMENT.md](ALIGNMENT_DEPLOYMENT.md) — Deployment checklist for alignments
- [RESEARCH_DIRECTIONS.md](RESEARCH_DIRECTIONS.md) — Gradual (SNP/indel) vs. modular ("repeat pangenome"/panconsensus) models of repeat divergence — an open research pathway, not yet implemented

## Related Repositories

- [SINEdb](https://github.com/Toki-bio/SINEdb) — SINE database (GitHub Pages)
- [SINE_consensus](https://github.com/Toki-bio/SINE_consensus) — Standalone consensus tools
- [SubFam](https://github.com/Toki-bio/SubFam) — Subfamily identification tool
