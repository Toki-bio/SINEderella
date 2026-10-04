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
| `SubFam` | Subfamily identification via chunk-sort-consensus-align |
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
| `tools/consensus_blocks.py`, `report_blocks.py` | Similarity blocks between consensuses (matrix + pair view in the report) |
| `tools/tsd_curve.py` | TSD share of a plate's copies against shuffled pairs, by minimum length |
| `canonicalize_consensus_bank.py`, `consensus_bank_lib.py` | Bank cleaning before step 1: RC duplicates merged (same-orientation pairs only at 98 %, lengths within 90 %), simple-repeat tail oriented 3′; length-version candidate pairs |
| `flankscan/` | Flank analysis in bash/awk: tandem repeats around copies, composites from junction peaks, candidate elements, chains, hierarchy, singletons/TSD, flank twins; `flankscan/HANDOFF.md` |
| `import_squamata_run.py` | Import SINEderella run results into SINEdb data format (requires sine-kb models) |

## Dependencies

- `ssearch36` (FASTA36 package)
- `mafft`
- `bedtools`
- `samtools`
- `seqkit`
- `cons` (EMBOSS)
- `trf` (Tandem Repeats Finder; the satellite stage is skipped with a log line without it), `gawk` (consensus audit, flankscan), `dustmasker` (flankscan)
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
