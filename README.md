# SINEderella

A reproducible Bash pipeline for genome-wide identification, classification, and characterization of SINE transposable elements.

## Pipeline Steps

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `step1_search_extract.sh` | Search genome for SINE hits (≥0.9 length, ≥0.65 similarity), merge overlaps, extract per-family FASTA/BED |
| 2 | `step2_asSINEment.sh` | Assign sequences to known subfamily consensuses: one ssearch36 search with fixed statistics, a clear-winner rule and a bitscore threshold |
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
| `SubFam` | Compression of many copies to chunk consensuses (k-mer guide-tree order, chunks, plurality consensus, final alignment); wrapper around SubFam 1.5.1 in `tools/vendor/SubFam.sh` with the old interface; the subfamily call itself is made on its output by a person (MANUAL §6.1) or the peel; the previous script is `tools/vendor/SubFam.old.sh` |
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
- `cons` and `seqret` (EMBOSS; only the old `tools/vendor/SubFam.old.sh` needs them, SubFam 1.2+ does not)
- `trf` (Tandem Repeats Finder; the satellite stage is skipped with a log line without it), `gawk` (consensus audit, flankscan), `dustmasker` (flankscan)
- `jellyfish` (conda `kmer-jellyfish`): 20-mer counts for the flank-twin check, made once per genome; without it genomes up to 300 Mb fall back to a slow awk count per family (about 0.3 GB of memory per Mb), and on larger genomes that check is skipped with a log line
- Python 3 with numpy and matplotlib (plots, report panels, consensus blocks)
- Tests: `python -m pytest` from the repo root (`pytest.ini`; the end-to-end tests need mafft, gawk, samtools and are skipped without them)
- `bc` (the consensus builder `sine_consensus.sh`: consensus audit, border loop); `minimap2` (`sear`, only for queries of 1000 bp or longer; without it `sear` falls back to ssearch36 with a warning); a C compiler (`cc`/`gcc`) is used by SubFam for its k-mer ordering program when present, otherwise the Python version gives the same order

All of these with the versions of the current server environment are in [environment.yml](environment.yml) (conda-forge + bioconda): `mamba env create -f environment.yml`.

## Run with a container

The image holds the tools of `environment.yml` and this repository in `/opt/SINEderella`; nothing has to be installed on the host.

```bash
# Apptainer (HPC, no root). Build from the repository root:
apptainer build --fakeroot --build-arg PIPELINE_COMMIT=$(git describe --always --dirty) sinederella.sif containers/sinederella.def
apptainer run  sinederella.sif genome.fa consensuses.fa            # = SINEderella genome.fa consensuses.fa
apptainer exec sinederella.sif publish_run.sh run_20260101_120000 mysp
apptainer exec --pwd /opt/SINEderella sinederella.sif python -m pytest -p no:cacheprovider   # the unit tests inside the image

# Docker: the same image
docker build -f containers/Dockerfile --build-arg PIPELINE_COMMIT=$(git describe --always --dirty) -t sinederella .
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work sinederella SINEderella genome.fa consensuses.fa
```

The run directory `run_<date>_<time>` is written to the working directory, as in a native run. Environment variables (`THREADS`, `SPECIES_CODE`, `PUBLISH`, `SKIP_*`, `TMPDIR`, ...) are read from the host shell by Apptainer; with Docker pass each with `-e`. Apptainer binds the home and the working directory; a genome elsewhere needs `--bind /data` (Docker: another `-v`). `PIPELINE_COMMIT` goes into `manifest.txt` of every run (there is no git in the image).

## Host-specific paths

No path of a particular server is needed. The ones that remain are defaults that are used only where they exist:

| Variable | Default | Used by | Meaning |
|---|---|---|---|
| `SINEDERELLA_TOOLCHAIN` | `/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin` (DRAGEN) | `publish/align_for_publish.sh`, `publish/verdict_chain.sh`, `publish/flank_border_iterate.py` | bin directories put first on `PATH` when the first one exists; set to an empty string to leave `PATH` alone (the container images do) |
| `DISC` (or `SINE_DISC`) | `disc/` of this repository; after it `/staging/tmp/sinedisc`, `$HOME/SINE_discriminator/site`, `$HOME/SINE_discriminator` | `publish_run.sh`, `publish/align_for_publish.sh` | SINE-discriminator modules; the later candidates are reached only if `disc/` is missing |
| `SINE_CONSENSUS_SH` | `/staging/tmp/SINEderella/sine_consensus.sh` | `publish/border_loop_subfam.py` run by hand (`--sine-script` overrides it; `verdict_chain.sh` always passes the repository's `sine_consensus.sh`) | consensus builder for the border loop |
| `PIPELINE_COMMIT` file | — | `SINEderella`, `publish/verdict_chain.sh` | read from the script directory for the manifest when `git describe` does not work there |

Scripts outside the main chain (`run_subfam_per_sf.sh`, `flank_border_consensus_test.sh/.py`, `tools/coseg_compare/`, `tests/bench/`, `tests/toy/run_toy.sh`) still carry the paths of the servers they were run on.

## Documentation

- [MANUAL.md](MANUAL.md) — Full user manual
- [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) — the same input gives the same output: why it did not (ssearch36 shuffle statistics), the fix (fixed statistics, a clear-winner rule), the two bugs that made the environment change the output (sort collation, ssearch36 and long paths; section 11) and the tests
- [docs/AIM.md](docs/AIM.md) — what the workflow is for; [docs/VERDICTS.md](docs/VERDICTS.md) — the verdict of every copy (SINE / grey zone / not SINE), `loci_verdict.tsv`, and how plates and page use it
- [docs/CONTINUATION.md](docs/CONTINUATION.md) — copies that continue past their ends into another consensus (a missing composite, or a copy labelled with a shorter consensus), checked in every run
- [disc/README.md](disc/README.md) — the SINE-discriminator modules the publish chain runs (kept in this repository since 2026-10-21)
- [docs/LENGTH_VARIANTS.md](docs/LENGTH_VARIANTS.md) — separate SINE or decayed 3′ end; [docs/CONSENSUS_AUDIT.md](docs/CONSENSUS_AUDIT.md) — consensus vs its own copies; [PUBLISH_WORKFLOW.md](PUBLISH_WORKFLOW.md) — publish flow, env vars, reproducible alignments (`mafft --threadit 0`)
- [QUALITY_FLAGGING_README.md](QUALITY_FLAGGING_README.md) — Consensus convergence QC system
- [PLAN_alignment_viewer.md](PLAN_alignment_viewer.md) — Alignment tier design rationale
- [ALIGNMENT_DEPLOYMENT.md](ALIGNMENT_DEPLOYMENT.md) — Deployment checklist for alignments
- [RESEARCH_DIRECTIONS.md](RESEARCH_DIRECTIONS.md) — Gradual (SNP/indel) vs. modular ("repeat pangenome"/panconsensus) models of repeat divergence — an open research pathway, not yet implemented

## Related Repositories

- [SINEdb](https://github.com/Toki-bio/SINEdb) — SINE database (GitHub Pages)
- [SINE_consensus](https://github.com/Toki-bio/SINE_consensus) — Standalone consensus tools
- [SubFam](https://github.com/Toki-bio/SubFam) — Subfamily identification tool
