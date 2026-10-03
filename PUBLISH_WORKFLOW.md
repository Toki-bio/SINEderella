# SINEderella publish workflow

One logical path from a finished run to a published HTML report (spline
divergence, gallery, PCA, publish alignments, optional discriminator
verdict columns).

## Repos

| Repo | Role |
|---|---|
| **SINEderella** (`Toki-bio/SINEderella`) | Pipeline + `publish_run.sh` + `step6_report.py` |
| **SINE-discriminator** | `verdict.py`, `boundary_justify.py`, `inject_disc_report.py` |
| **Your Pages repo** (e.g. Tal, SINE-discriminator site) | Where you push `report.html` + alignment FASTAs |

Set **`DISC`** to the SINE-discriminator tree (directory containing `verdict.py`).

## Operational order

```
1. SINEderella <genome> <consensus>     # steps 1–4 (+ basic step6 if not --publish)
        │
        ▼
2. publish_run.sh <RUN> <SPECIES>       # or: SINEderella … --publish with SPECIES_CODE set
        │
        ├─ step4 (if *_pctid.tsv missing)
        ├─ publish/align_for_publish.sh
        │     step7 boundary refine
        │     border loop (publish/*.py)
        │     step8a extract MSAs
        │     DISC: rebuild_consensus_row, boundary_justify, trim_display_flanks
        ├─ step6_report.py (--embed-images, pctid spline if TSVs exist)
        └─ DISC: inject_disc_report.py (verdict columns on alignment table)
        │
        ▼
3. results/report.html
        │
        ▼
4. Push HTML + alignments/ to your GitHub Pages tree
```

## Quick start (full publish)

```bash
export DISC=/path/to/SINE_discriminator/site   # or /staging/tmp/sinedisc
export SPECIES_CODE=mysp
export RAW_ALN_BASE=https://raw.githubusercontent.com/org/repo/main/mysp/alignments/
# Optional multi-species site nav:
# export PAGES_INDEX=https://your.pages.host/index.html

# New run with publish at the end:
SINEderella --publish genome.fa consensi.fa

# Or on an existing run:
./publish_run.sh /path/to/run_20260101_120000 mysp
```

## Environment variables

| Variable | Required | Meaning |
|---|---|---|
| `SPECIES_CODE` | with `--publish` | Filename prefix (e.g. `mysp`) |
| `DISC` | for publish alignments | SINE-discriminator root |
| `RAW_ALN_BASE` | for live MSA links | Raw URL through `…/alignments/` |
| `PAGES_INDEX` | optional | “All species” link in report header |
| `USE_DISC` | optional (default 1) | Run `inject_disc_report.py` |
| `SKIP_ALIGN` | optional | Skip align_for_publish if MSAs exist |
| `SKIP_STEP4` | optional | Skip step4 regen |
| `PEEL_FLAGS`, `PEEL_ALN_DIR` | optional | Border-loop hints from peel step4 |

| `SKIP_REBUILD_CONS` | optional | Skip copy-majority rebuild before step4 |
| `SKIP_CANONICALIZE` | optional | Skip RC merge on consensus bank |
| `CANON_MIN_ID` | optional (80) | RC merge threshold (%); same-orientation pairs need 98 % (`--direct-min-id`) and length ratio 0.9 (`--min-len-ratio`) |
| `SKIP_LENGTH_VARIANTS` | optional (0) | `1` skips the length-version test (shorter consensus = 5′ part of a longer one; see docs/LENGTH_VARIANTS.md) |
| `LENGTH_VARIANTS_MAX_PAIRS` | optional (12) | cap on candidate pairs tested per run |
| `SKIP_SATELLITES` | optional (0) | `1` skips the satellite screen inside step 1 (SINE-derived satellites and SINE-containing tandem arrays; docs/SATELLITES.md) |
| `SATELLITE_EXCLUDE` | optional (1) | `0` writes the satellite tables but removes no hits |
| `SATELLITE_EXCLUDE_B` | optional (flagged) | kind-B arrays removed for flagged consensuses only, `all`, or `none` |
| `SKIP_CONSENSUS_AUDIT` | optional (0) | `1` skips the consensus audit (each consensus rebuilt from its assigned copies; see docs/CONSENSUS_AUDIT.md) |
| `CONSENSUS_AUDIT_JOBS` | optional (8) | parallel rebuilds in the consensus audit |

## Consensus bank (RC merge + copy rebuild)

AnnoSINE can emit ± duplicates as separate seed names. Before step1,
`canonicalize_consensus_bank.py` merges clusters: reverse-complement pairs at ≥80% (oma ± pairs are ~83–84% RC) and
same-orientation pairs only at ≥98%, and only when the two lengths are within 90% of each other. A shorter variant is never an alias of a longer one
(2026-10-02: rsi r9 105 bp, r7 154 bp, MEG-RS 135 bp and MEG-RL 207 bp were merged at 89.5/97 % by a position-by-position identity, and the kept name carried the other's sequence)
and orients to AT-rich 3′. Before step4, `rebuild_consensus_bank.py` writes
`consensuses.rebuilt.fa` from assigned copies (no N ties). step4 pctid uses
`-3` (both strands), matching step2.

Oma repair procedure: [docs/OMA_CONSENSUS_REPAIR.md](docs/OMA_CONSENSUS_REPAIR.md).

### Tandem arrays

After assignment `tools/array_flag.py` writes `results/array_flag.tsv` (share of each family's copies in tandem arrays of regular spacing); step8a puts independent copies on the top-100 plate first and marks array rows `[array]`; the report shows "Tandem array" for a flagged family. See docs/ARRAYS.md.

### Consensus audit

After the length-version test, `tools/consensus_audit.py` rebuilds every consensus of the bank from its own assigned copies (bootstrap subsamples, two seeds; `tools/vendor/sine_consensus.sh`) and writes `results/consensus_audit/summary.tsv`: mismatches and gap columns against the bank, and a verdict `MATCH` / `SHORTER` / `LONGER` / `DIVERGED` / `UNSTABLE` / `SKIPPED`. The table is shown in the report. It is a decision input; the bank is never changed. Method and limits: docs/CONSENSUS_AUDIT.md.

### Length versions

A shorter consensus that is the 5' part of a longer one is never merged (see above). After assignment, `tools/length_variants_run.py` tests every such
pair on the copies (end modes, internal bases following the length, TSD after each end, residual similarity) and writes
`results/length_variants/summary.tsv`: `TWO_VERSIONS` / `UNLINKED_ENDS` / `SINGLE_MODE` / `UNRESOLVED`. Method, calibration on rle MEG-RS/MEG-RL and
rsi r9/r7/r8, limits: docs/LENGTH_VARIANTS.md. The verdict is a decision input; nothing is merged or renamed automatically.

## Divergence chart

When step4 `*_pctid.tsv` files exist, `step6_report.py` uses **variant 3**:
1% bins, Y = copy count, spline through bin midpoints (same metric as gallery PNGs).
Without pctid TSVs it falls back to bitscore KDE + violins.

## Discriminator overlay

`inject_disc_report.py` replaces the alignment table with copy counts + MSA links +
Flanks / Flank context / Element / Overall chips from `verdict.py`.

Skip with `USE_DISC=0` if you only want stock SINEderella links.

## Reproducible alignments (mafft `--threadit 0`)

mafft with iterative refinement (`--maxiterate`, and `--auto`, which picks it) gives a different alignment on every run when `--thread` > 1, because the refinement threads race. Tested on real rsi copies (top100 of r5, r7, r10): every such call differed between two identical runs; progressive calls (`--retree N --maxiterate 0`) did not. All iterative calls in the pipeline now pass `--threadit 0` (same result for 1, 3, 4, 16, 32 threads, no measurable slowdown). Alignments from runs made before this change (2026-10-02) are not reproducible bit for bit; the counts, assignments and length-version verdicts do not depend on them (ssearch36).
