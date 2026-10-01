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

## Consensus bank (RC merge + copy rebuild)

AnnoSINE can emit ± duplicates as separate seed names. Before step1,
`canonicalize_consensus_bank.py` merges clusters: reverse-complement pairs at ≥80% (oma ± pairs are ~83–84% RC) and
same-orientation pairs only at ≥98%, and only when the two lengths are within 90% of each other. A shorter variant is never an alias of a longer one
(2026-10-02: rsi r9 105 bp, r7 154 bp, MEG-RS 135 bp and MEG-RL 207 bp were merged at 89.5/97 % by a position-by-position identity, and the kept name carried the other's sequence)
and orients to AT-rich 3′. Before step4, `rebuild_consensus_bank.py` writes
`consensuses.rebuilt.fa` from assigned copies (no N ties). step4 pctid uses
`-3` (both strands), matching step2.

Oma repair procedure: [docs/OMA_CONSENSUS_REPAIR.md](docs/OMA_CONSENSUS_REPAIR.md).

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
