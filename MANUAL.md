# SINEderella — User Manual

**Version:** 2026-02
**Pipeline for genome-wide SINE family detection, subfamily assignment, and annotation**

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation & Dependencies](#2-installation--dependencies)
3. [Quick Start](#3-quick-start)
4. [Pipeline Architecture](#4-pipeline-architecture)
5. [Running Modes](#5-running-modes)
   - 5.1 [Full Mode](#51-full-mode-default)
   - 5.2 [Add Mode](#52-add-mode---add--a)
   - 5.3 [Exclude Mode](#53-exclude-mode---exclude--e)
6. [Step 1 — Search & Extract](#6-step-1--search--extract-step1_search_extractsh)
7. [Step 2 — Assignment (asSINEment)](#7-step-2--assignment-asssinement-step2_asssinementsh)
   - 7.1 [Full Assignment](#71-full-assignment)
   - 7.2 [Incremental Assignment](#72-incremental-assignment)
8. [Step 3 — Postprocessing](#8-step-3--postprocessing-step3_postprocesssh)
9. [Step 4 — Diagnostic Plots](#9-step-4--diagnostic-plots-step4_plotssh)
10. [Step 5 — Subfamily Alignments](#10-step-5--subfamily-alignments-step5_align_subfamiliessh)
11. [SINEderella_multi — Multi-Genome Mode](#11-sinederella_multi--multi-genome-mode)
12. [SINE-KB Report & Import](#12-sine-kb-report--import)
13. [GitHub Pages (SINEdb) Export & Re-Import Workflow](#13-github-pages-sinedb-export--re-import-workflow)
14. [Output Files Reference](#14-output-files-reference)
15. [Environment Variables](#15-environment-variables)
16. [Algorithm Details](#16-algorithm-details)
    - 16.1 [Sear Search](#161-sear-search)
    - 16.2 [10-Cycle Voting](#162-10-cycle-voting)
    - 16.3 [Bitscore Thresholds](#163-bitscore-thresholds)
    - 16.4 [Incremental BED-Overlap Strategy](#164-incremental-bed-overlap-strategy)
    - 16.5 [LEAK Detection](#165-leak-detection)
    - 16.6 [CONFLICT Annotation](#166-conflict-annotation)
    - 16.7 [Soft Assignment](#167-soft-assignment-of-unassigned-sequences)
    - 16.8 [Similarity Ratio (sim_ratio)](#168-similarity-ratio-sim_ratio)
17. [Directory Structure](#17-directory-structure)
18. [Troubleshooting](#18-troubleshooting)
19. [Examples](#19-examples)

---

## 1. Overview

SINEderella is a Bash-based bioinformatics pipeline that identifies SINE (Short Interspersed Nuclear Element) copies across a genome and assigns each copy to its correct subfamily consensus. The pipeline operates in three stages:

1. **Search & Extract** — Scans the genome with known SINE consensus sequences using `sear`, merges overlapping hits, and extracts candidate SINE sequences.
2. **Assignment** — Assigns each extracted sequence to a subfamily using 10 independent `ssearch36` voting cycles and strict bitscore thresholds.
3. **Postprocessing** — Computes per-copy similarity scores (single ssearch36 alignment vs. consensus), builds a comprehensive annotation table with LEAK detection, CONFLICT flagging, soft-assignment of unassigned sequences, and per-subfamily summary statistics.
4. **Diagnostic Plots** — Generates per-subfamily divergence histograms and nucleotide composition charts (PNG + PDF).

SINEderella supports three running modes:
- **Full** — Run the complete pipeline from scratch on a genome + consensus set.
- **Add** — Add new SINE consensus(es) to an existing run without repeating the full search for known families.
- **Exclude** — Remove a SINE family from an existing run and reassign the remaining sequences.

---

## 2. Installation & Dependencies

### Required Tools

| Tool | Purpose | Minimum Version |
|------|---------|-----------------|
| `sear` | Genome-wide SINE homology search | — |
| `ssearch36` | Smith-Waterman local alignment (FASTA suite) | 36.x |
| `seqkit` | FASTA/FASTQ manipulation | 2.x |
| `bedtools` | Genomic interval operations | 2.30+ |
| `samtools` | FASTA indexing | 1.10+ |
| `SubFam` | Sequence clustering | — |
| `mafft` | Multiple sequence alignment | 7.x |
| `python3` | Plotting (step 4) | 3.8+ |
| `matplotlib` | Python plotting library (step 4) | 3.x |
| `numpy` | Python numerical library (step 4) | 1.20+ |
| `awk`, `sort`, `grep`, `sed` | Standard Unix text processing | — |

### Installation

```bash
# Clone the repository
git clone https://github.com/Toki-bio/SINEderella-dev.git
cd SINEderella-dev

# Ensure the main script is executable
chmod +x SINEderella step1_search_extract.sh step2_asSINEment.sh step3_postprocess.sh step4_plots.sh step5_align_subfamilies.sh

# Verify all dependencies are available
for tool in sear ssearch36 seqkit bedtools samtools SubFam mafft python3; do
  command -v "$tool" >/dev/null && echo "OK: $tool" || echo "MISSING: $tool"
done
python3 -c "import matplotlib, numpy" && echo "OK: matplotlib+numpy" || echo "MISSING: python3 matplotlib/numpy"
```

All scripts (`SINEderella`, `step1_search_extract.sh`, `step2_asSINEment.sh`, `step3_postprocess.sh`, `step4_plots.sh`, `step5_align_subfamilies.sh`, `plot_subfamily.py`, `SubFam`) must be in the same directory.

---

## 3. Quick Start

### Full pipeline run
```bash
SINEderella genome.fa my_sine_consensuses.fa
```

### Add a new SINE family to an existing run
```bash
SINEderella --add new_family.fa
# or specify which run to extend:
SINEderella --add new_family.fa --run run_20260207_033350
```

### Remove a SINE family from an existing run
```bash
SINEderella --exclude Heno
# or specify which run:
SINEderella --exclude Heno --run run_20260207_033350
```

---

## 4. Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              SINEderella                                    │
│                         (master orchestrator)                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────┐  ┌─────────────────┐  ┌────────────────┐  ┌───────────┐ │
│  │  Step 1       │  │  Step 2          │  │  Step 3        │  │  Step 4   │ │
│  │  Search &     │─▶│  Assignment      │─▶│  Postprocess   │─▶│  Plots    │ │
│  │  Extract      │  │  (asSINEment)    │  │                │  │           │ │
│  │               │  │                  │  │                │  │           │ │
│  │  sear         │  │  ssearch36 ×10   │  │  sim_ratio     │  │ divergence│ │
│  │  bedtools     │  │  unanimous vote  │  │  ALL table     │  │ histogram │ │
│  │  SubFam       │  │  bitscore thr    │  │  LEAK detect   │  │ nucleotide│ │
│  │  MAFFT        │  │                  │  │  CONFLICT flag │  │ frequency │ │
│  │               │  │                  │  │  soft-assign   │  │ (png+pdf) │ │
│  └──────────────┘  └─────────────────┘  └────────────────┘  └───────────┘ │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 5. Running Modes

### 5.1 Full Mode (default)

```bash
SINEderella <GENOME_FASTA> <CONSENSUS_FASTA>
```

Runs the complete three-step pipeline from scratch:

1. **Input sanitization** — Cleans genome and consensus FASTA files (normalizes headers, removes whitespace/carriage returns, wraps at 60bp/line).
2. **Step 1** — Splits each consensus into individual query files, runs `sear` for each against the genome, merges overlapping BED intervals, extracts sequences, samples up to 30,000 for clustering, runs `SubFam` + `MAFFT`.
3. **Step 2** — Runs 10-cycle `ssearch36` voting and bitscore threshold assignment.
4. **Step 3** — Postprocesses into the final annotation table.

**Creates run directory:** `run_YYYYMMDD_HHMMSS/`

### 5.2 Add Mode (`--add` / `-a`)

```bash
SINEderella --add <newSINE.fas> [--run <run_dir>]
```

Adds one or more new SINE consensus sequences to an existing run. This mode is designed to be fast: it searches only for the new consensus(es) and uses **incremental assignment** to avoid redundant computation.

**Workflow:**

1. Identifies which consensus names are truly new (not already in the source run).
2. Copies step 1 data from the source run.
3. Runs `sear` only for the new consensus sequences.
4. Rebuilds the merged BED and re-extracts all sequences.
5. Runs **incremental step 2**: carries forward firm old assignments that don't overlap the new-consensus BED regions, and runs full 10-cycle ssearch36 only on sequences that are new, previously unassigned, or overlapping new BED regions.
6. Runs step 3 postprocessing.

**Creates run directory:** `run_add_YYYYMMDD_HHMMSS/`

**Incremental speedup:** In testing, ~75% of sequences are carried forward without ssearch36, reducing computation time proportionally.

If no previous assignment file is found, the mode falls back to full reassignment automatically.

### 5.3 Exclude Mode (`--exclude` / `-e`)

```bash
SINEderella --exclude <SINE_name> [--run <run_dir>]
```

Removes a SINE family from an existing run and reassigns remaining sequences.

**Workflow:**

1. Verifies the family exists in the source run's consensus file.
2. Creates a new consensus file without the excluded family.
3. Copies step 1 data and removes all BED/search files belonging to the excluded family.
4. Rebuilds the merged BED and re-extracts sequences from remaining search data.
5. Runs full step 2 assignment (no incremental mode — removal can shift all thresholds).
6. Runs step 3 postprocessing.

**Creates run directory:** `run_excl_YYYYMMDD_HHMMSS/`

---

## 6. Step 1 — Search & Extract (`step1_search_extract.sh`)

### Usage (standalone)
```bash
step1_search_extract.sh <consensus.fa> <genome.fa> [sample_size] [bin_size]
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `consensus.fa` | *required* | FASTA file of SINE consensus sequences |
| `genome.fa` | *required* | Genome FASTA file |
| `sample_size` | 30000 | Max sequences to sample for SubFam clustering |
| `bin_size` | 50 | SubFam bin size parameter |

### What It Does

1. **Index genome** — Creates `.fai` index via `samtools faidx` if not present.
2. **Split consensus** — Uses `seqkit split -i` to create one FASTA file per consensus sequence.
3. **Sear search** — For each consensus query, runs:
   ```
   sear -k <query> <genome> 0.8 65 50
   ```
   - `-k` flag: keeps genome split parts between queries for reuse.
   - `0.8`: minimum length fraction (SFL) — a hit must cover ≥80% of the query length.
   - `65`: minimum similarity percentage.
   - `50`: additional sear parameter.
4. **Build labeled BED** — Combines all per-query BED files into `all_hits.labeled.bed` (6-column: chr, start, end, subfamily, score, strand).
5. **Conflict regions** — Runs `bedtools merge` with `-c 4 -o distinct` to create `regions.by_subfam.bed`, identifying genomic regions where multiple subfamily hits overlap.
6. **Merge intervals** — Combines all BED files and merges overlapping intervals into `merged_hits.bed`.
7. **Extract sequences** — Uses `bedtools getfasta -s` (strand-aware) to extract SINE candidate sequences.
8. **Sample** — If more than `sample_size` sequences, randomly samples using `seqkit sample`.
9. **SubFam clustering** — Runs `SubFam` on sampled sequences to identify clusters.
10. **MAFFT alignment** — Aligns cluster representatives + consensus sequences using `mafft --localpair --maxiterate 1000`.

### Output Directory: `genome.clean_step1/`

| File/Directory | Description |
|----------------|-------------|
| `extracted.fasta` | All extracted SINE candidate sequences |
| `merged_hits.bed` | Merged genomic intervals |
| `sampled_30000.fasta` | Sampled subset for clustering |
| `genome.path` | Path to the genome file used |
| `consensus.path` | Path to the consensus file used |
| `run_info.txt` | Run metadata |
| `searches/` | Per-query BED files, labeled BED, conflict regions |
| `searches/all_hits.labeled.bed` | All hits with subfamily labels (6-col BED) |
| `searches/regions.by_subfam.bed` | Merged regions with comma-separated subfamily lists |
| `subfam_input/` | SubFam clustering output + MAFFT alignment |

---

## 7. Step 2 — Assignment (asSINEment) (`step2_asSINEment.sh`)

### Usage (standalone)
```bash
step2_asSINEment.sh [options] <consensus.fa> <sines.fa> [threads] [output_dir]
```

### Options

| Flag | Description |
|------|-------------|
| `--incremental` | Enable incremental mode (carry forward old assignments) |
| `--old-assignment <file>` | Path to previous `assignment_full.tsv` (required with `--incremental`) |
| `--new-beds <file>` | BED file of new consensus hit regions (used with `--incremental`) |
| `--help` | Show usage information |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `consensus.fa` | *required* | FASTA of subfamily consensus sequences |
| `sines.fa` | *required* | FASTA of SINE sequences to assign |
| `threads` | auto (`nproc`) | CPU threads for ssearch36 |
| `output_dir` | `step2_output` | Output directory name |

### Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `SPLIT_SIZE` | 20000 | Chunk size for splitting large SINE files |
| `MIN_REL_BITSCORE` | 0.45 | Relative bitscore threshold multiplier |

### 7.1 Full Assignment

The full assignment algorithm proceeds through these stages:

1. **Input preparation** — Clean headers, split large SINE files into chunks of 20,000 sequences.
2. **10-cycle ssearch36 voting** — Run `ssearch36` (Smith-Waterman local alignment) 10 independent times. Each cycle:
   - Aligns all consensus sequences against all SINE sequences.
   - Records the best-hit consensus for each SINE sequence.
   - The stochastic element of ssearch36 (Z-score normalization, shuffle seed) means each cycle may produce slightly different bitscores.
3. **Vote aggregation** — Count how many of the 10 cycles each sequence voted for each consensus.
4. **Unanimous filter** — Only sequences that voted for the same consensus in all 10/10 cycles pass.
5. **Bitscore threshold** — For each consensus:
   - Sort all unanimous sequences by bitscore (descending).
   - Take TopN = bitscore of the min(10, count)-th best sequence.
   - Threshold = `MIN_REL_BITSCORE × TopN` = `0.45 × TopN`.
   - Reject sequences below this threshold as `rejected_low_bitscore`.
6. **Output generation** — Create assigned/unassigned FASTA files, per-subfamily FASTA splits, statistics, and full assignment details.

### 7.2 Incremental Assignment

When called with `--incremental`, step 2 partitions sequences into two tiers:

**Tier 1 — Carry Forward** (no ssearch36 needed):
- Sequences that were firmly assigned (10/10 votes, status=assigned) in the old run.
- AND are still present in the current extracted set.
- AND do NOT overlap with any new-consensus BED regions.
- Their old assignment data is carried forward unchanged.

**Tier 2 — Re-evaluate** (full 10-cycle ssearch36):
- Sequences that are new (not in old run).
- Sequences that were previously unassigned/rejected.
- Sequences whose genomic region overlaps new-consensus BED hits.

After both tiers are resolved, the results are merged and thresholds are recalculated from the combined unanimous set to produce output files identical in format to full mode.

---

## 8. Step 3 — Postprocessing (`step3_postprocess.sh`)

### Usage (standalone)
```bash
step3_postprocess.sh <RUN_ROOT> [THREADS]
```

### What It Does

Step 3 transforms the raw assignment data into a comprehensive annotation table with quality flags and per-copy similarity scores.

#### Stage 0: Self-Bits Table
Creates `self_bits.tsv` — per-subfamily maximum threshold from `assignment_full.tsv`. Used as a reference for the highest-quality assignments.

#### Stage 0b: Single ssearch36 Similarity Scores
For each assigned copy, runs a single `ssearch36` alignment against its assigned consensus to compute a normalized similarity ratio:

1. **Real self-bits** (`self_bits_real.tsv`): Aligns each consensus against itself via `ssearch36` to obtain the maximum possible bitscore per subfamily.
2. **Per-copy alignment**: Groups assigned copies by subfamily, runs `ssearch36 consensus.fa copies.fa` per subfamily, extracts the best bitscore per copy.
3. **Similarity ratio**: `sim_ratio = per_copy_bitscore / self_bits_real[subfamily]`. A value of 1.0 means a perfect match to the consensus; lower values indicate divergence or truncation.
4. Outputs: `self_bits_real.tsv`, `sim_scores.tsv` (seqID, sim_bitscore, self_bits, sim_ratio).

#### Stage 1: Build ALL Table
Joins `extracted.fasta` headers with `assignment_full.tsv` and `sim_scores.tsv` to create the master 12-column table `all_sines.bedlike.ALL.tsv`:

| Column | Name | Description |
|--------|------|-------------|
| 1 | chr | Chromosome/contig name |
| 2 | start0 | 0-based start coordinate (BED format) |
| 3 | end | End coordinate |
| 4 | subfam_from_extracted | Original subfamily from extraction header |
| 5 | assigned_subfam | Subfamily from assignment (or extraction fallback) |
| 6 | strand | Strand (+/-/.) |
| 7 | best_bitscore | Best bitscore from assignment or extraction |
| 8 | leak_flag | `.` or `LEAK` |
| 9 | conflict_flag | `.` or `CONFLICT` |
| 10 | note | Semicolon-delimited metadata tags |
| 11 | sim_bitscore | Single ssearch36 bitscore vs assigned consensus (`.` if unassigned) |
| 12 | sim_ratio | sim_bitscore / self_bits_real (`.` if unassigned) |

#### Stage 2: LEAK Annotation
Scans column 10 for `runner_ratio=` tags. If the runner (second-best consensus) achieves ≥90% of the winner's bitscore, the sequence is flagged as `LEAK` in column 8. This indicates the assignment is close to ambiguous — the sequence may be a chimera or lie at a subfamily boundary.

#### Stage 3: CONFLICT Annotation
Uses `bedtools intersect` to check each SINE locus against `regions.by_subfam.bed`. If a locus falls within a genomic region where multiple subfamilies have overlapping search hits, it is flagged as `CONFLICT` in column 9. Additional tags are added:
- `conf_subfams=` — comma-separated list of overlapping subfamily names.
- `conf_alt=yes/no` — whether any conflicting subfamily differs from the assigned one.

#### Stage 3.5: Soft Assignment of Unassigned Sequences
For sequences that failed firm assignment (no unanimous vote, rejected low bitscore, or no ssearch data):
1. Intersects their genomic coordinates with `all_hits.labeled.bed`.
2. Identifies which original `sear` search query/queries produced the hit.
3. Picks the query with the highest search score as the "soft" subfamily assignment.
4. Cross-references with `all_votes.tsv` for ssearch36 data.
5. Records everything in `unassigned.tsv`.

#### Stage 4: Summary Statistics
Builds `summary.by_subfam.tsv` with per-subfamily counts:

| Column | Description |
|--------|-------------|
| subfam | Subfamily name |
| firm_assigned | Sequences with 10/10 unanimous assignment above threshold |
| soft_assigned | Unassigned sequences soft-assigned via search overlap |
| total_assigned | firm + soft |
| leak_n | Number of LEAK-flagged among firm assigned |
| conf_alt_n | Number with CONFLICT where conflicting subfamily ≠ assigned |
| firm_pct | Percentage of total extracted sequences |
| total_pct | Percentage of total extracted sequences |
| leak_pct | Percentage of firm assigned |
| conf_alt_pct | Percentage of firm assigned |
| sim_mean | Mean sim_ratio across firm assigned copies |
| sim_median | Median sim_ratio across firm assigned copies |

#### Stage 5: Sanity Checks
Logs counts of runner_ratio tags, LEAK rows, CONFLICT rows, conf_alt=yes rows, and sim_ratio scored counts for quick verification.

---

## 9. Step 4 — Diagnostic Plots (`step4_plots.sh`)

### Usage (standalone)
```bash
step4_plots.sh <RUN_ROOT> [THREADS] [MAX_SEQS]
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `RUN_ROOT` | *required* | Path to a SINEderella run directory |
| `THREADS` | auto (`nproc`) | CPU threads for ssearch36 and MAFFT |
| `MAX_SEQS` | 10000 | Maximum copies per subfamily to use for alignment (random subsample if exceeded) |

### What It Does

Step 4 generates two diagnostic plots per subfamily in both PNG and PDF formats:

#### Plot 1: Divergence Histogram

For each subfamily, computes `% divergence = 100 − %identity` between each assigned copy and its consensus via a single `ssearch36` alignment, then plots a histogram showing the distribution across all copies.

- X-axis: % divergence from consensus
- Y-axis: number of copies
- Mean and median divergence marked with dashed/dotted vertical lines
- Useful for assessing the age distribution and integrity of SINE copies within each subfamily

#### Plot 2: Nucleotide Composition (Per-Position Frequency)

Aligns all assigned copies to the consensus using MAFFT, then plots a stacked bar chart showing the frequency of each nucleotide (A, T, C, G, gap) at every position along the consensus.

- X-axis: consensus position (bp), labeled with the consensus nucleotide
- Y-axis: frequency (0–1.0)
- Colors: A=green, T=red, C=blue, G=orange, Gap=grey
- Useful for identifying conserved vs. variable positions, diagnostic SNPs between subfamilies, and degradation patterns

### Output

All plots are written to `<step2_output>/plots/`:

| File Pattern | Description |
|-------------|-------------|
| `<subfamily>_divergence.png` | Divergence histogram (PNG, 200 DPI) |
| `<subfamily>_divergence.pdf` | Divergence histogram (PDF, vector) |
| `<subfamily>_nucfreq.png` | Nucleotide frequency chart (PNG, 200 DPI) |
| `<subfamily>_nucfreq.pdf` | Nucleotide frequency chart (PDF, vector) |

### Dependencies

Step 4 requires `python3` with `matplotlib` and `numpy` installed. If these are not available, step 4 is skipped with a warning — it does not block the pipeline.

```bash
# Install if missing
pip3 install matplotlib numpy
# or via apt (Debian/Ubuntu):
sudo apt-get install python3-matplotlib python3-numpy
```

### Helper Script: `plot_subfamily.py`

The Python plotting logic is in `plot_subfamily.py`. It is called by `step4_plots.sh` for each subfamily and accepts:

```
python3 plot_subfamily.py \
  --subfamily <name> \
  --pctid <pctid.tsv> \
  --msa <msa.fasta> \
  --outdir <dir> \
  --consensus-name <name>
```

---

## 10. Step 5 — Subfamily Alignments (`step5_align_subfamilies.sh`)

### Usage (standalone)
```bash
step5_align_subfamilies.sh <consensi_bank.fa> <multi_run_dir>
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `consensi_bank.fa` | *required* | FASTA file of all consensus sequences (the common reference bank) |
| `multi_run_dir` | *required* | Path to a `SINEderella_multi` project directory (contains species subdirectories) |

### What It Does

Step 5 aligns subfamily copies from each species against a shared consensus bank. For every subfamily FASTA in each species' `subfamilies/` folder, the script applies a size-dependent strategy:

| Condition | Action |
|-----------|--------|
| **< 400 copies** | Sample 200 copies (or all if fewer), concatenate with consensus bank, align with `mafft --auto` |
| **≥ 400 copies** | Subsample 10,000 copies, run `SubFam` to generate chunk-consensus sequences (`.clw`), convert to aligned FASTA, then use `mafft --add` to insert the consensus bank into the existing SubFam alignment (avoids redundant full re-alignment) |

**Sampling** uses an `awk`+`shuf` approach (linearise FASTA records, shuffle, take the first N), not `seqkit sample`.

### Output

All alignment files are collected into a single aggregation directory:

```
<multi_run_dir>/alignments/
  <Species_1>/
    <subfamily_A>.aln.fa
    <subfamily_B>.aln.fa
    ...
  <Species_2>/
    ...
```

Each `.aln.fa` file is a MAFFT alignment in FASTA format containing the consensus bank sequences plus the subfamily copies (or SubFam chunk-consensuses for large subfamilies).

### Dependencies

| Tool | Purpose |
|------|---------|
| `mafft` | Multiple sequence alignment |
| `SubFam` | Sequence clustering for large subfamilies (must be in the same directory as the script) |
| `shuf` | Random sampling |

All intermediate files (SubFam work directories, temp FASTA samples) are cleaned up automatically.

---

## 11. SINEderella_multi — Multi-Genome Mode

`SINEderella_multi` is the multi-genome orchestration wrapper. It runs SINEderella on multiple genomes and produces a cross-species summary table plus an importable SINE-KB report.

### Usage

```bash
SINEderella_multi full    <genomes.tsv> <consensus.fa>  [options]
SINEderella_multi add     <genomes.tsv> <new_sine.fa>   [options]
SINEderella_multi exclude <genomes.tsv> <SINE_name>     [options]
SINEderella_multi plots   <genomes.tsv>                 [options]
SINEderella_multi summary                               --project <dir> [options]
```

### Subcommands

| Subcommand | Description |
|------------|-------------|
| `full` | Run the complete SINEderella pipeline on every genome in the list |
| `add` | Run `SINEderella --add` on every genome |
| `exclude` | Run `SINEderella --exclude` on every genome |
| `plots` | Run step 4 diagnostic plots on every genome |
| `summary` | (Re)generate the cross-species summary and SINE-KB report only |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--project <dir>` | auto-created | Project directory. Required for `summary` mode. Auto-created for `full`. |
| `--parallel <N>` | 1 | Run up to N genomes concurrently |
| `--threads <N>` | `nproc` | CPU threads per SINEderella invocation |

### genomes.tsv Format

Tab-separated file, one line per genome:

```
Species_Name    /path/to/genome.fa    [/consensus_override.fa]    [/workdir]
```

| Column | Required | Description |
|--------|----------|-------------|
| 1 | yes | Species name (used as directory name; must be unique) |
| 2 | yes | Path to genome FASTA |
| 3 | no | Per-species consensus override (uses shared consensus if omitted) |
| 4 | no | Working directory with existing `run_*` (for `add`/`exclude` without `--project`) |

Lines starting with `#` are comments; blank lines are ignored.

### Output

```
multi_YYYYMMDD_HHMMSS/
├── project_manifest.tsv               # Audit trail / provenance log
├── genomes.tsv                        # Copy of input genome list
├── <Species_Name>/                    # Per-species working directory
│   └── run_*/                         # SINEderella run directory
└── summary/
    ├── cross_species_summary.tsv      # Wide table (species × subfamily counts)
    ├── sinekb_report.json             # SINE-KB importable report
    ├── <Species>.summary.by_subfam.tsv      # (symlinks)
    ├── <Species>.assignment_full.tsv        # (symlinks)
    ├── <Species>.all_sines.bedlike.ALL.tsv  # (symlinks)
    └── <Species>.unassigned.tsv             # (symlinks)
```

### cross_species_summary.tsv

Wide-format table suitable for Google Sheets. One row per species:

| Column | Description |
|--------|-------------|
| Species | Species name |
| Total_Loci | Total extracted SINE loci |
| `<Subfamily>_firm` | Firmly assigned copies (10/10 unanimous + above threshold) |
| `<Subfamily>_total` | Firm + soft assigned copies |
| Unassigned_n | Number of unassigned loci |
| Unassigned_pct | Percentage unassigned |

### Example

```bash
# Run full pipeline on 4 genomes, 8 threads each, 2 in parallel
SINEderella_multi full genomes.tsv sine_consensuses.fa --parallel 2 --threads 8

# Re-generate summary for an existing project
SINEderella_multi summary --project multi_20260211_143000

# Add a new SINE family to all genomes
SINEderella_multi add genomes.tsv new_family.fa --project multi_20260211_143000
```

---

## 12. SINE-KB Report & Import

SINEderella_multi automatically generates a **`sinekb_report.json`** file that can be uploaded to the SINE Knowledge Base to create or update subfamily entries and per-genome taxon instances.

### Report Generation

The report is generated automatically at the end of every `SINEderella_multi` run (including `summary` mode). It can also be generated standalone:

```bash
python3 generate_sinekb_report.py <project_dir> [--out path/to/report.json]
```

Default output location: `<project_dir>/summary/sinekb_report.json`

### Report Contents

The JSON file contains:

| Field | Description |
|-------|-------------|
| `sinederella_report` | Always `true` — identifies the file as a valid report |
| `version` | Report format version |
| `generated` | Date of generation |
| `genomes[]` | Per-genome metadata: species name, genome path, consensus path |
| `subfamilies[]` | Per-subfamily data: name, consensus sequence, length, self_bitscore |
| `subfamilies[].per_genome{}` | Per-genome stats for each subfamily (keyed by species name) |

Per-genome stats included for each subfamily:

| Field | Description |
|-------|-------------|
| `firm_assigned` | Copies with 10/10 unanimous vote + above bitscore threshold |
| `soft_assigned` | Unassigned copies soft-assigned via search overlap |
| `total_assigned` | firm + soft |
| `leak_n` | Copies flagged as LEAK |
| `conf_alt_n` | Copies in CONFLICT regions (different subfamily) |
| `sim_mean` | Mean similarity ratio across firm copies |
| `sim_median` | Median similarity ratio across firm copies |

### Importing into SINE-KB

**Via the web interface:**

1. Open SINE-KB at `http://localhost:5555`
2. Go to **Tools** (top navigation)
3. Under **Import SINEderella Report**, click "Choose File" and select your `sinekb_report.json`
4. Click **Import Report**
5. A flash message will summarize what was created/updated

**Import behavior:**

- **Existing subfamilies** — matched by name. Their taxon instances are updated with new copy counts, sim_mean, and sim_median. User-edited fields (notes, NCBI taxid, aliases, module annotations, top100/rand100 identity) are preserved.
- **New subfamilies** — placed in an auto-created **"Unassigned"** family for manual triage. Move them to the correct family via the web UI.
- **Subfamilies with zero assignments** across all genomes are skipped.
- **Re-importing** the same report is safe — it updates rather than duplicates.

### Standalone Report Generation

If you have an existing `SINEderella_multi` project and want to (re)generate only the report:

```bash
# From an existing multi project
python3 generate_sinekb_report.py multi_20260211_143000/

# Custom output path
python3 generate_sinekb_report.py multi_20260211_143000/ --out ~/reports/eryx_report.json
```

The script requires only Python 3.8+ (no additional packages beyond the standard library).

---

## 13. GitHub Pages (SINEdb) Export & Re-Import Workflow

The SINE Knowledge Base can be exported to a static GitHub Pages site at **https://toki-bio.github.io/SINEdb/**. This allows public, browsable access to all SINE subfamily and instance data without requiring a running server.

### Architecture

The static site is a single-page JavaScript application in `docs/index.html` that reads from `docs/data.json`. It supports:

- **All tabs**: Subfamilies, Families, Instances, Taxonomy, Sequences, **Summary** (species × subfamily heatmap matrix)
- **Editing**: All data is editable via click/double-click; edits are stored in the browser's `localStorage`
- **Export/Import**: Users can download the edited `data.json` and re-import it
- **Summary table**: Full-width species × subfamily copy count matrix with heatmap coloring, sortable columns, filtering by SINE family, and per-cell editing

### Export to GitHub Pages

```bash
# Generate docs/data.json from the sine-kb YAML data
python3 sine-kb/export_to_github_pages.py

# Verify
python3 -c "import json; d=json.load(open('docs/data.json')); print(f'{len(d[\"families\"])} families exported')"

# Commit and push to deploy
git add docs/data.json docs/index.html
git commit -m "Update SINEdb data"
git push
```

GitHub Pages automatically serves the `docs/` folder. Changes appear at https://toki-bio.github.io/SINEdb/ within a few minutes of pushing.

### Complete Re-Import Workflow (from a SINEderella_multi run)

When a multi-genome run completes (or more genomes finish), use this workflow to import the results into the SINE-KB and push to the public site.

#### Step 1: Download results from server

```bash
# SSH to server and create tarball of completed results
ssh kit "cd /path/to/multi_run && tar -czh -f /tmp/sine_results.tar.gz \
  project_manifest.tsv genomes.tsv \
  */manifest.txt */results/summary.by_subfam.tsv */results/assignment_stats.tsv \
  */consensuses.clean.fa"

# SCP to local workspace
scp kit:/tmp/sine_results.tar.gz .
mkdir -p server_import && cd server_import
tar xzf ../sine_results.tar.gz
ssh kit "rm /tmp/sine_results.tar.gz"
```

> **Note:** The `-h` flag in `tar` follows symlinks, which is necessary because SINEderella stores most result files as symlinks in the `results/` directory.

#### Step 2: Generate SINE-KB report

```bash
# Set up the expected directory structure for the report generator
# (it looks for step2/step2_output/ under each species dir)
cd server_import
for sp in */; do
  sp="${sp%/}"
  [ -d "$sp/results" ] || continue
  mkdir -p "$sp/step2/step2_output"
  for f in "$sp/results/"*; do
    [ -f "$f" ] && ln -sf "../../results/$(basename $f)" "$sp/step2/step2_output/$(basename $f)"
  done
done

# Generate report
python3 ../generate_sinekb_report.py . --out sinekb_report.json
```

#### Step 3: Import into SINE-KB

**Option A: Via the web interface**

1. Start the SINE-KB server: `cd sine-kb && python3 server.py`
2. Open http://localhost:5555/tools
3. Upload `sinekb_report.json` under "Import SINEderella Report"
4. Review the results and manually move any new subfamilies from the "Unassigned" family to the correct family

**Option B: Via a custom import script**

For large imports with known subfamily-to-family mappings, write a Python import script (see `import_snake_run.py` as a template):

```python
#!/usr/bin/env python3
"""Import SINEderella multi-genome results into SINE-KB."""
import json, sys
sys.path.insert(0, 'sine-kb')
import models

models.init()

# Load report
with open('server_import/sinekb_report.json') as f:
    report = json.load(f)

# Define subfamily → family mapping
FAMILY_MAP = {
    'Sq3A_Eryx': 'Squam3',
    'sq1Pb1':    'Squam1',
    # ... add all subfamilies
}

# Import
result = models.import_sinekb_report(report, family_map=FAMILY_MAP)
print(f"Created: {result['subfamilies_created']} subfamilies, {result['instances_created']} instances")
print(f"Updated: {result['instances_updated']} instances")
```

#### Step 4: Update taxonomy (if needed)

If new species were added, verify the taxonomy tree in the SINE-KB web interface at http://localhost:5555/taxa. The import script adds species at the top level by default — rearrange them into the correct phylogenetic position.

#### Step 5: Export to GitHub Pages and deploy

```bash
# Export updated data
python3 sine-kb/export_to_github_pages.py

# Commit and push
git add docs/data.json
git commit -m "Import results from multi-genome run YYYY-MM-DD"
git push
```

### Summary of commands (quick reference)

```bash
# 1. Get data from server
scp kit:/path/to/multi_run/summary/sinekb_report.json server_import/

# 2. Import into SINE-KB (via web UI or script)
cd sine-kb && python3 server.py &
# → Upload at http://localhost:5555/tools

# 3. Export & deploy
python3 sine-kb/export_to_github_pages.py
git add docs/data.json && git commit -m "Update SINEdb" && git push
```

### Re-importing when more genomes complete

When additional genomes finish in an ongoing `SINEderella_multi` run:

1. **Re-generate the report** on the server (or locally after downloading):
   ```bash
   ssh kit "cd /path/to/multi_run && python3 generate_sinekb_report.py . --out summary/sinekb_report.json"
   ```
2. **Re-import** — the import is idempotent: existing instances are updated, new ones are created. No data is duplicated.
3. **Export and push** as above.

The report generator automatically picks up all species directories with completed `summary.by_subfam.tsv` files, so re-running it after more genomes finish will include the new species.

---

## 14. Output Files Reference

After a complete run, the `results/` directory contains symlinks to all key outputs:

### Primary Output Files

| File | Description |
|------|-------------|
| `summary.by_subfam.tsv` | Per-subfamily statistics (firm/soft assigned, LEAK, CONFLICT, sim_mean, sim_median) |
| `all_sines.bedlike.ALL.tsv` | Master annotation table (12 columns, BED-like, including sim_bitscore and sim_ratio) |
| `assignment_full.tsv` | Detailed per-sequence assignment (Sequence, Subfamily, Bitscore, Votes, Status, Threshold) |
| `assignment_stats.tsv` | Per-subfamily assignment statistics (Subfamily, Assigned, TopN_Bitscore, Threshold) |
| `unassigned.tsv` | Soft-assignment of sequences that failed firm assignment (12 columns) |
| `all_votes.tsv` | Per-sequence best vote from ssearch36 (all sequences, not just assigned) |
| `self_bits.tsv` | Per-subfamily maximum bitscore threshold (from assignment) |
| `self_bits_real.tsv` | Per-subfamily real self-alignment bitscore (from ssearch36 self-alignment) |
| `sim_scores.tsv` | Per-copy similarity scores (seqID, sim_bitscore, self_bits, sim_ratio) |

### Sequence Files

| File | Description |
|------|-------------|
| `assigned.fasta` | Assigned sequences (headers: `>seqID\|subfamily\|bitscore`) |
| `unassigned.fasta` | Sequences that failed assignment criteria |
| `extracted.fasta` | All candidate SINE sequences from step 1 |
| `subfamilies/` | Directory with per-subfamily FASTA files |
| `consensuses.fa` | Consensus sequences used in this run |

### Diagnostic Plots

| File | Description |
|------|-------------|
| `plots/<subfamily>_divergence.png` | Histogram of % divergence from consensus (PNG) |
| `plots/<subfamily>_divergence.pdf` | Histogram of % divergence from consensus (PDF) |
| `plots/<subfamily>_nucfreq.png` | Per-position nucleotide frequency chart (PNG) |
| `plots/<subfamily>_nucfreq.pdf` | Per-position nucleotide frequency chart (PDF) |

### Genomic Data

| File | Description |
|------|-------------|
| `all_hits.labeled.bed` | All search hits with subfamily labels (6-col BED) |
| `regions.by_subfam.bed` | Merged regions with comma-separated subfamily lists |

### Multi-Species Alignment Files (Step 5)

Generated by `step5_align_subfamilies.sh` in the multi-run project directory:

| File | Description |
|------|-------------|
| `alignments/<Species>/<subfamily>.aln.fa` | MAFFT alignment of consensus bank + subfamily copies (or SubFam chunk-consensuses) |

### SINE-KB Report (Multi-Genome Mode)

| File | Description |
|------|-------------|
| `sinekb_report.json` | JSON report importable into SINE-KB (generated by `generate_sinekb_report.py`) |

### Run Metadata

| File | Description |
|------|-------------|
| `manifest.txt` | Run parameters, paths, mode, date |
| `summary.txt` | Human-readable assignment summary report |

---

## 15. Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `THREADS` | `nproc` | Number of CPU threads for ssearch36 and MAFFT |
| `CHUNK_BP` | 30000 | Step 1 sample size (max sequences for SubFam clustering) |
| `FLANK` | 50 | Step 1 flank size parameter |

Example:
```bash
THREADS=32 CHUNK_BP=50000 SINEderella genome.fa consensuses.fa
```

---

## 16. Algorithm Details

### 16.1 Sear Search

`sear` performs a fast homology search of each SINE consensus against the genome. The parameters used:

```
sear -k <query> <genome> 0.8 65 50
```

- `-k` (keep): Retains genome split parts between queries, avoiding redundant genome fragmentation.
- `0.8` (SFL — Sequence Fraction Length): A match must cover at least 80% of the query consensus length.
- `65`: Minimum similarity percentage between query and hit.
- `50`: Additional parameter for match evaluation.

### 16.2 10-Cycle Voting

The core innovation of SINEderella's assignment is the **10-cycle unanimous voting** system:

```
For each of 10 independent cycles:
    Run ssearch36 (consensus vs. all SINE sequences)
    For each SINE sequence, record the best-matching consensus

A sequence is "unanimous" if and only if:
    The same consensus wins in all 10 out of 10 cycles
```

**Why 10 cycles?** `ssearch36` uses stochastic Z-score normalization (shuffling the database). This means bitscores fluctuate by ±1–3% between runs. By requiring 10/10 agreement, SINEderella ensures only robustly assignable sequences pass — eliminating borderline cases where stochastic variation could flip the winner.

**ssearch36 parameters used:**
```
ssearch36 -g -3 -T <threads> -Q -n -z 11 -E 2 -w 95 -W 70 -m 8
```

| Flag | Meaning |
|------|---------|
| `-g` | Global/gapped alignment |
| `-3` | Forward strand + reverse complement |
| `-T` | Thread count |
| `-Q` | Quiet mode |
| `-n` | Nucleotide comparison |
| `-z 11` | Z-score normalization method 11 (shuffle with composition) |
| `-E 2` | E-value threshold |
| `-w 95` | Alignment width |
| `-W 70` | Window size |
| `-m 8` | BLAST-tabular output format |

### 16.3 Bitscore Thresholds

After identifying unanimous sequences, per-subfamily thresholds prevent false assignments of divergent copies:

```
For each consensus subfamily C:
    Sort C's unanimous sequences by bitscore (descending)
    N = min(10, number_of_unanimous_for_C)
    TopN = bitscore of the N-th best sequence
    Threshold = floor(0.45 × TopN + 0.5)

A unanimous sequence is "assigned" if:
    its bitscore ≥ Threshold
Otherwise:
    status = "rejected_low_bitscore"
```

The `MIN_REL_BITSCORE = 0.45` factor was chosen to reject severely truncated or degraded copies while retaining most genuine subfamily members. The TopN reference uses up to the 10th-best sequence to be robust against outliers.

### 16.4 Incremental BED-Overlap Strategy

In `--add` mode, the pipeline uses a BED-overlap strategy to determine which existing assignments can be safely carried forward:

```
Tier 1 Criteria (carry forward — skip ssearch36):
    ✓ status=assigned in old run
    ✓ 10/10 unanimous votes in old run
    ✓ Still present in current extracted sequences
    ✗ Does NOT overlap with any new-consensus BED hit region

Tier 2 Criteria (re-evaluate — full 10-cycle ssearch36):
    - New sequences (not in old run)
    - Previously unassigned / rejected sequences
    - Sequences whose locus overlaps new-consensus BED regions
```

**Rationale:** Empirical testing showed that when a new consensus is added, no existing firm assignment (10/10 unanimous) ever flips subfamily. Bitscores drift ±1–26 stochastically but the 10/10 unanimity requirement prevents reassignment. The only cases where reassignment might matter are sequences whose genomic regions overlap new BED hits — these are conservatively re-evaluated.

**After merging tiers**, bitscore thresholds are recalculated from the combined unanimous set. This ensures thresholds reflect the full consensus landscape including the new subfamily.

### 16.5 LEAK Detection

A `LEAK` flag is applied when the **runner-up** consensus (second-best match) achieves ≥90% of the winner's bitscore:

```
If runner_ratio ≥ 0.90:
    mark as LEAK
```

LEAK-flagged sequences may represent:
- Chimeric insertions spanning two subfamily boundaries.
- Intermediate forms between subfamilies.
- Regions where subfamily differentiation is minimal.

LEAK sequences remain assigned but are flagged for review.

### 16.6 CONFLICT Annotation

A `CONFLICT` flag is applied when a SINE locus falls in a genomic region where multiple subfamilies had overlapping `sear` search hits:

```
If bedtools intersect shows the locus overlaps a region with comma-separated subfamily names:
    mark as CONFLICT
    If any conflicting subfamily ≠ assigned subfamily:
        conf_alt = yes
    Else:
        conf_alt = no
```

CONFLICT regions represent genomic hotspots where multiple SINE subfamilies co-occur — potentially indicating nested insertions, fragmented copies, or subfamily-dense regions.

### 16.7 Soft Assignment of Unassigned Sequences

Sequences that fail firm assignment are given a **soft assignment** based on which `sear` search query originally identified them:

1. The sequence's genomic coordinates are intersected with `all_hits.labeled.bed`.
2. All overlapping search hits are collected.
3. The hit with the highest search score determines the soft subfamily.
4. Cross-referenced with `all_votes.tsv` for ssearch36 voting data.

Soft-assigned sequences are included in `total_assigned` counts in the summary but NOT in `firm_assigned`. They are tracked separately in `unassigned.tsv` with columns documenting the reason for failure and all available evidence.

Possible failure reasons:
- `no_unanimous_votes` — No consensus won all 10/10 cycles.
- `rejected_low_bitscore` — Unanimous but below bitscore threshold.
- `no_ssearch_data` — No ssearch36 hit data (very short or degenerate sequence).

### 16.8 Similarity Ratio (sim_ratio)

Each firmly assigned copy receives a **normalized similarity score** computed in step 3:

```
1. Self-alignment: ssearch36 consensus.fa consensus.fa
   → per-subfamily self_bits = max bitscore where query == subject

2. Per-copy alignment: ssearch36 consensus_single.fa copies.fa
   → per-copy sim_bitscore = best bitscore vs its assigned consensus

3. sim_ratio = sim_bitscore / self_bits[assigned_subfamily]
```

**Interpretation:**
- `sim_ratio ≈ 1.0` — Copy is nearly identical to the consensus (young/intact insertion).
- `sim_ratio ≈ 0.5–0.7` — Typical for moderately diverged SINE copies.
- `sim_ratio < 0.3` — Highly degraded, truncated, or divergent copy.

**Why normalize?** Raw bitscores are not comparable across subfamilies of different lengths. Dividing by the self-alignment bitscore normalizes for consensus length and composition, making `sim_ratio` directly comparable across all subfamilies.

The summary table includes `sim_mean` and `sim_median` per subfamily. These provide a quick proxy for the "age" or conservation level of each subfamily within the genome.

---

## 17. Directory Structure

A complete full-mode run produces this directory tree:

```
run_YYYYMMDD_HHMMSS/
├── manifest.txt                         # Run metadata
├── genome.clean.fa                      # Sanitized genome
├── genome.clean.fa.fai                  # Genome index
├── consensuses.clean.fa                 # Sanitized consensuses
├── step1_search_extract.sh              # Frozen copy of step1
├── step2_asSINEment.sh                  # Frozen copy of step2
├── step3_postprocess.sh                 # Frozen copy of step3
├── step4_plots.sh                       # Frozen copy of step4
├── plot_subfamily.py                    # Frozen copy of plot helper
├── step1.stdout.log                     # Step1 stdout
├── step1.stderr.log                     # Step1 stderr
│
├── genome.clean_step1/                  # Step 1 outputs
│   ├── extracted.fasta                  # All extracted SINE sequences
│   ├── merged_hits.bed                  # Merged genomic intervals
│   ├── sampled_30000.fasta              # Sampled subset
│   ├── genome.path                      # Genome path reference
│   ├── consensus.path                   # Consensus path reference
│   ├── run_info.txt                     # Step1 run metadata
│   ├── searches/                        # Per-query search results
│   │   ├── gen-<query>.bed              # Per-query BED hits
│   │   ├── <query>.part_<name>.<ext>    # Split consensus query files
│   │   ├── all_hits.labeled.bed         # Combined labeled hits
│   │   └── regions.by_subfam.bed        # Conflict regions
│   └── subfam_input/                    # Clustering outputs
│       ├── input.fasta
│       ├── input.clw                    # SubFam output
│       ├── input_reps.fasta             # Cluster representatives
│       ├── combined_input.fasta         # Reps + consensuses
│       └── input.clw.al                 # MAFFT alignment
│
├── step2/                               # Step 2 outputs
│   └── step2_output/
│       ├── assigned.fasta               # Assigned sequences
│       ├── unassigned.fasta             # Unassigned sequences
│       ├── assignment_full.tsv          # Full assignment details
│       ├── assignment_stats.tsv         # Per-subfamily statistics
│       ├── all_votes.tsv                # All vote data
│       ├── all_sines.bedlike.ALL.tsv    # Master annotation table (12 cols)
│       ├── summary.by_subfam.tsv        # Summary statistics
│       ├── unassigned.tsv               # Soft assignment details
│       ├── self_bits.tsv                # Per-subfamily max threshold
│       ├── self_bits_real.tsv           # Per-subfamily self-alignment bitscore
│       ├── sim_scores.tsv              # Per-copy similarity scores
│       ├── summary.txt                  # Human-readable summary
│       ├── subfamilies/                 # Per-subfamily FASTA
│       │   ├── Sq3A.fasta
│       │   ├── Sq3C.fasta
│       │   └── ...
│       └── plots/                       # Diagnostic plots (step 4)
│           ├── Sq3A_divergence.png
│           ├── Sq3A_divergence.pdf
│           ├── Sq3A_nucfreq.png
│           ├── Sq3A_nucfreq.pdf
│           └── ...
│
└── results/                             # Symlinks to key outputs
    ├── summary.by_subfam.tsv
    ├── all_sines.bedlike.ALL.tsv
    ├── assignment_full.tsv
    ├── assignment_stats.tsv
    ├── unassigned.tsv
    ├── all_votes.tsv
    ├── self_bits.tsv
    ├── self_bits_real.tsv
    ├── sim_scores.tsv
    ├── assigned.fasta
    ├── unassigned.fasta
    ├── extracted.fasta
    ├── subfamilies/
    ├── plots/                           # Diagnostic plots (step 4)
    ├── all_hits.labeled.bed
    ├── regions.by_subfam.bed
    └── consensuses.fa
```

---

## 18. Troubleshooting

### Common Issues

**"Missing required tool: sear"**
Ensure all dependencies are installed and available on `$PATH`. Check with `which sear`.

**"No BED files generated by sear"**
Your consensus sequences may not have significant hits in the genome at the current sensitivity thresholds (SFL=0.8, similarity=65%). Check if the consensus is appropriate for this genome.

**"SubFam failed to produce input.clw"**
SubFam may fail on very small input or in certain environments. This is a step 1 clustering issue and does not affect the core assignment pipeline (steps 2–3 can be run independently if `extracted.fasta` exists).

**Step 2 produces 0 assigned sequences**
- Check that consensus sequences are in correct FASTA format.
- Verify that `extracted.fasta` contains sequences (not empty).
- Examine `ssearch.log` for alignment errors.
- Try reducing `MIN_REL_BITSCORE` if thresholds are too strict.

**"No existing run_* directory found"**
When using `--add` or `--exclude` without `--run`, SINEderella looks for the most recent `run_*` directory in the current working directory. Either `cd` to the correct directory or use `--run <path>`.

### Logs

| Log File | Content |
|----------|---------|
| `step1.stdout.log` | Step 1 standard output |
| `step1.stderr.log` | Step 1 error output (most informative) |
| `step2/step2_output/ssearch.log` | ssearch36 stderr from all cycles |
| `step2/step2_output/seqkit.log` | seqkit messages |
| `step2/step2_output/summary.txt` | Human-readable assignment summary |

---

## 19. Examples

### Example 1: Full analysis of a genome

```bash
# Run full pipeline with 16 threads
THREADS=16 ./SINEderella genome_assembly.fa sine_consensuses.fa

# Check results
cat run_*/results/summary.by_subfam.tsv
```

### Example 2: Add a newly discovered SINE family

```bash
cd /path/to/project
# Assuming run_20260207_033350/ exists from a previous full run

# Add new consensus (auto-detects latest run)
./SINEderella --add newly_discovered_SINE.fa

# Or specify the source run explicitly
./SINEderella --add newly_discovered_SINE.fa --run run_20260207_033350

# Compare results
diff <(cut -f1,2 run_20260207_033350/results/summary.by_subfam.tsv) \
     <(cut -f1,2 run_add_*/results/summary.by_subfam.tsv)
```

### Example 3: Remove a suspected false-positive family

```bash
# Exclude family "FalsePos1" from the latest run
./SINEderella --exclude FalsePos1

# Verify it's gone
grep FalsePos1 run_excl_*/results/summary.by_subfam.tsv
# (should produce no output)
```

### Example 4: Iterative refinement workflow

```bash
# 1. Initial full run
./SINEderella genome.fa initial_consensuses.fa

# 2. Add a candidate family
./SINEderella --add candidate_A.fa --run run_20260207_100000

# 3. Review results — candidate_A has too few assignments, remove it
./SINEderella --exclude CandidateA --run run_add_20260207_100500

# 4. Add a different candidate
./SINEderella --add candidate_B.fa --run run_excl_20260207_101000
```

### Example 5: Multi-genome run with SINE-KB import

```bash
# 1. Prepare genomes.tsv (tab-separated)
cat > genomes.tsv <<'EOF'
Eryx_jayakari	/data/genomes/eryx_jayakari.fa
Boa_constrictor	/data/genomes/boa_constrictor.fa
Python_bivittatus	/data/genomes/python_bivittatus.fa
EOF

# 2. Run SINEderella_multi on all genomes (4 parallel, 8 threads each)
./SINEderella_multi full genomes.tsv sine_consensuses.fa --parallel 4 --threads 8

# 3. The report is generated automatically:
#    multi_*/summary/sinekb_report.json
#    Upload it to SINE-KB via Tools → Import SINEderella Report

# 4. Or re-generate the report later:
python3 generate_sinekb_report.py multi_20260211_143000/
```

---

## License

See repository for license information.

## Citation

If you use SINEderella in your research, please cite the repository:
`https://github.com/Toki-bio/SINEderella-dev`
