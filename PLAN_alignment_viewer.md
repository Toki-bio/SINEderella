# Alignment & Viewer Integration Plan

Reference document for extraction, alignment, and MSA viewer integration into SINEdb.
**Created**: 2026-02-15
**Status**: In planning

---

## 1. Data Sources (per species × subfamily)

From SINEderella output (per run directory):
- **`.bnk` file**: FASTA of extracted SINE copies with flanks (default 50bp each side, configurable via `--slop`)
- **`.bed` file**: Genomic coordinates of SINE copies (6-col: chr, start, end, bitscore, homology, strand, ...)
- **`.best50.mafft.fa`**: Top-50-by-bitscore copies aligned with consensus (already produced by `sear --best`)
- **Genome FASTA** (`.fa` + `.fai`): For re-extraction with custom flanks

## 2. Alignment Tiers (per species × subfamily)

For each subfamily hit in a species, produce multiple alignment products:

### Tier A — Random 50
- **Up to 50 randomly sampled copies** (core SINE only, no flanks)
- Prepend consensus sequence as first entry
- Align with MAFFT (full params, see §4)
- Purpose: Unbiased view of subfamily diversity

### Tier B — Best 50
- **Up to 50 copies with highest bitscore to consensus**
- Already partially implemented by `sear --best 50`
- Prepend consensus, align with MAFFT
- Purpose: See best-preserved copies, estimate master sequence

### Tier C — SubFam consensuses (>400 copies only)
- If a subfamily has **>400 copies** in a species:
  - Subsample **2500 copies** randomly (50 chunks × 50 seqs/chunk = ~50 consensuses)
  - Run `SubFam` with default BnkSz=50
  - SubFam workflow: `mafft --retree 0 --reorder` → `seqkit split2 -s 50` → align each chunk → `cons -plurality 18` per chunk → final alignment of chunk consensuses → MSF
  - Output: FASTA of ~50 sub-subfamily consensuses
  - The subsample size (2500) is the variable, not BnkSz
- Purpose: Reveal internal subfamily structure / age layers

### Tier D — Zero/low-hit evidence
- For subfamilies with **0 hits** (or suspected false positives with very few hits):
  - Run raw ssearch36 against the genome with the subfamily consensus
  - Extract the **best insignificant bitscore** result
  - Store as proof of "not found" status (screenshot-equivalent)
- For **low copy number** subfamilies (1-few hits):
  - Flag for manual inspection to determine if hits are real or false positive
  - User will investigate manually
- Purpose: Document absence / justify "not found" status

### Minimum copies
- **Minimum = 1** (even a single copy could be the master copy and warrants investigation)

## 3. Flanking Sequences

For Tier A and Tier B alignments, each SINE copy is re-extracted from the genome with:
- **30 bp left flank**
- **70 bp right flank**

Flanks are included directly in the alignment (not as a separate file) to:
1. Prove both true ends of the SINE element exist
2. Inspect direct repeats at boundaries
3. Identify possible TSD (Target Site Duplication)

Implementation: Parse coordinates from subfamilies/*.fasta headers → BED →
`bedtools slop -l 30 -r 70 -s` → `bedtools getfasta` → prepend consensus → MAFFT.

## 4. MAFFT Parameters

All alignments use the same high-quality MAFFT invocation:

```bash
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) \
      --localpair --maxiterate 1000 --ep 0.123 \
      --nuc --reorder --preservecase --quiet \
      "$input" > "$output"
```

Notes:
- `--localpair --maxiterate 1000`: L-INS-i algorithm (most accurate for short sequences)
- `--ep 0.123`: Gap extension penalty (tuned for SINEs)
- `--preservecase`: Keep original case
- `--reorder`: Order by similarity in output

## 5. Output File Structure

```
docs/alignments/
  <Species>/
    <subfamily>.core.fa          # Tier A: random 50 with 30L+70R flanks, aligned
    <subfamily>.best50.fa        # Tier B: best 50 by bitscore with flanks, aligned
    <subfamily>.subfam.fa        # Tier C: SubFam consensuses (if >400 copies)
    <subfamily>.evidence.txt     # Tier D: ssearch36 best-hit for 0-hit subfams
```

Copies for Tier A and B are taken directly from step2 subfamilies/*.fasta
(no genome re-extraction for the copies themselves). Flanked versions are
re-extracted from the genome using coordinates parsed from FASTA headers.

## 6. Viewer Integration

### Phase 1: URL parameter for MSA viewer (start here)
- Add URL parameter support to `https://toki-bio.github.io/MSA-viewer/`
- Parameters: `?url=<path_to_fasta>` or `?data=<base64_encoded_fasta>`
- SINEdb heatmap cells get "View alignment" link → opens MSA viewer in new tab
- The MSA viewer fetches the alignment file from GitHub Pages

### Phase 2 (future): Embedded minimal viewer
- Extract core rendering functions from MSA viewer (~400 lines)
- Embed read-only display in SINEdb (panel or modal)
- Features: conservation shading, zoom, block/full mode, consensus line

## 7. Extraction Script Design

Script: `extract_alignments.sh`

```
Usage: extract_alignments.sh <consensus_bank.fa> <multi_run_dir> [output_dir]
```

For each species x subfamily:
1. Check copy count from step2 subfamilies/*.fasta
2. If 0 copies -> Tier D (run ssearch36, save evidence)
3. If 1-399 copies -> Tier A (random 50 + flanks) + Tier B (best 50 + flanks)
4. If >=400 copies -> Tier A + B + Tier C (SubFam)
5. For Tiers A & B, re-extract from genome with 30L+70R flanks
6. Align all with MAFFT (section 4 params)

## 8. Data Model Changes

In `data.json`, add to each instance:
```json
{
  "alignments": {
    "core": "alignments/Boa_constrictor/Squam1_1a.core.fa",
    "best50": "alignments/Boa_constrictor/Squam1_1a.best50.fa",
    "subfam": "alignments/Boa_constrictor/Squam1_1a.subfam.fa",
    "evidence": "alignments/Boa_constrictor/Squam1_1a.evidence.txt"
  }
}
```

Only populated fields are included (e.g., `subfam` only if >400 copies).

## 9. Open Questions

- [x] SubFam subsample size for >400 copies → **2500** (50 × BnkSz=50 = 50 consensuses)
- [x] Should flanked sequences be aligned or kept unaligned for manual inspection?
  -> Flanks are included in Tier A and B alignments directly (no separate file)
- [ ] How to handle species where genome is no longer accessible for re-extraction?
  - Fallback: use core sequences without flanks
- [ ] Tier D: exact ssearch36 parameters for "proof of absence" search?
- [ ] URL-parameter loading in MSA viewer: CORS considerations for GitHub Pages?
  - Same origin (toki-bio.github.io) should work fine for fetch()
- [ ] Maximum alignment file size for GitHub Pages? (100MB repo limit)

## 10. Dependencies

- `mafft` (with threading support)
- `SubFam` script (in SINEderella-dev)
- `bedtools` (slop, getfasta, merge)
- `samtools` (faidx)
- `seqkit` (split2, head, sample)
- `ssearch36` (for Tier D evidence)
- `cons` (EMBOSS, for SubFam consensus)
- `shuf` (GNU coreutils)

## 11. Implementation Order

1. [x] Create this plan document
2. [x] Add URL-parameter loading to MSA viewer (script.js)
3. [ ] Test URL loading works cross-origin on GitHub Pages
4. [x] Build `extract_alignments.sh` (Tiers A + B first)
5. [x] Add Tier C (SubFam integration)
6. [x] Add Tier D (zero-hit evidence)
7. [x] Update `export_to_github_pages.py` to include alignment paths in data.json
8. [x] Update SINEdb frontend: "View alignment" buttons in heatmap cells
9. [ ] Generate alignments for all 32 snake species
10. [ ] Push alignments to GitHub Pages, verify viewer links work
