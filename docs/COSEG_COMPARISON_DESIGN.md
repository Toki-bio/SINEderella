# SINEderella vs COSEG: test design (draft, 2026-10-08)

Status: **design only, nothing run. Metrics and decision rules below are fixed before any result is looked at; change them only by a dated edit to this file.**
Runs on KIT (tools: mafft, seqkit, cons, seqret, perl 5.22, hg38 at `/usr/local/genomes/hg38.mfa`; workdir `/data/V/toki/`, not `/data/W`). COSEG is not installed on KIT yet.

## 1. Question

Given the same set of copies, does SINEderella's route (SubFam chunks, then the peel, then asSINEment assignment of every copy) recover the subfamilies an expert recognises as well as, better than, or worse than COSEG (Price, Eskin, Pevzner 2004; maintained by Hubley/Smit/Siegel)? Both are tested as instruments handing an expert recalculated data, not as oracles.

What is **not** a question here: whether either is "better" in general. Only the families below.

## 2. The two sides are not the same kind of object

| | COSEG | SINEderella route |
|---|---|---|
| input | copies aligned to **one** reference (collinear), truncated copies hurt | copies, no reference; chunks of 50 by k-mer order |
| unit that is grouped | copy | chunk consensus (peel), then copy (step 2) |
| output | partition of copies + subfamily consensuses | peel groups -> consensuses -> flat 10-cycle assignment, firm / soft / unassigned |
| human step | review of the result | the subfamily call (MANUAL 6.1) |

A fair comparison therefore needs **copy-level partitions from both**: for SINEderella, peel groups -> consensus per group -> `step2_asSINEment.sh` on all copies. Arms:

* A. COSEG, default settings (`-k -d`, minimum sizes as in its README), reference = the family consensus, alignment built as in `benchmark/alu_konkel/to_coseg.py`.
* B. COSEG, parameters swept (min subfamily size, -m) and reported as a sweep, never one setting (SUBFAMILY_METHOD 7).
* C. SubFam chunks only (control: shows what chunking alone gives).
* D. SubFam + peel + asSINEment (the SINEderella route). Copies assigned soft or unassigned are kept and counted, not dropped.
* E. (optional) D with the old MAFFT ordering (`SUBFAM_ORDER=mafft`), to separate the effect of the new k-mer ordering.

## 3. Data (in order of independence of the truth)

1. **Konkel 2015, 343 loci KT305395-KT305737**: 295 usable copies; already run for A and C once (benchmark/alu_konkel). Labels are similarity-derived (best Price consensus) and therefore weak truth. Used as a smoke test and to find out why COSEG missed Ya5.
2. **hg38 chr21 Alu (10,521 copies, hs21 of FAMILY_SUBFAMILY_ASSIGNMENT.md)**: truth candidates = RepeatMasker subfamily calls (similarity-derived) and, for a sample, the expert's manual calls (section 5). Needs the copy set and its coordinates; the hs21 set exists on therioserver `~/tmp/votetest2/` (not checked from this session).
3. **Timema/Tal families where the owner has curated calls (`tim/`)**: no COSEG reference-collinearity assumption holds there; used only if COSEG can be run on a family consensus for it. Decide later.
4. **Simulation with known tree**: copies evolved along a chosen subfamily tree with gene conversion and indels, truncation added. The only data with true labels. Parameters to vary: subfamily size ratio, diagnostic positions per split (1, 2, 5), indel vs substitution diagnostics, fraction truncated.

## 4. Metrics (fixed now)

Primary, on copies present in both partitions, against the truth of each data set:

* adjusted Rand index and V-measure (homogeneity, completeness reported separately);
* number of groups at the truth's granularity, and the number of truth subfamilies with no group holding >=70 % of their copies (recovered / not recovered);
* the **50-best-copies separation test** (owner's criterion for near-identical subfamilies, `sine_his_separation_test_50best`): for each pair of called groups, do the 50 best copies of each separate by the diagnostic columns?
* stability: 10 random subsamples (50 % of copies); mean ARI between subsample partitions (Carey 2020 reports >10 % of replicate copies change RepeatMasker subfamily).
* coverage: fraction of copies unassigned / soft; for COSEG, fraction dropped.

Not used as evidence: raw bitscores, average pairwise identity, chunk-level separability (SUBFAMILY_METHOD sections 4c and 7).

Decision rule written before the run: a method is reported as better on a data set only if its ARI exceeds the other's by more than the spread over the 10 subsamples; otherwise "no difference detected". Wins on the simulation do not transfer to real families without the real-data result pointing the same way.

## 5. Manual review (the calibration loop)

For a random sample of **disagreements** (copies or groups where A and D place things differently) the owner sees the alignment in MSA-viewer without method labels and says what he sees. His calls are recorded verbatim and every variable is computed afterwards against them; nothing is predicted beforehand. Sample size: 30 disagreements per data set, drawn by a fixed seed.

## 6. Known risks

* Truth derived from similarity to consensuses favours whichever method also uses similarity to a consensus.
* COSEG may be run in a way its authors would not recommend on truncated copies; run it on a full-length-only subset as well as on all copies, and report both.
* The peel was calibrated on Timema (two cases), so Alu is out-of-sample for it; no parameter changes after seeing the Alu result.
* SubFam 1.2.0 k-mer ordering differs from the ordering the peel was calibrated on; arm E exists for this.
* Chunk names collide across runs (SUBFAMILY_METHOD 7a): compare by member sequence only.

## 7. Open decisions (owner)

1. Ground truth for data set 2: RepeatMasker calls only, or also a hand-called sample, and which sample.
2. Whether simulation (data set 4) is worth building before the real-data runs.
3. Whether the optional arm E is wanted.
4. Where COSEG gets installed on KIT (proposal: `/data/V/toki/coseg/`, from github.com/rmhubley/coseg, built by `make`).
