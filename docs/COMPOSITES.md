# Composite elements — design for SINEderella + SINE-discriminator (proposal, 2026-09-29)

Worked out by hand on rsi (Tal `rsi/REFINEMENT.md` §2–13); this generalizes it. Status: the two
tools exist (`tools/composite_scan.py`, `tools/build_composite.py`); **nothing below is wired into
the pipeline yet**.

## What a composite is, and what it is not

| case | what the copies show | example (rsi) | meaning |
|---|---|---|---|
| **structural dimer** | unit A then unit B, B starts at its consensus 5′ end, **fixed** spacer / fixed cut of A, recurring thousands of times | r1 + 39 bp + r3 (7 000+), r5 head cut at ~130 + r6 (4 400) | one element; the bank holds it as two pieces |
| **chance tandem** | A then B head to tail, **variable** cut and gap, rare | r8 rand100's 4 rows | two insertions; mark, don't rebuild |
| **tandem array** | ≥ 3 units at a periodic spacing | tbr MEG-RS, ~910 bp period | copies of a larger repeat (current `[array]`) |
| **piecewise match** | A then B where B starts **mid-consensus** | r10 + r3, r3 + r3 (internal ~50 bp repeat) | one element two consensuses each cover in part — a **consensus defect**, not two copies |

The `h` / `mid` test (does the downstream unit start at its consensus 5′ end?) separates the first
two from the last; the tightness of the spacer and cut separates structural from chance.

## 1. Detection — new step after assignment ("step3c")

`composite_scan.py RUN OUT --flank 400` on every assigned copy (rsi: 66 501 copies, ~25 s):

- units = consensus hits ≥ 40 bp in locus ± 400 bp, both strands, greedy non-overlapping by bits
  (ssearch36 `-z 11`: the library is all homologs; default statistics return nothing);
- each unit tagged `h` (starts ≤ 15 at its consensus 5′), `e` (reaches ≤ 15 of its 3′), else `mid`;
  links `+` ≤ 30 bp, `~` 30–200 bp, `..` further;
- per copy: its **chain** (e.g. `r1[he] ~ r3[he]`) and class single / left half / right half / middle;
- per subfamily: share of each chain, partners, gap and cut distributions, and the **random
  baseline** (copy density × (element length + 30 bp)) — rsi 0.4–1.3 % vs observed up to 40 %.

## 2. Classification rules (general, no per-species tuning)

A chain **L = A → B** is **structural** when all hold:
1. B is `h` (a new copy, not a piecewise continuation);
2. count ≥ max(50, 10 × random expectation for that pair);
3. tight geometry: IQR of the A–B gap ≤ 10 bp **or** IQR of A's cut position ≤ 10 bp
   (rsi r1 → r3: gap 39 / 39 / 42; r5 → r6: cut ~130).

Otherwise head-to-tail pairs are **chance tandems**. A → B with B `mid` (≥ 50 copies) is a
**piecewise match** → reported as a consensus problem of A / B (missing module, internal repeat,
two consensuses covering one element).

## 3. Actions

- **Structural dimer → candidate family.** `build_composite.py` (best 60 copies of that exact chain,
  cut first unit start → last unit end, L-INS-i, majority) gives the full consensus. Shown on the page
  as a *candidate*, never added silently.
- **Optional add-and-check loop** (his decision per candidate): `SINEderella --add`, re-assign,
  re-scan; accept when ≥ 70 % of the candidate's copies now read as **one full unit** (rsi: r5h_r6
  85 % → accept; r10_groupB 83 % → accept; r1_r3 50 % → needs a variant — its r3 part has an
  internal-repeat form, build the second chain `r1[he] ~ r3[h] + r3[e]`).
- **Piecewise matches → consensus repair queue** (feeds the border loop / consensus rebuild).

## 4. Page and plates

- Report: a **Composite** column per subfamily — "% single; main layouts: `r5 head + r6` 33 %,
  …" in plain words, plus the candidate consensuses with their check result.
- Plates: rows marked **`[dimer-L]` / `[dimer-R]`** (structural) or **`[tandem]`** (chance), like
  `[array]`; left out of the continuation decision (a dimer's partner reads as "sequence continues");
  top100 fills with single copies first. `[array]` tightened at the same time: ≥ 3 units **and** a
  consistent period (tbr: ~910 bp), not just "within 50 kb".

## 5. SINE-discriminator

- New verdict flag **`PART_OF_COMPOSITE`**: most copies are one half of a structural dimer → the
  Overall text says "this consensus is one part of a larger element (r1 + 39 bp + r3, 7 000 copies)"
  — the same kind of call as `FRAGMENT_OF_LONGER`, which reads the flank instead.
- The continuation measurement ignores `[dimer-*]` / `[tandem]` rows like `[array]` rows.

## 6. Origin of the parts (optional annotation)

Per consensus and per unit of a candidate: Dfam (`tRNA-*`, `5S`, known SINEs; curated, `--cut_ga`),
Rfam cmscan (7SL, 5S, U6, 7SK …), tRNAscan-SE. On rsi, Dfam was needed — Rfam / tRNAscan missed the
degenerate heads of r2–r4. Environment on therioserver: conda `rnatools` (Infernal, tRNAscan-SE),
`~/refs/smallrna/Rfam.cm`; Dfam via its web API (a local Dfam would make it offline).

## 7. Tests before any real run

- **Toy** (`tests/toy/`): plant into the toy genome (a) a structural dimer A + fixed 39 bp + B
  (60 copies), (b) chance tandems with random cut and gap, (c) a piecewise pair (B starting
  mid-consensus), (d) a periodic array; the scan must return each class, and build_composite must
  recover the planted dimer.
- **Real ground truth:** rsi — r5h_r6, r1_r3, r10_groupB structural; r8 rand100's 4 rows chance;
  r10 + r3 and r3 + r3 piecewise; tbr MEG-RS array.

## Known limits

- ssearch36 `-z 11` is non-deterministic: counts move ~1 % between runs.
- Weakly matching families (MEG in bats, sim_ratio 0.1–0.3) give too few hits to classify.
- Parts absent from the bank show up only as spacer: a long **fixed** spacer (r10_groupB's 105 bp
  middle) means a missing unit, not a linker — the builder captures it anyway because it cuts
  first unit → last unit.
