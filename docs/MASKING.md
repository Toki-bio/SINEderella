# Masking already-found SINE families before another ssearch36 pass

Written 2026-09-27. Implemented the same day (see **Usage**); `sear` and the scanner are unchanged.

## Recommendation

Mask by family, on a working copy of the genome, by rewriting the masked intervals to N. Keep ssearch36. Keep the original genome for extraction and coordinates.

- The mask is one interval set: every locus assigned to an accepted family, all of its subfamilies together. A copy the family consensus missed stays visible, because the mask is loci, not a classification level. Loci of a different family stay visible.
- The working copy has the same headers and the same lengths as the original, so every coordinate `sear` reports is already a coordinate on the original genome. Nothing maps back.
- Supply the mask as BED. A soft-masked genome can be the input too: its lowercase runs become the BED. Internally there is one path.
- Default use: the de novo rescan after known families are found, and tool comparisons ("what else is there once these families are gone"). Not the default for `--add`: `--add` re-votes loci where a new consensus overlaps an old family, and that overlap (MEG-RS inside MEG-RL, for example) is information the assignment and the CONFLICT flag need.
- Do not move the search to blastn, and do not use `ssearch36 -S` as the mask, not even a rebuilt one. The reasons are measured below.

## Measurements

Planted genome, 8 Mb, one sequence. DOM: 1,500 interspersed copies of a 120 bp family, the mask. RARE: 20 copies of a different family sharing DOM's 30 bp head. SUB: 20 copies of DOM with 8 substitutions, outside the mask. All runs `-T 1 -m 8C -E 2 -z 11`. Hits count when they cover half a locus; RARE and SUB were also checked for exact coordinates.

| binary | subject | DOM query | RARE query | SUB query | time per query |
|---|---|---|---|---|---|
| 36.3.8g conda, or current source stock, with or without `-S` | DOM lowercase | 1500 DOM, 20 SUB | 20 RARE | 20 SUB, 1500 DOM | ~37–46 s |
| any | DOM rewritten to N | 0 DOM, 20 SUB | 20 RARE | 20 SUB | ~4 s |
| any | DOM deleted (minus-bank) | 0 DOM | 20 RARE, shifted coordinates | — | ~3 s |
| current source built with `-DDNALIB_LC`, `-S` | DOM lowercase | 4 DOM, **0 SUB** | 20 RARE, exact | 20 SUB, exact; 2 DOM | ~3–4 s |

blastn 2.17.1+, a mutated 120 bp copy as its only subject, 20 replicates, E 10:

| identity | blastn (word 11) | blastn-short | ssearch36 |
|---|---|---|---|
| 100% | 20/20 | 20/20 | 20/20 |
| 85% | 4/20 | 20/20 | 20/20 |
| 75% and below | 0/20 | 17–20/20 down to 45% | 20/20 |

minimap2 (`-k 11 -w 5`) returned nothing for any family: repeat k-mers exceed `mid_occ`, and case is folded.

## What `-S` is in the source

`-S` sets `ext_sq_set`. The DNA alphabet doubles (`ntx` in `upam.h`: uppercase codes 1–17, lowercase 18–34). `init_pamx` (`initfa.c`) fills `pam2[1]` so a lowercase code scores as N, while `pam2[0]` scores the same code as the real base. `dropgsw2.c` builds both profiles.

In the shipped build the library never gets lowercase codes. `nascii` (`uascii.h`) maps `a/c/g/t` to the uppercase codes, and the only call that would rebuild the library map from `ntx` is inside `#ifdef DNALIB_LC` (`initfa.c`), a macro no Makefile defines. Protein libraries call `init_ascii` without that guard, which is why `-S` is a protein feature in practice. `dropfz2.c` says the DNA case was left unfinished. A query that is lowercase end to end is uppercased again with a warning (`comp_lib9.c`, `upper_seq`).

Rebuilt with `-DDNALIB_LC`, the mask engages, and it shows what `-S` was designed for: deciding which library sequences pass the scan. Lowercase scores as N while ssearch ranks library sequences. Each sequence that passes is then aligned again for display, and there the real bases are used ("treated as normal residues for the final alignment display", `fasta36.1`). For a protein database of many short sequences that is the intended behaviour. For a genome window it is the wrong one. In the planted genome, the one library sequence passed the scan because of the unmasked SUB copies. The display pass then chose the masked DOM copies (100% identity over the 93% SUB copies), printed four of them, and the E-values (computed on the masked statistics) cut the list before any SUB copy appeared. Any `sear` window that holds both a masked copy and an unmasked copy of a related sequence can lose the unmasked one the same way.

## Where the time goes

Masking 2% of the bases cut the search tenfold. The time is the known copies being aligned. A perfect 120 bp hit scores above the 8-bit limit, so that library sequence is redone in 16 bits (`dropgsw2.c`), and the repeat-alignment pass then walks every copy. The RARE query in the unmasked genome returned 1,544 lines, 20 of them RARE: its shared head aligned to every DOM copy. N removes that work and the alignments it would have produced; the dynamic programming still visits every base.

## Not on the default path

- BED exclusion after the search. Coordinates are right, but the search still aligns every known copy before the hits are discarded, so it saves nothing. A hit that straddles the edge of the mask is dropped with the known copy.
- Minus-bank (delete the intervals). Fast, but coordinates become a second system. `sine_scan.sh` maps them back only for headers of the form `>scaffold:start-end()`; a plain `getfasta` header is a different convention. `sear -m` stays what it is: several copies from one window, not a family mask.
- `ssearch36 -S`, stock or rebuilt: stock ignores DNA case; rebuilt, it masks library sequences, not alignments.
- blastn: default word 11 loses copies below ~85% identity; blastn-short keeps them only with a short-word seed, the class of miss full Smith–Waterman was chosen to avoid.
- minimap2: no SINE-scale sensitivity and no case.

## Implementation sketch

1. `tools/mask_genome.py <genome> <mask.bed>... -o <work.fa>`: copy of the genome with the BED intervals as N, headers and lengths identical; refuses to write if any length differs.
2. `SINEderella --mask-families <bed | run_dir>`: from a run directory, take the merged hit intervals of the listed accepted families (all subfamilies of each). Write `mask.bed` and `genome.masked.fa` into the new run, point `sear` and `sine_scan.sh` at `genome.masked.fa`, and keep `genome.clean.fa` for extraction.
3. Record in the run manifest: mask source, families masked, masked bp. The report says a mask was used, and which families it removed.
4. Self-test at startup: one known-family copy inside the mask and one outside; the search must return only the outside copy at its original coordinates.

## Implementation plan (proposed 2026-09-27, Claude, after reading the code)

The search and the extraction have to use different files: search on the N-masked copy, extract and report from
`genome.clean.fa`. Because the masked copy keeps every header and length, coordinates need no mapping.

1. **`tools/mask_genome.py GENOME MASK.bed... -o genome.masked.fa`** - writes the genome with the merged BED
   intervals as `N`, every other base (and its case) unchanged. Refuses to finish unless headers are identical, every
   length is identical, and the number of bases turned to N equals the merged BED length. Writes `mask.stats.tsv`
   (intervals, masked bp, fraction of genome). BED chromosome names are sanitised the way `sanitize_fasta` sanitises
   genome headers (`_` -> `@U@`), so a BED written against the original assembly still matches `genome.clean.fa`.
2. **Mask sources** (one code path, both end in `mask.bed`):
   - `--mask-bed FILE` (repeatable): any BED, e.g. RepeatMasker SINE rows, or a lowercase-run BED from a
     soft-masked assembly (`tools/lowercase_to_bed.py`).
   - `--mask-run RUN_DIR[:NAME,NAME...]`: the step1 hit intervals of the named consensuses from
     `genome.clean_step1/all_hits.labeled.bed` (column 4 = consensus), merged; no names = every consensus in that run.
     This is the family-level mask: pass every subfamily name of a family. Hits, not assignments, so a copy is masked
     even if it later failed the vote.
3. **`step1_search_extract.sh`**: one change - `sear -k "$query" "${SEARCH_GENOME:-$GENOME}"`. The later
   `bedtools getfasta -fi "$GENOME"` stays on the unmasked genome, so extracted copies and their flanks are real bases
   even where a neighbouring masked copy sits in the flank.
4. **Orchestrator (`mode_full`)**: after sanitising, if a mask was given, build `mask.bed` and `genome.masked.fa` in the
   run directory, export `SEARCH_GENOME`, and write `MASK_SOURCE`, `MASK_NAMES`, `MASK_BP`, `MASK_FRACTION` into
   `manifest.txt`. Step6 prints one line: which families were masked and how many bp. Off by default. Not offered in
   `--add` (the overlap of old and new families is what the re-vote and CONFLICT need).
5. **De novo scan**: no code change - `TARGET_GENOME=genome.clean.fa sine_scan.sh genome.masked.fa BANK`. The scanner
   already searches `SEARCH_DB` and extracts from `TARGET_GENOME`; its minus-bank check does not fire because the
   headers are unchanged. Document this as the "rescan after known families" recipe.
6. **Test** (`tests/test_mask.sh`, runs in seconds): a synthetic genome with planted copies of family A inside the mask
   and family B outside it; step1 on the masked genome must return every B copy at its planted coordinates and no A
   copy, and `mask_genome.py` must refuse a BED that runs off a sequence end.

What it buys, from the planted test above: a search for a different family went from ~46 s to ~4 s per query when the
masked family was 2 % of the bases, because the time was being spent aligning the known copies, not scanning bases.
What it does not do, by design: find a subfamily hiding inside the masked family - that is the peel's job.

## Usage (implemented 2026-09-27)

```bash
# mask every copy an earlier run's VES consensus hit, then search the genome for other families
SINEderella --mask-run /path/run_2026...:VES  genome.fa  other_families.fa

# several families (pass every subfamily name of each), or any BED, or both
SINEderella --mask-run RUN:r1_9seqs,r2_3seqs,VES --mask-bed extra.bed  genome.fa bank.fa

# a soft-masked assembly as the mask (note: covers every repeat class)
python3 tools/lowercase_to_bed.py genome.fa > lowercase.bed
SINEderella --mask-bed lowercase.bed genome.fa bank.fa

# de novo rescan after known families (no scanner change needed)
python3 tools/mask_genome.py genome.clean.fa mask.bed -o genome.masked.fa
TARGET_GENOME=genome.clean.fa bash sine_scan.sh genome.masked.fa BANK
```

The run directory gets `mask.bed`, `genome.masked.fa`, `mask.stats.tsv`, and `MASK_SOURCE`, `MASK_BP`,
`MASK_FRACTION`, `MASK_INTERVALS` in `manifest.txt`. `--resume` re-uses `genome.masked.fa`. `--add` and
`--exclude` refuse a mask. Test: `tests/test_mask.sh` (synthetic genome, about 20 s): unmasked finds 80/80
masked-family and 60/60 other-family copies; masked finds 0/80 and 60/60 at the planted coordinates.

### On a real genome (2026-09-27, *Tadarida brasiliensis*, 2.3 Gb)
`--mask-run tbr_run:VES` masked 640,519 loci, 136 Mb (6.0 %). Bank Rhin-1 + MEG-RL/RS/T2/TR, THREADS=12.

| query | hits unmasked | hits VES-masked | search time unmasked -> masked |
|---|---|---|---|
| Rhin-1 | 456 | 135 | 69 s -> 51 s |
| MEG-RS | 286 | 274 | ~55 s -> ~46 s |
| MEG-T2 | 90 | 49 | ~51 s -> ~51 s |
| MEG-TR | 40 | 29 | ~52 s -> ~49 s |
| MEG-RL | 12 | 8 | 64 s -> 76 s |

Speed: only queries that would align to many masked copies gain (Rhin-1 shares the tRNA-derived head with VES);
the Smith-Waterman scan over the genome costs the same, so an unrelated query is not faster. Expected larger gain:
a de novo fragment-bank rescan, where many tRNA-derived fragments hit the masked family. Content: 70 % of
Rhin-1 "hits" in *Tadarida* sit on VES copies and vanish; MEG-RS keeps 274 of 286 - not VES in disguise.
