# Fast mode: fasta36 first, ssearch36 on the minus bank

`SINEderella --fast` (or `SEAR_MODE=fast`, or `sear --mode fast`) changes only the
step 1 search engine. Steps 2 and later, and every other ssearch36 call in the
pipeline, are unchanged.

## What it does

`sear` already searches in rounds: a round finds `plus_to_add`, the minus bank is
the part with those loci cut out (`bedtools complement` + `getfasta`), and the
next round searches the minus bank (`--minus`). Fast mode uses that structure
with two engines:

| | regular (default) | fast |
|---|---|---|
| round 1 | `ssearch36 -g -3 -z 11 -E 2 ... 3` | `fasta36 -g -3 -z 11 -E 2 ... 6` (ktup 6) |
| length rule, round 1 | longest single HSP >= 0.8 x query | merged span of HSPs within 50 bp >= 0.8 x query |
| round 2 | only with `--minus` | always, once: `ssearch36 -z 1` on the minus bank |
| rounds 3+ | `--minus` only, ssearch36 `-z 11` | `--minus` only, ssearch36 `-z 1` |

The ssearch36 round is what makes fast mode safe: every locus fasta36 did not
claim is still in the minus bank and ssearch36 sees it with the same flags and
filters as in regular mode. What fasta36 finds, ssearch36 would also have found
(the same loci at 65 % identity and 0.8 length). If fasta36 finds nothing in a
part, the ssearch36 round runs on the whole part.

Why the span rule: fasta36 aligns inside a band of width 16 around the seed
diagonal, so a diverged copy often comes back as two or three HSPs instead of
one. With the single-HSP rule 18 % of loci would fail the length test although
they were found. Hits closer than 50 bp are one locus (`bedtools merge -d 50`),
and the merged span is tested.

Why `-z 1` in the ssearch36 round and not `-z 11`: on the depleted bank both
give the same loci (checked on 104 Mbp: identical) and `-z 11` costs 1.6x
because it shuffles every window. `-z 3` must not be used there: on a bank with
its homologs removed the precomputed statistics admit tens of thousands of
spurious HSPs (128,000 on 104 Mbp).

## Measured (D. valentini, dar_squam1 consensus, 1.65 Gbp, THREADS=4)

| | time | loci | note |
|---|---|---|---|
| `sear` regular | 524 s | 110,493 | |
| `sear --mode fast` | 257 s | 110,710 | covers 110,490 of the 110,493 regular loci (99.997 %); 248 new weak ones |

Times exclude the one-time `seqkit sliding` + `split2` step (25 s), which both modes
share. Fast mode's 257 s is 41 s of fasta36 and about 200 s of ssearch36 on the
minus bank (98 % of the genome); that scan is the floor for any scheme in which
ssearch36 must see every base.

On a 104 Mbp slice: regular 34 s, 6,450 loci; fast 21 s, 6,465 loci, all 6,450
regular loci covered, 15 new weak ones. The rescue round recovers about 1 % of
the loci, almost all at 60-70 % identity.

The gain is larger where the regular search is slow for a different reason:
on hit-rich families (VES in *Taphozous*, 620k copies) ssearch36's time is the
copies being re-aligned for display (docs/MASKING.md), not the scan. In fast
mode fasta36 claims those copies cheaply and the ssearch36 round sees a bank
with them removed, so only the scan term remains.

## What differs in the output

- `<prefix>-<query>.bed`: the same loci within 0.1 %; a few dozen extra weak
  loci per genome that `-z 1` statistics admit at the E <= 2 border.
- `<prefix>-<query>.rawhits.tsv` (read by the satellite stage): in fast mode it
  holds fasta36's round-1 hits plus the ssearch36 round's. fasta36 reports more
  short fragments (< 100 bp) than ssearch36: 2,569 vs 503 on the 104 Mbp slice.
  Of those 2,569, 2,477 lie within 100 bp of a hit the regular raw file also has
  (pieces of an alignment ssearch36 reports in one row), so the kind-A TRF
  windows they seed are the same windows; 92 are new short hits (median 52 bp,
  85 % identity). The >= 100 bp rows are the same in both modes.
- `manifest.txt` records `SEAR_MODE`; step 1's log prints `sear mode: fast`.
- `-z 11` is kept for fasta36 in round 1 and everywhere else in the pipeline.

## Not a mode: fasta36 alone

`SEAR_NO_RESCUE=1` skips the ssearch36 round: about 6x faster than regular
but loses ~1 % of loci, all old copies, which the voting and length-version
steps depend on. Undocumented in `--help` on purpose.

## Requirements

`fasta36` from the same FASTA36 package as `ssearch36`, in PATH. Checked at
start when the mode is fast.
