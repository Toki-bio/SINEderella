# Strand-aware merge of the hits (step 1) - 2026-10-13

**Problem.** Step 1 merged the hits of all consensuses with `bedtools merge -c 4,5,6 -o max,max,distinct`. The strand was only collected: a hit on "-" touching or
overlapping a hit of another family on "+" became one interval with strand "+,-", and `bedtools getfasta -s` (which flips only a pure "-") wrote it forward. Result:
the "-" element came out reverse-complemented, in a piece that holds two elements. Sicista (run_20261006_143837): 29,295 of 2,177,501 loci (1.3 %), 24,570 of 1,248,280 DIP-best-vote copies (2.0 %);
they filled the bottom chunks of the SubFam alignment (172 of the 300 copies in chunks 595-600 reverse, 98 % of them "+,-"). The other outliers of those chunks were not strand problems:
56 of 128 forward ones had no DIP hit at all (53 rejected_low_bitscore copies of RSINE1, IDL-Geo, B4, STRIDM, B2 ... whose best vote is DIP - an effect of taking "all copies" incl. rejected), 37 are diverged DIP
(< 80 % identity to the consensus), 30 are DIP plus another element on the same strand (350-380 bp pieces), 5 normal/truncated.

**Fix.** `tools/merge_hits.py` (default; `MERGE_STRAND_AWARE=0` = old merge) resolves each cluster of overlapping/book-ended hits:
single strand = the old result; both strands: hits of each strand are merged first; two opposite-strand intervals overlapping by >= 50 % of the shorter one are ONE locus seen from both
strands (one interval, strand of the highest-scoring hit); otherwise they are different elements and each keeps its interval and strand. The strand column is always "+" or "-".
Used in `step1_search_extract.sh` and in the orchestrator's `rebuild_from_searches` (add / exclude). The tool is copied into the run's `tools/`.

**Checked on the real Sicista hits (same 58 bed files):** cluster count identical to the old merge (2,177,501); all 2,148,206 pure-strand lines byte-identical; of the 29,295 mixed clusters 2,403 merged as one locus and
26,892 split into 54,440 pieces (2,205,049 intervals instead of 2,177,501). 600 resolved pieces containing a DIP hit of their own strand, extracted with `getfasta -s`: 601 of 601 forward against the DIP consensus; 600 of the
old "+,-" intervals with a DIP hit: 284 forward, 317 reverse. Toy (hs21 Alu, 10,537 clusters): 28 mixed, 0 one locus, 28 split into 56; the full pipeline with `--families` runs through. `tests/test_merge_hits.py` (7 tests); whole suite 96 passed.

**Not done.** The published Sicista alignments (firm and all-copies DIP, B1, B4) and assignments still come from the old merge; a rerun of step 1 from the existing searches plus step 2 (or an incremental step 2 for the changed loci) is needed to
get clean copies. Copies are still counted per piece, so a split locus counts as two copies (one per element), as it should.
