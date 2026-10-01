# Task C — CpG-adjusted divergence test across species

## 2026-10-01

- Why: docs/BORROWED_TRIAGE.md, section 5 item 2 + decision D2 — measure the
  CpG effect on divergence across species and SINEs before changing anything.
- tools/cpg_div/cpg_div.py (was empty, now implemented): per copy of a plate,
  raw_p, cpg_adj_p (RepeatMasker rule: C->T / G->A at consensus CpG sites
  counts 1/10), k2p_adj (Kimura 2-parameter after the same rule, NA when
  saturated), n_aligned, n_cpg_sites. CLI: cpg_div.py PLATE.aln.fa [--out TSV].
- tools/cpg_div/run_all.py (new): every plate under
  C:/work/Tal/*/alignments/*top100.aln.fa -> tools/cpg_div/summary.tsv
  (species, family, n_copies, consensus CpG share, median raw_p, median
  cpg_adj_p, Spearman raw vs adjusted, share of copies changing 5-bin rank
  bin); also echoes the TSV to stdout.
- docs/CPG_DIVERGENCE_TEST.md (new): method, decision rule and limitations
  written; ALL numbers PENDING until run_all.py has been run (no invented
  numbers).
- Pending: run the task C check (must print PASS); run run_all.py; fill the
  numbers in docs/CPG_DIVERGENCE_TEST.md from summary.tsv.
