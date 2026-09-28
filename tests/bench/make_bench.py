#!/usr/bin/env python3
"""make_bench.py OUT - tiny benchmark sets from real published plates (row 2 = consensus, copies'
header coordinates = extracted region incl. the base flanks). For each case: base.fa (as step8a round 0)
and long.fa (+600 bp both sides, as the last continuation round)."""
import os, re, subprocess, sys

OUT = sys.argv[1]
ONLY = sys.argv[2:]
CASES = [  # name, plate, genome
    ("cseRhin", "~/chiro/cse/run_add_20260927_172935/results/alignments/cse_Rhin-1_top100.aln.fa", "~/chiro/cse/run_add_20260927_172935/genome.clean.fa"),
    ("cseVES", "~/chiro/cse/run_add_20260927_172935/results/alignments/cse_VES_top100.aln.fa", "~/chiro/cse/run_add_20260927_172935/genome.clean.fa"),
    ("cseMEGRS", "~/chiro/cse/run_add_20260927_172935/results/alignments/cse_MEG-RS_top100.aln.fa", "~/chiro/cse/run_add_20260927_172935/genome.clean.fa"),
    ("cseMEGT2", "~/chiro/cse/run_add_20260927_172935/results/alignments/cse_MEG-T2_top100.aln.fa", "~/chiro/cse/run_add_20260927_172935/genome.clean.fa"),
    ("cseMEGRL", "~/chiro/cse/run_add_20260927_172935/results/alignments/cse_MEG-RL_top100.aln.fa", "~/chiro/cse/run_add_20260927_172935/genome.clean.fa"),
    ("rsiR1", "~/tmp/bw/new/rsi/rsi_r1_9seqs_top100.aln.fa", "~/rhin/rsi/run_add_20260927_180847/genome.clean.fa"),
]


def read_fa(p):
    n, s, cur = [], [], []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n:
                s.append("".join(cur))
            n.append(l[1:]); cur = []
        else:
            cur.append(l)
    s.append("".join(cur))
    return n, s


for name, plate, genome in CASES:
    if ONLY and name not in ONLY:
        continue
    plate, genome = os.path.expanduser(plate), os.path.expanduser(genome)
    if not (os.path.exists(plate) and os.path.exists(genome)):
        print("SKIP", name, "missing", plate if not os.path.exists(plate) else genome); continue
    d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
    n, s = read_fa(plate)
    cons = re.sub(r"[^A-Za-z]", "", s[1]).upper()          # row 2 = original
    bed = []
    for h in n[2:]:
        m = re.match(r"(_R_)?(\S+):(\d+)-(\d+)\(([+-])\)", h)
        if m:
            bed.append((m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)))
    for tag, ext in (("base", 0), ("long", 600)):
        with open(os.path.join(d, tag + ".bed"), "w") as fh:
            for i, (c, a, b, st) in enumerate(bed):
                fh.write("%s\t%d\t%d\tc%d\t0\t%s\n" % (c.replace("_", "@U@"), max(0, a - ext), b + ext, i, st))  # genome.clean.fa: _ -> @U@
        fa = subprocess.run(["bedtools", "getfasta", "-fi", genome, "-bed", os.path.join(d, tag + ".bed"), "-s"],
                            capture_output=True, text=True).stdout
        with open(os.path.join(d, tag + ".fa"), "w") as fh:
            fh.write(">CONSENSUS_%s\n%s\n" % (name, cons) + fa)
    print(name, "copies", len(bed), "cons_len", len(cons))
