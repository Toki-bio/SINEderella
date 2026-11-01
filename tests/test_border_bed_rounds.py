"""update_bed_5prime / update_bed_3prime over several rounds: from round 2 on the input bed is the output bed; before
2026-11-01 opening it for writing emptied it first (rsi MEG-RS: empty 3' bed, border loop crashed). python3 tests/test_border_bed_rounds.py"""
import os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "publish"))
import flank_border_iterate as F, border_loop_subfam as B
w = tempfile.mkdtemp()
open(os.path.join(w, "genome.sizes"), "w").write("c1\t10000\n")
open(os.path.join(w, "in.bed"), "w").write("c1\t1000\t1200\tx\t0\t+\nc1\t3000\t3200\ty\t0\t-\n")
b = F.update_bed_5prime(os.path.join(w, "in.bed"), 10, w); b = F.update_bed_5prime(b, 10, w)
c = B.update_bed_3prime(b, 5, w); c = B.update_bed_3prime(c, 5, w)
got = [l.split("\t")[:3] for l in open(c).read().splitlines()]
want = [["c1", "980", "1210"], ["c1", "2990", "3220"]]
assert got == want, (got, want)
print("test_border_bed_rounds: passed")
